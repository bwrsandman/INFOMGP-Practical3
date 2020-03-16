#include "scene.h"

#include <fstream>
#include <igl/readOFF.h>
#include <igl/readMESH.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>

Scene::Scene() = default;
Scene::~Scene() = default;

void Scene::global2Mesh() {
  for (int i = 0; i < meshes.size(); i++) {
    meshes[i].currPositions << globalPositions.segment(meshes[i].globalOffset, meshes[i].currPositions.size());
    meshes[i].currVelocities << globalVelocities.segment(meshes[i].globalOffset, meshes[i].currVelocities.size());
  }
}

void Scene::mesh2global() {
  for (int i = 0; i < meshes.size(); i++) {
    globalPositions.segment(meshes[i].globalOffset, meshes[i].currPositions.size()) << meshes[i].currPositions;
    globalVelocities.segment(meshes[i].globalOffset, meshes[i].currVelocities.size()) << meshes[i].currVelocities;
  }
}

void Scene::initScene(double timeStep, double alpha, double beta) {
  for (int i = 0; i < meshes.size(); i++) {
    if (!meshes[i].isFixed)
      meshes[i].createGlobalMatrices(timeStep, alpha, beta);
  }
  mesh2global();
}

void Scene::updateScene(double timeStep, double CRCoeff, double tolerance, int maxIterations) {
  /*******************1. Integrating velocity and position from external and internal forces************************************/
  for (int i = 0; i < meshes.size(); i++)
    meshes[i].integrate(timeStep);

  mesh2global();


  /*******************2. Creating and Aggregating constraints************************************/

  std::vector<Constraint> activeConstraints;

  //user constraints
  activeConstraints.insert(activeConstraints.end(), userConstraints.begin(), userConstraints.end());

  //barrier constraints
  //activeConstraints.insert(activeConstraints.end(), barrierConstraints.begin(), barrierConstraints.end());

  //collision constraints
  for (int i = 0; i < meshes.size(); i++)
    for (int j = i + 1; j < meshes.size(); j++)
      meshes[i].createCollisionConstraints(meshes[j], i == j, timeStep, CRCoeff, activeConstraints);


  /*******************3. Resolving velocity constraints iteratively until the velocities are valid************************************/
  int currIteration = 0;
  int zeroStreak =
      0;  //how many consecutive constraints are already below tolerance without any change; the algorithm stops if all are.
  int currConstIndex = 0;
  while (zeroStreak < activeConstraints.size() && (currIteration < maxIterations * globalPositions.size())) {

    Constraint currConstraint = activeConstraints[currConstIndex];

    Eigen::VectorXd constraintPositions(currConstraint.globalIndices.size());
    Eigen::VectorXd constraintVelocities(currConstraint.globalIndices.size());
    for (int i = 0; i < currConstraint.globalIndices.size(); i++) {
      constraintPositions(i) = globalPositions(currConstraint.globalIndices(i));
      constraintVelocities(i) = globalVelocities(currConstraint.globalIndices(i));
    }

    //generating impulses
    Eigen::VectorXd generatedImpulses;
    if (currConstraint.resolveVelocityConstraint(constraintPositions,
                                                 constraintVelocities,
                                                 generatedImpulses,
                                                 tolerance))
      zeroStreak++;
    else
      zeroStreak = 0;

    //cout<<"zeroStreak: "<<zeroStreak;

    //correcting velocities according to impulses
    for (int i = 0; i < currConstraint.globalIndices.size(); i++) {
      int currIndex = currConstraint.globalIndices(i);
      globalVelocities(currIndex) += globalInvMasses(currIndex) * generatedImpulses(i);
      //TODO: velocity bias
      //if (timeStep>tolerance)
      //  rawImpulses(fullConstraints[currConstraint].particleIndices(i))+=CRCoeff*currDiff(i)/timeStep;
    }
    currIteration++;
    currConstIndex = (currConstIndex + 1) % (activeConstraints.size());
  }

  global2Mesh();

  /*******************4. Solving for position drift************************************/

  mesh2global();

  currIteration = 0;
  zeroStreak =
      0;  //how many consecutive constraints are already below tolerance without any change; the algorithm stops if all are.
  currConstIndex = 0;
  while (zeroStreak < activeConstraints.size() && (currIteration < maxIterations * globalPositions.size())) {

    Constraint currConstraint = activeConstraints[currConstIndex];

    Eigen::VectorXd constraintPositions(currConstraint.globalIndices.size());
    Eigen::VectorXd constraintVelocities(currConstraint.globalIndices.size());
    for (int i = 0; i < currConstraint.globalIndices.size(); i++) {
      constraintPositions(i) = globalPositions(currConstraint.globalIndices(i));
      constraintVelocities(i) = globalVelocities(currConstraint.globalIndices(i));
    }

    //generating impulses
    Eigen::VectorXd generatedPosDiffs;
    if (currConstraint.resolvePositionConstraint(constraintPositions,
                                                 constraintVelocities,
                                                 generatedPosDiffs,
                                                 tolerance))
      zeroStreak++;
    else
      zeroStreak = 0;

    //cout<<"zeroStreak: "<<zeroStreak<<endl;

    //correcting velocities according to impulses
    for (int i = 0; i < currConstraint.globalIndices.size(); i++) {
      int currIndex = currConstraint.globalIndices(i);
      globalPositions(currIndex) += generatedPosDiffs(i);
    }
    currIteration++;
    currConstIndex = (currConstIndex + 1) % (activeConstraints.size());
  }
  global2Mesh();
}

#if 0
void Scene::setPlatformBarriers(const MatrixXd& platV, const double CRCoeff) {

  RowVector3d minPlatform = platV.colwise().minCoeff();
  RowVector3d maxPlatform = platV.colwise().maxCoeff();

  //y value of maxPlatform is lower bound
  for (int i = 1; i < globalPositions.size(); i += 3) {
    VectorXi coordIndices(1);
    coordIndices(0) = i;
    VectorXd constraintInvMasses(1);
    constraintInvMasses(0) = globalInvMasses(i);
    barrierConstraints.push_back(Constraint(BARRIER,
                                            INEQUALITY,
                                            coordIndices,
                                            constraintInvMasses,
                                            MatrixXd::Zero(1, 1),
                                            maxPlatform(1),
                                            CRCoeff));
  }
}
#endif

void Scene::addMesh(const Eigen::MatrixXd& V, const Eigen::MatrixXi& boundF, const Eigen::MatrixXi& T, double youngModulus, double PoissonRatio, double density, bool isFixed, const Eigen::RowVector3d& userCOM, const Eigen::RowVector4d& userOrientation) {
  Eigen::VectorXd Vxyz(3 * V.rows());
  for (int i = 0; i < V.rows(); i++)
    Vxyz.segment(3 * i, 3) = V.row(i).transpose();

  // std::cout << "Vxyz: " << Vxyz << std::endl;
  Mesh m(Vxyz, boundF, T, globalPositions.size(), youngModulus, PoissonRatio, density, isFixed, userCOM, userOrientation);
  meshes.push_back(m);
  int oldTsize = globalT.rows();
  globalT.conservativeResize(globalT.rows() + T.rows(), 4);
  globalT.block(oldTsize, 0, T.rows(), 4) = T.array() + globalPositions.size() / 3;  //to offset T to global index
  globalPositions.conservativeResize(globalPositions.size() + Vxyz.size());
  globalVelocities.conservativeResize(globalPositions.size());
  int oldIMsize = globalInvMasses.size();
  globalInvMasses.conservativeResize(globalPositions.size());
  for (int i = 0; i < m.invMasses.size(); i++)
    globalInvMasses.segment(oldIMsize + 3 * i, 3) = Eigen::Vector3d::Constant(m.invMasses(i));

  mesh2global();
}

bool Scene::loadScene(const std::string& dataFolder, const std::string& sceneFileName, const std::string& constraintFileName) {
  std::ifstream sceneFileHandle;
  std::ifstream constraintFileHandle;
  sceneFileHandle.open(dataFolder + std::string("/") + sceneFileName);
  if (!sceneFileHandle.is_open())
    return false;

  constraintFileHandle.open(dataFolder + std::string("/") + constraintFileName);
  if (!constraintFileHandle.is_open())
    return false;
  int numofObjects, numofConstraints;

  currTime = 0;
  sceneFileHandle >> numofObjects;
  for (int i = 0; i < numofObjects; i++) {
    Eigen::MatrixXi objT, objF;
    Eigen::MatrixXd objV;
    std::string MESHFileName;
    bool isFixed;
    double youngModulus, poissonRatio, density;
    Eigen::RowVector3d userCOM;
    Eigen::RowVector4d userOrientation;
    sceneFileHandle >> MESHFileName >> density >> youngModulus >> poissonRatio >> isFixed >> userCOM(0) >> userCOM(1)
                    >> userCOM(2) >> userOrientation(0) >> userOrientation(1) >> userOrientation(2)
                    >> userOrientation(3);
    userOrientation.normalize();
    //if the mesh is an OFF file, tetrahedralize it
    if (MESHFileName.find(".off") != std::string::npos) {
      Eigen::MatrixXd VOFF;
      Eigen::MatrixXi FOFF;
      igl::readOFF(dataFolder + std::string("/") + MESHFileName, VOFF, FOFF);
      Eigen::RowVectorXd mins = VOFF.colwise().minCoeff();
      Eigen::RowVectorXd maxs = VOFF.colwise().maxCoeff();
      for (int k = 0; k < VOFF.rows(); k++)
        VOFF.row(k) << 25.0 * (VOFF.row(k) - mins).array() / (maxs - mins).array();

      if (!isFixed)
        igl::copyleft::tetgen::tetrahedralize(VOFF, FOFF, "pq1.1", objV, objT, objF);
      else
        igl::copyleft::tetgen::tetrahedralize(VOFF, FOFF, "pq1.414Y", objV, objT, objF);
    } else {
      igl::readMESH(dataFolder + std::string("/") + MESHFileName, objV, objT, objF);
    }

    //fixing weird orientation problem
    Eigen::MatrixXi tempF(objF.rows(), 3);
    tempF << objF.col(2), objF.col(1), objF.col(0);
    objF = tempF;

    //cout<<"objF: "<<objF<<endl;
    //cout<<"viewerF: "<<viewerF<<endl;
    addMesh(objV, objF, objT, youngModulus, poissonRatio, density, isFixed, userCOM, userOrientation);
  }

  //reading intra-mesh attachment constraints
  constraintFileHandle >> numofConstraints;
  for (int i = 0; i < numofConstraints; i++) {
    int attachM1, attachM2, attachV1, attachV2;
    constraintFileHandle >> attachM1 >> attachV1 >> attachM2 >> attachV2;

    Eigen::VectorXi coordIndices(6);
    coordIndices << meshes[attachM1].globalOffset + 3 * attachV1,
        meshes[attachM1].globalOffset + 3 * attachV1 + 1,
        meshes[attachM1].globalOffset + 3 * attachV1 + 2,
        meshes[attachM2].globalOffset + 3 * attachV2,
        meshes[attachM2].globalOffset + 3 * attachV2 + 1,
        meshes[attachM2].globalOffset + 3 * attachV2 + 2;

    Eigen::VectorXd constraintInvMasses(6);
    constraintInvMasses << meshes[attachM1].invMasses(attachV1),
        meshes[attachM1].invMasses(attachV1),
        meshes[attachM1].invMasses(attachV1),
        meshes[attachM2].invMasses(attachV2),
        meshes[attachM2].invMasses(attachV2),
        meshes[attachM2].invMasses(attachV2);
    double refValue = (meshes[attachM1].currPositions.segment(3 * attachV1, 3)
        - meshes[attachM2].currPositions.segment(3 * attachV2, 3)).norm();
    userConstraints.push_back(Constraint(DISTANCE,
                                         EQUALITY,
                                         coordIndices,
                                         constraintInvMasses,
                                         Eigen::VectorXd::Zero(0),
                                         refValue,
                                         0.0));
  }
  return true;
}
