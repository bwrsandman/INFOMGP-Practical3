#include "mesh.h"

#include <list>
#include <iostream>
#include <Eigen/Eigen>

#include "auxfunctions.h"
#include "ccd.h"

/*****************************Auxiliary functions for collision detection. Do not need updating********************************/

/** Support function for libccd*/
void support(const void* _obj, const ccd_vec3_t *_d, ccd_vec3_t *_p) {
  // assume that obj_t is user-defined structure that holds info about
  // object (in this case box: x, y, z, pos, quat - dimensions of box,
  // position and rotation)
  //std::cout<<"calling support"<<std::endl;
  auto obj = reinterpret_cast<const Eigen::MatrixXd*>(_obj);
  Eigen::RowVector3d p;
  Eigen::RowVector3d d;
  for (int i = 0; i < 3; i++)
    d(i) = _d->v[i]; //p(i)=_p->v[i];


  d.normalize();
  //std::cout<<"d: "<<d<<std::endl;

  Eigen::RowVector3d objCOM = obj->colwise().mean();
  int maxVertex = -1;
  int maxDotProd = -32767.0;
  for (int i = 0; i < obj->rows(); i++) {
    double currDotProd = d.dot(obj->row(i) - objCOM);
    if (maxDotProd < currDotProd) {
      maxDotProd = currDotProd;
      //std::cout<<"maxDotProd: "<<maxDotProd<<std::endl;
      maxVertex = i;
    }

  }
  //std::cout<<"maxVertex: "<<maxVertex<<std::endl;

  for (int i = 0; i < 3; i++)
    _p->v[i] = (*obj)(maxVertex, i);

  //std::cout<<"end support"<<std::endl;
}

void stub_dir(const void *obj1, const void *obj2, ccd_vec3_t *dir) {
  dir->v[0] = 1.0;
  dir->v[1] = 0.0;
  dir->v[2] = 0.0;
}

void center(const void* obj, ccd_vec3_t *center) {
  Eigen::RowVector3d objCOM = reinterpret_cast<const Eigen::MatrixXd *>(obj)->colwise().mean();
  for (int i = 0; i < 3; i++)
    center->v[i] = objCOM(i);
}

Mesh::Mesh(const Eigen::VectorXd& _origPositions, const Eigen::MatrixXi& boundF, const Eigen::MatrixXi& _T, const int _globalOffset, const double _youngModulus, const double _poissonRatio, const double _density, const bool _isFixed, const Eigen::RowVector3d& userCOM, const Eigen::RowVector4d& userOrientation) {
  origPositions = _origPositions;
  //cout<<"original origPositions: "<<origPositions<<endl;
  T = _T;
  F = boundF;
  isFixed = _isFixed;
  globalOffset = _globalOffset;
  density = _density;
  poissonRatio = _poissonRatio;
  youngModulus = _youngModulus;
  currVelocities = Eigen::VectorXd::Zero(origPositions.rows());

  Eigen::VectorXd naturalCOM = initializeVolumesAndMasses();
  //cout<<"naturalCOM: "<<naturalCOM<<endl;


  origPositions -= naturalCOM.replicate(origPositions.rows() / 3,
                                        1);  //removing the natural COM of the OFF file (natural COM is never used again)
  //cout<<"after natrualCOM origPositions: "<<origPositions<<endl;

  for (int i = 0; i < origPositions.size(); i += 3)
    origPositions.segment(i, 3)
        << (QRot(origPositions.segment(i, 3).transpose(), userOrientation) + userCOM).transpose();

  currPositions = origPositions;

  if (isFixed)
    invMasses.setZero();

  //finding boundary tets
  Eigen::VectorXi boundVMask(origPositions.rows() / 3);
  boundVMask.setZero();
  for (int i = 0; i < boundF.rows(); i++)
    for (int j = 0; j < 3; j++)
      boundVMask(boundF(i, j)) = 1;

  std::cout << "boundVMask.sum(): " << boundVMask.sum() << std::endl;

  std::vector<int> boundTList;
  for (int i = 0; i < T.rows(); i++) {
    int incidence = 0;
    for (int j = 0; j < 4; j++)
      incidence += boundVMask(T(i, j));
    if (incidence > 2)
      boundTList.push_back(i);
  }

  boundTets.resize(boundTList.size());
  for (int i = 0; i < boundTets.size(); i++)
    boundTets(i) = boundTList[i];

  ASolver = NULL;
}

Mesh::~Mesh() {
  if (ASolver != NULL)
    delete ASolver;
}

bool Mesh::isBoxCollide(const Mesh& m2) {
  Eigen::RowVector3d XMin1 = Eigen::RowVector3d::Constant(3276700.0);
  Eigen::RowVector3d XMax1 = Eigen::RowVector3d::Constant(-3276700.0);
  Eigen::RowVector3d XMin2 = Eigen::RowVector3d::Constant(3276700.0);
  Eigen::RowVector3d XMax2 = Eigen::RowVector3d::Constant(-3276700.0);
  for (int i = 0; i < origPositions.size(); i += 3) {
    XMin1 = XMin1.array().min(currPositions.segment(i, 3).array().transpose());
    XMax1 = XMax1.array().max(currPositions.segment(i, 3).array().transpose());
  }
  for (int i = 0; i < m2.origPositions.size(); i += 3) {
    XMin2 = XMin2.array().min(m2.currPositions.segment(i, 3).array().transpose());
    XMax2 = XMax2.array().max(m2.currPositions.segment(i, 3).array().transpose());
  }

#if 0
  double rmax1 = vertexSphereRadii.maxCoeff();
  double rmax2 = m2.vertexSphereRadii.maxCoeff();
  XMin1.array() -= rmax1;
  XMax1.array() += rmax1;
  XMin2.array() -= rmax2;
  XMax2.array() += rmax2;
#endif

  // checking all axes for non-intersection of the dimensional interval
  for (int i = 0; i < 3; i++)
    if ((XMax1(i) < XMin2(i)) || (XMax2(i) < XMin1(i)))
      return false;

  return true;  // all dimensional intervals are overlapping = possible intersection
}

bool Mesh::isNeighborTets(const Eigen::RowVector4i& tet1, const Eigen::RowVector4i& tet2) {
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      if (tet1(i) == tet2(j)) //shared vertex
        return true;

  return false;
}

void Mesh::createCollisionConstraints(const Mesh& m, const bool sameMesh, const double timeStep, const double CRCoeff, std::vector<Constraint>& activeConstraints) {
  //collision between bounding boxes
  if (!isBoxCollide(m))
    return;

  if ((isFixed && m.isFixed))  //collision does nothing
    return;

#if 0
  // creating tet spheres
  MatrixXd c1(T.rows(), 3);
  MatrixXd c2(m.T.rows(), 3);
  VectorXd r1(T.rows());
  VectorXd r2(m.T.rows());
#endif

  Eigen::MatrixXd maxs1(boundTets.rows(), 3);
  Eigen::MatrixXd mins1(boundTets.rows(), 3);
  Eigen::MatrixXd maxs2(m.boundTets.rows(), 3);
  Eigen::MatrixXd mins2(m.boundTets.rows(), 3);

  for (int i = 0; i < boundTets.size(); i++) {
    Eigen::MatrixXd tet1(4, 3);
    tet1 << currPositions.segment(3 * T(boundTets(i), 0), 3).transpose(),
        currPositions.segment(3 * T(boundTets(i), 1), 3).transpose(),
        currPositions.segment(3 * T(boundTets(i), 2), 3).transpose(),
        currPositions.segment(3 * T(boundTets(i), 3), 3).transpose();

    //c1.row(i) = tet1.colwise().mean();
    //r1(i) = ((c1.row(i).replicate(4, 1) - tet1).rowwise().norm()).maxCoeff();
    mins1.row(i) = tet1.colwise().minCoeff();
    maxs1.row(i) = tet1.colwise().maxCoeff();

  }

  for (int i = 0; i < m.boundTets.size(); i++) {
    Eigen::MatrixXd tet2(4, 3);
    tet2 << m.currPositions.segment(3 * m.T(m.boundTets(i), 0), 3).transpose(),
        m.currPositions.segment(3 * m.T(m.boundTets(i), 1), 3).transpose(),
        m.currPositions.segment(3 * m.T(m.boundTets(i), 2), 3).transpose(),
        m.currPositions.segment(3 * m.T(m.boundTets(i), 3), 3).transpose();

    // c2.row(i) = tet2.colwise().mean();
    // r2(i) = ((c2.row(i).replicate(4, 1) - tet2).rowwise().norm()).maxCoeff();
    mins2.row(i) = tet2.colwise().minCoeff();
    maxs2.row(i) = tet2.colwise().maxCoeff();
  }

  //checking collision between every tetrahedrons
  std::list<Constraint> collisionConstraints;
  for (int i = 0; i < boundTets.size(); i++) {
    for (int j = 0; j < m.boundTets.size(); j++) {

      //not checking for collisions between tetrahedra neighboring to the same vertices
      if (sameMesh)
        if (isNeighborTets(T.row(boundTets(i)), m.T.row(m.boundTets(j))))
          continue;  //not creating collisions between neighboring tets

      bool overlap = true;
      for (int k = 0; k < 3; k++)
        if ((maxs1(i, k) < mins2(j, k)) || (maxs2(j, k) < mins1(i, k)))
          overlap = false;

      if (!overlap)
        continue;

      Eigen::VectorXi globalCollisionIndices(24);
      Eigen::VectorXd globalInvMasses(24);
      for (int t = 0; t < 4; t++) {
        globalCollisionIndices.segment(3 * t, 3) << globalOffset + 3 * (T(boundTets(i), t)), globalOffset
            + 3 * (T(boundTets(i), t)) + 1, globalOffset + 3 * (T(boundTets(i), t)) + 2;
        globalInvMasses.segment(3 * t, 3) << invMasses(T(boundTets(i), t)), invMasses(T(boundTets(i), t)), invMasses(T(
            boundTets(i),
            t));
        globalCollisionIndices.segment(12 + 3 * t, 3) << m.globalOffset + 3 * m.T(m.boundTets(j), t), m.globalOffset
            + 3 * m.T(m.boundTets(j), t) + 1, m.globalOffset + 3 * m.T(m.boundTets(j), t) + 2;
        globalInvMasses.segment(12 + 3 * t, 3) << m.invMasses(m.T(m.boundTets(j), t)), m.invMasses(m.T(m.boundTets(j),
                                                                                                       t)), m.invMasses(
            m.T(m.boundTets(j), t));
      }

      ccd_t ccd;
      CCD_INIT(&ccd);
      ccd.support1 = support; // support function for first object
      ccd.support2 = support; // support function for second object
      ccd.center1 = center;
      ccd.center2 = center;

      ccd.first_dir = stub_dir;
      ccd.max_iterations = 100;     // maximal number of iterations

      Eigen::MatrixXd tet1(4, 3);
      tet1 << currPositions.segment(3 * T(boundTets(i), 0), 3).transpose(),
          currPositions.segment(3 * T(boundTets(i), 1), 3).transpose(),
          currPositions.segment(3 * T(boundTets(i), 2), 3).transpose(),
          currPositions.segment(3 * T(boundTets(i), 3), 3).transpose();

      Eigen::MatrixXd tet2(4, 3);
      tet2 << m.currPositions.segment(3 * m.T(m.boundTets(j), 0), 3).transpose(),
          m.currPositions.segment(3 * m.T(m.boundTets(j), 1), 3).transpose(),
          m.currPositions.segment(3 * m.T(m.boundTets(j), 2), 3).transpose(),
          m.currPositions.segment(3 * m.T(m.boundTets(j), 3), 3).transpose();

      void *obj1 = (void *) &tet1;
      void *obj2 = (void *) &tet2;

      ccd_real_t _depth;
      ccd_vec3_t dir, pos;

      int nonintersect = ccdMPRPenetration(obj1, obj2, &ccd, &_depth, &dir, &pos);

      if (nonintersect)
        continue;

      Eigen::Vector3d intNormal, intPosition;
      double depth;
      for (int k = 0; k < 3; k++) {
        intNormal(k) = dir.v[k];
        intPosition(k) = pos.v[k];
      }

      depth = _depth;
      intPosition -= depth * intNormal / 2.0;

      Eigen::Vector3d p1 = intPosition + depth * intNormal;
      Eigen::Vector3d p2 = intPosition;

      //getting barycentric coordinates of each point

      Eigen::MatrixXd PMat1(4, 4);
      PMat1 << 1.0, currPositions.segment(3 * T(boundTets(i), 0), 3).transpose(),
          1.0, currPositions.segment(3 * T(boundTets(i), 1), 3).transpose(),
          1.0, currPositions.segment(3 * T(boundTets(i), 2), 3).transpose(),
          1.0, currPositions.segment(3 * T(boundTets(i), 3), 3).transpose();
      PMat1.transposeInPlace();

      Eigen::Vector4d rhs1;
      rhs1 << 1.0, p1;

      Eigen::Vector4d B1 = PMat1.inverse() * rhs1;

      Eigen::MatrixXd PMat2(4, 4);
      PMat2 << 1.0, m.currPositions.segment(3 * m.T(m.boundTets(j), 0), 3).transpose(),
          1.0, m.currPositions.segment(3 * m.T(m.boundTets(j), 1), 3).transpose(),
          1.0, m.currPositions.segment(3 * m.T(m.boundTets(j), 2), 3).transpose(),
          1.0, m.currPositions.segment(3 * m.T(m.boundTets(j), 3), 3).transpose();
      PMat2.transposeInPlace();

      Eigen::Vector4d rhs2;
      rhs2 << 1.0, p2;

      Eigen::Vector4d B2 = PMat2.inverse() * rhs2;

      //cout<<"B1: "<<B1<<endl;
      //cout<<"B2: "<<B2<<endl;

      //Matrix that encodes the vector between interpenetration points by the c
      Eigen::MatrixXd v2cMat1(3, 12);
      v2cMat1.setZero();
      for (int k = 0; k < 3; k++) {
        v2cMat1(k, k) = B1(0);
        v2cMat1(k, 3 + k) = B1(1);
        v2cMat1(k, 6 + k) = B1(2);
        v2cMat1(k, 9 + k) = B1(3);
      }

      Eigen::MatrixXd v2cMat2(3, 12);
      v2cMat2.setZero();
      for (int k = 0; k < 3; k++) {
        v2cMat2(k, k) = B2(0);
        v2cMat2(k, 3 + k) = B2(1);
        v2cMat2(k, 6 + k) = B2(2);
        v2cMat2(k, 9 + k) = B2(3);
      }

      Eigen::MatrixXd v2dMat(3, 24);
      v2dMat << -v2cMat1, v2cMat2;
      Eigen::VectorXd constVector = intNormal.transpose() * v2dMat;

      //cout<<"intNormal: "<<intNormal<<endl;
      //cout<<"n*(p2-p1): "<<intNormal.dot(p2-p1)<<endl;
      collisionConstraints.emplace_back(COLLISION,
                                        INEQUALITY,
                                        globalCollisionIndices,
                                        globalInvMasses,
                                        constVector,
                                        0,
                                        CRCoeff);

      //i=10000000;
      //break;

    }
  }

  activeConstraints.insert(activeConstraints.end(), collisionConstraints.begin(), collisionConstraints.end());
}

void Mesh::createGlobalMatrices(const double timeStep, const double _alpha, const double _beta) {
  /*
   * TODO: create the M, D, K matrices from alpha, beta, poisson ratio, and
   *       Young's modulus as learnt in class. Afterward create the matrix "A"
   *       with the given timeStep that is the left hand side of the entire
   *       system.
   */

  A = M + D * timeStep + K * (timeStep * timeStep);

  //Should currently fail since A is empty
  if (ASolver == nullptr)
    ASolver = new Eigen::SimplicialLLT<Eigen::SparseMatrix<double>>();
  ASolver->analyzePattern(A);
  ASolver->factorize(A);

}

Eigen::Vector3d Mesh::initializeVolumesAndMasses() {
  // TODO: compute tet volumes and allocate to vertices
  tetVolumes.conservativeResize(T.rows());
  voronoiVolumes.conservativeResize(origPositions.size() / 3);
  voronoiVolumes.setZero();
  invMasses.conservativeResize(origPositions.size() / 3);
  Eigen::Vector3d COM;
  COM.setZero();
  for (int i = 0; i < T.rows(); i++) {
    Eigen::Vector3d e01 = origPositions.segment(3 * T(i, 1), 3) - origPositions.segment(3 * T(i, 0), 3);
    Eigen::Vector3d e02 = origPositions.segment(3 * T(i, 2), 3) - origPositions.segment(3 * T(i, 0), 3);
    Eigen::Vector3d e03 = origPositions.segment(3 * T(i, 3), 3) - origPositions.segment(3 * T(i, 0), 3);
    Eigen::Vector3d tetCentroid = (origPositions.segment(3 * T(i, 0), 3) + origPositions.segment(3 * T(i, 1), 3)
        + origPositions.segment(3 * T(i, 2), 3) + origPositions.segment(3 * T(i, 3), 3)) / 4.0;
    tetVolumes(i) = std::abs(e01.dot(e02.cross(e03))) / 6.0;
    for (int j = 0; j < 4; j++)
      voronoiVolumes(T(i, j)) += tetVolumes(i) / 4.0;

    COM += tetVolumes(i) * tetCentroid;
  }

  COM.array() /= tetVolumes.sum();
  for (int i = 0; i < origPositions.size() / 3; i++)
    invMasses(i) = 1.0 / (voronoiVolumes(i) * density);

  return COM;
}

void Mesh::integrateVelocity(double timeStep) {
  if (isFixed)
    return;

  /*
   * TODO: construct rhs (right-hand side) and use ASolver->solve(rhs)
   *       to solve for velocities
   */

  Eigen::VectorXd rhs = Eigen::VectorXd::Zero(currVelocities.size());  //REMOVE THIS! it's a stub
  currVelocities = ASolver->solve(rhs);
}

void Mesh::integratePosition(double timeStep) {
  if (isFixed)
    return;  // a fixed object is immobile

  currPositions += currVelocities * timeStep;
  // std::cout << "currPositions: " << currPositions << std::endl;
}

void Mesh::integrate(double timeStep) {
  integrateVelocity(timeStep);
  integratePosition(timeStep);
}
