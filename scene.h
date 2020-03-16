#ifndef SCENE_HEADER_FILE
#define SCENE_HEADER_FILE

#include <vector>
#include <Eigen/Core>

#include "mesh.h"

/// This class contains the entire scene operations, and the engine time loop.
class Scene {
public:
  double currTime;

  Eigen::VectorXd globalPositions;   ///< 3*|V| all positions
  Eigen::VectorXd globalVelocities;  ///< 3*|V| all velocities
  Eigen::VectorXd globalInvMasses;   ///< 3*|V| all inverse masses  (NOTE: the invMasses in the Mesh class is |v| (one per vertex)!
  Eigen::MatrixXi globalT;           ///< |T|x4 tetraheda in global index

  std::vector<Mesh> meshes;

  std::vector<Constraint> userConstraints;     ///< provided from the scene
  std::vector<Constraint> barrierConstraints;  ///< provided by the platform

  Scene();
  ~Scene();
  /// Updates from global values back into mesh values
  void global2Mesh();
  /// Update from mesh current values into global values
  void mesh2global();
  /// This should be called whenever the timestep changes
  void initScene(double timeStep, double alpha, double beta);
  /*********************************************************************
   This function handles a single time step
   1. Integrating velocities and position from forces and previous impulses
   2. detecting collisions and generating collision constraints, alongside with given user constraints
   3. Resolving constraints iteratively by updating velocities until the system is valid (or maxIterations has passed)
   *********************************************************************/
  void updateScene(double timeStep, double CRCoeff, double tolerance, int maxIterations);
  // void setPlatformBarriers(const MatrixXd& platV, const double CRCoeff)
  /// Adding an object.
  void addMesh(const Eigen::MatrixXd& V, const Eigen::MatrixXi& boundF, const Eigen::MatrixXi& T, double youngModulus, double PoissonRatio,  double density, bool isFixed, const Eigen::RowVector3d& userCOM, const Eigen::RowVector4d& userOrientation);
  /// loading a scene from the scene .txt files
  /// you do not need to update this function
  bool loadScene(const std::string& dataFolder, const std::string& sceneFileName, const std::string& constraintFileName);
};

#endif
