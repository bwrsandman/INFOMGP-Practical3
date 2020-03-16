#ifndef MESH_HEADER_FILE
#define MESH_HEADER_FILE

#include <Eigen/Core>
#include <Eigen/Sparse>

#include "constraints.h"

/// the class the contains each individual rigid objects and their functionality
class Mesh {
public:
  Mesh(const Eigen::VectorXd& _origPositions, const Eigen::MatrixXi& boundF, const Eigen::MatrixXi& _T, const int _globalOffset, const double _youngModulus, const double _poissonRatio, const double _density, const bool _isFixed, const Eigen::RowVector3d& userCOM, const Eigen::RowVector4d& userOrientation);
  virtual ~Mesh();
  /// Quick-reject checking collision between mesh bounding boxes.
  bool isBoxCollide(const Mesh& m2);
  static bool isNeighborTets(const Eigen::RowVector4i& tet1, const Eigen::RowVector4i& tet2);
  /// this function creates all collision constraints between vertices of the two meshes
  void createCollisionConstraints(const Mesh& m, const bool sameMesh, const double timeStep, const double CRCoeff, std::vector<Constraint>& activeConstraints);
  void createGlobalMatrices(const double timeStep, const double _alpha, const double _beta);
  /// returns center of mass
  Eigen::Vector3d initializeVolumesAndMasses();
  /// performing the integration step of the soft body.
  void integrateVelocity(double timeStep);
  /// Update the current position with the integrated velocity
  void integratePosition(double timeStep);
  /// the full integration for the time step (velocity + position)
  void integrate(double timeStep);

  // position
  Eigen::VectorXd origPositions;     ///< 3|V|x1 original vertex positions in xyzxyz format - never change this!
  Eigen::VectorXd currPositions;     ///< 3|V|x1 current vertex positions in xyzxyz format
  // kinematics
  bool isFixed;                      ///< is the object immobile (infinite mass)
  Eigen::VectorXd currVelocities;    ///< 3|V|x1 velocities per coordinate in xyzxyz format.
  Eigen::MatrixXi T;                 ///< |T|x4 tetrahdra
  Eigen::MatrixXi F;                 ///< |F|x3 boundary faces
  Eigen::VectorXd invMasses;         ///< |V|x1 inverse masses of vertices, computed in the beginning as 1.0/(density * vertex voronoi area)
  Eigen::VectorXd voronoiVolumes;    ///< |V|x1 the voronoi volume of vertices
  Eigen::VectorXd tetVolumes;        ///< |T|x1 tetrahedra volumes
  int globalOffset;                  ///< the global index offset of the of opositions/velocities/impulses from the beginning of the global coordinates array in the containing scene class
  Eigen::VectorXi boundTets;         ///< just the boundary tets, for collision
  double youngModulus, poissonRatio, density, alpha, beta;
  Eigen::SparseMatrix<double> A, K, M, D;   ///< The soft-body matrices
  /// the solver for the left-hand side matrix constructed for FEM
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double>>* ASolver;
};

#endif  // MESH_HEADER_FILE
