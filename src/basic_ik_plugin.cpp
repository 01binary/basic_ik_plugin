#include <Eigen/Geometry>

using namespace Eigen;

MatrixXd inverseKinematics(Matrix4d pose)
{
  MatrixXd angles(3, 1);
  angles << 0.0, 0.0, 0.0;
  return angles;
}
