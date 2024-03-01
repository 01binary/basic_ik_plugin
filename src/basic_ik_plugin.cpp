#include <cmath>
#include <pluginlib/class_list_macros.h>
#include <eigen_conversions/eigen_msg.h>
#include <tf2_eigen/tf2_eigen.h>
#include "basic_ik_plugin.h"

using namespace std;
using namespace kinematics;
using namespace moveit::core;
using namespace robot_state;
using namespace moveit_msgs;
using namespace Eigen;
using namespace basic;

IKPlugin::IKPlugin() :
    m_pPlanningGroup(NULL),
    m_pBaseJoint(NULL),
    m_pShoulderJoint(NULL),
    m_pElbowJoint(NULL),
    m_pWristJoint(NULL)
{
}

MatrixXd IKPlugin::solve(Matrix4d pose) const
{
  // Normal
  double nx = pose(0, 0);
  double ny = pose(1, 0);
  double nz = pose(2, 0);

  // Orientation
  double ox = pose(0, 1);
  double oy = pose(1, 1);
  double oz = pose(2, 1);

  // Approach
  double ax = pose(0, 2);
  double ay = pose(1, 2);
  double az = pose(2, 2);

  // Position
  double px = pose(0, 3);
  double py = pose(1, 3);
  double pz = pose(2, 3);

  // Solve base
  double base = atan2(-ox,oy);

  // Solve elbow
  double elbow = asin(
    (
      (ax*ny-ay*nx)*(nx*oy-ny*ox-nx*py*2.0E+1+ny*px*2.0E+1)*3.5E+2
    )
    /
    (
      (ax*ax)*(ny*ny)*7.744E+3 +
      (ay*ay)*(nx*nx)*7.744E+3 +
      (ax*ax)*(oy*oy)*2.5E+1 +
      (ay*ay)*(ox*ox)*2.5E+1 +
      (ax*ax)*(py*py)*1.0E+4 +
      (ay*ay)*(px*px)*1.0E+4 +
      (nx*nx)*(oy*oy)*2.5E+1 +
      (ny*ny)*(ox*ox)*2.5E+1 +
      (nx*nx)*(py*py)*1.0E+4 +
      (ny*ny)*(px*px)*1.0E+4 -
      (ay*ay)*ox*px*1.0E+3 -
      (ax*ax)*oy*py*1.0E+3 -
      (ny*ny)*ox*px*1.0E+3 -
      (nx*nx)*oy*py*1.0E+3 +
      (ay*ay)*nx*ox*8.8E+2 +
      (ax*ax)*ny*oy*8.8E+2 -
      (ay*ay)*nx*px*1.76E+4 -
      (ax*ax)*ny*py*1.76E+4 -
      ax*ay*nx*oy*8.8E+2 -
      ax*ay*ny*ox*8.8E+2 +
      ax*ay*nx*py*1.76E+4 +
      ax*ay*ny*px*1.76E+4 -
      ax*ay*ox*oy*5.0E+1 +
      ax*ay*ox*py*1.0E+3 +
      ax*ay*oy*px*1.0E+3 -
      ax*ay*px*py*2.0E+4 -
      nx*ny*ox*oy*5.0E+1 +
      nx*ny*ox*py*1.0E+3 +
      nx*ny*oy*px*1.0E+3 -
      nx*ny*px*py*2.0E+4 -
      ax*ay*nx*ny*1.5488E+4
    )
  );

  // Solve shoulder
  double shoulder = atan2(
    (
      (
        ax *
        sqrt(
          pow(ax*ny-ay*nx, 2.0) *
          pow(nx*oy-ny*ox-nx*py*2.0E+1+ny*px*2.0E+1, 2.0) *
          1.0 / pow(
            (ax*ax)*(ny*ny)*7.744E+3 +
            (ay*ay)*(nx*nx)*7.744E+3 +
            (ax*ax)*(oy*oy)*2.5E+1 +
            (ay*ay)*(ox*ox)*2.5E+1 +
            (ax*ax)*(py*py)*1.0E+4 +
            (ay*ay)*(px*px)*1.0E+4 +
            (nx*nx)*(oy*oy)*2.5E+1 +
            (ny*ny)*(ox*ox)*2.5E+1 +
            (nx*nx)*(py*py)*1.0E+4 +
            (ny*ny)*(px*px)*1.0E+4 -
            (ay*ay)*ox*px*1.0E+3 -
            (ax*ax)*oy*py*1.0E+3 -
            (ny*ny)*ox*px*1.0E+3 -
            (nx*nx)*oy*py*1.0E+3 +
            (ay*ay)*nx*ox*8.8E+2 +
            (ax*ax)*ny*oy*8.8E+2 -
            (ay*ay)*nx*px*1.76E+4 -
            (ax*ax)*ny*py*1.76E+4 -
            ax*ay*nx*oy*8.8E+2 -
            ax*ay*ny*ox*8.8E+2 +
            ax*ay*nx*py*1.76E+4 +
            ax*ay*ny*px*1.76E+4 -
            ax*ay*ox*oy*5.0E+1 +
            ax*ay*ox*py*1.0E+3 +
            ax*ay*oy*px*1.0E+3 -
            ax*ay*px*py*2.0E+4 -
            nx*ny*ox*oy*5.0E+1 +
            nx*ny*ox*py*1.0E+3 +
            nx*ny*oy*px*1.0E+3 -
            nx*ny*px*py*2.0E+4 -
            ax*ay*nx*ny*1.5488E+4,
            2.0
          )
          *
          -1.225E+5+1.0
        )
        -
        (nx*(ax*ny-ay*nx)*(nx*oy-ny*ox-nx*py*2.0E+1+ny*px*2.0E+1)*3.5E+2)
        /
        (
          (ax*ax)*(ny*ny)*7.744E+3 +
          (ay*ay)*(nx*nx)*7.744E+3 +
          (ax*ax)*(oy*oy)*2.5E+1 +
          (ay*ay)*(ox*ox)*2.5E+1 +
          (ax*ax)*(py*py)*1.0E+4 +
          (ay*ay)*(px*px)*1.0E+4 +
          (nx*nx)*(oy*oy)*2.5E+1 +
          (ny*ny)*(ox*ox)*2.5E+1 +
          (nx*nx)*(py*py)*1.0E+4 +
          (ny*ny)*(px*px)*1.0E+4 -
          (ay*ay)*ox*px*1.0E+3 -
          (ax*ax)*oy*py*1.0E+3 -
          (ny*ny)*ox*px*1.0E+3 -
          (nx*nx)*oy*py*1.0E+3 +
          (ay*ay)*nx*ox*8.8E+2 +
          (ax*ax)*ny*oy*8.8E+2 -
          (ay*ay)*nx*px*1.76E+4 -
          (ax*ax)*ny*py*1.76E+4 -
          ax*ay*nx*oy*8.8E+2 -
          ax*ay*ny*ox*8.8E+2 +
          ax*ay*nx*py*1.76E+4 +
          ax*ay*ny*px*1.76E+4 -
          ax*ay*ox*oy*5.0E+1 +
          ax*ay*ox*py*1.0E+3 +
          ax*ay*oy*px*1.0E+3 -
          ax*ay*px*py*2.0E+4 -
          nx*ny*ox*oy*5.0E+1 +
          nx*ny*ox*py*1.0E+3 +
          nx*ny*oy*px*1.0E+3 -
          nx*ny*px*py*2.0E+4 -
          ax*ay*nx*ny*1.5488E+4
        )
      )
      *
      sqrt(ox*ox+oy*oy)
    )
    /
    oy,
    (
      (
        nx *
        sqrt(
          pow(ax*ny-ay*nx, 2.0) *
          pow(nx*oy-ny*ox-nx*py*2.0E+1+ny*px*2.0E+1, 2.0) *
          1.0 / pow(
            (ax*ax)*(ny*ny)*7.744E+3 +
            (ay*ay)*(nx*nx)*7.744E+3 +
            (ax*ax)*(oy*oy)*2.5E+1 +
            (ay*ay)*(ox*ox)*2.5E+1 +
            (ax*ax)*(py*py)*1.0E+4 +
            (ay*ay)*(px*px)*1.0E+4 +
            (nx*nx)*(oy*oy)*2.5E+1 +
            (ny*ny)*(ox*ox)*2.5E+1 +
            (nx*nx)*(py*py)*1.0E+4 +
            (ny*ny)*(px*px)*1.0E+4 -
            (ay*ay)*ox*px*1.0E+3 -
            (ax*ax)*oy*py*1.0E+3 -
            (ny*ny)*ox*px*1.0E+3 -
            (nx*nx)*oy*py*1.0E+3 +
            (ay*ay)*nx*ox*8.8E+2 +
            (ax*ax)*ny*oy*8.8E+2 -
            (ay*ay)*nx*px*1.76E+4 -
            (ax*ax)*ny*py*1.76E+4 -
            ax*ay*nx*oy*8.8E+2 -
            ax*ay*ny*ox*8.8E+2 +
            ax*ay*nx*py*1.76E+4 +
            ax*ay*ny*px*1.76E+4 -
            ax*ay*ox*oy*5.0E+1 +
            ax*ay*ox*py*1.0E+3 +
            ax*ay*oy*px*1.0E+3 -
            ax*ay*px*py*2.0E+4 -
            nx*ny*ox*oy*5.0E+1 +
            nx*ny*ox*py*1.0E+3 +
            nx*ny*oy*px*1.0E+3 -
            nx*ny*px*py*2.0E+4 -
            ax*ay*nx*ny*1.5488E+4,
            2.0
          )
          *
          -1.225E+5+1.0
        ) +
        (ax*(ax*ny-ay*nx)*(nx*oy-ny*ox-nx*py*2.0E+1+ny*px*2.0E+1)*3.5E+2)
        /
        (
          (ax*ax)*(ny*ny)*7.744E+3 +
          (ay*ay)*(nx*nx)*7.744E+3 +
          (ax*ax)*(oy*oy)*2.5E+1 +
          (ay*ay)*(ox*ox)*2.5E+1 +
          (ax*ax)*(py*py)*1.0E+4 +
          (ay*ay)*(px*px)*1.0E+4 +
          (nx*nx)*(oy*oy)*2.5E+1 +
          (ny*ny)*(ox*ox)*2.5E+1 +
          (nx*nx)*(py*py)*1.0E+4 +
          (ny*ny)*(px*px)*1.0E+4 -
          (ay*ay)*ox*px*1.0E+3 -
          (ax*ax)*oy*py*1.0E+3 -
          (ny*ny)*ox*px*1.0E+3 -
          (nx*nx)*oy*py*1.0E+3 +
          (ay*ay)*nx*ox*8.8E+2 +
          (ax*ax)*ny*oy*8.8E+2 -
          (ay*ay)*nx*px*1.76E+4 -
          (ax*ax)*ny*py*1.76E+4 -
          ax*ay*nx*oy*8.8E+2 -
          ax*ay*ny*ox*8.8E+2 +
          ax*ay*nx*py*1.76E+4 +
          ax*ay*ny*px*1.76E+4 -
          ax*ay*ox*oy*5.0E+1 +
          ax*ay*ox*py*1.0E+3 +
          ax*ay*oy*px*1.0E+3 -
          ax*ay*px*py*2.0E+4 -
          nx*ny*ox*oy*5.0E+1 +
          nx*ny*ox*py*1.0E+3 +
          nx*ny*oy*px*1.0E+3 -
          nx*ny*px*py*2.0E+4 -
          ax*ay*nx*ny*1.5488E+4
        )
      )
      *
      sqrt(ox*ox+oy*oy)
    )/oy
  );

  double wrist = atan((az-sqrt(az*az+oz*oz))/oz)*-2.0;

  MatrixXd angles(4, 1);
  angles << base, shoulder, elbow, wrist;
  return angles;
}

bool IKPlugin::initialize(
  const RobotModel &robot_model,
  const string &group_name,
  const string &base_frame,
  const vector<string> &tip_frames,
  double search_discretization)
{
  ROS_INFO("Basic IK Plugin Initializing");

  m_pPlanningGroup = robot_model.getJointModelGroup(group_name);
  auto joints = m_pPlanningGroup->getJointModels();
  auto jointNames = m_pPlanningGroup->getJointModelNames();

  for (auto pJoint: joints)
  {
    if (pJoint->getType() == JointModel::REVOLUTE)
    {
      auto limits = pJoint->getVariableBoundsMsg().front();
      m_joints.push_back(pJoint);
    }
    else if (pJoint->getType() == JointModel::PRISMATIC)
    {
      auto limits = pJoint->getVariableBoundsMsg().front();
      m_joints.push_back(pJoint);
    }
  }

  vector<string> chainTips;
  auto chains = m_pPlanningGroup->getConfig().chains_;
  auto chainTip = chains.front().second;
  chainTips.push_back(chainTip);

  KinematicsBase::storeValues(
    robot_model,
    group_name,
    base_frame,
    chainTips,
    search_discretization);

  m_pBaseJoint = getJoint(JointModel::REVOLUTE);
  m_pShoulderJoint = getJoint(JointModel::REVOLUTE, m_pBaseJoint);
  m_pElbowJoint = getJoint(JointModel::REVOLUTE, m_pShoulderJoint);
  m_pWristJoint = getJoint(JointModel::REVOLUTE, m_pElbowJoint);

  // Initialize state
  m_pState.reset(new robot_state::RobotState(robot_model_));
  m_pState->setToDefaultValues();

  return true;
}

bool IKPlugin::supportsGroup(
  const JointModelGroup *jmg, string *error_text_out) const
{
  return true;
}

const vector<string> &IKPlugin::getJointNames() const
{
  return m_pPlanningGroup->getJointModelNames();
}

const vector<string> &IKPlugin::getLinkNames() const
{
  return m_pPlanningGroup->getLinkModelNames();
}

bool IKPlugin::getPositionFK(
  const vector<string> &link_names,
  const vector<double> &joint_angles,
  vector<geometry_msgs::Pose> &poses) const
{
  return false;
}

bool IKPlugin::getPositionIK(
  const geometry_msgs::Pose &ik_pose,
  const vector<double> &ik_seed_state,
  vector<double> &solution,
  MoveItErrorCodes &error_code,
  const KinematicsQueryOptions &options) const
{
  return searchPositionIK(
    ik_pose,
    ik_seed_state,
    1.0,
    solution,
    error_code,
    options);
}

bool IKPlugin::searchPositionIK(
  const geometry_msgs::Pose &ik_pose,
  const vector<double> &ik_seed_state,
  double timeout,
  vector<double> &solution,
  MoveItErrorCodes &error_code,
  const KinematicsQueryOptions &options) const
{
  const IKCallbackFn solution_callback = NULL;
  vector<double> consistency_limits;
  vector<geometry_msgs::Pose> poses;
  poses.push_back(ik_pose);

  return searchPositionIK(
    poses,
    ik_seed_state,
    timeout,
    consistency_limits,
    solution,
    solution_callback,
    error_code,
    options);
}

bool IKPlugin::searchPositionIK(
  const geometry_msgs::Pose &ik_pose,
  const vector<double> &ik_seed_state,
  double timeout,
  const vector<double> &consistency_limits,
  vector<double> &solution,
  MoveItErrorCodes &error_code,
  const KinematicsQueryOptions &options) const
{
  const IKCallbackFn solution_callback = NULL;

  return searchPositionIK(
    ik_pose,
    ik_seed_state,
    timeout,
    consistency_limits,
    solution,
    solution_callback,
    error_code,
    options);
}

bool IKPlugin::searchPositionIK(
  const geometry_msgs::Pose &ik_pose,
  const vector<double> &ik_seed_state,
  double timeout,
  vector<double> &solution,
  const IKCallbackFn &solution_callback,
  MoveItErrorCodes &error_code,
  const KinematicsQueryOptions &options) const
{
  vector<double> consistency_limits;

  return searchPositionIK(
    ik_pose,
    ik_seed_state,
    timeout,
    consistency_limits,
    solution,
    solution_callback,
    error_code,
    options);
}

bool IKPlugin::searchPositionIK(
  const geometry_msgs::Pose &ik_pose,
  const vector<double> &ik_seed_state,
  double timeout,
  const vector<double> &consistency_limits,
  vector<double> &solution,
  const IKCallbackFn &solution_callback,
  MoveItErrorCodes &error_code,
  const KinematicsQueryOptions &options) const
{
  vector<geometry_msgs::Pose> poses;
  poses.push_back(ik_pose);

  return searchPositionIK(
    poses,
    ik_seed_state,
    timeout,
    consistency_limits,
    solution,
    solution_callback,
    error_code,
    options);
}

bool IKPlugin::searchPositionIK(
  const vector<geometry_msgs::Pose> &ik_poses,
  const vector<double> &ik_seed_state,
  double timeout,
  const vector<double> &consistency_limits,
  vector<double> &solution,
  const IKCallbackFn &solution_callback,
  MoveItErrorCodes &error_code,
  const KinematicsQueryOptions &options,
  const robot_state::RobotState *context_state) const
{
  // Initialize solution
  solution.resize(m_joints.size());

  // Solve inverse kinematics
  Vector3d origin = getOrigin();
  Matrix4d goal = getGoal(ik_poses, origin);
  MatrixXd angles = solve(goal);

  // Return solution
  setJointState(m_pBaseJoint, angles(BASE, 0), solution);
  setJointState(m_pShoulderJoint, angles(SHOULDER, 0), solution);
  setJointState(m_pElbowJoint, angles(ELBOW, 0), solution);
  setJointState(m_pWristJoint, angles(WRIST, 0), solution);

  error_code.val = error_code.SUCCESS;

  if (!solution_callback.empty())
  {
      solution_callback(ik_poses.front(), solution, error_code);
  }

  return true;
}

const JointModel* IKPlugin::getJoint(
  JointModel::JointType type, const JointModel* parent) const
{
  for (auto joint : m_joints)
  {
    if (parent)
    {
      if (joint->getType() == type && joint->getParentLinkModel() == parent->getChildLinkModel())
        return joint;
    }
    else
    {
      if (joint->getType() == type)
        return joint;
    }
  }

  return nullptr;
}

Matrix4d IKPlugin::getGoal(
  const std::vector<geometry_msgs::Pose>& ik_poses,
  const Eigen::Vector3d& origin) const
{
  Isometry3d goal;
  auto goalPose = ik_poses.back();
  tf2::fromMsg(goalPose, goal);

  Matrix4d goalMatrix = goal.matrix();
  goalMatrix(0, 3) = goalMatrix(0, 3) - origin.x();
  goalMatrix(1, 3) = goalMatrix(1, 3) - origin.y();
  goalMatrix(2, 3) = goalMatrix(2, 3) - origin.z();

  return goalMatrix;
}

Vector3d IKPlugin::getGoalPosition(
  const vector<geometry_msgs::Pose>& ik_poses,
  const Eigen::Vector3d& origin) const
{
  Isometry3d goal;
  auto goalPose = ik_poses.back();
  tf2::fromMsg(goalPose, goal);

  Vector3d pos = goal.translation();
  pos.x() = pos.x() - origin.x();
  pos.y() = pos.y() - origin.y();
  pos.z() = pos.z() - origin.z();

  return pos;
}

Vector3d IKPlugin::getOrigin() const
{
  return m_pState->getGlobalLinkTransform(m_pBaseJoint->getChildLinkModel())
    .translation();
}

Isometry3d IKPlugin::setJointState(
  const JointModel* pJoint,
  double angle,
  std::vector<double>& states) const
{
  const Vector3d& axis = getJointAxis(pJoint);
  const JointLimits& limits = pJoint->getVariableBoundsMsg().front();

  double jointState = isnan(angle)
    ? limits.max_position
    : clamp(angle, limits.min_position, limits.max_position);

  size_t index = find(m_joints.begin(), m_joints.end(), pJoint) - m_joints.begin();
  states[index] = jointState;

  m_pState->setJointPositions(pJoint, &jointState);
  m_pState->enforceBounds(pJoint);

  return m_pState->getJointTransform(pJoint);
}

const Vector3d& IKPlugin::getJointAxis(const JointModel* pJoint)
{
  if (pJoint->getType() == JointModel::REVOLUTE)
  {
    return dynamic_cast<const RevoluteJointModel*>(pJoint)->getAxis();
  }
  else
  {
    return dynamic_cast<const PrismaticJointModel*>(pJoint)->getAxis();
  }
}

PLUGINLIB_EXPORT_CLASS(basic::IKPlugin, kinematics::KinematicsBase);
