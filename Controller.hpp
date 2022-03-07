#pragma once

#include <Eigen/Core>
#include "dart/dart.hpp"
#include "MPCSolver.hpp"
#include "Utility.hpp"
#include  <fstream>
#include "utils.cpp"
#include "types.hpp"
//#include "FootstepPlan.hpp"
#include "StateFiltering.hpp"
#include <memory>
#include <random>

class Controller
{
public:

  // fundamental methods
  Controller(dart::dynamics::SkeletonPtr _robot, dart::simulation::WorldPtr _world);
  virtual ~Controller();
  void update();
  void setInitialConfiguration();
  
  // various methods of the controller class
  Eigen::Vector3d getZmpFromExternalForces();
  Eigen::MatrixXd getTorsoAndSwfJacobian();
  
  void ArmSwing();
  void FixWaistToChestJoints();

  Eigen::VectorXd getJointVelocitiesStacked(Eigen::VectorXd, Eigen::VectorXd,
		  Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd);
  Eigen::VectorXd getJointVelocitiesStacked_worldframe(Eigen::VectorXd, Eigen::VectorXd,
		  Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd);

  Eigen::Vector3d getRPY(dart::dynamics::BodyNode*, dart::dynamics::BodyNode*);
  Eigen::Vector3d getRPY(dart::dynamics::BodyNode*);

  void storeData();
  void drawArrow(Eigen::Vector3d);
  void createForceLine();

private:

  // dart related pointers
  dart::dynamics::SkeletonPtr mRobot;
  dart::dynamics::BodyNode* mTorso;
  dart::simulation::WorldPtr mWorld;
  dart::dynamics::BodyNode* mLeftFoot;
  dart::dynamics::BodyNode* mRightFoot;
  dart::dynamics::BodyNode* mSupportFoot;
  dart::dynamics::BodyNode* mSwingFoot;
  dart::dynamics::BodyNode* mBase;
  dart::dynamics::LineSegmentShapePtr mForceLine;
    
  // some flages
  bool balancePoint;
  bool TORSO = false;
  bool COM = true;

  // some useful parameters
  int wait = 300;

  // matrices and vectors
  Eigen::MatrixXd ftsp_and_time;
  Eigen::VectorXd dq;
  Eigen::MatrixXd A_obs;
  Eigen::VectorXd B_obs;
  Eigen::MatrixXd G_obs;
  Eigen::MatrixXd C_obs;
  Eigen::VectorXd measured_state_x;
  Eigen::VectorXd measured_state_y;   
  Eigen::Vector3d previousCoM, currentCoM;
  Eigen::MatrixXd _taskGain, _R_CoM, _R_swgf;
  Eigen::MatrixXd _Rot_buffer;
  Eigen::VectorXd desired_pos;  
  Eigen::VectorXd actual_pos;  
  Eigen::VectorXd ComVref;
  Eigen::MatrixXd Jacobian_tot, PseudoJacobian_tot;
  const Eigen::Vector4d DefaultForceLineColor
    = Eigen::Vector4d(1.0, 0.0, 0.0, 1.0);
      
  // MPC class
  MPCSolver* solver;

  // data structures to store the state
  State desired;
  State current;
  WalkState walkState;

};
