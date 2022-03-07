#include "Controller.hpp"
#include "parameters.cpp"

using namespace std;
using namespace dart::dynamics;
using namespace dart::simulation;

static std::ofstream foutDebug(realpath("../data/debug.txt", NULL), std::ofstream::trunc);

Controller::Controller(dart::dynamics::SkeletonPtr _robot, dart::simulation::WorldPtr _world)
	: mRobot(_robot), mWorld(_world)
{
    // some useful pointers to robot limbs
    mLeftFoot = mRobot->getBodyNode("l_sole");
    mRightFoot = mRobot->getBodyNode("r_sole");
    mBase = mRobot->getBodyNode("base_link");
    mTorso = mRobot->getBodyNode("body");

    // Initialize walk state
    walkState.mpcIter = 0;
    walkState.controlIter = 0;
    walkState.footstepCounter = 0;
    walkState.supportFoot = true;

    balancePoint = TORSO;

    Eigen::Vector3d comInitialPosition;
    comInitialPosition = mBase->getCOM();
    
    // Desired footsteps and timings in world frame

    int N_footsteps = 40;
    ftsp_and_time = Eigen::MatrixXd::Zero(N_footsteps,4);
    int start = 1;
    for (int i = start; i < N_footsteps; i++) {
         ftsp_and_time(i,0) = (i-start)*0.2; //0.2
         ftsp_and_time(i,1) = pow(-(double)1,(i-start))*0.09;  //0.08
         ftsp_and_time(i,2) = 0.0;
         ftsp_and_time(i,3) = (double)((double)mpcTimeStep/(double)controlTimeStep)*(S+F)*i;  // time: BAD IMPLEMENTATION
    }
 
    std::cout << "Footstep plan - x, y, z, time" << std::endl;
    std::cout << ftsp_and_time << std::endl;

    bool firstSupportFootIsLeft = true;
    Eigen::VectorXd leftFootPose(6), rightFootPose(6);
    leftFootPose << getRPY(mLeftFoot), mLeftFoot->getCOM();
    rightFootPose << getRPY(mRightFoot), mRightFoot->getCOM();
    
    // Instantiate MPC solver
    const Eigen::MatrixXd& ftsp_and_time_ref = ftsp_and_time;
    solver = new MPCSolver(ftsp_and_time_ref);
    current.comPos = mBase->getCOM();
    current.comVel = mBase->getCOMLinearVelocity();
    current.comAcc = mBase->getCOMLinearAcceleration();
    desired.comPos = Eigen::Vector3d(current.comPos(0), current.comPos(1), comTargetHeight);
    desired.comPos << -midrange_x/(eta*eta),-midrange_y/(eta*eta), comTargetHeight;
    desired.comVel = Eigen::Vector3d::Zero();
    desired.zmpPos = Eigen::Vector3d(current.comPos(0), current.comPos(1),0.0);
    desired.leftFootPos = Eigen::Vector3d(0.0, 0.08, 0.0);
    desired.rightFootPos = Eigen::Vector3d(0.0, -0.08, 0.0);
    desired.torsoOrient = Eigen::Vector3d(0.0, 0.0, 0.0);
    desired.leftFootOrient = Eigen::Vector3d(0.0, 0.0, 0.0);
    desired.rightFootOrient = Eigen::Vector3d(0.0, 0.0, 0.0);
    desired.comAcc = eta * eta * (desired.comPos - desired.zmpPos);
    current.disturbance << midrange_x, midrange_y, 0.0;
    previousCoM = current.comPos;


    double ch = cosh(eta*mpcTimeStep);
    double sh = sinh(eta*mpcTimeStep);
    A_obs = Eigen::MatrixXd::Zero(4,4);
    B_obs = Eigen::VectorXd::Zero(4);
    G_obs = Eigen::MatrixXd::Zero(4,2);
    C_obs = Eigen::MatrixXd::Zero(2,4);
    C_obs(0,0) = 1.0;
    C_obs(1,2) = 1.0;
    A_obs<<ch,sh/eta,1-ch,0.0, 
           eta*sh,ch,-eta*sh,mpcTimeStep,
           0.0,0.0,1.0,0.0,
           0.0,0.0,0.0,1.0;
    B_obs << mpcTimeStep-sh/eta,1-ch,mpcTimeStep,0.0;
    _R_CoM = Eigen::MatrixXd::Identity(6,6);
    _R_swgf = Eigen::MatrixXd::Identity(6,6); 
    _Rot_buffer = Eigen::MatrixXd::Identity(3,3);  
    desired.zmpDot = Eigen::VectorXd::Zero(2);
    desired.disturbance << midrange_x, midrange_y, 0.0;
    desired.obsState_x << 0.0, 0.0, 0.0, desired.disturbance(0);
    desired.obsState_y << 0.0, 0.0, 0.0, desired.disturbance(1);
    G_obs <<      0.90131,     0.00000,
                 26.18387,     0.03266,
                 -0.00000,     0.50000,
                239.94774,     6.53262;
                
    std::cout << "disturbance range is: " << std::endl;
    std::cout << "SIMULATION IS READY, PRESS THE SPACEBAR TO START!" << std::endl;
    std::cout << "   " << std::endl;    
}

Controller::~Controller() {}



void Controller::update() {

    // stop the simulation at last_simulation_frame
    int last_simulation_frame = 3000;
    if (mWorld->getSimFrames() == last_simulation_frame) exit(1);

    // start measuring the control routine execution time (just for check)
    auto start = std::chrono::high_resolution_clock::now();
    
    /*
    if (mWorld->getSimFrames() > wait-10 ){
       mBase->addExtForce(Eigen::Vector3d(0.5*50.0,0,0)); //(0.25*50,13,0)
       drawArrow(Eigen::Vector3d(0.65,0,0));
    }

    if (mWorld->getSimFrames() > 800 && mWorld->getSimFrames() < 810){
       mBase->addExtForce(Eigen::Vector3d(0,180.0,0)); //(0.25*50,13,0)
    }

    if (mWorld->getSimFrames() > 800 && mWorld->getSimFrames() < 850){
       drawArrow(Eigen::Vector3d(0,0.65,0));
    }

    if (mWorld->getSimFrames() > 1000 && mWorld->getSimFrames() < 1010){
       mBase->addExtForce(Eigen::Vector3d(180.0,180.0,0)); //(0.25*50,13,0)
    }

    if (mWorld->getSimFrames() > 1000 && mWorld->getSimFrames() < 1050){
       drawArrow(Eigen::Vector3d(0.65,0.65,0));
    }
    /**/

    if (((int) walkState.simulationTime) % 5 == 0 ) {

        previousCoM = current.comPos;

    }

    if (walkState.mpcIter >= S+F || walkState.footstepCounter == 0) { // && walkState.footstepCounter > 3
    	walkState.controlIter = 0;
    	walkState.mpcIter = 0;
        walkState.footstepCounter++;
	walkState.supportFoot = !walkState.supportFoot;
        //std::cout << "Iteration " << walkState.controlIter << " Footstep " << walkState.footstepCounter << std::endl;
        std::cout <<  "new step starts at frame " << mWorld->getSimFrames() << std::endl;
    }
                  
    // simulation time, wait some time before starting the gait
    if (mWorld->getSimFrames() < wait) { walkState.simulationTime = 0; walkState.controlIter = 0; walkState.mpcIter = 0;}
    else walkState.simulationTime = mWorld->getSimFrames()-wait;
    

    // get measurements
    current.comPos = mBase->getCOM();
    current.comVel = mBase->getCOMLinearVelocity();
    current.comAcc = mBase->getCOMLinearAcceleration();
    current.zmpPos = current.comPos - current.comAcc / (eta*eta);
    current.leftFootPos = mLeftFoot->getCOM();
    current.rightFootPos = mRightFoot->getCOM();
    current.torsoOrient = getRPY(mBase);
    current.leftFootOrient = getRPY(mLeftFoot);
    current.rightFootOrient = getRPY(mRightFoot);

    // luenberger observer
    measured_state_x = Eigen::MatrixXd::Zero(2,1);
    measured_state_y = Eigen::MatrixXd::Zero(2,1); 

    measured_state_x << current.comPos(0) - 0.05, desired.zmpPos(0);
    measured_state_y << current.comPos(1), desired.zmpPos(1);

    if (mWorld->getSimFrames() > 10) {
       desired.obsState_x << A_obs*desired.obsState_x + B_obs*desired.zmpDot(0) - G_obs*(C_obs*desired.obsState_x - measured_state_x);
       desired.obsState_y << A_obs*desired.obsState_y + B_obs*desired.zmpDot(1) - G_obs*(C_obs*desired.obsState_y - measured_state_y);
    }
    
    // store observed disturbance in the data structure that is passed to the MPC class
    desired.disturbance << desired.obsState_x(3), desired.obsState_y(3), 0.0;          

    // combine the measurements with the LIP model for quasi-closed loop behavior
    // do this after the robot starts walking
    if (((int) walkState.simulationTime) % 50 == 0 && walkState.footstepCounter > 2) {

        Eigen::Vector3d pos = mTorso->getCOM();
        Eigen::Vector3d vel = mTorso->getCOMLinearVelocity();
        double gain = 0.075/2.0;
        desired.comPos(0) += gain * (current.comPos(0) - desired.comPos(0));
        desired.comVel(0) += gain * (current.comVel(0) - desired.comVel(0));
        desired.comPos(1) += gain * (current.comPos(1) - desired.comPos(1));
        desired.comVel(1) += gain * (current.comVel(1) - desired.comVel(1));

        desired.leftFootPos += 0*0.02 * (current.leftFootPos - desired.leftFootPos);
        desired.rightFootPos += 0*0.02 * (current.rightFootPos - desired.rightFootPos);

    }
  
    desired.comPosmeas = current.comPos;
    desired.comVelmeas = current.comVel;

    // compute the new desired state using MPC
    if (mWorld->getSimFrames() >= wait) desired = solver->solve(desired, walkState, ftsp_and_time); 
    // get the new plan, in case of modifications
    if (mWorld->getSimFrames() >= wait) ftsp_and_time = solver->getFtstpAndTimeMatrix();
    // get the new timing, in case of modifications    
    if (mWorld->getSimFrames() >= wait) walkState.controlIter = solver->getModifiedControlIter();

 
    // Compute inverse kinematics
    Eigen::VectorXd qDot =  getJointVelocitiesStacked(desired.getRelComPose(walkState.supportFoot),  current.getRelComPose(walkState.supportFoot),
  		 desired.getRelSwingFootPose(walkState.supportFoot),  current.getRelSwingFootPose(walkState.supportFoot),  desired.comVel);

    // Set the velocity of the floating base to zero
    for (int i = 0; i < 6; ++i){
        mRobot->setCommand(i, 0);
    }
    // set joint commands
    for (int i = 0; i < 50; ++i){
        mRobot->setCommand(i+6,qDot(i)); //velocity joint control
    }
    
    // Store the results in files (for plotting)
    //storeData();
    // Arm swing
    ArmSwing();
    FixWaistToChestJoints();

    // update the iteration counters, if a step is finished reset and change support foot
    ++walkState.controlIter;
    walkState.mpcIter = floor(walkState.controlIter*controlTimeStep/mpcTimeStep);
 
    // footstep plan update, to maintain consistency in the anticipative tail
    if (false && mWorld->getSimFrames()>= wait && walkState.footstepCounter > 2) {

       // in practice translate the remaining steps in the plan
       double offset_x = current.getSupportFootPose(walkState.supportFoot)(3) - ftsp_and_time(walkState.footstepCounter-1,0);
       double offset_y = current.getSupportFootPose(walkState.supportFoot)(4) - ftsp_and_time(walkState.footstepCounter-1,1);
       std::cout << offset_x << " " << offset_y << std::endl; 

       for (int i = walkState.footstepCounter ; i < ftsp_and_time.rows(); i++) {
           ftsp_and_time(i,0) = ftsp_and_time(i,0) + offset_x;
           ftsp_and_time(i,1) = ftsp_and_time(i,1) + offset_y;
       }
    }

    auto finish = std::chrono::high_resolution_clock::now();
    auto interval = std::chrono::time_point_cast<std::chrono::microseconds>(finish) - std::chrono::time_point_cast<std::chrono::microseconds>(start);
    //std::cout << "elapsed time for a control cycle: " << interval.count() << std::endl;

}

Eigen::MatrixXd Controller::getTorsoAndSwfJacobian() { 

    if (walkState.supportFoot == true) {
        mSupportFoot = mRobot->getBodyNode("r_sole");
        mSwingFoot = mRobot->getBodyNode("l_sole");
    } else {
        mSupportFoot = mRobot->getBodyNode("l_sole");
        mSwingFoot = mRobot->getBodyNode("r_sole");
    }

    Eigen::MatrixXd Jacobian_supportToBase;
    Jacobian_supportToBase =  mRobot->getJacobian(mTorso,mSupportFoot) - mRobot->getJacobian(mSupportFoot,mSupportFoot); //(mBase,mSupportFoot)
    _Rot_buffer << 1.0, 0.0, sin(current.torsoOrient(1)),
                   0.0, cos(current.torsoOrient(0)), -cos(current.torsoOrient(1)) * sin(current.torsoOrient(0)),
                   0.0, sin(current.torsoOrient(0)), -sin(current.torsoOrient(1)) * cos(current.torsoOrient(0)); 
    _R_CoM.block(0,0,3,3) << _Rot_buffer.inverse();
    Jacobian_supportToBase = _R_CoM * Jacobian_supportToBase;

    for (unsigned int i=0; i<44; i++){ 
    if (i!=5 || i!=6)  Jacobian_supportToBase.col(i).setZero();
    }

    Eigen::MatrixXd Jacobian_SupportToSwing =  mRobot->getJacobian(mSwingFoot,mSupportFoot) - mRobot->getJacobian(mSupportFoot,mSupportFoot);
    _Rot_buffer << 1.0, 0.0, sin(getRPY(mSwingFoot)(1)),
                   0.0, cos(getRPY(mSwingFoot)(0)), -cos(getRPY(mSwingFoot)(1)) * sin(getRPY(mSwingFoot)(0)),
                   0.0, sin(getRPY(mSwingFoot)(0)), -sin(getRPY(mSwingFoot)(1)) * cos(getRPY(mSwingFoot)(0)); 
    _R_swgf.block(0,0,3,3) << _Rot_buffer.inverse();
    Jacobian_SupportToSwing = _R_swgf * Jacobian_SupportToSwing;
    
    Eigen::MatrixXd Jacobian_tot_(12, 56);

    Jacobian_tot_ << Jacobian_supportToBase, Jacobian_SupportToSwing;

    Eigen::MatrixXd Jacobian_tot(12, 50);

    // Remove the floating base columns
    Jacobian_tot = Jacobian_tot_.block<12,50>(0, 6);

    return Jacobian_tot;
    
}

Eigen::VectorXd Controller::getJointVelocitiesStacked(Eigen::VectorXd desider_pos_base, Eigen::VectorXd actPosBase,
		Eigen::VectorXd desider_pos_SwingFoot, Eigen::VectorXd actPosSwingFoot, Eigen::VectorXd desider_com_vel){

        ComVref = Eigen::VectorXd::Zero(12);
        ComVref << 0.0, 0.0, 0.0,
                   desider_com_vel(0), desider_com_vel(1), desider_com_vel(2), 
                   0.0,0.0,0.0,
                   0.0,0.0,0.0;
        desired_pos = Eigen::VectorXd::Zero(12);           
	desired_pos << desider_pos_base, desider_pos_SwingFoot;

        // adding a pitch forward inclination can help!
        desired_pos(1) = desired_pos(1) + 0.01;
        
	// Assemble actual positions and orientations
	actual_pos = Eigen::VectorXd::Zero(12);
	actual_pos << actPosBase, actPosSwingFoot;

	// Get the proper jacobian and pseudoinvert it
	Jacobian_tot = getTorsoAndSwfJacobian();
	PseudoJacobian_tot = (Jacobian_tot.transpose())*(Jacobian_tot*Jacobian_tot.transpose()).inverse();

        // position error gain
	_taskGain = Eigen::MatrixXd::Identity(12,12);

	// Torso Orientation
	_taskGain(0,0) = 0.1;//0.1
	_taskGain(1,1) = 3.5;//0.1;
	_taskGain(2,2) = 0.1;//0.1;

	// CoM Position
	_taskGain(3,3) = 5;  //0.1  10
	_taskGain(4,4) = 5;
	_taskGain(5,5) = 2; //1 IT WAS ONE BEFORE !!!!

	// Swing Foot Orientation
	_taskGain(6,6) = 8.5; //5.5
	_taskGain(7,7) = 8.0; //7.0
	_taskGain(8,8) = 1;

	// Swing Foot Position
	_taskGain(9,9) = 5;
	_taskGain(10,10) = 5;
	_taskGain(11,11) = 5;	
	
        double ikGain = 9;

	Eigen::VectorXd qDot(50);
	qDot = PseudoJacobian_tot*(ComVref+ikGain*_taskGain*(desired_pos - actual_pos));
	
        dq = mRobot->getVelocities();  
        
	return (qDot - dq.segment(6,50))/controlTimeStep;
	
}

Eigen::Vector3d Controller::getRPY(dart::dynamics::BodyNode* body, dart::dynamics::BodyNode* referenceFrame) {
    Eigen::MatrixXd rotMatrix = body->getTransform(referenceFrame).rotation();

    Eigen::Vector3d RPY;
    RPY << atan2(rotMatrix(2,1),rotMatrix(2,2)),
        atan2(-rotMatrix(2,0),sqrt(rotMatrix(2,1)*rotMatrix(2,1)+rotMatrix(2,2)*rotMatrix(2,2))),
        atan2(rotMatrix(1,0),rotMatrix(0,0));

    return RPY;
}

Eigen::Vector3d Controller::getRPY(dart::dynamics::BodyNode* body) {
    Eigen::MatrixXd rotMatrix = body->getTransform().rotation();

    Eigen::Vector3d RPY;
    RPY << atan2(rotMatrix(2,1),rotMatrix(2,2)),
        atan2(-rotMatrix(2,0),sqrt(rotMatrix(2,1)*rotMatrix(2,1)+rotMatrix(2,2)*rotMatrix(2,2))),
        atan2(rotMatrix(1,0),rotMatrix(0,0));

    return RPY;
}

Eigen::Vector3d Controller::getZmpFromExternalForces()
{
    Eigen::Vector3d zmp_v;
    bool left_contact = false;
    bool right_contact = false;

    Eigen::Vector3d left_cop;
    if(abs(mLeftFoot->getConstraintImpulse()[5]) > 0.01){
        left_cop << -mLeftFoot->getConstraintImpulse()(1)/mLeftFoot->getConstraintImpulse()(5), mLeftFoot->getConstraintImpulse()(0)/mLeftFoot->getConstraintImpulse()(5), 0.0;
        Eigen::Matrix3d iRotation = mLeftFoot->getWorldTransform().rotation();
        Eigen::Vector3d iTransl   = mLeftFoot->getWorldTransform().translation();
        left_cop = iTransl + iRotation*left_cop;
        left_contact = true;
    }

    Eigen::Vector3d right_cop;
    if(abs(mRightFoot->getConstraintImpulse()[5]) > 0.01){
        right_cop << -mRightFoot->getConstraintImpulse()(1)/mRightFoot->getConstraintImpulse()(5), mRightFoot->getConstraintImpulse()(0)/mRightFoot->getConstraintImpulse()(5), 0.0;
        Eigen::Matrix3d iRotation = mRightFoot->getWorldTransform().rotation();
        Eigen::Vector3d iTransl   = mRightFoot->getWorldTransform().translation();
        right_cop = iTransl + iRotation*right_cop;
        right_contact = true;
    }

    if(left_contact && right_contact){
        zmp_v << (left_cop(0)*mLeftFoot->getConstraintImpulse()[5] + right_cop(0)*mRightFoot->getConstraintImpulse()[5])/(mLeftFoot->getConstraintImpulse()[5] + mRightFoot->getConstraintImpulse()[5]),
                 (left_cop(1)*mLeftFoot->getConstraintImpulse()[5] + right_cop(1)*mRightFoot->getConstraintImpulse()[5])/(mLeftFoot->getConstraintImpulse()[5] + mRightFoot->getConstraintImpulse()[5]),
		 0.0;
    }else if(left_contact){
        zmp_v << left_cop(0), left_cop(1), 0.0;
    }else if(right_contact){
        zmp_v << right_cop(0), right_cop(1), 0.0;
    }else{
        // No contact detected
        zmp_v << 0.0, 0.0, 0.0;
    }

    return zmp_v;
}

void Controller::setInitialConfiguration() {
  Eigen::VectorXd q = mRobot->getPositions();


    // Floating Base
    q[0] = 0.0;
    q[1] = 2*M_PI/180;
    q[2] = 0.0;
    q[3] = 0.03;
    q[4] = 0.00;
    q[5] = 0.753;

    // Right Leg
    q[44] = 0.0;            // hip yaw
    q[45] = 3*M_PI/180;    // hip roll
    q[46] = -25*M_PI/180;   // hip pitch
    q[47] = 50*M_PI/180;   // knee pitch
    q[48] = -30*M_PI/180; // ankle pitch
    q[49] = -4*M_PI/180;   // ankle roll         
    // Left Leg
    q[50] = 0.0;            // hip yaw
    q[51] = -3*M_PI/180;    // hip roll
    q[52] = -25*M_PI/180;   // hip pitch
    q[53] = 50*M_PI/180;   // knee pitch
    q[54] = -30*M_PI/180; // ankle pitch
    q[55] = 4*M_PI/180;   // ankle roll       

    mRobot->setPositions(q);

    // Additional arm position setting
    mRobot->setPosition(mRobot->getDof("R_SHOULDER_P")->getIndexInSkeleton(), (4)*M_PI/180 );
    mRobot->setPosition(mRobot->getDof("R_SHOULDER_R")->getIndexInSkeleton(), -8*M_PI/180  );
    mRobot->setPosition(mRobot->getDof("R_SHOULDER_Y")->getIndexInSkeleton(), 0 );

    mRobot->setPosition(mRobot->getDof("R_ELBOW_P")->getIndexInSkeleton(), -25*M_PI/180 );
    mRobot->setPosition(mRobot->getDof("R_ELBOW_P")->getIndexInSkeleton(), -25*M_PI/180 );

    mRobot->setPosition(mRobot->getDof("L_SHOULDER_P")->getIndexInSkeleton(), (4)*M_PI/180  );
    mRobot->setPosition(mRobot->getDof("L_SHOULDER_R")->getIndexInSkeleton(), 8*M_PI/180  );
    mRobot->setPosition(mRobot->getDof("L_SHOULDER_Y")->getIndexInSkeleton(), 0 );

    mRobot->setPosition(mRobot->getDof("L_ELBOW_P")->getIndexInSkeleton(), -25*M_PI/180 ); 
    mRobot->setPosition(mRobot->getDof("L_ELBOW_P")->getIndexInSkeleton(), -25*M_PI/180 ); 

    // Set initial chest position 
    mRobot->setPosition(mRobot->getDof("CHEST_P")->getIndexInSkeleton(), 0.025 );
    mRobot->setPosition(mRobot->getDof("CHEST_Y")->getIndexInSkeleton(), 0.025 );

}


void Controller::ArmSwing() {

  if (mWorld->getSimFrames() > wait) {
    mRobot->setPosition(mRobot->getDof("R_SHOULDER_P")->getIndexInSkeleton(), (4-6*sin(2*M_PI*0.01*(mWorld->getSimFrames())/(2.0*(singleSupportDuration+doubleSupportDuration))))*M_PI/180 );   // 2.5 amplitude
    mRobot->setPosition(mRobot->getDof("L_SHOULDER_P")->getIndexInSkeleton(), (4+6*sin(2*M_PI*0.01*(mWorld->getSimFrames())/(2.0*(singleSupportDuration+doubleSupportDuration))))*M_PI/180 );

    mRobot->setPosition(mRobot->getDof("R_ELBOW_P")->getIndexInSkeleton(), (-30-7*sin(2*M_PI*0.01*(mWorld->getSimFrames())/(2.0*(singleSupportDuration+doubleSupportDuration))))*M_PI/180 );   // 2.5 amplitude
    mRobot->setPosition(mRobot->getDof("L_ELBOW_P")->getIndexInSkeleton(), (-30+7*sin(2*M_PI*0.01*(mWorld->getSimFrames())/(2.0*(singleSupportDuration+doubleSupportDuration))))*M_PI/180 );
}

}

void Controller::FixWaistToChestJoints() {

     double inclination = 0.55;

     mRobot->setPosition(mRobot->getDof("CHEST_P")->getIndexInSkeleton(), mRobot->getPosition(mRobot->getDof("CHEST_P")->getIndexInSkeleton()) - 0.05*(mRobot->getPosition(mRobot->getDof("CHEST_P")->getIndexInSkeleton())) + inclination*0.009);

}


void Controller::storeData() {

     Eigen::VectorXd COMPOS_meas  = Eigen::VectorXd::Zero(3);
     COMPOS_meas = mRobot->getCOM();
     Eigen::VectorXd COMVEL_meas  = Eigen::VectorXd::Zero(3);
     COMVEL_meas = mTorso->getCOMLinearVelocity();
     Eigen::VectorXd FOOT_meas  = Eigen::VectorXd::Zero(3);
     FOOT_meas = mSupportFoot->getCOM();

     Eigen::VectorXd ZMPPOS_meas_cop =  Eigen::VectorXd::Zero(3);
     ZMPPOS_meas_cop = COMPOS_meas - mRobot->getCOMLinearAcceleration()/(eta*eta);
     ZMPPOS_meas_cop = getZmpFromExternalForces();
     

     ofstream myfile;

     myfile.open ("./Data/x_RF.txt",ios::app);
     myfile << mSupportFoot->getCOM()(0) <<endl; 
     myfile.close();
     myfile.open ("./Data/y_RF.txt",ios::app);
     myfile << mSupportFoot->getCOM()(1) <<endl; 
     myfile.close();

     myfile.open ("./Data/x_m.txt",ios::app);
     myfile << COMPOS_meas(0) <<endl; 
     myfile.close();
     myfile.open ("./Data/y_m.txt",ios::app);
     myfile << COMPOS_meas(1) <<endl; 
     myfile.close();
     myfile.open ("./Data/xz_m_cop.txt",ios::app);
     myfile << ZMPPOS_meas_cop(0) <<endl; //current.zmpPos
     myfile.close();
     myfile.open ("./Data/yz_m_cop.txt",ios::app);
     myfile << ZMPPOS_meas_cop(1) <<endl; 
     myfile.close();
     myfile.open ("./Data/x.txt",ios::app);
     myfile << desired.comPos(0) <<endl; 
     myfile.close();
     myfile.open ("./Data/y.txt",ios::app);
     myfile << desired.comPos(1) <<endl; 
     myfile.close();
     myfile.open ("./Data/xz.txt",ios::app);
     myfile << desired.zmpPos(0) <<endl; 
     myfile.close();
     myfile.open ("./Data/yz.txt",ios::app);
     myfile << desired.zmpPos(1) <<endl; 
     myfile.close();

}

void Controller::createForceLine()
  {

    dart::dynamics::SimpleFramePtr lineFrame
        = std::make_shared<dart::dynamics::SimpleFrame>(
            dart::dynamics::Frame::World());

    mForceLine = std::make_shared<dart::dynamics::LineSegmentShape>( current.comPos, previousCoM, 3.0);
    mForceLine->addDataVariance(dart::dynamics::Shape::DYNAMIC_VERTICES);
    mForceLine->addConnection(0,1) ;		


    lineFrame->setShape(mForceLine);
    lineFrame->createVisualAspect();
    lineFrame->getVisualAspect()->setColor(DefaultForceLineColor);

    mWorld->addSimpleFrame(lineFrame);

  }

void Controller::drawArrow(Eigen::Vector3d force){

    auto visualShapeNodes = mTorso->getShapeNodesWith<VisualAspect>();

    if(visualShapeNodes.size() == 3u )
    {
    //assert(visualShapeNodes[2]->getShape() == mArrow);
        visualShapeNodes[2]->remove();
    }

    if(visualShapeNodes.size() == 4u )
    {
    //assert(visualShapeNodes[2]->getShape() == mArrow);
        visualShapeNodes[2]->remove();
        visualShapeNodes[3]->remove();
    }

    std::shared_ptr<ArrowShape> mArrow;
    ArrowShape::Properties arrow_properties;
    arrow_properties.mRadius = 0.05;
    Eigen::Vector3d tail_pos = Eigen::Vector3d(-0.1, -0.1, 0.15) ;
    Eigen::Vector3d tail_offset = force;
    Eigen::Vector3d head_pos = tail_pos - tail_offset;
    mArrow = std::shared_ptr<ArrowShape>(new ArrowShape(
             Eigen::Vector3d(0.0, 0.0, 0.0),
             Eigen::Vector3d(0.2, 0.05, 0.05),
             arrow_properties, dart::Color::Red(1.0)));
    mArrow->setPositions(
            head_pos,
            tail_pos);
    mTorso->createShapeNodeWith<VisualAspect>(mArrow);
   
}
