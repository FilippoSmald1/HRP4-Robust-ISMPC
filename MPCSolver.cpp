#include "MPCSolver.hpp"

MPCSolver::MPCSolver(const Eigen::MatrixXd& ftsp_and_timings) {

    // Matrices for cost function 
    costFunctionH = Eigen::MatrixXd::Zero(2*(N+M),2*(N+M));
    costFunctionF = Eigen::VectorXd::Zero(2*(N+M));

    // Matrices for stability constraint
    Aeq = Eigen::MatrixXd::Zero(2,(N*2)+(M*2));
    beq = Eigen::VectorXd::Zero(2);

    // Matrices for ZMP constraint
    AZmp = Eigen::MatrixXd::Zero(2*N,2*(N+M));
    bZmpMax = Eigen::VectorXd::Zero(2*N);
    bZmpMin = Eigen::VectorXd::Zero(2*N);

    // Matrices for kinematic constraint
    AFootsteps = Eigen::MatrixXd::Zero(2*M,2*(M+N));
    bFootstepsMax = Eigen::VectorXd::Zero(2*M);
    bFootstepsMin = Eigen::VectorXd::Zero(2*M);

    // Matrices for double support constraint
    ADoubleSupport = Eigen::MatrixXd::Zero(2,2*(M+N));
    bDoubleSupport = Eigen::VectorXd::Zero(2);

    // Matrices for fixed footsteps constraint
    AFixedFootsteps = Eigen::MatrixXd::Zero(2*M,2*(M+N));
    bFixedFootsteps = Eigen::VectorXd::Zero(2*M);

    // Matrices for all constraints stacked
    AConstraint = Eigen::MatrixXd::Zero(2*(N+M)+2,2*(N+M));
    bConstraintMin = Eigen::VectorXd::Zero(2*(N+M)+2);
    bConstraintMax = Eigen::VectorXd::Zero(2*(N+M)+2);

    // Matrices for ZMP prediction (might be useful)
    p = Eigen::VectorXd::Ones(N);
    P = Eigen::MatrixXd::Ones(N,N)*mpcTimeStep;
    for(int i = 0; i < N; ++i) {
    	for(int j = 0; j < N; ++j) {
            if (j > i) P(i, j)=0;
        }
    }

    halfAZmp = Eigen::MatrixXd::Zero(N,N+M);
    CcFull = Eigen::MatrixXd::Zero(N,M+1);
    Ic = Eigen::MatrixXd::Identity(N,N);
    differenceMatrix = Eigen::MatrixXd::Identity(M,M);
    for (int i = 0; i < M-1; ++i) {
        differenceMatrix(i+1,i) = -1;
    }
    pFr = Eigen::VectorXd::Zero(M);
    pFl = Eigen::VectorXd::Zero(M);
    pF = Eigen::VectorXd::Ones(M);
    currentFootstepKinematic = Eigen::VectorXd::Zero(M);
    halfH = Eigen::MatrixXd::Zero(N+M, N+M);
    xCandidateFootsteps = Eigen::VectorXd::Zero(M);
    yCandidateFootsteps = Eigen::VectorXd::Zero(M);
    nConstraints = Aeq.rows() + AFootsteps.rows() + AZmp.rows() + ADoubleSupport.rows() + AFixedFootsteps.rows();
    footstepsOptimalX = Eigen::VectorXd::Zero(M);
    footstepsOptimalY = Eigen::VectorXd::Zero(M);
    decision_variables = Eigen::VectorXd::Zero(2*(N+M));
    fix_in_double_support = Eigen::VectorXd::Zero(2);

    cost_dist_x = Eigen::MatrixXd::Constant(1,1,1.0);
    cost_dist_y = Eigen::MatrixXd::Constant(1,1,1.0);
    A_dist_x = Eigen::MatrixXd::Constant(2,1,1.0);
    A_dist_y = Eigen::MatrixXd::Constant(2,1,1.0); 
    f_dist_x = Eigen::VectorXd::Zero(1);
    f_dist_y = Eigen::VectorXd::Zero(1);
    b_dist_min_x = Eigen::VectorXd::Zero(2);
    b_dist_max_x = Eigen::VectorXd::Zero(2);
    b_dist_min_y = Eigen::VectorXd::Zero(2);
    b_dist_max_y = Eigen::VectorXd::Zero(2);
    w_bar = Eigen::Vector2d(midrange_x,midrange_y);
    footstepPredicted = Eigen::Vector4d(ftsp_and_timings(1,0),ftsp_and_timings(1,1),0.0,0.0);
    exponential_weight = Eigen::MatrixXd::Constant(1,N,0.0);
    
    double lambda_ = exp(-eta*mpcTimeStep);

    for (double i=0.0; i < N; i++) {
        Aeq(0,i) = ((1.0/eta)*(1.0-lambda_)/(1.0-pow(lambda_,N)))*exp(-eta*mpcTimeStep*i)-mpcTimeStep*exp(-eta*mpcTimeStep*N);
        Aeq(1,N+M+i) = ((1.0/eta)*(1.0-lambda_)/(1.0-pow(lambda_,N)))*exp(-eta*mpcTimeStep*i)-mpcTimeStep*exp(-eta*mpcTimeStep*N);
        exponential_weight(0,i) = exp(-eta*mpcTimeStep*i);
    }

    double ch = cosh(eta*mpcTimeStep);
    double sh = sinh(eta*mpcTimeStep);
    A_upd = Eigen::MatrixXd::Zero(3,3);
    B_upd = Eigen::VectorXd::Zero(3);
    A_upd<<ch,sh/eta,1-ch,eta*sh,ch,-eta*sh,0,0,1;
    B_upd<<mpcTimeStep-sh/eta,1-ch,mpcTimeStep;

    old_fsCount = 0;
    ct = 0;
    ds_samples = F;
    itr = 0;

    xz_dot = 0.0;
    yz_dot = 0.0;
 
    // build midpoint sequence
    ftsp_midpoint = Eigen::MatrixXd::Zero(ftsp_and_timings.rows()*(S+F),3);
    Eigen::MatrixXd double_support_transition = Eigen::MatrixXd::Zero(F,1);
    Eigen::MatrixXd p_dst = Eigen::MatrixXd::Ones(F,1);
    for (int i = 0; i < F; i++) double_support_transition(i,0) = (double)i/(double)F;
    for (int i = 0; i < ftsp_and_timings.rows() - 1 ; i++) {
        // single support phases
        ftsp_midpoint.block(i*(S+F),0,S,1) = ftsp_and_timings(i,0)*Eigen::MatrixXd::Ones(S,1);
        ftsp_midpoint.block(i*(S+F),1,S,1) = ftsp_and_timings(i,1)*Eigen::MatrixXd::Ones(S,1);
        ftsp_midpoint.block(i*(S+F),2,S,1) = ftsp_and_timings(i,2)*Eigen::MatrixXd::Ones(S,1);
        // double support phases
        ftsp_midpoint.block(i*(S+F)+S,0,F,1) = ftsp_and_timings(i,0)*p_dst + (ftsp_and_timings(i+1,0)-ftsp_and_timings(i,0))*double_support_transition;
        ftsp_midpoint.block(i*(S+F)+S,1,F,1) =  ftsp_and_timings(i,1)*p_dst + (ftsp_and_timings(i+1,1)-ftsp_and_timings(i,1))*double_support_transition;
        ftsp_midpoint.block(i*(S+F)+S,2,F,1) =  ftsp_and_timings(i,2)*p_dst + (ftsp_and_timings(i+1,2)-ftsp_and_timings(i,2))*double_support_transition;
    }
    
    // useful multiplier
    deltas = Eigen::MatrixXd::Zero(1,N);
    for (int i = 0; i<N; i++) deltas(0,i) = exp(-mpcTimeStep*eta*i);

    restriction = Eigen::VectorXd::Zero(2*N);
    res = Eigen::VectorXd::Zero(N);
    for (double i = 0.0; i < N; i++) res(i) = (i / (double) N) * footConstraintSquareWidth/2.0;
    restriction << Ic*res, Ic*res;
        
}

MPCSolver::~MPCSolver() {}

State MPCSolver::solve(State current, WalkState walkState, const Eigen::MatrixXd& ftsp_and_timings) {

    itr = walkState.controlIter;
    fsCount = walkState.footstepCounter;
    ftsp_and_timings_ = ftsp_and_timings;
    
    // Output struct
    State next = current;

    if (walkState.mpcIter <= 1){

       activate_feasibility_recovery = false;
       sta_counter = 0;
       
    } 


    Eigen::MatrixXd state_z = Eigen::MatrixXd::Zero(2,1);
    state_z << current.comPos(2), current.comVel(2);

    int iter_step = walkState.mpcIter;

    // DISTURBANCE SATURATION
    f_dist_x << -(current.disturbance(0)+midrange_x);
    f_dist_y << -(current.disturbance(1)+midrange_y);
    b_dist_min_x << midrange_x-dist_range_x/2.0, (w_bar(0)-dist_range_x*(exp(eta*mpcTimeStep)-1.0)*(1.0+alpha)/(eta*eta));
    b_dist_max_x << midrange_x+dist_range_x/2.0, (w_bar(0)+dist_range_x*(exp(eta*mpcTimeStep)-1.0)*(1.0+alpha)/(eta*eta));
    b_dist_min_y << (midrange_y-dist_range_y/2.0), (w_bar(1)-dist_range_y*(exp(eta*mpcTimeStep)-1.0)*(1.0+alpha)/(eta*eta));
    b_dist_max_y << (midrange_y+dist_range_y/2.0), (w_bar(1)+dist_range_y*(exp(eta*mpcTimeStep)-1.0)*(1.0+alpha)/(eta*eta));

    w_bar << solveQP_hpipm(cost_dist_x,f_dist_x,A_dist_x,b_dist_min_x,b_dist_max_x,0), 
             solveQP_hpipm(cost_dist_y,f_dist_y,A_dist_y,b_dist_min_y,b_dist_max_y,0);
             
    //std::cout << "current.disturbance " << current.disturbance.transpose() << std::endl;            

    //////// IMPLEMENT HERE IS-MPC //////////////
    //////////////////////////
    /////////////////////////
    /////////////////////////

    // Reset constraint matrices
    AZmp.setZero();
    AFootsteps.setZero();

    // Get the pose of the support foot in the world frame

    Eigen::VectorXd supportFootPose = Eigen::VectorXd::Zero(6);
    if (walkState.footstepCounter > 2 && false) {
      supportFootPose << current.getSupportFootPose(walkState.supportFoot);
      supportFootPose(5) = 0.0;
    } else {
      supportFootPose << Eigen::VectorXd::Zero(3), ftsp_and_timings(fsCount-1,0), ftsp_and_timings(fsCount-1,1), 0.0;
    }

    // STABILITY CONSTRAINT
    anticipative_x = 0.0;
    anticipative_y = 0.0;

    for (double i=0.0; i<Prev; i++) {
        anticipative_x = anticipative_x + 
          		  eta*mpcTimeStep*exp(-eta*mpcTimeStep*(N+i))*(ftsp_midpoint((int)(walkState.simulationTime/(mpcTimeStep/controlTimeStep))+N+i,0));
        anticipative_y = anticipative_y + 
        	         eta*mpcTimeStep*exp(-eta*mpcTimeStep*(N+i))*(ftsp_midpoint((int)(walkState.simulationTime/(mpcTimeStep/controlTimeStep))+N+i,1));
    }

    anticipative_x = anticipative_x + 
     		     exp(-eta*mpcTimeStep*(N+Prev))*(ftsp_midpoint((int)(walkState.simulationTime/(mpcTimeStep/controlTimeStep))+N+Prev,0)-ftsp_and_timings(fsCount-1,0));
    anticipative_y = anticipative_y + 
    		     exp(-eta*mpcTimeStep*(N+Prev))*(ftsp_midpoint((int)(walkState.simulationTime/(mpcTimeStep/controlTimeStep))+N+Prev,0)-ftsp_and_timings(fsCount-1,1));

    beq << current.comPos(0) + current.comVel(0)/eta - current.zmpPos(0)*(1.0-exp(-eta*mpcTimeStep*N)) - 
    	   anticipative_x + (activate_observer_x*(w_bar(0)-midrange_x) + activate_midrange_x*midrange_x)/(eta*eta),
          current.comPos(1) + current.comVel(1)/eta - current.zmpPos(1)*(1.0-exp(-eta*mpcTimeStep*N)) - 
          anticipative_y + (activate_observer_y*(w_bar(1)-midrange_y) + activate_midrange_y*midrange_y)/(eta*eta);

    
    // FEASIBILITY REGION WITH FIXED FOOTSTEPS
    Eigen::VectorXd product = exponential_weight*(restriction.segment(0,N)-p*footConstraintSquareWidth/2.0 + ftsp_midpoint.block((int)(walkState.simulationTime/(mpcTimeStep/controlTimeStep)),0,N,1));
    x_u_m = eta*mpcTimeStep*product(0) + anticipative_x - (activate_observer_x*(w_bar(0)-midrange_x) + activate_midrange_x*midrange_x)/(eta*eta);
    product = exponential_weight*(-restriction.segment(0,N)+p*footConstraintSquareWidth/2.0+ftsp_midpoint.block((int)(walkState.simulationTime/(mpcTimeStep/controlTimeStep)),0,N,1));
    x_u_M = eta*mpcTimeStep*product(0) + anticipative_x - (activate_observer_x*(w_bar(0)-midrange_x) + activate_midrange_x*midrange_x)/(eta*eta);

    product = exponential_weight*(restriction.segment(0,N)-p*footConstraintSquareWidth/2.0+ftsp_midpoint.block((int)(walkState.simulationTime/(mpcTimeStep/controlTimeStep)),1,N,1));
    y_u_m = eta*mpcTimeStep*product(0) + anticipative_y - (activate_observer_x*(w_bar(1)-midrange_y) + activate_midrange_y*midrange_y)/(eta*eta);
    product = exponential_weight*(-restriction.segment(0,N)+p*footConstraintSquareWidth/2.0+ftsp_midpoint.block((int)(walkState.simulationTime/(mpcTimeStep/controlTimeStep)),1,N,1));
    y_u_M = eta*mpcTimeStep*product(0) + anticipative_y - (activate_observer_y*(w_bar(1)-midrange_y) + activate_midrange_y*midrange_y)/(eta*eta);

    // Construct the Cc matrix, which maps iterations to predicted footsteps
    CcFull.setZero();

    int fsAhead = 0;
    for (int i = 0; i < N; ++i) {
        // if with the prediction index i we have reach the next footstep increase fsAhead
        if (ftsp_and_timings(walkState.footstepCounter-1,3) + walkState.mpcIter + i + 1 >= (int)ftsp_and_timings(walkState.footstepCounter + fsAhead ,3) ) 
fsAhead++;    
        // how many samples are left till the end of the current footstep?
        int samplesTillNextFootstep = (int)ftsp_and_timings(walkState.footstepCounter + fsAhead ,3) - ((int)ftsp_and_timings(walkState.footstepCounter-1,3) + walkState.mpcIter + i + 1);
        // If it is the current footstep, it does not go in Cc, but subsequent ones do
        if (samplesTillNextFootstep > ds_samples) {
            CcFull(i, fsAhead) = 1;
        } else {
            CcFull(i, fsAhead) = (double)samplesTillNextFootstep / (double)ds_samples;
            CcFull(i, fsAhead + 1) = 1.0 - (double)samplesTillNextFootstep / (double)ds_samples;
        }
    }

    Eigen::VectorXd currentFootstepZmp = CcFull.col(0);
    Cc = CcFull.block(0,1,N,M);

    // ZMP CONSTRAINTS
    // Construct the ZMP constraint
    // ****************************
    halfAZmp << Ic*P, -Ic*Cc;
    AZmp.block(0,0,N,N+M) = halfAZmp;
    AZmp.block(N,N+M,N,N+M) = halfAZmp;

    // Construct the b vector of the ZMP constraint
    Eigen::VectorXd bZmpSizeTerm = Eigen::VectorXd::Zero(2*N);
    Eigen::VectorXd bZmpStateTerm = Eigen::VectorXd::Zero(2*N);

    bZmpSizeTerm << Ic*p*footConstraintSquareWidth/2.0, Ic*p*footConstraintSquareWidth/2.0;
 
    bZmpStateTerm << Ic*(-p*current.zmpPos(0)+currentFootstepZmp*supportFootPose(3)), Ic*(-p*current.zmpPos(1)+currentFootstepZmp*supportFootPose(4));

    if (activate_restriction == 1.0) {
       bZmpMin = - bZmpSizeTerm + bZmpStateTerm + activate_restriction*restriction;
       bZmpMax =   bZmpSizeTerm + bZmpStateTerm - activate_restriction*restriction;
    } else {
       bZmpMin = - bZmpSizeTerm + bZmpStateTerm;
       bZmpMax =   bZmpSizeTerm + bZmpStateTerm;
    }


    // KINEMATIC CONSTRAINTS

    // Construct the matrices that activate the left or right constraint
    pFr.setZero();
    pFl.setZero();
 
     // Vector for the current footstep position
    currentFootstepKinematic.setZero();
    currentFootstepKinematic(0) = 1;

    // Assemble the A matrix for the kinematic constraint, and rotate
    AFootsteps.block(0,N,M,M) = differenceMatrix;
    AFootsteps.block(M,2*N+M,M,M) = differenceMatrix;

    for (int i = 0; i < M; i++) pFl(i) = pow(-1.0, walkState.supportFoot+i);

    bFootstepsMin << -pF*wkx/2.0 + supportFootPose(3)*currentFootstepKinematic, -pF*wky/2.0 + pFl*ell + supportFootPose(4)*currentFootstepKinematic;
    bFootstepsMax <<  pF*wkx/2.0 + supportFootPose(3)*currentFootstepKinematic, pF*wky/2.0 + pFl*ell + supportFootPose(4)*currentFootstepKinematic;

    if (walkState.footstepCounter <= 1) {
       AFootsteps.setZero();
       bFootstepsMin.setZero();
       bFootstepsMax.setZero();
    }

    // DOUBLE SUPPORT CONSTRAINTS
    ADoubleSupport.setZero();
    bDoubleSupport.setZero();
    if (walkState.mpcIter > S ) {
       ADoubleSupport(0,N) = 1.0;
       ADoubleSupport(1,2*N+M) = 1.0;
       bDoubleSupport(0) = fix_in_double_support(0);
       bDoubleSupport(1) = fix_in_double_support(1);
    }

    // FIXED FOOTSTEP CONSTRAINTS
    AFixedFootsteps.setZero();
    bFixedFootsteps.setZero();
    if (activate_fixed_footstep_constraint>0) {
       AFixedFootsteps.block(0,N,M,M) = Eigen::MatrixXd::Identity(M,M);
       AFixedFootsteps.block(M,2*N+M,M,M) = Eigen::MatrixXd::Identity(M,M);
       bFixedFootsteps << ftsp_and_timings.block(fsCount,0,M,1) , ftsp_and_timings.block(fsCount,1,M,1); //// CAREFUL, THERE MIGHT NEED TO BE VECTOR!
    }

    // COST FUNCTION
    halfH.setZero();
    halfH << qZd*Eigen::MatrixXd::Identity(N,N) + qZ*P.transpose()*P,
             -qZ*P.transpose()*Cc, -qZ*Cc.transpose()*P, qZ*Cc.transpose()*Cc + qF*Eigen::MatrixXd::Identity(M,M);
    costFunctionH.block(0,0,N+M,N+M) = halfH;
    costFunctionH.block(N+M,N+M,N+M,N+M) = halfH;

    // PROBLEM HERE
    // Contruct candidate footstep vectors
    for (int i = 0; i < M; i++) {
 
         xCandidateFootsteps(i) = ftsp_and_timings(fsCount+i,0);
         yCandidateFootsteps(i) = ftsp_and_timings(fsCount+i,1);

    } 

    costFunctionF << qZ*P.transpose()*(p*current.zmpPos(0) - currentFootstepZmp*supportFootPose(3)),
		     -qZ*Cc.transpose()*(p*current.zmpPos(0) - currentFootstepZmp*supportFootPose(3)) - qF*xCandidateFootsteps,
		     qZ*P.transpose()*(p*current.zmpPos(1) - currentFootstepZmp*supportFootPose(4)),
		     -qZ*Cc.transpose()*(p*current.zmpPos(1) - currentFootstepZmp*supportFootPose(4)) - qF*yCandidateFootsteps;

    // STACK CONSTRAINT MATRICES
    AConstraint.resize(nConstraints, 2*(N+M));
    bConstraintMin.resize(nConstraints);
    bConstraintMax.resize(nConstraints);
  
    if (activate_fixed_footstep_constraint){

       AConstraint    << Aeq, 0.0*ADoubleSupport, activate_fixed_footstep_constraint*AFixedFootsteps, 0.0*AFootsteps, AZmp;
       bConstraintMin << beq, 0.0*bDoubleSupport, activate_fixed_footstep_constraint*bFixedFootsteps, 0.0*bFootstepsMin, bZmpMin;
       bConstraintMax << beq, 0.0*bDoubleSupport, activate_fixed_footstep_constraint*bFixedFootsteps, 0.0*bFootstepsMax, bZmpMax;

    } else {
   
      double ds_constr = 0.0;
      //if (walkState.mpcIter >= S) ds_constr = 1.0;
      AConstraint    << Aeq, ds_constr*ADoubleSupport, activate_fixed_footstep_constraint*AFixedFootsteps, AFootsteps, AZmp;
      bConstraintMin << beq, ds_constr*bDoubleSupport, activate_fixed_footstep_constraint*bFixedFootsteps, bFootstepsMin, bZmpMin;
      bConstraintMax << beq, ds_constr*bDoubleSupport, activate_fixed_footstep_constraint*bFixedFootsteps, bFootstepsMax, bZmpMax;

    }

    if (activate_feasibility_recovery && walkState.footstepCounter > 2 ) {
      double ds_constr = 0.0;
      //if (walkState.mpcIter >= S) ds_constr = 1.0; // TO FIX (if ds constraints a re activated it becomes unfeasible (?))
      AConstraint    << Aeq, ds_constr*ADoubleSupport, activate_fixed_footstep_constraint*AFixedFootsteps, AFootsteps, AZmp;
      bConstraintMin << beq, ds_constr*bDoubleSupport, activate_fixed_footstep_constraint*bFixedFootsteps, bFootstepsMin, bZmpMin;
      bConstraintMax << beq, ds_constr*bDoubleSupport, activate_fixed_footstep_constraint*bFixedFootsteps, bFootstepsMax, bZmpMax;
    }

    // SOLVE QP
    Eigen::VectorXd decision_variables_ = solveQP(costFunctionH, costFunctionF, AConstraint, bConstraintMin, bConstraintMax);

    // UPDATE TAIL AND FOOTSTEP PLAN IF NEEDED
    Eigen::VectorXd zDotOptimalX(N);
    Eigen::VectorXd zDotOptimalY(N);
    Eigen::VectorXd footstepsOptimalX(M);
    Eigen::VectorXd footstepsOptimalY(M);

    zDotOptimalX = (decision_variables_.head(N));
    zDotOptimalY = (decision_variables_.segment(N+M,N));
    footstepsOptimalX = decision_variables_.segment(N,M);
    footstepsOptimalY = decision_variables_.segment(2*N+M,M);

    next.zmpDot << zDotOptimalX(0), zDotOptimalY(0);

    if (walkState.mpcIter >= S-4 && walkState.mpcIter < S){ 
       fix_in_double_support << footstepsOptimalX(0), footstepsOptimalY(1);
    }

    Eigen::Vector2d push = Eigen::Vector2d(0.0,0.0);
  
    Eigen::Vector3d currentStateX = Eigen::Vector3d(current.comPos(0), current.comVel(0), current.zmpPos(0));
    Eigen::Vector3d currentStateY = Eigen::Vector3d(current.comPos(1), current.comVel(1), current.zmpPos(1));
    Eigen::Vector3d nextStateX = A_upd*currentStateX + B_upd*zDotOptimalX(0) + Eigen::Vector3d(0.0,controlTimeStep,0.0)*(activate_observer_x*(w_bar(0)-midrange_x) + activate_midrange_x*midrange_x + push(0));
    Eigen::Vector3d nextStateY = A_upd*currentStateY + B_upd*zDotOptimalY(0) + Eigen::Vector3d(0.0,controlTimeStep,0.0)*(activate_observer_y*(w_bar(1)-midrange_y) + activate_midrange_y*midrange_y + push(1));  // 1.6

    next.comPos = Eigen::Vector3d(nextStateX(0), nextStateY(0), comTargetHeight);
    next.comVel = Eigen::Vector3d(nextStateX(1), nextStateY(1), 0.0);
    next.zmpPos = Eigen::Vector3d(nextStateX(2), nextStateY(2), 0.0);
    next.comAcc = eta*eta * (next.comPos - next.zmpPos);
    next.torsoOrient = Eigen::Vector3d(0.0, 0.0, 0.0);

    // IF FEASIBILITY RECOVERY == TRUE : ADAPT PLAN
    if (activate_feasibility_recovery) {

    // MODIFY COMPUTED STEPS IN THE PLAN

       for (int i = 0; i < 1; i++) { // up to M in theory
           ftsp_and_timings_(fsCount+i,0) = footstepsOptimalX(i);
           ftsp_and_timings_(fsCount+i,1) = footstepsOptimalY(i);
       }

       // TRANSLATE THE REMAINING STEPS 
       double offset_x = footstepsOptimalX(0) - ftsp_and_timings_(fsCount,0);
       double offset_y = footstepsOptimalY(0) - ftsp_and_timings_(fsCount,1);
       //std::cout << "offset_x " << offset_x << " offset_y " << offset_y << std::endl;
       //int stop = getchar();
       for (int i = fsCount+1 ; i < ftsp_and_timings_.rows(); i++) {
           ftsp_and_timings_(i,0) = ftsp_and_timings_(i,0) + offset_x;
           ftsp_and_timings_(i,1) = ftsp_and_timings_(i,1) + offset_y;
       }

       // RECOMPUTE CENTERLINE SEQUENCE FOR TAIL
       Eigen::MatrixXd double_support_transition = Eigen::MatrixXd::Zero(F,1);
       Eigen::MatrixXd p_dst = Eigen::MatrixXd::Ones(F,1);
       for (int i = 0; i < F; i++) double_support_transition(i,0) = (double)i/(double)F;
       for (int i = 0; i < ftsp_and_timings.rows() - 1 ; i++) {
           // single support phases
           ftsp_midpoint.block(i*(S+F),0,S,1) = ftsp_and_timings(i,0)*Eigen::MatrixXd::Ones(S,1);
           ftsp_midpoint.block(i*(S+F),1,S,1) = ftsp_and_timings(i,1)*Eigen::MatrixXd::Ones(S,1);
           ftsp_midpoint.block(i*(S+F),2,S,1) = ftsp_and_timings(i,2)*Eigen::MatrixXd::Ones(S,1);
           // double support phases
           ftsp_midpoint.block(i*(S+F)+S,0,F,1) = ftsp_and_timings(i,0)*p_dst + (ftsp_and_timings(i+1,0)-ftsp_and_timings(i,0))*double_support_transition;
           ftsp_midpoint.block(i*(S+F)+S,1,F,1) =  ftsp_and_timings(i,1)*p_dst + (ftsp_and_timings(i+1,1)-ftsp_and_timings(i,1))*double_support_transition;
           ftsp_midpoint.block(i*(S+F)+S,2,F,1) =  ftsp_and_timings(i,2)*p_dst + (ftsp_and_timings(i+1,2)-ftsp_and_timings(i,2))*double_support_transition;
       }
    }


    // Swing foot trajectory generator
    footstepPredicted << footstepsOptimalX(0), footstepsOptimalY(0), 0.0, 0.0; //ftsp_and_timings(fsCount+0,0), ftsp_and_timings(fsCount+0,1), 0.0, 0.0;
    double timeSinceFootstepStart = (walkState.controlIter) * controlTimeStep;
    double singleSupportEnd = ((double)(S)*(mpcTimeStep/controlTimeStep)) * controlTimeStep;
    double swingFootHeight = -(4*stepHeight/pow(singleSupportEnd,2)) * timeSinceFootstepStart * (timeSinceFootstepStart - singleSupportEnd);
    int samplesTillNextFootstep = (S+F)*(mpcTimeStep/controlTimeStep) - walkState.controlIter;

    if (walkState.footstepCounter > 1) {

      if (walkState.footstepCounter <= 1 || samplesTillNextFootstep <= F) swingFootHeight = 0.00;

      // If support foot is left, move right foot
      if (walkState.supportFoot == 1) {

         if (samplesTillNextFootstep > F) {
	   next.rightFootPos += kSwingFoot * (footstepPredicted.head(3) - current.rightFootPos);
         }

         next.rightFootPos(2) = swingFootHeight;
         next.leftFootPos(2) = 0.0;
         if (samplesTillNextFootstep < F*(mpcTimeStep/controlTimeStep)) next.leftFootPos(2) = 0.0;
         else next.leftFootPos(2) = 0.0;
         next.rightFootVel = Eigen::Vector3d::Zero();
         next.rightFootAcc = Eigen::Vector3d::Zero();
         next.rightFootOrient(2) += 5 * timeSinceFootstepStart * kSwingFoot * wrapToPi(footstepPredicted(3) - current.rightFootOrient(2));
         next.rightFootOrient(1) = 0*0.0523599 ;

      } else {

         if (samplesTillNextFootstep > F) {
	   next.leftFootPos += kSwingFoot * (footstepPredicted.head(3) - current.leftFootPos);
         }

         next.leftFootPos(2) = swingFootHeight;
         next.rightFootPos(2) = 0.0;
         if (samplesTillNextFootstep < F*(mpcTimeStep/controlTimeStep)) next.rightFootPos(2) = 0.0;
         else next.rightFootPos(2) = 0.0;
         next.leftFootVel = Eigen::Vector3d::Zero();
         next.leftFootAcc = Eigen::Vector3d::Zero();
         next.leftFootOrient(2) += 5 * timeSinceFootstepStart * kSwingFoot * wrapToPi(footstepPredicted(3) - current.leftFootOrient(2));

      }

    }

     next.torsoOrient(2) = wrapToPi((next.leftFootOrient(2) + next.rightFootOrient(2)) / 2.0);

    // Return next robot state
    return next;
}

Eigen::MatrixXf MPCSolver::getSwingFootTrajectoryBezier(
    const Eigen::VectorXf& optimalFootstep,
    const Eigen::VectorXf& swingFootStartingPosition,
    const Eigen::VectorXf& swingFootActualPosition,
    float stepHeight, float singleSupportDuration,
    float timeSim, float ti) 
{
 
  Eigen::VectorXf swingFootPosition(6);
  Eigen::VectorXf swingFootVelocity(6);
  swingFootPosition.setZero();
  swingFootVelocity.setZero();
   
  float t = timeSim - ti;
  float s, s_dot;

  if (timeSim < ti) {
    // First DS configuration, feet do not move.
    swingFootPosition.head<3>() = swingFootStartingPosition.head<3>();
    swingFootPosition.tail<3>().setZero();

    swingFootVelocity.setZero();
  } else if (timeSim >= ti + singleSupportDuration) {
    // DS configuration, swing foot stays still.
    swingFootPosition = swingFootActualPosition;
    swingFootPosition(3) = 0;
    swingFootPosition(4) = 0;

    swingFootVelocity.setZero();
  } else {
    // SS configuration, swing foot moves.
    // Compute control points:
    Eigen::Vector3f P0 = swingFootStartingPosition.head(3);
    Eigen::Vector3f P2 = optimalFootstep.head(3);
    Eigen::Vector3f P1;
    float diffz = 0.01; 
    if ((P2(2)-P0(2)) > diffz) {
      // Climbing stairs:
      P1 << P0(0) + (P2(0) - P0(0)) * -0.2f, (P0(1)+P2(1))/2.0, stepHeight + P0(2);
    } else if ((P2(2)-P0(2)) < -diffz) {
      // Going down the stairs:
      P1 << P0(0) + (P2(0) - P0(0)) * 1.2f, (P0(1)+P2(1))/2.0, stepHeight + P0(2);
    } else {
      // Even terrain walking:
      P1 << (P0(0)+P2(0))/2.0, (P0(1)+P2(1))/2.0, 0.04 + stepHeight + P0(2);
    }
   
    // NOTE: s \in [0, 1]
    auto quadratic_bezier = [&P0, &P1, &P2](float s) {
      return std::pow(1.0f - s, 2.0f) * P0 +
          2.0f * (1.0 - s) * s * P1 +
          std::pow(s, 2.0f) * P2;
    };

    auto quadratic_bezier_derivative = [&P0, &P1, &P2](float s) {
      return 2.0f * (1.0f - s) * (P1 - P0) + 2.0f * s * (P2 - P1);
    };

    // Trapezoidal speed:
    float alpha = 0.45f;
    auto ta = ti + alpha * singleSupportDuration;
    auto tb = ti + (1.0f - alpha) * singleSupportDuration;
    auto tf = ti + singleSupportDuration;

    auto a_max = std::pow(
        (std::pow(ta, 2.0f) / 2.0f - ti * ta + std::pow(ti, 2.0f) / 2.0f) + 
        (ta - ti) * (tf - ta) -
        (std::pow(tf, 2.0f) / 2.0f - tb * tf + std::pow(tb, 2.0f) / 2.0f),
        -1.0f);
    // auto v_max = a_max * (ta - ti);

    // Timing law with trapezoidal speed:
    // s(t): [ti, tf] -> [0, 1]
    // NOTE: ti < ta < tb < tf
    auto timing_law = [ti, ta, tb, a_max](float t) {
      auto s_0 = [ti, a_max](float t) {
        return a_max *
            (std::pow(t, 2.0f) / 2.0f - ti * t + std::pow(ti, 2.0f) / 2.0f);
      };
      auto s_1 = [ti, ta, a_max, &s_0](float t) {
        return s_0(ta) + a_max * (ta - ti) * (t - ta);
      };
      auto s_2 = [tb, a_max, &s_1](float t) {
        return s_1(t) - a_max *
            (std::pow(t, 2.0f) / 2.0f - tb * t + std::pow(tb, 2.0f) / 2.0f);
      };
      if (t <= ta) {
        return s_0(t);
      } else if(t <= tb) {
        return s_1(t);
      }
      return s_2(t);
    };

    auto timing_law_derivative = [ti, ta, tb, a_max](float t) {
      auto sdot_0 = [ti, a_max](float t) {
        return a_max * (t - ti);
      };
      auto sdot_1 = [ta, a_max, &sdot_0]() {
        return sdot_0(ta);
      };
      auto sdot_2 = [tb, a_max, &sdot_1](float t) {
        return sdot_1() - a_max * (t - tb);
      };
      if (t <= ta) {
        return sdot_0(t);
      } else if (t <= tb) {
        return sdot_1();
      }
      return sdot_2(t);
    };

    // Compute timing law:
    //s = t / singleSupportDuration;
    //s_dot = 1.0f / singleSupportDuration;
    s = timing_law(ti + t);
    s_dot = timing_law_derivative(ti + t);

    // Compute geometric path:
    Eigen::Vector3f b = quadratic_bezier(s);
    Eigen::Vector3f b_dot = quadratic_bezier_derivative(s) * s_dot;

    swingFootPosition << b, 0.0, 0.0, 0.0; //swingFootStartingPosition(5) + s * argos::angle_difference(optimalFootstep(3), swingFootStartingPosition(5))
    swingFootVelocity << b_dot, 0.0, 0.0, 0.0; //argos::angle_difference(optimalFootstep(3), swingFootStartingPosition(5)) * s_dot
  }

  Eigen::MatrixXf ret(6,2);
  ret.col(0) = swingFootPosition;
  ret.col(1) = swingFootVelocity;

  return ret;
}

Eigen::MatrixXd MPCSolver::getFtstpAndTimeMatrix(){

    Eigen::MatrixXd modified_matrix = ftsp_and_timings_;
    return modified_matrix;

}

int MPCSolver::getModifiedControlIter(){

    return itr;

}
