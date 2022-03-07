#pragma once

#include <Eigen/Core>
#include "qpOASES/qpOASES.hpp"
#include "types.hpp"
#include "parameters.cpp"
#include "utils.cpp"
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>     
#include <chrono>

class MPCSolver{
    public:
    MPCSolver(const Eigen::MatrixXd& ftsp_and_timings);
    ~MPCSolver();

    // Compute the next desired state starting from the current state
    State solve(State current, WalkState walkState, const Eigen::MatrixXd& ftsp_and_timings);
    Eigen::MatrixXd getFtstpAndTimeMatrix();
    Eigen::MatrixXf getSwingFootTrajectoryBezier(
    const Eigen::VectorXf& optimalFootstep,
    const Eigen::VectorXf& swingFootStartingPosition,
    const Eigen::VectorXf& swingFootActualPosition,
    float stepHeight, float singleSupportDuration,
    float timeSim, float ti);

    int getModifiedControlIter();
    Eigen::MatrixXd ftsp_and_timings_;

    private:
    // Matrices for prediction
    Eigen::VectorXd p;
    Eigen::MatrixXd P;

    // Matrices for cost function
    Eigen::MatrixXd costFunctionH;
    Eigen::VectorXd costFunctionF;

    // Matrices for stability constraint
    Eigen::MatrixXd Aeq;
    Eigen::VectorXd beq;

    //Matrices for balance constraint
    Eigen::MatrixXd AZmp;
    Eigen::VectorXd bZmpMax;
    Eigen::VectorXd bZmpMin;

    // Matrices for kinematic constraints
    Eigen::MatrixXd AFootsteps;
    Eigen::VectorXd bFootstepsMax;
    Eigen::VectorXd bFootstepsMin;

    // Matrices for double support constraints
    Eigen::MatrixXd ADoubleSupport;
    Eigen::VectorXd bDoubleSupport;

    // Matrices for fixed footstep constraints
    Eigen::MatrixXd AFixedFootsteps;
    Eigen::VectorXd bFixedFootsteps;

    // Matrices for the stacked constraints
    Eigen::MatrixXd AConstraint;
    Eigen::VectorXd bConstraintMax;
    Eigen::VectorXd bConstraintMin;

    // Matrices for Z QP and constraints
    Eigen::VectorXd lambda, etas;
    Eigen::MatrixXd deltas;
    Eigen::VectorXd optimal_restriction_x, optimal_restriction_y;
    Eigen::MatrixXd A_virtual_ZMP;
    Eigen::VectorXd b_virtual_ZMP_x, b_virtual_ZMP_y;
    Eigen::Vector4d footstepPredicted;

    // Matrices
    Eigen::MatrixXd halfAZmp, CcFull, Cc, Ic, differenceMatrix, halfH;
    Eigen::VectorXd pFr, pFl, pF, currentFootstepKinematic, xCandidateFootsteps, yCandidateFootsteps, footstepsOptimalX, footstepsOptimalY;
    Eigen::VectorXd decision_variables, fix_in_double_support;
    Eigen::MatrixXd cost_dist_x, A_dist_x, cost_dist_y, A_dist_y;
    Eigen::VectorXd b_dist_min_x, b_dist_min_y, b_dist_max_x, b_dist_max_y, f_dist_x, f_dist_y;
    Eigen::VectorXd w_bar;
    Eigen::MatrixXd exponential_weight;
    Eigen::VectorXd restriction, res;

    //Matrices State update
    Eigen::MatrixXd A_upd = Eigen::MatrixXd::Zero(3,3);
    Eigen::VectorXd B_upd = Eigen::VectorXd::Zero(3);
    
    // Midpoint of ZMP constraint
    Eigen::MatrixXd ftsp_midpoint; 

    // useful variables
    double anticipative_x, anticipative_y;
    double alpha = 2.0;
    double xu_offset = 0.011;
    int sta_counter = 0;
    double x_u_m, x_u_M, y_u_m, y_u_M;
    int nConstraints;
    int itr;
    int fsCount, old_fsCount, adaptation_memo, ds_samples, ct;
    double xz_dot, yz_dot;
    
    // flags
    double activate_fixed_footstep_constraint = 1.0;
    double activate_restriction = 1.0;
    double activate_observer_x = 0.0;
    double activate_midrange_x = 1.0;
    double activate_observer_y = 0.0;
    double activate_midrange_y = 0.0;
    bool recovery_sim = true;
    bool activate_feasibility_recovery = true;
    int sta_limit = 3;
    int measured_state_update = 0;

};
