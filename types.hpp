#pragma once

#include <Eigen/Core>
#include "utils.cpp"

// Contains the state of the LIP robot
struct State {
    Eigen::Vector3d comPos;
    Eigen::Vector3d comVel;
    Eigen::Vector3d comPosmeas;
    Eigen::Vector3d comVelmeas;
    Eigen::Vector3d comAcc;
    Eigen::Vector3d zmpPos;
    Eigen::Vector3d disturbance;
    Eigen::Vector2d zmpDot;
    Eigen::Vector4d obsState_x;
    Eigen::Vector4d obsState_y;

    Eigen::Vector3d leftFootPos;
    Eigen::Vector3d leftFootVel;
    Eigen::Vector3d leftFootAcc;

    Eigen::Vector3d rightFootPos;
    Eigen::Vector3d rightFootVel;
    Eigen::Vector3d rightFootAcc;

    Eigen::Vector3d torsoOrient;
    Eigen::Vector3d leftFootOrient;
    Eigen::Vector3d rightFootOrient;

    inline Eigen::VectorXd getComPose() {
	Eigen::VectorXd comPose(6);
        comPose << torsoOrient, comPos;
        return comPose;
    }

    inline Eigen::VectorXd getSupportFootPose(bool supportFoot) {
	Eigen::VectorXd sfPose(6);
        if (supportFoot == 0) sfPose << leftFootOrient, leftFootPos;
        else sfPose << rightFootOrient, rightFootPos;
        return sfPose;
    }

    inline Eigen::VectorXd getSwingFootPose(bool supportFoot) {
	Eigen::VectorXd sfPose(6);
        if (supportFoot == 1) sfPose << leftFootOrient, leftFootPos;
        else sfPose << rightFootOrient, rightFootPos;
        return sfPose;
    }

    inline Eigen::VectorXd getRelComPose(bool supportFoot) {
	return vvRel(getComPose(), getSupportFootPose(supportFoot));
    }

    inline Eigen::VectorXd getRelSwingFootPose(bool supportFoot) {
	return vvRel(getSwingFootPose(supportFoot), getSupportFootPose(supportFoot));
    }
};

struct WalkState {
    bool supportFoot;
    double simulationTime;
    int mpcIter, controlIter, footstepCounter, indInitial;
};

struct Vref {
    Vref(double _x, double _y, double _omega) : x(_x), y(_y), omega(_omega) {}

    double x = 0;
    double y = 0;
    double omega = 0;
}; 
