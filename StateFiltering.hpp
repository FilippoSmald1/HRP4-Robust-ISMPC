#pragma once


#include <Eigen/Core>
#include "types.hpp"
#include "parameters.cpp"
#include "utils.cpp"
#include <vector>
#include <iostream>
#include <fstream>
#include <chrono>

class StateFiltering {  

      public:

      //ATTENTION: the state is composed by (position, velocity, acceleration, ext force, ext force derivative)
      //ATTENTION: the measurement covariance is on joint encoders (direct kinematics), IMU, fsr (noisy)
      
      StateFiltering(Eigen::Vector3f& state_x0, Eigen::Matrix2f& q_x_process, Eigen::Matrix3f& q_x_measurement, 
                     Eigen::Vector3f& state_y0, Eigen::Matrix2f& q_y_process, Eigen::Matrix3f& q_y_measurement,
                     Eigen::Vector3f& state_z0, Eigen::Matrix2f& q_z_process, Eigen::Matrix3f& q_z_measurement, 
                     float& h_com, float& mass, float& sampling_time); 
      ~StateFiltering(); // if the distructor is explicitly defined ---> define also:
      StateFiltering(StateFiltering&& var) = default; // move constructor - set to default
      StateFiltering(const StateFiltering& var) = default;  // copy constructor - set to default
      StateFiltering& operator = (StateFiltering&& var) = default; // move assignememnt operator - set to default
      StateFiltering& operator = (const StateFiltering& var) = default; // copy assignement operator - set to default
                                                                                                       
      const float GetTestVariable();                                                                   
      void FilterWithKalman(Eigen::Vector3f& x_measurements, float& input_x,
                            Eigen::Vector3f& y_measurements, float& input_y,
                            Eigen::Vector3f& z_measurements, float& input_z);

      // return position, velocity, acceleration, ext force, covariances ...
      Eigen::VectorXf GetStateX();
      Eigen::VectorXf GetStateY();
      Eigen::VectorXf GetStateZ();

      Eigen::MatrixXf GetCovarianceX();
      Eigen::MatrixXf GetCovarianceY();
      Eigen::MatrixXf GetCovarianceZ();

      Eigen::MatrixXf GetZMP();

      private:

      float test_variable = 3.0f;
      float h_com, mass, sampling_time, input_x, input_y, input_z;
      float g = 9.81f;
      int time_count = 0;
      Eigen::VectorXf state_x, state_y, state_z;
      Eigen::MatrixXf sigma_x, sigma_y, sigma_z;
      Eigen::MatrixXf A, B, C_z, C_xy;
      Eigen::Matrix2f q_x_process_in, q_y_process_in, q_z_process_in; // process input noise covariance
      Eigen::Matrix3f q_x_measurement, q_y_measurement, q_z_measurement; // measurement noise covariance
      Eigen::Vector3f x_measurements, y_measurements, z_measurements;
      void predict_z();
      void update_z();
      void predict_xy();
      void update_xy();   




};


