#include "StateFiltering.hpp"

StateFiltering::StateFiltering(Eigen::Vector3f& state_x0, Eigen::Matrix2f& q_x_process_in_, Eigen::Matrix3f& q_x_measurement_, 
                               Eigen::Vector3f& state_y0, Eigen::Matrix2f& q_y_process_in_, Eigen::Matrix3f& q_y_measurement_,
                               Eigen::Vector3f& state_z0, Eigen::Matrix2f& q_z_process_in_, Eigen::Matrix3f& q_z_measurement_, 
                               float& h_com_, float& mass_, float& sampling_time_) 
{

                h_com = h_com_;
                mass = mass_;
                sampling_time = sampling_time_;
                q_x_process_in = q_x_process_in_;
                q_y_process_in = q_y_process_in_;
                q_z_process_in = q_z_process_in_;
                q_x_measurement = q_x_measurement_;
                q_y_measurement = q_y_measurement_;  
                q_z_measurement = q_z_measurement_;  

                state_x = Eigen::VectorXf::Zero(5);
                state_y = Eigen::VectorXf::Zero(5);
                state_z = Eigen::VectorXf::Zero(5);
                state_x << state_x0, 0.0f, 0.0f;
                state_y << state_y0, 0.0f, 0.0f;
                state_z << state_z0, 0.0f, 0.0f;

                A = Eigen::MatrixXf::Zero(5,5);
                B = Eigen::MatrixXf::Zero(5,2);
                C_z = Eigen::MatrixXf::Zero(3,5);
                C_xy = Eigen::MatrixXf::Zero(3,5);
                sigma_x = Eigen::MatrixXf::Identity(5,5); // set to zero for initialization...
                sigma_y = Eigen::MatrixXf::Identity(5,5); // set to zero for initialization...
                sigma_z = (Eigen::MatrixXf::Identity(5,5)); // set to zero for initialization...



                A << 1.0f, sampling_time, sampling_time*sampling_time/2,          0.0f,          0.0f,
                     0.0f,          1.0f,                 sampling_time, sampling_time,          0.0f,
                     0.0f,          0.0f,                          1.0f,          0.0f,          0.0f, 
                     0.0f,          0.0f,                          0.0f,          1.0f, sampling_time, 
                     0.0f,          0.0f,                          0.0f,          0.0f,          1.0f;

                B <<  sampling_time*sampling_time*sampling_time/6,                          0.0f,
                                    sampling_time*sampling_time/2,                          0.0f,
                                                    sampling_time,                          0.0f,
                                                             0.0f, sampling_time*sampling_time/2,
                                                             0.0f,                 sampling_time;

               C_z << 1.0f, 0.0f,  0.0f, 0.0f, 0.0f,
                      0.0f, 0.0f,  1.0f, 0.0f, 0.0f,  
                      0.0f, 0.0f, -mass, 1.0f, 0.0f;  

               C_xy << 1.0f, 0.0f, 0.0f, 0.0f, 0.0f,
                       0.0f, 0.0f, 1.0f, 0.0f, 0.0f,  
                       1.0f, 0.0f, 0.0f, 0.0f, 0.0f; 

}

StateFiltering::~StateFiltering()
{

  std::cout << "Filter object deleted. Total steps = " << time_count << std::endl;


}

void StateFiltering::FilterWithKalman(Eigen::Vector3f& x_measurements_, float& input_x_,
                                      Eigen::Vector3f& y_measurements_, float& input_y_,
                                      Eigen::Vector3f& z_measurements_, float& input_z_)
{


                x_measurements = x_measurements_;
                y_measurements = y_measurements_;
                z_measurements = z_measurements_;
                input_x = input_x_;
                input_y = input_y_;
                input_z = input_z_;
                
                predict_z();
                update_z();
                
                predict_xy();
                update_xy();

                time_count++;

}



void StateFiltering::predict_z()
{

                state_z  = A*state_z + B*Eigen::Vector2f(input_z,0.0f);
                sigma_z  = A*sigma_z*A.transpose() + B*q_z_process_in*B.transpose();

}
void StateFiltering::update_z()
{
                //std::cout << C_z*state_z << std::endl;                  
                Eigen::MatrixXf K_gain_z = sigma_z*C_z.transpose()*((q_z_measurement+C_z*sigma_z*C_z.transpose()).inverse());  
                state_z = state_z + K_gain_z*(z_measurements - (C_z*state_z + Eigen::Vector3f(0.0f,0.0f,-g*mass)));
                sigma_z = sigma_z - K_gain_z*C_z*sigma_z;

}              
               
 
void StateFiltering::predict_xy()
{
                 
                state_x  = A*state_x + B*Eigen::Vector2f(input_x,0.0f);
                sigma_x  = A*sigma_x*A.transpose() + B*q_x_process_in*B.transpose();

                state_y  = A*state_y + B*Eigen::Vector2f(input_y,0.0f);
                sigma_y  = A*sigma_y*A.transpose() + B*q_y_process_in*B.transpose();
              
}
void StateFiltering::update_xy()
{

                float f_n = - mass*g - mass*state_z(2) + state_z(3);
                C_xy(2,2) =  mass*state_z(0)/f_n;
                C_xy(2,3) = - state_z(0)/f_n;

                Eigen::MatrixXf K_gain_x = sigma_x*C_xy.transpose()*((q_x_measurement+C_xy*sigma_x*C_xy.transpose()).inverse());              
                state_x = state_x + K_gain_x*(x_measurements - C_xy*state_x);
                sigma_x = sigma_x - K_gain_x*C_xy*sigma_x;

                Eigen::MatrixXf K_gain_y = sigma_y*C_xy.transpose()*((q_y_measurement+C_xy*sigma_y*C_xy.transpose()).inverse());              
                state_y = state_y + K_gain_y*(y_measurements - C_xy*state_y);
                sigma_y = sigma_y - K_gain_y*C_xy*sigma_y;

}

// GETTERS

const float StateFiltering::GetTestVariable() {

     return test_variable;                

}

Eigen::VectorXf StateFiltering::GetStateX(){

     return state_x;                

}
Eigen::MatrixXf StateFiltering::GetCovarianceX(){

     return sigma_x;                

}

Eigen::VectorXf StateFiltering::GetStateY(){

     return state_y;                

}
Eigen::MatrixXf StateFiltering::GetCovarianceY(){

     return sigma_y;                

}

Eigen::VectorXf StateFiltering::GetStateZ(){

     return state_z;                

}
Eigen::MatrixXf StateFiltering::GetCovarianceZ(){

     return sigma_z;                

}

Eigen::MatrixXf StateFiltering::GetZMP(){

     Eigen::MatrixXf result = Eigen::MatrixXf::Zero(2,1);
     result << C_xy.block(2,0,1,5)*state_x,C_xy.block(2,0,1,5)*state_y;
     return result;                

}
