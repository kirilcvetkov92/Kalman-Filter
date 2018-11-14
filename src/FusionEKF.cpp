#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
    is_initialized_ = false;
    
    previous_timestamp_ = 0;
    
    // initializing matrices
    R_laser_ = MatrixXd(2, 2);
    R_radar_ = MatrixXd(3, 3);
    H_laser_ = MatrixXd(2, 4);
    Hj_ = MatrixXd(3, 4);
    
    //measurement covariance matrix - laser
    R_laser_ << 0.0225, 0,
    0, 0.0225;
    
    //measurement covariance matrix - radar
    R_radar_ << 0.09, 0, 0,
                0, 0.0009, 0,
                0, 0, 0.09;
    
    

    
    /**
     * Set the matrices which are the same for radar and lidar
     */
    
    //Set the acceleration noise components
    noise_ax=5;
    noise_ay=5;
    
    // Set the process and noises matrices
    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ << 1, 0, 0, 0,
               0, 1, 0, 0,
               0, 0, 1000, 0,
               0, 0, 0, 1000;
    
    H_laser_ << 1, 0, 0, 0,
                0, 1, 0, 0;
    
    ekf_.F_ = MatrixXd(4, 4);
    ekf_.F_ << 1, 0, 1, 0,
               0, 1, 0, 1,
               0, 0, 1, 0,
               0, 0, 0, 1;
    
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

float FusionEKF::ProcessNoise(const float &dt, const float &sigma, const float &power)
{
    return sigma*pow(dt,power);
}


void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
    
    
    /*****************************************************************************
     *  Initialization
     ****************************************************************************/
    if (!is_initialized_) {
        /**
         TODO:
         * Initialize the state ekf_.x_ with the first measurement.
         * Create the covariance matrix.
         * Remember: you'll need to convert radar from polar to cartesian coordinates.
         */
        // first measurement
        
        cout << "EKF: " << endl;
        ekf_.x_ = VectorXd(4);
        ekf_.x_ << 1, 1, 1, 1;
        
        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            /**
             Convert radar from polar to cartesian coordinates and initialize state.
             */
            float radius = measurement_pack.raw_measurements_[0];
            float angle = measurement_pack.raw_measurements_[1];
            float range_rate = measurement_pack.raw_measurements_[2]; //radial velocity
            
            float px = radius * cos(angle);
            float py = radius * sin(angle);
            
            float vx = range_rate * cos(angle);
            float vy = range_rate * sin(angle);
            
            ekf_.x_ << px, py, vx, vy;
            
        }
        else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
            
            /**
             * Initialize state.
             */
            
            //set the state with the initial location and zero velocity
            ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
        }
        // done initializing, no need to predict or update
        is_initialized_ = true;
        previous_timestamp_ = measurement_pack.timestamp_;
        return;
    }
    
    /*****************************************************************************
     *  Prediction - It's the same for Radar and Lidar
     ****************************************************************************/
    
    /**
     * Update the state transition matrix F according to the new elapsed time.
     - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
     */
    
    //1. Modify the F matrix so that the time is integrated
    
    float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;    //dt - expressed in seconds
    previous_timestamp_ = measurement_pack.timestamp_;
    
    ekf_.F_ << 1, 0, dt, 0,
    0, 1, 0, dt,
    0, 0, 1, 0,
    0, 0, 0, 1;
    
    //2. Set the process covariance matrix Q
    
    float x4 = ProcessNoise(dt, noise_ax, 4)/4;
    float x3 = ProcessNoise(dt, noise_ax, 3)/2;
    float x2 = ProcessNoise(dt, noise_ax, 2);
    float y4 = ProcessNoise(dt, noise_ay, 4)/4;
    float y3 = ProcessNoise(dt, noise_ay, 3)/2 ;
    float y2 = ProcessNoise(dt, noise_ay, 2);
    
    ekf_.Q_ = MatrixXd(4, 4);
    ekf_.Q_ << x4, 0, x3, 0,
    0, y4, 0, y3,
    x3, 0, x2, 0,
    0, y3, 0, y2;
    
    //3. Call the Extended Kalman Filter predict() function
    ekf_.Predict();
    
    /*****************************************************************************
     *  Update
     ****************************************************************************/
    
    /**
     TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
     */
    
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        // Radar updates
        
        // calculate Jacobian based on converted x measurements from polar to eucledian coordinate system
        Hj_ = tools.CalculateJacobian(ekf_.x_);
        // adjust the H and R matrices that suits for Radar measurements
        ekf_.H_ = Hj_;
        ekf_.R_ = R_radar_;
        // update the kalman radar filter parameters based on new measurements
        ekf_.UpdateEKF(measurement_pack.raw_measurements_);

    } else {
        // Laser updates
        //compute the time elapsed between the current and previous measurements
        
        // adjust the H and R matrices that suits for Radar measurements
        ekf_.H_ = H_laser_;
        ekf_.R_ = R_laser_;
        // update the calman laser filter paramenters based on new measurements
        ekf_.Update(measurement_pack.raw_measurements_);
    }
    
    // print the output
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;
}
