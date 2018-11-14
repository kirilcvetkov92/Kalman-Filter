#ifndef FusionEKF_H_
#define FusionEKF_H_

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include "kalman_filter.h"
#include "tools.h"

class FusionEKF {
public:
    /**
     * Constructor.
     */
    FusionEKF();
    
    /**
     * Destructor.
     */
    virtual ~FusionEKF();
    
    /**
     * Run the whole flow of the Kalman Filter from here.
     */
    void ProcessMeasurement(const MeasurementPackage &measurement_pack);
    
    /**
     * Kalman Filter update and prediction math lives in here.
     */
    
    
    /**
     * Laser function that help to calculate Q matrix parameters
     */
    float ProcessNoise(const float &dt, const float &sigma, const float &power);
    
    KalmanFilter ekf_;
    
private:
    // check whether the tracking toolbox was initialized or not (first measurement)
    bool is_initialized_;
    
    float noise_ax;
    float noise_ay;
    
    // previous timestamp
    long long previous_timestamp_;
    
    // tool object used to compute Jacobian and RMSE
    Tools tools;
    Eigen::MatrixXd R_laser_;
    Eigen::MatrixXd R_radar_;
    Eigen::MatrixXd H_laser_;
    Eigen::MatrixXd Hj_;
};

#endif /* FusionEKF_H_ */
