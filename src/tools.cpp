#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
    * Calculate the RMSE
  */
    
    VectorXd rmse(4);
    rmse << 0,0,0,0;
    
    // check the validity of the following inputs:
    // the estimation vector size should not be zero

    if(estimations.size() != ground_truth.size()
       || estimations.size() == 0){
        cout << "Invalid estimation or ground_truth data" << endl;
        return rmse;
    }
    //the estimation vector size should equal ground truth vector size
    //accumulate squared residuals
    
    for(int i=0; i < estimations.size(); ++i){
        VectorXd error = (estimations[i] - ground_truth[i]);
        error =  (error.array() * error.array());
        rmse += error;
    }
    
    rmse  = (rmse.array()/estimations.size()).sqrt();
    return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
    
    
    MatrixXd Hj(3,4);
    //recover state parameters
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);
    
    float divisor = px*px+py*py;
    float d = sqrt(divisor);
    
    //check division by zero
    if (fabs(divisor)<1e-8)
    {
        cout << "CalculateJacobian () - Error - Division by Zero" << endl;
        return Hj;
    }
    
    
    Hj << px/d, py/d, 0, 0,
         -py/(divisor), px/(divisor), 0, 0,
         (py*(vx*py-vy*px))/pow(divisor,3/2.0f), (px*(vy*px-vx*py))/pow((divisor),3/2.0f), px/d, py/d;
    
    
    return Hj;
}
