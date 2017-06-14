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
  MatrixXd F_ini = MatrixXd(4, 4);
  MatrixXd P_ini = MatrixXd(4, 4);
  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
     * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */

  H_laser_ << 1, 0, 0, 0,
		  0, 1, 0, 0;

  noise_ax = 9;
  noise_ay = 9;


  VectorXd x_ini = VectorXd(4);
  x_ini << 0, 0, 0, 0;


  F_ini<< 1, 0, 1, 0,
			  0, 1, 0, 1,
			  0, 0, 1, 0,
			  0, 0, 0, 1;

  P_ini << 1, 0, 0, 0,
			  0, 1, 0, 0,
			  0, 0, 1000, 0,
			  0, 0, 0, 1000;


  /* Initialize the Kalman filter, Q matrix has same dimensions as P_ini*/
  ekf_.Init(x_ini, P_ini,F_ini,H_laser_,R_laser_,P_ini);


}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**

      * Initialize the state ekf_.x_ with the first measurement.

    */
    // first measurement

	cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */

    	// px = rho* cos( phi)
    	// py = rho* sin( phi)

    	ekf_.x_(0) = measurement_pack.raw_measurements_(0) * cos(float(measurement_pack.raw_measurements_(1)));
    	ekf_.x_(1) = measurement_pack.raw_measurements_(0) * sin(float(measurement_pack.raw_measurements_(1)));
    	ekf_.x_(2) = measurement_pack.raw_measurements_(2) * sin(float(measurement_pack.raw_measurements_(1)));
    	ekf_.x_(3) = measurement_pack.raw_measurements_(2) * sin(float(measurement_pack.raw_measurements_(1)));

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */

    	ekf_.x_(0) = measurement_pack.raw_measurements_(0);

    	ekf_.x_(1) = measurement_pack.raw_measurements_(1);


    }
    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

	float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
	previous_timestamp_ = measurement_pack.timestamp_;


	//1. Modify the F matrix so that the time is integrated
	ekf_.F_ << 1, 0, dt, 0,
			  0, 1, 0, dt,
			  0, 0, 1, 0,
			  0, 0, 0, 1;

	//2. Set the process covariance matrix Q

	float dt_s = dt * dt;
	float dt_c = dt_s * dt ;
	float dt_f = dt_s * dt;

	/*Setting thew new vals for Q matrix*/
	ekf_.Q_ << noise_ax* dt_f/4,  0,              noise_ax* dt_c/4,        0 ,
			  0,              noise_ay* dt_f/4,   0,         noise_ax* dt_c/4,
	         noise_ax* dt_c/4,     0,                   noise_ax* dt_s ,   0,
           	0,              noise_ay* dt_c/4,   0,              noise_ax* dt_s      ;



  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**

     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */
  const VectorXd z = measurement_pack.raw_measurements_;

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    ekf_.R_  = R_radar_;
    //z = measurement_pack.raw_measurements_;
    ekf_.UpdateEKF(z);


  }
  else {
    // Laser updates

	  ekf_.H_  = H_laser_;
	  ekf_.R_  = R_laser_;
	  ekf_.Update(z);

  }





  // print the output
 // cout << "x_ = " << ekf_.x_ << endl;
 // cout << "P_ = " << ekf_.P_ << endl;
}
