#include "kalman_filter.h"

#include <math.h>
using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
     * predict the state
  */

	x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;


}

void KalmanFilter::Update(const VectorXd &z) {
  /**
      * update the state by using Kalman Filter equations; taken from lecture notes
  */

	VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;









}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
    * update the state by using Extended Kalman Filter equations
  */


	VectorXd Hx = VectorXd(3);
	double pi = 3.14159;
	float  px = x_(0);
	float  py = x_(1);
	float  vx = x_(2);
	float  vy = x_(3);

	/*Conversion from cartesian to polar coordinates*/
	float c1 = px*px + py*py;

	float range = sqrt(c1);         //magnitude
	float angle_rad = atan2(py,px);	 //phi in radians
	float  range_dot = (px*vx+py*vy)/range;   //range rate

	Hx << range, angle_rad, range_dot;


	VectorXd y = z - Hx;      // kalman filter equation


	/*Normalization of y to keep between -pi and  pi */
	float andle_deg = y(1)*180/pi;    // Conversion to degrees to keep angle between +-180
	andle_deg = fmod(andle_deg + 180,360);
	if (andle_deg < 0)
	    	andle_deg += 360;
	andle_deg = andle_deg - 180; // limit angle between -180 and 180
	y(1) = andle_deg*pi/180; //convert back into radians


	/*KF equations */
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;



}
