#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;


/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
  
  // Default to false, this flag is being set once the ProcessMeasurement() has been called for the fist time 
  is_initialized_ = false; 
  
  // Time stamp is set to zero by default
  time_us_ = 0; 
  
  // Set Sigma point creation parameters
  n_x_ = 5; 
  n_aug_ = 7; 
  lambda_ = 3 - n_aug_; 
  
  // Create empty sigma point matrix (for augemented states)
  Xsig_pred_ = MatrixXd(n_x_ , 2*n_aug_ + 1);
  
  // Create sigma points weight vector
  weights_ = VectorXd(2*n_aug_ + 1); 
  double w0 = lambda_ / (lambda_ + n_aug_);
  weights_(0) = w0; 
  for (int i = 1; i < 2*n_aug_ + 1; ++i){
  	double w = 0.5 / (lambda_ + n_aug_); 
    weights_(i) = w; 
  } 
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  
  // I N I T 
  // Set is_initialized flag to true once this method is called for the first time (i.e the first measurement was received)
  if (! is_initialized_){
	// Initialize the position (x and y values) based on the first received measurement
    // Distinguish RADAR and LIDAR measurments 
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      float rho = meas_package.raw_measurements_[0];  
      float phi = meas_package.raw_measurements_[1]; 
      float phi_dot = meas_package.raw_measurements_[2]; 
      float px = rho * cos(phi); 
      float py = rho * sin(phi);
      x_(0) = px;  
      x_(1) = py;
      x_(2) = 0;  
      x_(3) = 0;
      x_(4) = 0;
      //x_(2) = phi_dot;  
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      float px = meas_package.raw_measurements_[0]; 
      float py = meas_package.raw_measurements_[1]; 
      x_(0) = px;  
      x_(1) = py; 
      x_(2) = 0;  
      x_(3) = 0;
      x_(4) = 0;  
    }
    
    // Initialize process (state) covariance matrix P
    P_.fill(0.0); 
    P_(0,0) = 1;  // px 
    P_(1,1) = 1;  // py
    P_(2,2) = 800;  // v
    P_(3,3) = 800; // phi
    P_(4,4) = 1; // phi dot 
    
    // Store the time stamp of the first measurment
    time_us_ = meas_package.timestamp_;  
    
    // Make sure that the position initialization based on the raw measurements is only done once  
    is_initialized_ = true; 
    return;  
  }
  
  // P R E D I C T I O N
  float delta_t = (meas_package.timestamp_ - time_us_) * 1e-6; 
  Prediction(delta_t);
  
  // U P D A T E
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_){
  	UpdateRadar(meas_package);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
    UpdateLidar(meas_package);
  }
  
  // Store current timestamp for delta t calculation in the next measurement process iteration
  time_us_ = meas_package.timestamp_;   
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
  
  // 1.  C r e a t e   S i g m a   P o i n t s
  // create augmented mean vector
  VectorXd x_aug = VectorXd(7);
  x_aug.head(5) = x_; 
  x_aug(5) = 0; 
  x_aug(6) = 0; 
  // Create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;
  // create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  // create square root matrix
  MatrixXd L = P_aug.llt().matrixL();
  // create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; ++i) {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }
  
  // 2.  P r e d i c t   S i g m a   P o i n t s   
  // (from time k to time k+1 using the CTVR model)
  for (int i = 0; i< 2*n_aug_+1; ++i) {
    // extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);
    // predicted state values
    double px_p, py_p;
    // avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    } else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }
    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;
    // add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;
    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;
    // write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
  
  // 3.  R e c r e a t e   m e a n   a n d   c o v a r i a n c e 
  // predicted state mean
  // create vector for predicted state
  VectorXd x = VectorXd(n_x_);
  x.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  
    x = x + weights_(i) * Xsig_pred_.col(i);
  }
  // predicted state covariance matrix
  MatrixXd P = MatrixXd(n_x_, n_x_); 
  P.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
	// Sum up P
    P = P + weights_(i) * x_diff * x_diff.transpose() ;
  } 
  // Writing results back to attributes
  x_ = x; 
  P_ = P; 
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
  // 4.  P r e d i c t   M e a s u r e m e n t
  // Define the measurment space dimension
  int z_dim = 2; 
  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(z_dim, 2 * n_aug_ + 1);
  Zsig.fill(0.0); 
  // mean predicted measurement
  VectorXd z_pred = VectorXd(z_dim);
  z_pred.fill(0.0); 
  // measurement covariance matrix S
  MatrixXd S = MatrixXd(z_dim, z_dim);
  S.fill(0.0);
  // transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  
    // extract values for better readability
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    // measurement model
    Zsig(0,i) = p_x;                       // x position is directly measured by lidar
    Zsig(1,i) = p_y;                       // y position is directly measured by lidar
  }
  // Recreate mean from predicted measurement sigma points
  for (int i=0; i < 2*n_aug_ + 1; ++i) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }
  // Recreate innovation covariance matrix S from predicted measurment sigma points
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  // add measurement noise covariance matrix
  MatrixXd R = MatrixXd(z_dim, z_dim);
  R <<  std_laspx_*std_laspx_, 0, 
        0, std_laspy_*std_laspy_;
  S = S + R;
  
  // 5.  U P D A T E 
  // Get lidar measurments
  VectorXd z = VectorXd(z_dim); 
  z(0) = meas_package.raw_measurements_[0]; // x position in m
  z(1) = meas_package.raw_measurements_[1]; // y position in m 
  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, z_dim);
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
	// Sum up Tc
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  // Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  // residual
  VectorXd z_diff = z - z_pred;
  
  // update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();  
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
  // 4.  P r e d i c t   M e a s u r e m e n t 
  // Define the measurement space dimension
  int z_dim = 3; 
  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(z_dim, 2 * n_aug_ + 1);
  // mean predicted measurement
  VectorXd z_pred = VectorXd(z_dim);
  z_pred.fill(0.0);
  // measurement covariance matrix S
  MatrixXd S = MatrixXd(z_dim, z_dim);
  S.fill(0.0); 
  // transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  
    // extract values for better readability
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);
    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;
    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                       // r
    Zsig(1,i) = atan2(p_y,p_x);                                // phi
    Zsig(2,i) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y);   // r_dot
  }
  // Recreate mean from predicted measurement sigma points
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_ + 1; ++i) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }
  // Recreate innovation covariance matrix S from predicted measurment sigma points
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
	// Calculate S
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  // add measurement noise covariance matrix
  MatrixXd R = MatrixXd(z_dim, z_dim);
  R <<  std_radr_*std_radr_, 0, 0,
        0, std_radphi_*std_radphi_, 0,
        0, 0,std_radrd_*std_radrd_;
  S = S + R;
  
  // 5.  U P D A T E 
  // Get radar measurments
  VectorXd z = VectorXd(z_dim); 
  z(0) = meas_package.raw_measurements_[0]; // rho in m
  z(1) = meas_package.raw_measurements_[1]; // phi in rad
  z(2) = meas_package.raw_measurements_[2]; // phi dot in m/s  
  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, z_dim);
  // calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
	// Sum up Tc
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  // Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  // residual
  VectorXd z_diff = z - z_pred;
  // angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
  // update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose(); 
}