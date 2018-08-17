#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
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
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  // Setting up the initialization flag
  is_initialized_ = false;

  // Time when state is true
  time_us_ = 0.0;

  // State Dimension
  n_x_ = 5;

  // Agumented state dimension
  n_aug_ = 7;

  // Sigma Point, lambda parameter
   lambda_ = 3 - n_x_;

  // Predicted Sigma points
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // Create Vector for Weights
  weights_ = VectorXd(2 * n_aug_ + 1);

  NIS_radar_ = 0.0;
  NIS_laser_ = 0.0;
  //Complete the initialization. See ukf.h for other member properties.

  //Hint: one or more values initialized above might be wildly off...
  
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if((meas_package.sensor_type_ ==  MeasurementPackage::RADAR && use_radar_) ||
     (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)) {

       /* Initilization of State Parameters and Covariance Matrix */
       if( !is_initialized_) {
         // Initial Measurment
         x_ << 1, 1, 1, 1, 0.1;

         // Covariane Matrix
         P_ << 0.15, 0, 0, 0, 0,
               0, 0.15, 0, 0, 0,
               0, 0, 1, 0, 0,
               0, 0, 0, 1, 0,
               0, 0, 0, 0, 1;

         // Init timestamp
         time_us_ = meas_package.timestamp_;

         // Reading initial position of the vehicle
         if(meas_package.sensor_type_ ==  MeasurementPackage::LASER && use_laser_) {
            x_(0) = meas_package.raw_measurements_(0);
            x_(1) = meas_package.raw_measurements_(1);
         } else if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
           // Polar to Cartesian coordinates conversion
           float ro = meas_package.raw_measurements_(0);
           float phi = meas_package.raw_measurements_(1);
           float ro_dot = meas_package.raw_measurements_(2);
           x_(0) = ro * cos(phi);
           x_(1) = ro * sin(phi);
         }
         is_initialized_ = true;
         return;
    }

    // Prediction Step
    // Time Elapsed between the current and previous measurments
    float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
    time_us_ = meas_package.timestamp_;

    Prediction(dt);

    // Update Step
     if(meas_package.sensor_type_ == MeasurementPackage::LASER) {
       UpdateLidar(meas_package);
     } else if(meas_package.sensor_type_ == MeasurementPackage::RADAR) {
       UpdateRadar(meas_package);
     }
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
 
  // Generating Sigma Points matrix
  MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);
  
  // Calculate Sqrt of P
  MatrixXd A = P_.llt().matrixL();

  // Lambda for non-agumented Sigma Points
  lambda_ = 3 - n_x_;
  
  //Set first Col of Sigma Points
  Xsig.col(0) = x_;

  // Remaing Sigma Points
  for(int i = 0; i < n_x_; i++)
  {
    Xsig.col(i + 1) = x_ +  sqrt(lambda_ + n_x_) * A.col(i);
    Xsig.col(i + 1 + n_x_) = x_ - sqrt(lambda_ + n_x_) * A.col(i);
  }

  // Agument Sigma Points
  // Mean Vector
  VectorXd x_aug = VectorXd(n_aug_);
 
  // State Covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  // Sigma Point Matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  // Lambda for Agumented Sigma Points
  lambda_ = 3 - n_aug_;

  // Mean State
  x_aug.head(5) = x_;
  x_aug(5) = 0;  
  x_aug(6) = 0;

  // Covariance Matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5, 5) = P_;
  P_aug(5, 5) = std_a_ * std_a_;
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;

  // Sqrt Matrix
  MatrixXd L = P_aug.llt().matrixL();

  // Agumented Sigma Points
  Xsig_aug.col(0) = x_aug;
  for(int i = 0; i < n_aug_; i++)
  {
    Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ +  n_aug_) * L.col(i);
  }

  // Predict Sigma Points
  for(int i = 0; i < 2 * n_aug_ + 1; i ++)
  {
    double p_x  = Xsig_aug(0, i);
    double p_y  = Xsig_aug(1, i);
    double v    = Xsig_aug(2, i);
    double yaw  = Xsig_aug(3, i);
    double yawd = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);

    // State Values
    double px_pred, py_pred;
    if(fabs(yawd) > 0.001) {
     px_pred = p_x + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
     py_pred = p_y + v / yawd * (cos(yaw + yawd * delta_t));
    } else {
     px_pred = p_x + v * delta_t * cos(yaw);
     py_pred = p_y + v * delta_t * sin(yaw);
   }

   double v_pred, yaw_pred, yawd_pred;
   v_pred = v;
   yaw_pred - yaw + yawd * delta_t;
   yawd_pred = yawd;

   // Adding noise parameter
   px_pred = px_pred + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
   py_pred = py_pred + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
   v_pred = v_pred + nu_a * delta_t;
   
   yaw_pred = yaw_pred + 0.5 * nu_yawdd * delta_t * delta_t;
   yawd_pred = yawd_pred + nu_yawdd * delta_t;

   // Predicted Sigma Points
   Xsig_pred_(0, i) = px_pred;
   Xsig_pred_(1, i) = py_pred;
   Xsig_pred_(2, i) = py_pred;
   Xsig_pred_(3, i) = yaw_pred;
   Xsig_pred_(4, i) = yawd_pred;
  }
  
  // Predicted Sigma Points to Mean and Convariance
  //Weights
  double w_0 = lambda_ / (lambda_ + n_aug_);
  weights_(0) = w_0;
  for(int i = 1; i < 2 * n_aug_ + 1; i++) {
    double w = 0.5 / (n_aug_ + lambda_);
    weights_(i) =  w;
  }

  // Stea Mean
  x_.fill(0.0);
  for(int i = 0; i < 2 * n_aug_ + 1; i++) {
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  // Covariance Matrix
  P_.fill(0.0);
  for(int i = 0; i < 2 * n_aug_ + 1; i++) {
   VectorXd x_diff = Xsig_pred_.col(i)  - x_;
   while( x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
   while( x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;
   P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  
  //Extracting measurements
  VectorXd z = meas_package.raw_measurements_;

  //Measurement dimension -  p_x and p_y
  int n_z = 2;

  //Sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //Transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) { 

    double pred_x = Xsig_pred_(0, i);
    double pred_y = Xsig_pred_(1, i);

    Zsig(0, i) = pred_x;
    Zsig(1, i) = pred_y;
  }

  //Mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //Measurement covariance matrix
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) { 
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //Adding measurement noise 
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_laspx_*std_laspx_, 0,
       0, std_laspy_*std_laspy_;
  S = S + R;

  //Matrix for cross correlation
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  // UKF Update for Lidar
  //Calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //calculate NIS
  NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;

  // Updated state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  //Extracting measurement
  VectorXd z = meas_package.raw_measurements_;

  //Masurement dimension - r, phi, and r_dot
  int n_z = 3;

  //Matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //Transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  

    double pred_x = Xsig_pred_(0, i);
    double pred_y = Xsig_pred_(1, i);
    double v   = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model - r, phi and r_dot
    Zsig(0, i) = sqrt(pred_x*pred_x + pred_y*pred_y);                       
    Zsig(1, i) = atan2(pred_y, pred_x);                                 
    Zsig(2, i) = (pred_x*v1 + pred_y*v2) / sqrt(pred_x*pred_x + pred_y*pred_y);  
  }

  //Mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //Measurement covariance matrix
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //Normalization of angle
    while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //Measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_radr_*std_radr_, 0, 0,
       0, std_radphi_*std_radphi_, 0,
       0,   0, std_radrd_*std_radrd_;
  S = S + R;

  //Mtrix for cross correlation
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  // UKF Update for Radar
  // Cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  

    VectorXd z_diff = Zsig.col(i) - z_pred;
    while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  VectorXd z_diff = z - z_pred;

  while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

  //Calculating NIS
  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;

  //Update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
}
