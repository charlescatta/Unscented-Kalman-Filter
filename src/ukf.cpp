#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include "tools.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

#define EPSILON 0.0001

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  is_initialized_ = false;
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

  n_x_ = x_.size();
  n_aug_ = n_x_ + 2;
  n_sig_ = 2 * n_aug_ + 1;
  lambda_ = 3 - n_aug_;

  Xsig_pred_ = MatrixXd(n_x_, n_sig_);
  weights_ = VectorXd(n_sig_);

  double std_radr_2 = std_radr_ * std_radr_;
  double std_radphi_2 = std_radphi_ * std_radphi_;
  double std_radrd_2 = std_radrd_ * std_radrd_;
  double std_laspx_2 = std_laspx_ * std_laspx_;
  double std_laspy_2 = std_laspy_ * std_laspy_;

  R_radar_ = MatrixXd(3, 3);
  R_radar_ <<   std_radr_2,            0,           0,
                         0, std_radphi_2,           0,
                         0,            0, std_radrd_2;
  
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << std_laspx_2,            0,
                        0,  std_laspy_2;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  if (!is_initialized_) {
    P_ = MatrixXd::Identity(5, 5);

    // Initialize weights
    weights_(0) = lambda_ / (n_aug_ + lambda_);
    for (int i = 1; i < weights_.size(); i++) {
      weights_(i) = (1/2) / (n_aug_ + lambda_); 
    }

    // Radar Init
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      double rho = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      double drho_dt = meas_package.raw_measurements_[2];

      // Conversion
      double px = rho * cos(phi);
      double py = rho * sin(phi);
      double vx = drho_dt * cos(phi);
      double vy = drho_dt * sin(phi);
      double vel = sqrt((vx * vx) + (vy * vy));
      
      x_ << px,
            py,
            vel,
            0,
            0;
    }

    // Laser Init
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      // Velocity unknown from first measurement
      double px = meas_package.raw_measurements_[0];
      double py = meas_package.raw_measurements_[1];
      if (fabs(px) < EPSILON) {
        px = EPSILON;
      }
      if (fabs(py) < EPSILON) {
        py = EPSILON;
      }
      x_ << px, 
            py, 
            0, 
            0, 
            0;
    }
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
  }
  time_us_ = meas_package.timestamp_;
  double dt = (meas_package.timestamp_ - time_us_) / 1000000;
  Prediction(dt);
  if (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
  }
  if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  // Augmented Mean
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.fill(0);
  x_aug.head(n_x_) = x_;

  // Sigma Points of Augmented State
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sig_);
  Xsig_aug.col(0) = x_aug;

  // Augmented State Covarience
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0);
  P_aug.topLeftCorner(n_x_, n_x_);
  P_aug(5, 5) = pow(std_a_, 2);
  P_aug(6, 6) = pow(std_yawdd_, 2);

  // Square root of Augmented State Covarience
  MatrixXd P_aug_sqrt = P_aug.llt().matrixL();

  // Sigma Point Creation
  VectorXd sqrt_n_aug_lambda;
  for (int i = 1; i <= n_sig_; i++) {
    sqrt_n_aug_lambda = sqrt(lambda_ + n_aug_) * P_aug_sqrt.col(i - 1);
    Xsig_aug.col(i) = x_aug + sqrt_n_aug_lambda;
    Xsig_aug.col(i + n_aug_) = x_aug - sqrt_n_aug_lambda;
  }

  // Sigma Point Prediction 
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
}

void UKF::CommonUpdate(MeasurementPackage meas_package) {

}
