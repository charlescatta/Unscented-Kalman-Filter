#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include "tools.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

#define EPSILON 0.0001

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
  for (int i = 0; i < n_sig_; i++) {
    MatrixXd x = MatrixXd(n_aug_, 1);
    MatrixXd CTRV_det = MatrixXd(n_aug_, 1);
    MatrixXd CTRV_nondet = MatrixXd(n_aug_, 1);

    double px       = Xsig_aug(0, i);
    double py       = Xsig_aug(1, i);
    double v        = Xsig_aug(2, i);
    double psi      = Xsig_aug(3, i);
    double psi_dot  = Xsig_aug(4, i);
    double nu_aug   = Xsig_aug(5, i);
    double nu_2_dot = Xsig_aug(6, i);
    double half_dt_2 = 0.5 * pow(delta_t, 2);
    
    x << px,
       py,
       v,
       psi,
       psi_dot;
           
    if (psi_dot == 0) {
      CTRV_det << v * cos(psi) * delta_t,
                  v * sin(psi) * delta_t,
                  0,
                  psi_dot * delta_t,
                  0;
    } else {
        double v_over_psi_dot = v / psi_dot;
        double psi_dot_dt = psi_dot * delta_t;
        // Deterministic part of the CTRV model
        CTRV_det << v_over_psi_dot * (sin(psi + psi_dot_dt) - sin(psi)),
                    v_over_psi_dot * (-cos(psi + psi_dot_dt) + cos(psi)),
                    0,
                    psi_dot_dt,
                    0;
    }
    // Non deterministic part of the CTRV model
    CTRV_nondet <<  half_dt_2 * cos(psi) * nu_aug,
                    half_dt_2 * sin(psi) * nu_aug,
                    delta_t * nu_aug,
                    half_dt_2 * nu_2_dot,
                    delta_t * nu_2_dot;
    
    Xsig_pred_.col(i) = x + CTRV_det + CTRV_nondet; 
  }
  x_ = Xsig_pred_ * weights_;
  P_.fill(0);

  for (int i = 0; i < n_sig_; i++) {
    VectorXd x_delta = Xsig_pred_.col(i) = x_;
    Tools::NormalizeAngle(&x_delta(3));
    P_ += weights_(i) * x_delta * x_delta.transpose();
  } 
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  int n = 2;
  MatrixXd sig = Xsig_pred_.block(0, 0, n, n_sig_);
  CommonUpdate(meas_package, n, sig);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  int n = 3;
  MatrixXd sig = MatrixXd(n, n_sig_);

  for (int i = 0; i < n_sig_; i++) {
    double px       = Xsig_pred_(0, i);
    double py       = Xsig_pred_(1, i);
    double v        = Xsig_pred_(2, i);
    double psi      = Xsig_pred_(3, i);
    sig(0, i) = sqrt(px * px + py * py);
    sig(1, i) = atan2(py, px);
    sig(2, i) = (px * cos(psi) * v) + (py * sin(psi) * v);
  }
  CommonUpdate(meas_package, n, sig);
}

/**
 * Updates the state and the state covariance matrix using the common KF update procedures.
 * @param {MeasurementPackage} meas_package
 * @param {n} input measurement dimension
 * @param {sig} sigma point matrix
 */
void UKF::CommonUpdate(MeasurementPackage meas_package, int n, MatrixXd sig) {
  // Predicted Mean
  VectorXd z_pred = VectorXd(n);

  // Cross Correlation Matrix
  MatrixXd TC = MatrixXd(n_x_, n);
  TC.fill(0);

  // Measurement Noise Covariance Matrix
  MatrixXd R = MatrixXd(n, n);
  
  if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    R = R_laser_;
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    R = R_radar_;
  }

  // Measurement Covariance Matrix
  MatrixXd S = MatrixXd(n, n);
  S.fill(0);

  for (int i = 0; i < n_sig_; i++) {
    VectorXd z_delta = sig.col(i) - z_pred;
    Tools::NormalizeAngle(&z_delta(1));
    S += weights_(i) * z_delta * z_delta.transpose();
  }
  S += R;


  // Compute Cross Correlation
  for (int i = 0; i < n_sig_; i++) {
    VectorXd z_delta = sig.col(i) - z_pred;
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      Tools::NormalizeAngle(&z_delta(3));
    }
    VectorXd x_delta = Xsig_pred_.col(i) - z_pred;
    TC += weights_(i) * x_delta * z_delta.transpose();
  }
  // Compute Kalman Gain
  MatrixXd K = TC * S.inverse();
  VectorXd z_residual = meas_package.raw_measurements_ - z_pred;
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    Tools::NormalizeAngle(&z_residual(1));
  }

  // Update
  x_ += K * z_residual;
  P_ += K * S * K.transpose();
}
