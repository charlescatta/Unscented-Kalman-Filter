#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd RMSE(4);
  RMSE.fill(0);
  if (estimations.size() != ground_truth.size()) {
    std::cout << "Estimation and ground truth data are of different size" << std::endl;
    return RMSE;
  }
  for(unsigned int i=0; i < estimations.size(); ++i){
    VectorXd residual = estimations[i] - ground_truth[i];

    // Coefficient-wise multiplication
    residual = residual.array() * residual.array();
    RMSE += residual;
  }
  // Calculate the mean
  RMSE = RMSE / estimations.size();
  RMSE = RMSE.array().sqrt();
  return RMSE;

}

void Tools::NormalizeAngle(double *angle) {
  while (*angle < -M_PI) {
    *angle += M_PI * 2;
  }
  while (*angle > M_PI) {
    *angle -= M_PI * 2;
  }
}