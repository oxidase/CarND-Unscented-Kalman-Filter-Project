#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF()
{
    // initially set to false, set to true in first call of ProcessMeasurement
    is_initialized_ = false;

    // State dimension
    n_x_ = 5;

    // Augmented state dimension
    n_aug_ = 7;

    // Sigma point spreading parameter
    lambda_ = 3 - n_aug_;

    // if this is false, laser measurements will be ignored (except during init)
    use_laser_ = true;

    // if this is false, radar measurements will be ignored (except during init)
    use_radar_ = true;

    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 9;

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 6;

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

    // initial state vector
    x_ = VectorXd(n_x_);

    // initial covariance matrix
    P_ = MatrixXd(n_x_, n_x_);
    P_ <<
        1, 0,  0,  0,  0,
        0,  1, 0,  0,  0,
        0,  0,  1, 0,  0,
        0,  0,  0,  1, 0,
        0,  0,  0,  0, 1;

    // predicted sigma points matrix
    Xsig_pred_ = MatrixXd::Zero(n_x_, 2 * n_aug_ + 1);

    // Weights of sigma points
    weights_ = VectorXd::Zero(2 * n_aug_ + 1);
    weights_(0) = lambda_ / (lambda_ + n_aug_);
    for (int i= 1; i < 2 * n_aug_ + 1; i++)
    {
        weights_(i) = 1. / (2. * (lambda_ + n_aug_));
    }
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package)
{
    if (!is_initialized_)
    {
        if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
        {
            double rho = meas_package.raw_measurements_[0];
            double phi = meas_package.raw_measurements_[1];
            double rho_dot = meas_package.raw_measurements_[2];
            x_ << rho * std::cos(phi), rho * std::sin(phi), rho_dot, 0., 0.;
        }
        else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
        {
            double x = meas_package.raw_measurements_[0];
            double y = meas_package.raw_measurements_[1];
            x_ << x, y, 0., 0., 0.;
        }

        is_initialized_ = true;
    }
    else
    {
        Prediction((meas_package.timestamp_ - previous_timestamp_) / 1000000.0);

        if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_)
        {
            UpdateRadar(meas_package);
        }
        else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)
        {
            UpdateLidar(meas_package);
        }
    }

    previous_timestamp_ = meas_package.timestamp_;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t)
{
    /**
       TODO:

       Complete this function! Estimate the object's location. Modify the state
       vector, x_. Predict sigma points, the state, and the state covariance matrix.
    */
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package)
{
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
void UKF::UpdateRadar(MeasurementPackage meas_package)
{
    /**
       TODO:

       Complete this function! Use radar data to update the belief about the object's
       position. Modify the state vector, x_, and covariance, P_.

       You'll also need to calculate the radar NIS.
    */
}
