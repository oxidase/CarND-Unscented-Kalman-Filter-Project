#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

// https://stackoverflow.com/a/11126083
// Bring the 'difference' between two angles into [-pi; pi]
template <int K, typename T>
T normalize(T rad) {
  // Copy the sign of the value in radians to the value of pi.
  T signed_pi = std::copysign(M_PI, rad);
  // Set the value of difference to the appropriate signed value between pi and -pi.
  rad = std::fmod(rad + K * signed_pi,(2 * M_PI)) - K * signed_pi;
  return rad;
}

VectorXd ctrv_model(double delta_t, const VectorXd &x)
{
    double px = x(0);
    double py = x(1);
    double v = x(2);
    double psi = x(3);
    double psid = x(4);
    double nu_a = x(5);
    double nu_psidd = x(6);

    // Predicted state values
    double px_p, py_p, v_p, psi_p, psid_p;

    if (std::fabs(psid) > 1e-6)
    {
        px_p = px + v / psid * ( std::sin(psi + psid * delta_t) - std::sin(psi));
        py_p = py + v / psid * (-std::cos(psi + psid * delta_t) + std::cos(psi));
    }
    else
    {
        px_p = px + v * delta_t * std::cos(psi);
        py_p = py + v * delta_t * std::sin(psi);
    }

    v_p = v;
    psi_p = psi + psid * delta_t;
    psid_p = psid;

    // Add noise
    double delta_t2 = delta_t * delta_t / 2.;
    px_p   += delta_t2 * std::cos(psi) * nu_a;
    py_p   += delta_t2 * std::sin(psi) * nu_a;
    v_p    += delta_t * nu_a;
    psi_p  += delta_t2 * nu_psidd;
    psid_p += delta_t * nu_psidd;

    VectorXd x_p(5);
    x_p << px_p, py_p, v_p, psi_p, psid_p;
    return x_p;
}

VectorXd radar_model(const VectorXd &x)
{
    double px = x(0);
    double py = x(1);
    double v = x(2);
    double psi = x(3);

    double v1 = cos(psi)*v;
    double v2 = sin(psi)*v;
    double r = sqrt(px * px + py * py);

    // measurement model
    VectorXd z_p(3);
    z_p(0) = r;
    z_p(1) = atan2(py, px);                                        // phi
    z_p(2) = (std::fabs(r) > 1e-6) ? (px * v1 + py * v2 ) / r : v; // r_dot
    return z_p;
}

VectorXd lidar_model(const VectorXd &x)
{
    double px = x(0);
    double py = x(1);

    // measurement model
    VectorXd z_p(2);
    z_p(0) = px;
    z_p(1) = py;
    return z_p;
}

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF()
{
    // initially set to false, set to true in first call of ProcessMeasurement
    is_initialized = false;

    // State dimension
    n_x = 5;

    // Augmented state dimension
    n_aug = 7;

    // Sigma point spreading parameter
    lambda = 3 - n_aug;

    // if this is false, laser measurements will be ignored (except during init)
    use_laser = true;

    // if this is false, radar measurements will be ignored (except during init)
    use_radar = true;

    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a = 2;

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd = 1;

    // Laser measurement noise standard deviation position1 in m
    std_laspx = 0.15;

    // Laser measurement noise standard deviation position2 in m
    std_laspy = 0.15;

    // Radar measurement noise standard deviation radius in m
    std_radr = 0.3;

    // Radar measurement noise standard deviation angle in rad
    std_radphi = 0.03;

    // Radar measurement noise standard deviation radius change in m/s
    std_radrd = 0.3;

    // initial state vector
    x = VectorXd(n_x);

    // initial covariance matrix
    P = MatrixXd(n_x, n_x);
    P <<
        1, 0,  0,  0,  0,
        0,  1, 0,  0,  0,
        0,  0,  1, 0,  0,
        0,  0,  0,  1, 0,
        0,  0,  0,  0, 1;

    R_lidar = MatrixXd(2, 2);
    R_lidar <<
        std_laspx * std_laspx, 0,
        0, std_laspy * std_laspy;

    R_radar = MatrixXd(3, 3);
    R_radar <<
        std_radr * std_radr, 0, 0,
        0, std_radphi * std_radphi, 0,
        0, 0, std_radrd * std_radrd;

    // predicted sigma points matrix
    Xsig_pred = MatrixXd::Zero(n_x, 2 * n_aug + 1);

    // Weights of sigma points
    weights = VectorXd::Zero(2 * n_aug + 1);
    weights(0) = lambda / (lambda + n_aug);
    for (int i = 1; i < 2 * n_aug + 1; i++)
    {
        weights(i) = 1. / (2. * (lambda + n_aug));
    }
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package)
{
    double delta_t = (meas_package.timestamp_ - previous_timestamp) / 1000000.0;
    if (!is_initialized || std::abs(delta_t) > 100.)
    {
        if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
        {
            double rho = meas_package.raw_measurements_[0];
            double phi = meas_package.raw_measurements_[1];
            double rho_dot = meas_package.raw_measurements_[2];
            x << rho * std::cos(phi), rho * std::sin(phi), rho_dot, 0., 0.;
        }
        else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
        {
            x << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0., 0., 0.;
        }

        is_initialized = true;
    }
    else
    {
        Prediction(delta_t);

        if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar)
        {
            UpdateRadar(meas_package);
        }
        else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser)
        {
            UpdateLidar(meas_package);
        }
    }

    previous_timestamp = meas_package.timestamp_;
}

/**
 * Predicts sigma points Xsig_pred, the state x, and the state covariance matrix P.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t)
{
    // Create augmented mean state x_{a,k|k}
    VectorXd x_aug = VectorXd(n_aug);
    x_aug.head(5) = x;
    x_aug(5) = 0.;
    x_aug(6) = 0.;

    // Create augmented covariance matrix P_{a,k|k}
    MatrixXd P_aug = MatrixXd(n_aug, n_aug);
    P_aug.fill(0.0);
    P_aug.topLeftCorner(5,5) = P;
    P_aug(5,5) = std_a * std_a;
    P_aug(6,6) = std_yawdd * std_yawdd;

    // Create square root matrix \sqrt{P_{a,k|k}}
    MatrixXd L = P_aug.llt().matrixL();

    // Create augmented sigma points X_{a,k|k}
    MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);
    Xsig_aug.col(0) = x_aug;   // x_{a,k|k}
    for (int i = 0; i < n_aug; ++i)
    {
        Xsig_aug.col(i+1)       = x_aug + sqrt(lambda + n_aug) * L.col(i); // x_{a,k|k} + \sqrt{(\lambda + n_a) P_{a,k|k}}
        Xsig_aug.col(i+1+n_aug) = x_aug - sqrt(lambda + n_aug) * L.col(i); // x_{a,k|k} - \sqrt{(\lambda + n_a) P_{a,k|k}}
    }

    //std::cout << Xsig_aug << "\n\n";

    // Compute process model prediction X_{k+1|k}
    for (int i = 0; i < Xsig_aug.cols(); ++i)
    {
        Xsig_pred.col(i) = ctrv_model(delta_t, Xsig_aug.col(i));
    }
    //std::cout << Xsig_pred << "\n\n";

    // Predicted mean x_{k+1|k}
    x.fill(0);
    for (int i = 0; i < Xsig_pred.cols(); ++i)
    {
        x += weights(i) * Xsig_pred.col(i);
    }

    //std::cout << x << "\n\n";

    // Predicted covariance P_{k+1|k}
    P.fill(0);
    for (int i = 0; i < Xsig_pred.cols(); ++i)
    {
        VectorXd x_diff = Xsig_pred.col(i) - x;
        x_diff(3) = normalize<1>(x_diff(3)); // angle normalization
        P += weights(i) * x_diff * x_diff.transpose();
    }
    //std::cout << P << "\n\n";
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

    // Measurement Prediction \mathcal{Z}_{k+1|k,i}
    MatrixXd Z_k1k = MatrixXd(2, Xsig_pred.cols());
    for (int i = 0; i < Xsig_pred.cols(); i++)
    {
        Z_k1k.col(i) = lidar_model(Xsig_pred.col(i));
    }

    // Predicted Measurement Mean z_{k+1|k}
    VectorXd z_k1k = VectorXd::Zero(2);
    for (int i = 0; i < Xsig_pred.cols(); ++i)
    {
        z_k1k += weights(i) * Z_k1k.col(i);
    }

    // Predicted Covariance S_{k+1|k}
    MatrixXd S_k1k = MatrixXd::Zero(2, 2);
    for (int i = 0; i < Xsig_pred.cols(); ++i)
    {
        VectorXd z_diff = Z_k1k.col(i) - z_k1k;
        S_k1k += weights(i) * z_diff * z_diff.transpose();
    }

    S_k1k += R_lidar;

    double nis = Update(Z_k1k, z_k1k, S_k1k, meas_package.raw_measurements_ - z_k1k);

    std::cout << "L NIS = " << nis << "\n";
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

    // Measurement Prediction \mathcal{Z}_{k+1|k,i}
    MatrixXd Z_k1k = MatrixXd(3, Xsig_pred.cols());
    for (int i = 0; i < Xsig_pred.cols(); i++)
    {
        Z_k1k.col(i) = radar_model(Xsig_pred.col(i));
    }

    // Predicted Measurement Mean z_{k+1|k}
    VectorXd z_k1k = VectorXd::Zero(3);
    for (int i = 0; i < Xsig_pred.cols(); ++i)
    {
        z_k1k += weights(i) * Z_k1k.col(i);
    }

    // Predicted Covariance S_{k+1|k}
    MatrixXd S_k1k = MatrixXd::Zero(3, 3);
    for (int i = 0; i < Xsig_pred.cols(); ++i)
    {
        VectorXd z_diff = Z_k1k.col(i) - z_k1k;
        z_diff(1) = normalize<1>(z_diff(1)); // angle normalization
        S_k1k += weights(i) * z_diff * z_diff.transpose();
    }

    S_k1k += R_radar;

    VectorXd z_diff = meas_package.raw_measurements_ - z_k1k;
    z_diff(1) = normalize<1>(z_diff(1));
    double nis = Update(Z_k1k, z_k1k, S_k1k, z_diff);

    std::cout << "R NIS = " << nis << "\n";
}


double UKF::Update(const MatrixXd& Z, const VectorXd& z, const MatrixXd& S, const VectorXd& z_res)
{
    // Cross-correlation between sigma points in state space and measurement space T_{k+1|k}
    MatrixXd T = MatrixXd::Zero(n_x, z.size());
    for (int i = 0; i < Xsig_pred.cols(); i++)
    {
        T += weights(i) * (Xsig_pred.col(i) - x) * (Z.col(i) - z).transpose();
    }

    //std::cout << "T = " << T << "\n\n";

    // Kalman gain K_{k+1|k}
    MatrixXd K = T * S.inverse();

    // State update x_{k+1|k+1}
    x = x + K * z_res;

    // Covariance matrix update P_{k+1|k+1}
    P = P - K * S * K.transpose();

    // Compute Normalized Innovation Squared \varepsilon
    return z_res.transpose() * S.inverse() * z_res;
}
