#include <uWS/uWS.h>
#include <iostream>
#include "json.hpp"
#include <math.h>
#include "ukf.h"
#include "tools.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
std::string hasData(std::string s) {
    auto found_null = s.find("null");
    auto b1 = s.find_first_of("[");
    auto b2 = s.find_first_of("]");
    if (found_null != std::string::npos) {
        return "";
    }
    else if (b1 != std::string::npos && b2 != std::string::npos) {
        return s.substr(b1, b2 - b1 + 1);
    }
    return "";
}

void process_file(char *input_file)
{
    std::ifstream in(input_file);
    std::vector<MeasurementPackage> measurement_pack_list;
    std::vector<Eigen::VectorXd> true_values;
    std::string line;

    // prep the measurement packages (each line represents a measurement at a
    // timestamp)
    while (getline(in, line))
    {
        std::string sensor_type;
        std::istringstream iss(line);
        MeasurementPackage meas_package;
        long long timestamp;

        // reads first element from the current line
        iss >> sensor_type;

        if (sensor_type.compare("L") == 0) {
            // laser measurement

            // read measurements at this timestamp
            meas_package.sensor_type_ = MeasurementPackage::LASER;
            meas_package.raw_measurements_ = VectorXd(2);
            float px;
            float py;
            iss >> px;
            iss >> py;
            meas_package.raw_measurements_ << px, py;
            iss >> timestamp;
            meas_package.timestamp_ = timestamp;
            measurement_pack_list.push_back(meas_package);
        } else if (sensor_type.compare("R") == 0) {
            // radar measurement
            meas_package.sensor_type_ = MeasurementPackage::RADAR;
            meas_package.raw_measurements_ = VectorXd(3);
            float rho;
            float phi;
            float rho_dot;
            iss >> rho;
            iss >> phi;
            iss >> rho_dot;
            meas_package.raw_measurements_ << rho, phi, rho_dot;
            iss >> timestamp;
            meas_package.timestamp_ = timestamp;
            measurement_pack_list.push_back(meas_package);
        }

        // read ground truth data to compare later
        Eigen::VectorXd gt_package(4);
        float x_gt, y_gt, vx_gt, vy_gt, a, b;
        iss >> x_gt >> y_gt >> vx_gt  >> vy_gt >> a >> b;
        gt_package << x_gt, y_gt, vx_gt, vy_gt;
        true_values.push_back(gt_package);

        static int xxx = 0;
        if (++xxx == 2) break;
    }

    // Create a UKF instance
    UKF ukf;

    // used to compute the RMSE later
    vector<VectorXd> estimations;
    vector<VectorXd> ground_truth;

    size_t number_of_measurements = measurement_pack_list.size();

    // column names for output file
    std::cout << "px" << "\t";
    std::cout << "py" << "\t";
    std::cout << "v" << "\t";
    std::cout << "yaw_angle" << "\t";
    std::cout << "yaw_rate" << "\t";
    std::cout << "px_measured" << "\t";
    std::cout << "py_measured" << "\t";
    std::cout << "px_true" << "\t";
    std::cout << "py_true" << "\t";
    std::cout << "vx_true" << "\t";
    std::cout << "vy_true" << "\t";
    std::cout << "NIS" << "\n";


    for (size_t k = 0; k < number_of_measurements; ++k)
    {
        // Call the UKF-based fusion
        ukf.ProcessMeasurement(measurement_pack_list[k]);

        // output the estimation
        std::cout << ukf.x_(0) << "\t"; // pos1 - est
        std::cout << ukf.x_(1) << "\t"; // pos2 - est
        std::cout << ukf.x_(2) << "\t"; // vel_abs -est
        std::cout << ukf.x_(3) << "\t"; // yaw_angle -est
        std::cout << ukf.x_(4) << "\t"; // yaw_rate -est

        // output the measurements
        if (measurement_pack_list[k].sensor_type_ == MeasurementPackage::LASER)
        {
            // output the estimation

            // p1 - meas
            std::cout << measurement_pack_list[k].raw_measurements_(0) << "\t";

            // p2 - meas
            std::cout << measurement_pack_list[k].raw_measurements_(1) << "\t";
        }
        else if (measurement_pack_list[k].sensor_type_ == MeasurementPackage::RADAR)
        {
            // output the estimation in the cartesian coordinates
            float rho = measurement_pack_list[k].raw_measurements_(0);
            float phi = measurement_pack_list[k].raw_measurements_(1);
            std::cout << rho * cos(phi) << "\t"; // p1_meas
            std::cout << rho * sin(phi) << "\t"; // p2_meas
        }

        // output the ground truth packages
        std::cout << true_values[k](0) << "\t";
        std::cout << true_values[k](1) << "\t";
        std::cout << true_values[k](2) << "\t";
        std::cout << true_values[k](3) << "\t";

        // output the NIS values

        if (measurement_pack_list[k].sensor_type_ == MeasurementPackage::LASER)
        {
            //std::cout << ukf.NIS_laser_ << "\n";
        }
        else if (measurement_pack_list[k].sensor_type_ == MeasurementPackage::RADAR)
        {
            //std::cout << ukf.NIS_radar_ << "\n";
        }
        std::cout << "\n";


        // convert ukf x vector to cartesian to compare to ground truth
        VectorXd ukf_x_cartesian_ = VectorXd(4);

        float x_estimate_ = ukf.x_(0);
        float y_estimate_ = ukf.x_(1);
        float vx_estimate_ = ukf.x_(2) * cos(ukf.x_(3));
        float vy_estimate_ = ukf.x_(2) * sin(ukf.x_(3));

        ukf_x_cartesian_ << x_estimate_, y_estimate_, vx_estimate_, vy_estimate_;

        estimations.push_back(ukf_x_cartesian_);
        ground_truth.push_back(true_values[k]);
    }

    // compute the accuracy (RMSE)
    Tools tools;
    std::cout << "Accuracy - RMSE:" << endl << tools.CalculateRMSE(estimations, ground_truth) << endl;
}

int main(int argc, char* argv[])
{
    if (argc > 1)
    {
        process_file(argv[1]);
        return 0;
    }

    uWS::Hub h;

    // Create a Kalman Filter instance
    UKF ukf;

    // used to compute the RMSE later
    Tools tools;
    vector<VectorXd> estimations;
    vector<VectorXd> ground_truth;

    h.onMessage([&ukf,&tools,&estimations,&ground_truth](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length, uWS::OpCode opCode) {
                    // "42" at the start of the message means there's a websocket message event.
                    // The 4 signifies a websocket message
                    // The 2 signifies a websocket event

                    if (length && length > 2 && data[0] == '4' && data[1] == '2')
                    {

                        auto s = hasData(std::string(data));
                        if (s != "") {

                            auto j = json::parse(s);

                            std::string event = j[0].get<std::string>();

                            if (event == "telemetry") {
                                // j[1] is the data JSON object

                                string sensor_measurment = j[1]["sensor_measurement"];

                                MeasurementPackage meas_package;
                                istringstream iss(sensor_measurment);
                                long long timestamp;

                                // reads first element from the current line
                                string sensor_type;
                                iss >> sensor_type;

                                if (sensor_type.compare("L") == 0) {
                                    meas_package.sensor_type_ = MeasurementPackage::LASER;
                                    meas_package.raw_measurements_ = VectorXd(2);
                                    float px;
                                    float py;
                                    iss >> px;
                                    iss >> py;
                                    meas_package.raw_measurements_ << px, py;
                                    iss >> timestamp;
                                    meas_package.timestamp_ = timestamp;
                                } else if (sensor_type.compare("R") == 0) {

                                    meas_package.sensor_type_ = MeasurementPackage::RADAR;
                                    meas_package.raw_measurements_ = VectorXd(3);
                                    float ro;
                                    float theta;
                                    float ro_dot;
                                    iss >> ro;
                                    iss >> theta;
                                    iss >> ro_dot;
                                    meas_package.raw_measurements_ << ro,theta, ro_dot;
                                    iss >> timestamp;
                                    meas_package.timestamp_ = timestamp;
                                }
                                float x_gt;
                                float y_gt;
                                float vx_gt;
                                float vy_gt;
                                iss >> x_gt;
                                iss >> y_gt;
                                iss >> vx_gt;
                                iss >> vy_gt;
                                VectorXd gt_values(4);
                                gt_values(0) = x_gt;
                                gt_values(1) = y_gt;
                                gt_values(2) = vx_gt;
                                gt_values(3) = vy_gt;
                                ground_truth.push_back(gt_values);

                                //Call ProcessMeasurment(meas_package) for Kalman filter
                                ukf.ProcessMeasurement(meas_package);

                                //Push the current estimated x,y positon from the Kalman filter's state vector

                                VectorXd estimate(4);

                                double p_x = ukf.x_(0);
                                double p_y = ukf.x_(1);
                                double v  = ukf.x_(2);
                                double yaw = ukf.x_(3);

                                double v1 = cos(yaw)*v;
                                double v2 = sin(yaw)*v;

                                estimate(0) = p_x;
                                estimate(1) = p_y;
                                estimate(2) = v1;
                                estimate(3) = v2;

                                estimations.push_back(estimate);

                                VectorXd RMSE = tools.CalculateRMSE(estimations, ground_truth);

                                json msgJson;
                                msgJson["estimate_x"] = p_x;
                                msgJson["estimate_y"] = p_y;
                                msgJson["rmse_x"] =  RMSE(0);
                                msgJson["rmse_y"] =  RMSE(1);
                                msgJson["rmse_vx"] = RMSE(2);
                                msgJson["rmse_vy"] = RMSE(3);
                                auto msg = "42[\"estimate_marker\"," + msgJson.dump() + "]";
                                // std::cout << msg << std::endl;
                                ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);

                            }
                        } else {

                            std::string msg = "42[\"manual\",{}]";
                            ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
                        }
                    }

                });

    // We don't need this since we're not using HTTP but if it's removed the program
    // doesn't compile :-(
    h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data, size_t, size_t) {
                        const std::string s = "<h1>Hello world!</h1>";
                        if (req.getUrl().valueLength == 1)
                        {
                            res->end(s.data(), s.length());
                        }
                        else
                        {
                            // i guess this should be done more gracefully?
                            res->end(nullptr, 0);
                        }
                    });

    h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
                       std::cout << "Connected!!!" << std::endl;
                   });

    h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code, char *message, size_t length) {
                          ws.close();
                          std::cout << "Disconnected" << std::endl;
                      });

    int port = 4567;
    if (h.listen(port))
    {
        std::cout << "Listening to port " << port << std::endl;
    }
    else
    {
        std::cerr << "Failed to listen to port" << std::endl;
        return -1;
    }
    h.run();
}
