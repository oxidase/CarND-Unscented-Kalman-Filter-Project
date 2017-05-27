#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    VectorXd rmse(4);
    rmse << 0,0,0,0;

    //accumulate squared residuals
    for (int i=0; i < estimations.size(); ++i)
    {
        VectorXd res = estimations[i]-ground_truth[i];
        res = res.array()*res.array();
        rmse += res;
    }

    //calculate the mean
    rmse /= estimations.size();

    //calculate the squared root
    return rmse.array().sqrt();
}
