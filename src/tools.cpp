#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using namespace std;

Tools::Tools() {
}

Tools::~Tools() {
}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
		const vector<VectorXd> &ground_truth) {
	/**
	 * Calculate the RMSE here.
	 */
	//Sum(sqrt(1/n(a-b)^2))
	std::cout <<"Calculating RMSE" <<endl;
	VectorXd rmse(4);
	rmse << 0,0,0,0;

	long estimations_length = estimations.size();

	if((estimations_length==0) || (estimations_length!=ground_truth.size())){
		std::cout<<"length of the estimation is either 0 or does not match with ground_truth!!!"<<endl;
		return rmse;
	}

	for(unsigned long i=0;i<estimations_length;i++){
		VectorXd diff = (estimations[i]-ground_truth[i]);
		diff = diff.array()*diff.array();
		rmse+= diff;
	}
	//calculating mean
	rmse = rmse/estimations_length;
	//calculating sqare root
	rmse = rmse.array().sqrt();
	std::cout<<"RMSE END"<<endl;
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
	/**
	 * Calculate a Jacobian here.
	 */
	std::cout<<"Calculating Jacobian"<<endl;
	MatrixXd Hj = MatrixXd::Zero(3, 4);
	float min_val_threshold = 0.0001;
	double px = x_state(0);
	double py = x_state(1);
	double vx = x_state(2);
	double vy = x_state(3);

	double c1 = px*px+py*py;
	double c2 = sqrt(c1);
	double c3 = (c1*c2);

	if(fabs(c1)>=min_val_threshold){
		Hj << (px/c2),(py/c2),0.0,0.0,
				-(py/c1), (px/c1), 0.0, 0.0,
				py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;
	}
	std::cout<<"Jacobian End"<<endl;
	return Hj;
}
