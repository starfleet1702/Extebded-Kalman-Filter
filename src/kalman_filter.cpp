#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace std;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {
	MatrixXd I_= MatrixXd::Identity(2,2);
}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;

}

void KalmanFilter::Predict() {
  /**
    * predict the state
  */
	std::cout<<"Predict Start"<<endl;
	x_=F_*x_;
	std::cout<<"1 Done"<<endl;
	P_=F_*P_*F_.transpose()+Q_;
	std::cout<<"Predict End"<<endl;

}

void KalmanFilter::Update(const VectorXd &z) {
  /**
    * update the state by using Kalman Filter equations
  */
	std::cout<<"Update Start"<<endl;
	VectorXd y = z - H_*x_;
	cout<<"1"<<endl;
	MatrixXd S = H_*P_*H_.transpose()+R_;
	cout<<"2"<<endl;
	MatrixXd K = P_*H_.transpose()*S.inverse();
	cout<<"3"<<endl;
	x_ = x_+K*y;
	cout<<"4"<<endl;
	int x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size,x_size);
	P_ = (I-K*H_)*P_;
	std::cout<<"Update End"<<endl;

}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
    * update the state by using Extended Kalman Filter equations
  */
	std::cout<<"Update EKF Start"<<endl;
	VectorXd x_polar(3);
	double px = x_[0];
	double py = x_[1];
	double vx = x_[2];
	double vy = x_[3];
	
	double r = sqrt(px*px+py*py);
	double phi = r==0?0.0:atan2(py,px);
	double drdt = r>0.0001 ?(px*vx+py*vy)/r:0.0;
	x_polar<< r,phi,drdt;
	
	cout<<"x_polar calculated"<<endl;
	
  	VectorXd y = z - x_polar; //
	//Normalizing phi in y
	y(1) = atan2(sin(y(1)), cos(y(1)));
	MatrixXd S = H_*P_*H_.transpose()+R_;
	MatrixXd K = P_*H_.transpose()*S.inverse();
	cout<<y.size()<<endl;
	
	x_ = x_+K*y;
	int x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size,x_size);
	P_ = (I-K*H_)*P_;
	std::cout<<"Update EKF End"<<endl;
}
