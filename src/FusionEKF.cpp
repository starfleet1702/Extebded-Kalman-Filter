#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

#define NOISE_AX 9;
#define NOISE_AY 9;
/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
	is_initialized_ = false;

	previous_timestamp_ = 0;

	// initializing matrices
	R_laser_ = MatrixXd(2, 2);
	R_radar_ = MatrixXd(3, 3);
	H_laser_ = MatrixXd(2, 4);
	H_laser_ << 1, 0, 0, 0,
				0, 1, 0, 0;
	Hj_ = MatrixXd(3, 4);
	
	Hj_ << 1,1,0,0,
		   1,1,0,0,
		   1,1,1,1;
	
	/**
	 * Finish initializing the FusionEKF.
	 * Set the process and measurement noises
	 */

	
	//measurement covariance matrix - laser
	R_laser_ << 0.0225, 0, 0, 0.0225;

	//measurement covariance matrix - radar
	R_radar_ << 0.09, 0, 0, 0, 0.0009, 0, 0, 0, 0.09;
	
	ekf_.P_ = MatrixXd(4, 4);
	ekf_.P_ << 1,0,0,0,
				0,1,0,0,
				0,0,1000,0,
				0,0,0,1000;
				
	ekf_.F_ = MatrixXd(4, 4);
	ekf_.F_ << 1,0,1,0,
				0,1,0,1,
				0,0,1,0,
				0,0,0,1;
	
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {
}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

	/*****************************************************************************
	 *  Initialization
	 ****************************************************************************/
	if (!is_initialized_) {
		/**
		 * Initialize the state ekf_.x_ with the first measurement.
		 * Create the covariance matrix.
		 * Remember: you'll need to convert radar from polar to cartesian coordinates.
		 */
		// first measurement
		cout << "FusionEKF Initialization: " << endl;
		ekf_.x_ = VectorXd(4);
		


		if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
			/**
			 Convert radar from polar to cartesian coordinates and initialize state.
			 */
			float r = measurement_pack.raw_measurements_[0];
			float phi = measurement_pack.raw_measurements_[1];
			double drdt = measurement_pack.raw_measurements_[2];
			ekf_.x_ << r * cos(phi), r * sin(phi), 0, 0;
		} else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
			/**
			 Initialize state.
			 */
			float px = measurement_pack.raw_measurements_[0];
			float py = measurement_pack.raw_measurements_[1];
			ekf_.x_ << px, py, 0, 0;
		}
		previous_timestamp_ = measurement_pack.timestamp_;
		// done initializing, no need to predict or update
		is_initialized_ = true;
		return;
	}

	/*****************************************************************************
	 *  Prediction
	 ****************************************************************************/

	/**
	 * Update the state transition matrix F according to the new elapsed time.
	 - Time is measured in seconds.
	 * Update the process noise covariance matrix.
	 * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
	 */
	//Calculating the time difference
	float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;

	//Updating state transition matrix F
	ekf_.F_(0, 2) = dt;
	ekf_.F_(1, 3) = dt;

	double dt2 = dt * dt;
	double dt3 = dt2 * dt;
	double dt4 = dt3 * dt;
	float noise_ax = 9.0;
	float noise_ay = 9.0;

	//updating noise covariance matrix
	ekf_.Q_ = MatrixXd(4, 4);
	ekf_.Q_ << (dt4 * noise_ax / 4.0), 0, (dt3 * noise_ax / 2.0), 0, 0, (dt4
			* noise_ay / 4.0), 0, (dt3 * noise_ay / 2.0), (dt3 * noise_ax / 2.0), 0, dt2
			* noise_ax, 0, 0, (dt3 * noise_ay / 2.0), 0, (dt2 * noise_ay);
	std::cout<<"Predicting..."<<endl;
	ekf_.Predict();

	/*****************************************************************************
	 *  Update
	 ****************************************************************************/

	/**
	 * Use the sensor type to perform the update step.
	 * Update the state and covariance matrices.
	 */

	if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
		//RADAR Update
		ekf_.R_ = R_radar_;
		Hj_ = tools.CalculateJacobian(ekf_.x_);
		ekf_.H_ = Hj_;
		ekf_.UpdateEKF(measurement_pack.raw_measurements_);
	} else {
		// Laser updates
		ekf_.R_ = R_laser_;
		ekf_.H_ = H_laser_;
		ekf_.Update(measurement_pack.raw_measurements_);
		
	}

	previous_timestamp_ = measurement_pack.timestamp_;

	// print the output
	cout << "x_ = " << ekf_.x_ << endl;
	cout << "P_ = " << ekf_.P_ << endl;
}
