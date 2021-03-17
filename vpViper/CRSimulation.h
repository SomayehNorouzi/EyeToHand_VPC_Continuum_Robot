#ifndef CRSimulation_h
#define CRSimulation_h

#include <visp3/core/vpCameraParameters.h>
#include <visp3/core/vpHomogeneousMatrix.h>
#include <visp3/core/vpImage.h>
#include <visp3/core/vpRGBa.h>
#include <visp3/core/vpVelocityTwistMatrix.h>
#include <visp3/robot/vpRobotException.h>

class CRSimulation
{
public:
	CRSimulation();
	virtual ~CRSimulation() {};


	vpHomogeneousMatrix get_fMe(const vpColVector& q);

	void setfMc(vpHomogeneousMatrix& fMc);
	vpHomogeneousMatrix get_fMc();
	void get_fMw(const vpColVector& q, vpHomogeneousMatrix& fMw);
	void get_wMe(vpHomogeneousMatrix& wMe);
	//void get_eMc(vpHomogeneousMatrix &eMc);
	void get_eMs(vpHomogeneousMatrix& eMs);

	void get_fJw(const vpColVector& q, vpMatrix& fJw);
	void get_fJe(const vpColVector& q, vpMatrix& fJe);
	void get_eJe(const vpColVector& q, vpMatrix& eJe);
	void get_q(const vpColVector& CamVel, vpColVector& q);
	vpColVector get_qdot(const vpColVector& CamVel, const vpColVector& q);
	
	vpColVector getJointMin();
	vpColVector getJointMax();

	void setSamplingTime(double SamplingTime);
	double getSamplingTime();

	friend VISP_EXPORT std::ostream& operator<<(std::ostream& os, const NewVpViper& viper);

private:
	bool convertJointPositionInLimits(unsigned int joint, const double& q, double& q_mod,
		const bool& verbose = false);

public:
	static const unsigned int njoint; ///< Number of joint.

protected:

	// Denavit-Hartenberg parameters
	double L;
	double d7;
	double DK;
	double RC;
	double K;

	double _SamplingTime;
	vpHomogeneousMatrix _fMc;

	// Software joint limits in radians
	vpColVector joint_max; // Maximal value of the joints
	vpColVector joint_min; // Minimal value of the joints
};

#endif#pragma once


