/****************************************************************************
*
* ViSP, open source Visual Servoing Platform software.
* Copyright (C) 2005 - 2019 by Inria. All rights reserved.
*
* This software is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or
* (at your option) any later version.
* See the file LICENSE.txt at the root directory of this source
* distribution for additional information about the GNU GPL.
*
* For using ViSP with software that can not be combined with the GNU
* GPL, please contact Inria about acquiring a ViSP Professional
* Edition License.
*
* See http://visp.inria.fr for more information.
*
* This software was developed at:
* Inria Rennes - Bretagne Atlantique
* Campus Universitaire de Beaulieu
* 35042 Rennes Cedex
* France
*
* If you have questions regarding the use of this file, please contact
* Inria at visp@inria.fr
*
* This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
* WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
*
* Description:
* Interface for a  generic ADEPT Viper (either 650 or 850) robot.
*
* Authors:
* Fabien Spindler
*
*****************************************************************************/

#ifndef NewVpViper_h
#define NewVpViper_h

#include <visp3/core/vpCameraParameters.h>
#include <visp3/core/vpHomogeneousMatrix.h>
#include <visp3/core/vpImage.h>
#include <visp3/core/vpRGBa.h>
#include <visp3/core/vpVelocityTwistMatrix.h>
#include <visp3/robot/vpRobotException.h>

class NewVpViper
{
public:
	NewVpViper();
	virtual ~NewVpViper() {};

	/** @name Inherited functionalities from NewVpViper */
	//@{
	vpHomogeneousMatrix getForwardKinematics(const vpColVector &q);
	//unsigned int getInverseKinematicsWrist(const vpHomogeneousMatrix &fMw, vpColVector &q,
		const bool &verbose = false;
	//unsigned int getInverseKinematics(const vpHomogeneousMatrix &fMc, vpColVector &q, const bool &verbose = false);
	vpHomogeneousMatrix get_fMc(const vpColVector &q);
	void get_fMw(const vpColVector &q, vpHomogeneousMatrix &fMw);
	void get_wMe(vpHomogeneousMatrix &wMe);
	void get_eMc(vpHomogeneousMatrix &eMc);
	void get_eMs(vpHomogeneousMatrix &eMs);
	void get_fMe(const vpColVector &q, vpHomogeneousMatrix &fMe);
	void get_fMc(const vpColVector &q, vpHomogeneousMatrix &fMc);

	void get_cMe(vpHomogeneousMatrix &cMe);
	void get_cVe(vpVelocityTwistMatrix &cVe);
	void get_fJw(const vpColVector &q, vpMatrix &fJw);
	void get_fJe(const vpColVector &q, vpMatrix &fJe);
	void get_eJe(const vpColVector &q, vpMatrix &eJe);
	void get_q(const vpColVector &CamVel, vpColVector &q);
	vpColVector getJointMin();
	vpColVector getJointMax();
	double getCoupl56();
	//@}
	void setSamplingTime(double SamplingTime);
	double getSamplingTime();
		friend VISP_EXPORT std::ostream &operator<<(std::ostream &os, const NewVpViper &viper);

private:
	bool convertJointPositionInLimits(unsigned int joint, const double &q, double &q_mod,
		const bool &verbose = false);

public:
	static const unsigned int njoint; ///< Number of joint.

protected:
	vpHomogeneousMatrix eMc; //!< End effector to camera transformation
							 // Minimal representation of eMc
	vpTranslationVector etc; // meters
	vpRxyzVector erc;        // radian

	// Denavit-Hartenberg parameters
	double a1, d1; 
	double a2;     
	double a3;    
	double d4;     
	double d6;     
	double d7;     
	double c56;    
	double _SamplingTime;
	// Software joint limits in radians
	vpColVector joint_max; // Maximal value of the joints
	vpColVector joint_min; // Minimal value of the joints
};

#endif