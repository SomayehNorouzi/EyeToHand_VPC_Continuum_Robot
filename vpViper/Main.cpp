#include "stdafx.h"
#include <visp3/gui/vpDisplayGDI.h>
#include <visp3/gui/vpDisplayOpenCV.h>
#include <visp3/gui/vpDisplayX.h>
#include <visp3/vs/vpServoDisplay.h>
#include <visp3/vs/vpServo.h>
#include <visp3/gui/vpProjectionDisplay.h>			
#include <visp3/visual_features/vpFeatureBuilder.h> 
#include "ControlAlgDisplay.h"
#include <visp/vpViper.h>                          
#include <stdio.h>                           
#include "CRSimulation.h"
#include <nlopt.hpp>
#include "ControlAlg.h"
#include <cmath>    


void display_trajectory(const vpImage<unsigned char> &I, std::vector<vpPoint> &point, const vpHomogeneousMatrix &cMo,
	const vpCameraParameters &cam)
{
	static std::vector<vpImagePoint> traj[4];// 2D point
	vpImagePoint cog;

	for (unsigned int i = 0; i < 4; i++) {
		point[i].project(cMo);
		vpMeterPixelConversion::convertPoint(cam, point[i].get_x(), point[i].get_y(), cog);
		traj[i].push_back(cog);
	}
	for (unsigned int i = 0; i < 4; i++) {
		for (unsigned int j = 1; j < traj[i].size(); j++) {
			vpDisplay::displayLine(I, traj[i][j - 1], traj[i][j], vpColor::green);
		}
	}
}

double set_min_obj_wrapper(const std::vector<double>& x, std::vector<double>& grad, void* data) {

	ControlAlg *obj = static_cast<ControlAlg*>(data);
	return obj->costcont(x, grad, NULL);
}

double set_vision_Constraint_wrapper(const std::vector<double>& x, std::vector<double>& grad, void* data) {

	ControlAlg* obj = static_cast<ControlAlg*>(data);
	return obj->vision_const(x, grad, NULL);
}

int main()
{
	try {
		CRSimulation robot;
		ControlAlg task;
	    nlopt_opt opt;

		nlopt::opt myOptContLaw(nlopt::algorithm::GN_ISRES, 3);
		void* data_wrapper = &task;
		myOptContLaw.set_min_objective(set_min_obj_wrapper, data_wrapper); 
		std::vector <double> lb = { vpMath::rad(-10),vpMath::rad(-10), -0.05 };
		std::vector <double> ub = { vpMath::rad(10),vpMath::rad(10), 0.05 };
		myOptContLaw.set_lower_bounds(lb);
		myOptContLaw.set_upper_bounds(ub);

		vpCameraParameters cam;  // Create a camera parameter container

		double px = 700; double py = 700; double u0 = 320; double v0 = 240; double zestimated = 0.24;
		cam.initPersProjWithoutDistortion(px, py, u0, v0); // Camera initialization with a perspective projection without distortion model
		
		task.set_cam_intrinsic(px, py, u0, v0);
		task.set_vision_boundaries(600, 50, 50, 400);//set_vision_boundaries(double RightBound, double LeftBound, double UpperBound, double DownBound)
		
		void* vision_data_wrapper = &task;
		myOptContLaw.add_inequality_constraint(set_vision_Constraint_wrapper, vision_data_wrapper);
	
		myOptContLaw.set_xtol_rel(1e-4);
		
		vpHomogeneousMatrix eMc;
		vpTranslationVector _etc;
		vpRxyzVector _erc;

		_etc[0] = 0;
		_etc[1] = 0.22;
		_etc[3] = 0;//0
		_erc[0] = 0.5 * M_PI;
		_erc[1] = 0;
		_erc[2] = -0.5 * M_PI;//
		vpRotationMatrix eRc(_erc);
		eMc.buildFrom(_etc, eRc); //this matrix is subject to change in ETH config. as the robot moves

		vpColVector q(3);
		vpColVector qd(3);

		q[0] = vpMath::rad(0);
		q[1] = vpMath::rad(0);
		q[2] = 0.00;
	
		qd[0] = vpMath::rad(10);
		qd[1] = vpMath::rad(70);
		qd[2] = 0.00;

		vpHomogeneousMatrix cMo = eMc.inverse();
		vpHomogeneousMatrix fMe = robot.get_fMe(q);
		vpHomogeneousMatrix fMc = fMe * eMc;

		robot.setfMc(fMc);

		vpHomogeneousMatrix fMed = robot.get_fMe(qd);
		vpHomogeneousMatrix edMc = fMed.inverse() * fMc;
		vpHomogeneousMatrix cMod = edMc.inverse();

		vpVelocityTwistMatrix cVe;   
		vpHomogeneousMatrix cMe;     
		vpMatrix eJe;                
		eJe.resize(6,3);

		std::vector<vpPoint> point(4);
		point[0] = vpPoint(-0.0035, 0, -0.0035);// meter
		point[1] = vpPoint(0.0035, 0, -0.0035);
		point[2] = vpPoint(0.0035, 0, 0.0035);
		point[3] = vpPoint(-0.0035, 0, 0.0035);

		task.set_cVe(cMo);     
		robot.get_eJe(q, eJe);  
		task.set_eJe(eJe);     

		vpFeaturePoint p[4], pd[4], Sm[4];
		vpColVector sEstim(8);
		vpColVector epsil(8);
		vpColVector sDesi(8);

		for (unsigned int i = 0; i < 4; i++) {
			point[i].track(cMod); 
			vpFeatureBuilder::create(pd[i], point[i]);
			point[i].track(cMo);
			vpFeatureBuilder::create(p[i], point[i]);
			task.addFeature(p[i], pd[i]);

			sEstim[2 * i] = p[i].get_x();
			sEstim[2 * i +1] = p[i].get_y();

		}
			
		vpImage<unsigned char> Iint(480, 640, 0);
		vpDisplayGDI displayInt(Iint, 0, 0, "Internal view");
		
		double Landa = 1;
		int Np = 3;
		double Ts = 0.4;
		double minf;
		std::vector<double>  jointvel{ 0.001, 0.001, 0.001 };
		vpColVector qdot(3);
		qdot[0] = 0.001;
		qdot[1] = 0.001; 
		qdot[2] = 0.001; 
		double NormError = 50;
		double error_value = 10;

		task.Set_landa(Landa);
		task.SetNp(Np);
		robot.setSamplingTime(Ts);
		task.set_delta_t(Ts);
		task.set_fMc(fMc);
		
	while (1) 
	{
		fMe = robot.get_fMe(q);
		eMc = fMe.inverse() * fMc;
		cMe = eMc.inverse();
		cMo = eMc.inverse();
		vpMatrix L_s(2, 6);
		vpMatrix L_total;
		for (unsigned int i = 0; i < 4; i++)
		{
			point[i].track(cMo);
			vpFeatureBuilder::create(p[i], point[i]);
			L_s = p[i].interaction(vpBasicFeature::FEATURE_ALL);
			L_total.stack(L_s);
		}
		task.SetInteractionMatrix(L_total);
		task.Setq(q);
		robot.get_eJe(q, eJe);
		task.set_eJe(eJe);
		
		vpMatrix Maxi;
		Maxi.resize(8,3);
		vpColVector Ma;
		Ma.resize(8);

		vpVelocityTwistMatrix cVe;
		cVe.buildFrom(cMe);

		Maxi = -1 * L_total * cVe * eJe;
		Ma = Ts * Landa * Maxi * qdot;
		
		for (unsigned int i = 0; i < 8; i++) {
			sEstim[i] += Ma[i];
		}
		epsil = task.computeEpsilon(sEstim);
		sDesi = task.computeSDes(epsil);
		task.Set_sDesi(sDesi);
		task.Set_sEstim(sEstim);

		myOptContLaw.optimize(jointvel, minf);	

		qdot = vpColVector(jointvel);
		q += Ts * qdot;
		error_value = task.computeError();


		vpDisplay::display(Iint);
		display_trajectory(Iint, point, cMo, cam);
		ControlAlgDisplay::display(task, cam, Iint, vpColor::green, vpColor::red);
		vpDisplay::flush(Iint);

		if (error_value < 1e-6)
		{
			std::cout << "Control task is successfully completed!" << std::endl;
			return 0;
		}
    }
	task.kill();
 }
 catch (const vpException &e) {
	    std::cout << "Catch an exception: " << e << std::endl;
		std::cin.get();
 }
}