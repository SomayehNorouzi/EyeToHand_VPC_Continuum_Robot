#include <sstream>
#include "stdafx.h"
#include "ControlAlg.h"
#include <visp3/core/vpException.h>
#include <visp3/core/vpDebug.h>
#include "CRSimulation.h"
#include <visp3/core/vpMatrix.h>
#include <visp/vpTranslationVector.h>
#include <cmath>   
#include <visp3/visual_features/vpFeatureBuilder.h>
ControlAlg::ControlAlg()
	: L(), error(), epsilon(), J1(), J1p(), s(), sDes(), sStar(), sEstimated(), e1(), e(), q_dot(), v(), servoType(), featureList(),
	desiredFeatureList(), DesFeatureList(), estimatedFeatureList(), featureSelectionList(), lambda(), signInteractionMatrix(1),
	interactionMatrixType(DESIRED), inversionType(PSEUDO_INVERSE), q_axi(), _Ls(), _Np(), sDesi_axi(), sEstim_axi(), cVe(),
	init_cVe(false), cVf(), init_cVf(false), fVe(), init_fVe(false), _landa(0.5), eJe(), init_eJe(false), fJe(),
	init_fJe(false), errorComputed(false), epsilonComputed(false), sDesComputed(false), interactionMatrixComputed(false),
	dim_task(0), taskWasKilled(false), forceInteractionMatrixComputation(false), WpW(), I_WpW(), P(), sv(), mu(4), e1_initial(),
	iscJcIdentity(false), _ts(0), cJc(6,6)
{
	cJc.eye();
}

ControlAlg::ControlAlg(vpServoType servo_type)
	: L(), error(), epsilon(), J1(), J1p(), s(), sDes(), sStar(), sEstimated(), e1(), e(), q_dot(), v(), servoType(), featureList(),
	desiredFeatureList(), DesFeatureList(), estimatedFeatureList(), featureSelectionList(), lambda(), signInteractionMatrix(1),
	interactionMatrixType(DESIRED), inversionType(PSEUDO_INVERSE), q_axi(), _Ls(), _Np(), sDesi_axi(), sEstim_axi(), cVe(),
	init_cVe(false), cVf(), init_cVf(false), fVe(), init_fVe(false), _landa(0.5), eJe(), init_eJe(false), fJe(),
	init_fJe(false), errorComputed(false), epsilonComputed(false), sDesComputed(false), interactionMatrixComputed(false),
	dim_task(0), taskWasKilled(false), forceInteractionMatrixComputation(false), WpW(), I_WpW(), P(), sv(), mu(4), e1_initial(),
	iscJcIdentity(false), _ts(0), cJc(6,6)
{
	cJc.eye();
}

ControlAlg::~ControlAlg()
{
	if (taskWasKilled == false) {
		vpTRACE("--- Begin Warning Warning Warning Warning Warning ---");
		vpTRACE("--- You should explicitly call ControlAlg.kill()...  ---");
		vpTRACE("--- End Warning Warning Warning Warning Warning   ---");
	}
}

void ControlAlg::init()
{
	servoType = ControlAlg::NONE;
	init_cVe = false;// Twist transformation matrix
	init_cVf = false;
	init_fVe = false;
	init_eJe = false;
	init_fJe = false;
	dim_task = 0;
	featureList.clear();
	desiredFeatureList.clear();
	DesFeatureList.clear();
	estimatedFeatureList.clear();
	featureSelectionList.clear();
	signInteractionMatrix = 1;
	interactionMatrixType = DESIRED;
	inversionType = PSEUDO_INVERSE;
	interactionMatrixComputed = false;
	errorComputed = false;
	taskWasKilled = false;
	forceInteractionMatrixComputation = false;
	rankJ1 = 0;
}

void ControlAlg::kill()
{
	if (taskWasKilled == false) {
		// kill the current and desired feature lists

		// current list
		for (std::list<vpBasicFeature *>::iterator it = featureList.begin(); it != featureList.end(); ++it) {
			if ((*it)->getDeallocate() == vpBasicFeature::vpServo) {
				delete (*it);
				(*it) = NULL;
			}
		}

		// estimated list
		for (std::list<vpBasicFeature*>::iterator it = estimatedFeatureList.begin(); it != estimatedFeatureList.end(); ++it) {
			if ((*it)->getDeallocate() == vpBasicFeature::vpServo) {
				delete (*it);
				(*it) = NULL;
			}
		}

		// desired list
		for (std::list<vpBasicFeature *>::iterator it = desiredFeatureList.begin(); it != desiredFeatureList.end(); ++it) {
			if ((*it)->getDeallocate() == vpBasicFeature::vpServo) {
				delete (*it);
				(*it) = NULL;
			}
		}
		// VPC desired list (Des list)
		for (std::list<vpBasicFeature*>::iterator it = DesFeatureList.begin(); it != DesFeatureList.end(); ++it) {
			if ((*it)->getDeallocate() == vpBasicFeature::vpServo) {
				delete (*it);
				(*it) = NULL;
			}
		}

		featureList.clear();
		estimatedFeatureList.clear();
		desiredFeatureList.clear();
		DesFeatureList.clear();
		taskWasKilled = true;
	}
}


void ControlAlg::setServo(const vpServoType &servo_type)
{
	this->servoType = servo_type;

	if ((servoType == EYEINHAND_CAMERA) || (servoType == EYEINHAND_L_cVe_eJe))
		signInteractionMatrix = 1;
	else
		signInteractionMatrix = -1;

	if (servoType == EYEINHAND_CAMERA) {
		vpVelocityTwistMatrix _cVe;
		set_cVe(_cVe);

		vpMatrix _eJe;
		_eJe.eye(6);
		set_eJe(_eJe);
	};
}

void ControlAlg::addFeature(vpBasicFeature &s_cur, vpBasicFeature &s_star, unsigned int select)
{
	featureList.push_back(&s_cur);
	desiredFeatureList.push_back(&s_star);
	
	featureSelectionList.push_back(select);

}

void ControlAlg::setInteractionMatrixType(const vpServoIteractionMatrixType &interactionMatrix_type,
	const vpServoInversionType &interactionMatrixInversion)
{
	this->interactionMatrixType = interactionMatrix_type;
	this->inversionType = interactionMatrixInversion;
}

static void computeInteractionMatrixFromList(const std::list<vpBasicFeature *> &featureList,
	const std::list<unsigned int> &featureSelectionList, vpMatrix &L)
{
	if (featureList.empty()) {
		vpERROR_TRACE("feature list empty, cannot compute Ls");
		throw(vpServoException(vpServoException::noFeatureError, "feature list empty, cannot compute Ls"));
	}

	unsigned int rowL = L.getRows();
	const unsigned int colL = 6;
	if (0 == rowL) {
		rowL = 1;
		L.resize(rowL, colL);
	}

	vpMatrix matrixTmp;

	unsigned int cursorL = 0;

	std::list<vpBasicFeature *>::const_iterator it;
	std::list<unsigned int>::const_iterator it_select;

	for (it = featureList.begin(), it_select = featureSelectionList.begin(); it != featureList.end(); ++it, ++it_select) {
		
		matrixTmp = (*it)->interaction(*it_select);
		unsigned int rowMatrixTmp = matrixTmp.getRows();
		unsigned int colMatrixTmp = matrixTmp.getCols();

		while (rowMatrixTmp + cursorL > rowL) {
			rowL *= 2;
			L.resize(rowL, colL, false);
			vpDEBUG_TRACE(15, "Realloc!");
		}

		/* Copy the temporarily matrix into L. */
		for (unsigned int k = 0; k < rowMatrixTmp; ++k, ++cursorL) {
			for (unsigned int j = 0; j < colMatrixTmp; ++j) {
				L[cursorL][j] = matrixTmp[k][j];
			}
		}
	}
	L.resize(cursorL, colL, false);
	return;
}

vpMatrix ControlAlg::computeInteractionMatrix()
{
	try {

		switch (interactionMatrixType) {
		case CURRENT: {
			try {
				computeInteractionMatrixFromList(this->featureList, this->featureSelectionList, L);
				dim_task = L.getRows();
				interactionMatrixComputed = true;
			}

			catch (...) {
				throw;
			}
		} break;
		case DESIRED: {
			try {
				if (interactionMatrixComputed == false || forceInteractionMatrixComputation == true) {
					computeInteractionMatrixFromList(this->desiredFeatureList, this->featureSelectionList, L);

					dim_task = L.getRows();
					interactionMatrixComputed = true;
				}
			}
			catch (...) {
				throw;
			}
		} break;
		case MEAN: {
			vpMatrix Lstar(L.getRows(), L.getCols());
			try {
				computeInteractionMatrixFromList(this->featureList, this->featureSelectionList, L);
				computeInteractionMatrixFromList(this->desiredFeatureList, this->featureSelectionList, Lstar);
			}
			catch (...) {
				throw;
			}
			L = (L + Lstar) / 2;

			dim_task = L.getRows();
			interactionMatrixComputed = true;
		} break;
		case USER_DEFINED:
			// dim_task = L.getRows() ;
			interactionMatrixComputed = false;
			break;
		}
	}
	catch (...) {
		throw;
	}
	return L;
}


vpColVector 
ControlAlg::computeErrorVPC(const vpColVector& sDes, const vpColVector& sEstimated)
{
	try {
		error.resize(8);
		for (int i = 0; i < 8; i++) {
			error[i] = sDes[i] - sEstimated[i];
		}
	}
	catch (...) {
		throw;
	}
	return error;
}

double
ControlAlg::computeError()
{
	vpBasicFeature* current_s;
	vpBasicFeature* star_s;

	unsigned int dimS = s.getRows();
	unsigned int dimSStar = sStar.getRows();

	if (0 == dimS) {
		dimS = 1;
		s.resize(dimS);
	}
	if (0 == dimSStar) {
		dimSStar = 1;
		sStar.resize(dimSStar);
	}
	vpColVector vectTmp;
	unsigned int cursorS = 0;

	vpColVector vectTmpStar;
	unsigned int cursorSStar = 0;

	std::list<vpBasicFeature*>::const_iterator it_s;
	std::list<unsigned int>::const_iterator it_select;

	std::list<vpBasicFeature*>::const_iterator it_sStar;

	for (it_s = featureList.begin(), it_select = featureSelectionList.begin();
		it_s != featureList.end(); ++it_s, ++it_select) {
		current_s = (*it_s);// 

		unsigned int select = (*it_select);

		vectTmp = current_s->get_s(select); // a column of 2*1
		unsigned int dimVectTmp = vectTmp.getRows();
		while (dimVectTmp + cursorS > dimS) {
			dimS *= 2;
			s.resize(dimS, false);
			vpDEBUG_TRACE(15, "Realloc!");
		}
		for (unsigned int k = 0; k < dimVectTmp; ++k) {
			s[cursorS++] = vectTmp[k];
		}
	}
	
	for (it_sStar = desiredFeatureList.begin(), it_select = featureSelectionList.begin();
		it_sStar != desiredFeatureList.end(); ++it_sStar, ++it_select) {
		star_s = (*it_sStar);// 
		unsigned int select = (*it_select);

		vectTmpStar = star_s->get_s(select); // a column of 2*1
		unsigned int dimVectTmp = vectTmpStar.getRows();

		while (dimVectTmp + cursorSStar > dimSStar) {
			dimSStar *= 2;
			sStar.resize(dimSStar, false);
			vpDEBUG_TRACE(15, "Realloc!");
		}
		for (unsigned int k = 0; k < dimVectTmp; ++k) {
			sStar[cursorSStar++] = vectTmpStar[k];
		}
	}
	vpColVector sStar_sError;
	sStar_sError.resize(8);
	for (int i = 0; i < 8; i++) {

		sStar_sError[i] = sStar[i] - s[i];
	}
	
	double ErrorValue = 0;
	for (unsigned int i = 1; i < 8; i++) {
		ErrorValue += sStar_sError[i] * sStar_sError[i];
	}
	
	return ErrorValue;
}

vpColVector
ControlAlg::computeEpsilon(const vpColVector& sEstimated)
{
	try {
		vpBasicFeature* current_s;

		unsigned int dimS = s.getRows();

		if (0 == dimS) {
			dimS = 1;
			s.resize(dimS);
		}
		vpColVector vectTmp;
		unsigned int cursorS = 0;

		std::list<vpBasicFeature*>::const_iterator it_s;
		std::list<unsigned int>::const_iterator it_select;

		for (it_s = featureList.begin(), it_select = featureSelectionList.begin();
			it_s != featureList.end(); ++it_s, ++it_select) {
			current_s = (*it_s);// 
			unsigned int select = (*it_select);
			vectTmp = current_s->get_s(select); // a column of 2*1
			unsigned int dimVectTmp = vectTmp.getRows();
			//unsigned int dimvectTmpSEst = vectTmpSEst.getRows();
			while (dimVectTmp + cursorS > dimS) {
				dimS *= 2;
				s.resize(dimS, false);
				vpDEBUG_TRACE(15, "Realloc!");
			}
			for (unsigned int k = 0; k < dimVectTmp; ++k) {
				s[cursorS++] = vectTmp[k];
			}
		}
		epsilon.resize(8);
		for (int i = 0; i < 8; i++) {
			epsilon[i] = s[i] - sEstimated[i];
		}
	}
	catch (...) {
		throw;
	}
	return epsilon;
}
vpColVector ControlAlg::computeSDes(const vpColVector& epsilon)
{
	   vpBasicFeature* star_s;

		unsigned int dimSStar = sStar.getRows();


		if (0 == dimSStar) {
			dimSStar = 1;
			sStar.resize(dimSStar);
		}
		vpColVector vectTmp;
		unsigned int cursorSStar = 0;

		std::list<vpBasicFeature*>::const_iterator it_sStar;
		std::list<unsigned int>::const_iterator it_select;

		for (it_sStar = desiredFeatureList.begin(), it_select = featureSelectionList.begin();
			it_sStar != desiredFeatureList.end(); ++it_sStar, ++it_select) {
			star_s = (*it_sStar);// 
			unsigned int select = (*it_select);

			vectTmp = star_s->get_s(select); // a column of 2*1
			unsigned int dimVectTmp = vectTmp.getRows();

			while (dimVectTmp + cursorSStar > dimSStar) {
				dimSStar *= 2;
				sStar.resize(dimSStar, false);
				vpDEBUG_TRACE(15, "Realloc!");
			}
			for (unsigned int k = 0; k < dimVectTmp; ++k) {
				sStar[cursorSStar++] = vectTmp[k];
			}
		}
		sDes.resize(8);
		for (int i = 0; i < 8; i++) {
			sDes[i] = sStar[i] - epsilon[i];
		}
	return sDes;
}

void
ControlAlg::Setq(vpColVector& q)
{
	q_axi.resize(3, 1);

	q_axi[0] = q[0];
	q_axi[1] = q[1];
	q_axi[2] = q[2];
	return;
}

void
ControlAlg::Set_landa(double landa) {

	_landa = landa;
	return;
}

void 
ControlAlg::SetInteractionMatrix(vpMatrix& Lint) {

	_Ls.resize(8, 6);

	for (unsigned int i = 0; i < 8; i++) {
		for (unsigned int j = 0; j < 6; j++) {
			_Ls[i][j] = Lint[i][j];
		}
	}
	return;
}

void 
ControlAlg::Set_sDesi(const vpColVector sDesi) {
	sDesi_axi.resize(8);

	for (unsigned int i = 0; i < 8; i++) {
			sDesi_axi[i] = sDesi[i];
		}
	return;
}

void
ControlAlg::Set_sEstim(const vpColVector sEstim) {

	sEstim_axi.resize(8);

	for (unsigned int i = 0; i < 8; i++) {
		sEstim_axi[i] = sEstim[i];
	}
	return;
}


void 
ControlAlg::SetNp(int Np) {
	_Np = Np;
}
void
ControlAlg::set_delta_t(double SamplingTime)
{
	_ts = SamplingTime;
}

double ControlAlg::costcont(const std::vector<double> &qdot_v, std::vector<double> &grad, void *data)
{
	CRSimulation robot_axi;
	vpColVector qdot;
	qdot.resize(3);
	qdot = vpColVector(qdot_v);

	vpColVector _q(3);
	_q[0] = q_axi[0]; _q[1] = q_axi[1]; _q[2] = q_axi[2];

	vpColVector _sEstim(8);
	for (int i = 0; i < 8; i++) {
		_sEstim[i] = sEstim_axi[i];  //_sEstim[1] = sEstim_axi[1];  _sEstim[2] = sEstim_axi[2]; _sEstim[3] = sEstim_axi[3]; _sEstim[4] = sEstim_axi[4]; _sEstim[5] = sEstim_axi[5]; _sEstim[6] = sEstim_axi[6]; _sEstim[7] = sEstim_axi[7];
	}
	vpHomogeneousMatrix fMe_axi;
	vpHomogeneousMatrix eMc_axi;
	vpHomogeneousMatrix cMe_axi;
	vpHomogeneousMatrix cMo_axi;

	vpColVector ErrorS;
	ErrorS.resize(8);

	double cost_val = 0;

	vpVelocityTwistMatrix cVe_axi;
	vpMatrix Maxi_axi;
	Maxi_axi.resize(8, 3);
	vpColVector Ma_axi;
	Ma_axi.resize(8);

	for (int i = 0; i < _Np; i++){
		fMe_axi = robot_axi.get_fMe(_q);
		eMc_axi = fMe_axi.inverse() * fMc;
		cMe_axi = eMc_axi.inverse();
		cMo_axi = eMc_axi.inverse();
		cVe_axi.buildFrom(cMe_axi);
		Maxi_axi = -1 * _Ls * cVe_axi* eJe;
		Ma_axi = _ts * _landa * Maxi_axi * qdot;

		for (unsigned int i = 0; i < 8; i++) {
			_sEstim[i] += Ma_axi[i];
	    }
	_q += _ts * qdot;
	}

	ErrorS = computeErrorVPC(sDesi_axi, _sEstim);
	for (unsigned int i = 0; i < 8; i++) {
		cost_val += ErrorS[i] * ErrorS[i];
	}
	
	return cost_val;
}

double 
ControlAlg::compute_vision_error_vpc(const vpColVector sEstim)
{
	vpColVector upConst(4), downConst(4), leftConst(4), rightConst(4); 
	
	for (unsigned int i = 0; i < 4; i++) {
		upConst[i] = _UpperBound - ((_py * sEstim[2*i+1]) + _v0);
		downConst[i] = ((_py * sEstim[2*i+1]) + _v0) - _DownBound;
		rightConst[i] = ((_px * sEstim[2*i]) + _u0) - _RightBound;
		leftConst[i] = _LeftBound - ((_px * sEstim[2*i]) + _u0);
	}
	
	double max_val;
	for (unsigned int i = 0; i < 4; i++) {
	max_val = vpMath::maximum(max_val, upConst[i]);
	max_val = vpMath::maximum(max_val, downConst[i]);
	max_val = vpMath::maximum(max_val, rightConst[i]);
	max_val = vpMath::maximum(max_val, leftConst[i]);
	}
	return max_val;
}

double ControlAlg::vision_const(const std::vector<double>& x, std::vector<double>& grad, void* data)

{
	CRSimulation robot_axi;
	vpColVector qdot(3);
	qdot = vpColVector(x);
	
	vpColVector _q(3);
	_q[0] = q_axi[0]; _q[1] = q_axi[1]; _q[2] = q_axi[2];
	_q += _ts * qdot;

	vpHomogeneousMatrix fMe_axi = robot_axi.get_fMe(_q);
	vpHomogeneousMatrix eMc_axi = fMe_axi.inverse() * fMc;
	vpHomogeneousMatrix cMo_axi = eMc_axi.inverse();

	std::vector<vpPoint> point_axi(4);
	point_axi[0] = vpPoint(-0.0035, 0, -0.0035);// meter
	point_axi[1] = vpPoint(0.0035, 0, -0.0035);
	point_axi[2] = vpPoint(0.0035, 0, 0.0035);
	point_axi[3] = vpPoint(-0.0035, 0, 0.0035);

	vpFeaturePoint p_axi[4];
	vpColVector sEstim(8);
		for (unsigned int i = 0; i < 4; i++) {

			point_axi[i].track(cMo_axi);
			vpFeatureBuilder::create(p_axi[i], point_axi[i]);
			
			sEstim[2 * i] = p_axi[i].get_x(); 	sEstim[2 * i + 1] = p_axi[i].get_y();
		}
	
	double vision_cost_val = 0;
	double ErrorVision=5;
	ErrorVision = compute_vision_error_vpc(sEstim);

	return ErrorVision;
}
