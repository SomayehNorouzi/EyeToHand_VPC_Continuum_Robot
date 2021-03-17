#ifndef ControlAlg_H
#define ControlAlg_H
#include <list>
#include <visp3/core/vpMatrix.h>
#include <visp3/core/vpVelocityTwistMatrix.h>
#include <visp3/visual_features/vpBasicFeature.h>
#include <visp3/vs/vpAdaptiveGain.h>
#include <visp3/vs/vpServoException.h>
#include <visp/vpTranslationVector.h>

class ControlAlg
{
public:
	typedef enum {
		NONE,
		EYEINHAND_CAMERA,
		EYEINHAND_L_cVe_eJe,
		EYETOHAND_L_cVe_eJe,
		EYETOHAND_L_cVf_fVe_eJe,
		EYETOHAND_L_cVf_fJe
		} vpServoType;

	typedef enum {
		CURRENT,
		DESIRED,
		MEAN,
		USER_DEFINED
		} vpServoIteractionMatrixType;

	typedef enum {
		TRANSPOSE,     /*!< In the control law (see vpServo::vpServoType), uses the transpose instead of the pseudo inverse. */
		PSEUDO_INVERSE /*!< In the control law (see vpServo::vpServoType), uses the pseudo inverse. */
	} vpServoInversionType;

	typedef enum {
		ALL,                /*!< Print all the task information. */
		CONTROLLER,         /*!< Print the type of controller law. */
		ERROR_VECTOR,       /*!< Print the error vector \f$\bf e = (s-s^*)\f$. */
		FEATURE_CURRENT,    /*!< Print the current features \f$\bf s\f$. */
		FEATURE_DESIRED,    /*!< Print the desired features \f${\bf s}^*\f$. */
		GAIN,               /*!< Print the gain \f$\lambda\f$. */
		INTERACTION_MATRIX, /*!< Print the interaction matrix. */
		MINIMUM             /*!< Same as vpServo::vpServoPrintType::ERROR_VECTOR. */
	} vpServoPrintType;

public:
	
	ControlAlg();// default constructor
	explicit ControlAlg(vpServoType servoType);// constructor with Choice of the visual servoing control law
	virtual ~ControlAlg();// destructor
	void addFeature(vpBasicFeature &s, vpBasicFeature &s_star, const unsigned int select = vpBasicFeature::FEATURE_ALL);// create a new ste of  two visual features
	vpColVector computeErrorVPC(const vpColVector& sDes, const vpColVector& sEstimated);
	double computeError();
	vpColVector computeEpsilon(const vpColVector& sEstimated);
	vpColVector computeSDes(const vpColVector& epsilon);
	vpMatrix computeInteractionMatrix();
	inline vpColVector getEpsilon() const { return epsilon; }
	inline vpMatrix getInteractionMatrix() const { return L; }
	
	double costcont(const std::vector<double>& qdot_v, std::vector<double>& grad, void* data);
	//double vision_const(const std::vector<double>& q, void* data, double tol);
	inline vpColVector getTaskSingularValues() const { return sv; }
	
	vpVelocityTwistMatrix get_cVe() const { return cVe; }
	vpVelocityTwistMatrix get_cVf() const { return cVf; }
	vpVelocityTwistMatrix get_fVe() const { return fVe; }
	vpMatrix get_eJe() const { return eJe; }
	vpMatrix get_fJe() const { return fJe; }
	double vision_const(const std::vector<double>& x, std::vector<double>& grad, void* data);
	double compute_vision_error_vpc(const vpColVector sEstim);

	void kill();
   	void setForceInteractionMatrixComputation(bool force_computation)
	{
		this->forceInteractionMatrixComputation = force_computation;
	}

	void setInteractionMatrixType(const vpServoIteractionMatrixType &interactionMatrixType,
		const vpServoInversionType &interactionMatrixInversion = PSEUDO_INVERSE);
	void setLambda(double c) { lambda.initFromConstant(c); }
	void setLambda(const double gain_at_zero, const double gain_at_infinity, const double slope_at_zero)
	{
		lambda.initStandard(gain_at_zero, gain_at_infinity, slope_at_zero);
	}
	void setLambda(const vpAdaptiveGain &l) { lambda = l; }
	void setMu(double mu_) { this->mu = mu_; }
	void setServo(const vpServoType &servo_type);
	void set_cVe(const vpVelocityTwistMatrix &cVe_)
	{
		this->cVe = cVe_;
		init_cVe = true;
	}
	
	void set_cVe(const vpHomogeneousMatrix &cMe)
	{
		cVe.buildFrom(cMe);
		init_cVe = true;
	}
	
	void set_cVf(const vpVelocityTwistMatrix &cVf_)
	{
		this->cVf = cVf_;
		init_cVf = true;
	}
	
	void set_cVf(const vpHomogeneousMatrix &cMf)
	{
		cVf.buildFrom(cMf);
		init_cVf = true;
	}
	
	void set_fVe(const vpVelocityTwistMatrix &fVe_)
	{
		this->fVe = fVe_;
		init_fVe = true;
	}
	
	void set_fVe(const vpHomogeneousMatrix &fMe)
	{
		fVe.buildFrom(fMe);
		init_fVe = true;
	}

	void set_eJe(const vpMatrix &eJe_)
	{
		this->eJe = eJe_;
		init_eJe = true;
	}
	
	void set_fMc(const vpHomogeneousMatrix& fMc_){this->fMc = fMc_;}


	void set_fJe(const vpMatrix &fJe_)
	{
		this->fJe = fJe_;
		init_fJe = true;
	}

	void set_vision_boundaries(double RightBound, double LeftBound, double UpperBound, double DownBound) {
		_UpperBound = UpperBound;
		_DownBound  = DownBound;
		_RightBound = RightBound;
		_LeftBound  = LeftBound;
	}
	
	void set_cam_intrinsic(double px, double py, double u0, double v0) {
		_px = px;
		_py = py;
		_u0 = u0;
		_v0 = v0;

	}

	void Setq(vpColVector& q);
	void Set_landa(double landa);
	void SetInteractionMatrix(vpMatrix& Lint);
	void SetNp(int Np);
	void Set_sDesi(const vpColVector sDesi);
	void Set_sEstim(const vpColVector sEstim);
	void set_delta_t(double SamplingTime);

protected:
	void init();
public:
	vpMatrix L;
	vpColVector error;
	vpColVector epsilon;
	vpMatrix J1;
	vpMatrix J1p;
	vpColVector s;
	vpColVector sDes;
	vpColVector sStar;
	vpColVector sEstimated;
	vpColVector e1;
	vpColVector e;
	vpColVector q_dot;
	vpColVector v;
	vpServoType servoType;
	unsigned int rankJ1;
	std::list<vpBasicFeature *> featureList;
	std::list<vpBasicFeature *> desiredFeatureList;
	std::list<vpBasicFeature*> DesFeatureList;
	std::list<vpBasicFeature*> estimatedFeatureList;

	std::list<unsigned int> featureSelectionList;

	vpAdaptiveGain lambda;
	int signInteractionMatrix;
	vpServoIteractionMatrixType interactionMatrixType;
	vpServoInversionType inversionType;
	
protected:
	vpColVector q_axi;
	vpColVector sDesi_axi;
	vpColVector sEstim_axi;
	vpColVector sv;
	vpColVector e1_initial;

	vpVelocityTwistMatrix cVe;
	vpVelocityTwistMatrix cVf;
	vpVelocityTwistMatrix fVe;

	vpHomogeneousMatrix fMc;

	vpMatrix eJe;
	vpMatrix fJe;
	vpMatrix _Ls;
	vpMatrix WpW;
	vpMatrix I_WpW;
	vpMatrix P;
	vpMatrix cJc;

	int _Np;
	unsigned int dim_task;
	double _landa;
	double mu;
	double _ts;

	double _UpperBound;
	double _DownBound;
	double _RightBound;
	double _LeftBound;

	double _px;
	double _py;
	double _u0;
	double _v0;

	bool init_cVe;
	bool init_cVf;
	bool init_fVe;
	bool init_eJe;
	bool init_fJe;
	bool errorComputed;
	bool epsilonComputed;
	bool sDesComputed;
	bool interactionMatrixComputed;
	bool taskWasKilled;
	bool forceInteractionMatrixComputation;
	bool iscJcIdentity;
};

#endif
