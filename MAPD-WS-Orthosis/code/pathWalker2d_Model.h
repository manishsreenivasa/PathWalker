#ifndef _PATHWALKER2D_MODEL_H
#define _PATHWALKER2D_MODEL_H

#include <vector>
#include <rbdl/rbdl.h>
#include <rbdl/addons/customforces/customforces.h>
#include <rbdl/addons/geometry/geometry.h>
#include <limits>

#include "pathWalker2d_Enums.h"
#include "pathWalker2d_EnumMap.h"

struct Point {
	Point() :
		name ("unknown"),
		body_id (-1),
		body_name (""),
		point_local (
				std::numeric_limits<double>::signaling_NaN(),
				std::numeric_limits<double>::signaling_NaN(),
				std::numeric_limits<double>::signaling_NaN()
				)
	{ }

	std::string name;
	unsigned int body_id;
	std::string body_name;
	RigidBodyDynamics::Math::Vector3d point_local;
};

/** Data of a single constraint */
struct ConstraintInfo {
	ConstraintInfo() :
		point_id (PointNameLast),
		point_name (""),
		normal (
				std::numeric_limits<double>::signaling_NaN(),
				std::numeric_limits<double>::signaling_NaN(),
				std::numeric_limits<double>::signaling_NaN()
				) {
	}
	unsigned int point_id;
	std::string point_name;
	RigidBodyDynamics::Math::Vector3d normal;
};

/** Structure that holds data of a complete constraint set */
struct ConstraintSetInfo {
	ConstraintSetInfo() :
		name ("undefined") {
	}
	std::vector<ConstraintInfo> constraints;
	std::string name;
};

struct OrthosisInfo {
  std::string orthosis_name;
  double angleHardStopInPlantarFlexion, angleHardStopInDorsiFlexion;
  double angleCenterOffset;
  double stiffnessEndRangeInDorsiFlexion, stiffnessEndRangeInPlantarFlexion;
  double torquePreload, anglePreloadWindow;
  double stiffnessInDorsiFlexion, stiffnessInPlantarFlexion;
  bool loadParamsFromFile(const char* filename);
  double jointOffsetAngle, jointAngleSign, jointTorqueSign;
};

struct WalkerModel {
	WalkerModel();

	unsigned int activeConstraintSet;
	unsigned int nDof;
	unsigned int nActuatedDof;

	bool dynamicsComputed;
	bool kinematicsUpdated;

	RigidBodyDynamics::Math::VectorNd q;
	RigidBodyDynamics::Math::VectorNd q_max;
	RigidBodyDynamics::Math::VectorNd q_min;
	RigidBodyDynamics::Math::VectorNd qdot_max;
	RigidBodyDynamics::Math::VectorNd qdot;
	RigidBodyDynamics::Math::VectorNd qddot;
	RigidBodyDynamics::Math::VectorNd tau_condensed;
	RigidBodyDynamics::Math::VectorNd tau_expanded;
	RigidBodyDynamics::Math::VectorNd a;
	RigidBodyDynamics::Math::VectorNd adot;
	RigidBodyDynamics::Math::VectorNd e;
	RigidBodyDynamics::Math::VectorNd p;

	RigidBodyDynamics::Model model;

	RigidBodyDynamics::Addons::CustomForces::Anderson2007TorqueMuscle RightHipExtension_TorqueMuscle;
	RigidBodyDynamics::Addons::CustomForces::Anderson2007TorqueMuscle RightHipFlexion_TorqueMuscle;
	RigidBodyDynamics::Addons::CustomForces::Anderson2007TorqueMuscle LeftHipExtension_TorqueMuscle;
	RigidBodyDynamics::Addons::CustomForces::Anderson2007TorqueMuscle LeftHipFlexion_TorqueMuscle;	
	RigidBodyDynamics::Addons::CustomForces::Anderson2007TorqueMuscle RightKneeExtension_TorqueMuscle;
	RigidBodyDynamics::Addons::CustomForces::Anderson2007TorqueMuscle RightKneeFlexion_TorqueMuscle;
	RigidBodyDynamics::Addons::CustomForces::Anderson2007TorqueMuscle LeftKneeExtension_TorqueMuscle;
	RigidBodyDynamics::Addons::CustomForces::Anderson2007TorqueMuscle LeftKneeFlexion_TorqueMuscle;
	RigidBodyDynamics::Addons::CustomForces::Anderson2007TorqueMuscle RightAnkleExtension_TorqueMuscle;
	RigidBodyDynamics::Addons::CustomForces::Anderson2007TorqueMuscle RightAnkleFlexion_TorqueMuscle;
	RigidBodyDynamics::Addons::CustomForces::Anderson2007TorqueMuscle LeftAnkleExtension_TorqueMuscle;
	RigidBodyDynamics::Addons::CustomForces::Anderson2007TorqueMuscle LeftAnkleFlexion_TorqueMuscle;	
	RigidBodyDynamics::Addons::CustomForces::Anderson2007TorqueMuscle TorsoExtension_TorqueMuscle;
	RigidBodyDynamics::Addons::CustomForces::Anderson2007TorqueMuscle TorsoFlexion_TorqueMuscle;

	// Modeling orthosis spring characteristics
	RigidBodyDynamics::Addons::CustomForces::NeuroSwingTorsionSpring orthosis_r;
	RigidBodyDynamics::Addons::CustomForces::NeuroSwingTorsionSpring orthosis_l;
	OrthosisInfo orthosisParams_r, orthosisParams_l;
	
	Point pointInfos[PointNameLast];
	/// Information of the constraint sets (mostly used when parsing Lua file)
	ConstraintSetInfo constraintSetInfos[ConstraintSetNameLast];

	/// RDBL constraint sets that are used during forward dynamics and collision computations
	std::vector<RigidBodyDynamics::ConstraintSet> constraintSets;

	/// Copies state information from MUSCOD to the model and switches to the given constraint set.
	void updateState (const double *sd, const double *u, const double *p, const ConstraintSetName cs_name);

	/// Compute torque from q, qdot and activation
	void computeTorqueFromActivation ();
	
	/// Compute torque from q, qdot and activation
	void computeActivationDynamics ();
	
	/// Compute excitation dot from control
	void computeExcitation(const double *u);

	/// Updates the kinematics of the RDBL model (is called automatically when required)
	void updateKinematics ();

	/// Computes the forward dynamics for the model and active constraint set
	void calcForwardDynamicsRhs (double *res);
	
	RigidBodyDynamics::Math::Vector3d getPointPosition (const PointName &point_name);
	RigidBodyDynamics::Math::Vector3d getPointVelocity (const PointName &point_name);
	RigidBodyDynamics::Math::Vector3d getPointForce (const PointName &point_name);

	bool loadFromFile (const char* filename, std::string datfilename, bool verbose=false);
	void setupTorqueMuscles (std::string datfilename, bool verbose=false);
	bool loadPoints (const char* filename, bool verbose=false);
	bool loadConstraintSets (const char* filename, bool verbose=false);
	bool loadOrthosisParams (const char* filename, bool verbose=false);
};

/* _WALKER_MODEL_H */
#endif
