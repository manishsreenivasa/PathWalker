#include <sys/types.h>
#include <sys/wait.h>
#include <iomanip>
#include "pathWalker2d_Model.h"

#include <rbdl/addons/luamodel/luamodel.h>
#include <rbdl/addons/luamodel/luatables.h>
#include <rbdl/addons/customforces/customforces.h>
#include <rbdl/addons/geometry/geometry.h>

#include "LuaTypes.h"
#include <boost/config/no_tr1/complex.hpp>

#include <string>

using namespace std;
using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;
using namespace RigidBodyDynamics::Addons::CustomForces;
using namespace RigidBodyDynamics::Addons::Geometry;
#include "datfileutils.h"

WalkerModel::WalkerModel() :
		activeConstraintSet (-1),
		nDof (-1),
		nActuatedDof(0) {
	// Create invalid values for point infos
	for (unsigned int i = 0; i < PointNameLast; i++) {
		pointInfos[i] = Point();
	}

	// Create invalid values for constraint set infos
	for (unsigned int i = 0; i < ConstraintSetNameLast; i++) {
		constraintSetInfos[i] = ConstraintSetInfo();
	}

	p = VectorNd::Zero (ParamNameLast);
}

void WalkerModel::updateState (const double *sd, const double *u, const double *p_in, const ConstraintSetName cs_name) {
	activeConstraintSet = cs_name;

	dynamicsComputed = false;
	kinematicsUpdated = false;
	
	unsigned int k = 0;
	for (unsigned int i = 0; i < PosNameLast; i++) {
	  q[i] = sd[i];
	  qdot[i] = sd[i + PosNameLast];
	  k = i + PosNameLast + 1;
	}
	for (unsigned int j = 0; j < ControlNameLast; j++) {
	  a[j] = sd[k];
	  k = k + 1;
	}
	
	// First apply parameters so that torques are computed with new muscle torque values
	for (unsigned int i = 0; i < ParamNameLast; i++) {
	  p[i] = p_in[i];
	}

	computeExcitation(u);
	computeActivationDynamics();
	computeTorqueFromActivation();
	
	unsigned int ctrl_idx = 0;
	for (unsigned int i = 0; i < nActuatedDof; i++) {
	  tau_condensed[i + nDof - nActuatedDof] = tau_expanded[ctrl_idx] + tau_expanded[ctrl_idx+1];
	  ctrl_idx+=2;
	}
	
	// Update orthosis parameters
	orthosis_r.setStiffnessInDorsiFlexion (p[ParamRightOrthosisStiffnessInDorsiFlexion]);
	orthosis_r.setStiffnessInPlantarFlexion (p[ParamRightOrthosisStiffnessInPlantarFlexion]);
	orthosis_l.setStiffnessInDorsiFlexion (p[ParamLeftOrthosisStiffnessInDorsiFlexion]);
	orthosis_l.setStiffnessInPlantarFlexion (p[ParamLeftOrthosisStiffnessInPlantarFlexion]);

	// Compute torque contribution from orthosis springs
	tau_condensed[PosRightAnkleRotY] = tau_condensed[PosRightAnkleRotY] + orthosis_r.calcJointTorque (sd[StateRightAnkleRotY]);
	tau_condensed[PosLeftAnkleRotY]  = tau_condensed[PosLeftAnkleRotY] + orthosis_l.calcJointTorque (sd[StateLeftAnkleRotY]);
}

// From current activations compute agonist/antagonist muscle torques
void WalkerModel::computeTorqueFromActivation () {
  
  // Compute torque muscle dynamics based on Anderson model
  tau_expanded[ControlRightHipExtensionRotY] 	= RightHipExtension_TorqueMuscle.calcJointTorque(q[3], qdot[3], a[ControlRightHipExtensionRotY]);  
  tau_expanded[ControlRightHipFlexionRotY] 	= RightHipFlexion_TorqueMuscle.calcJointTorque(q[3], qdot[3], a[ControlRightHipFlexionRotY]);
  tau_expanded[ControlRightKneeExtensionRotY] 	= RightKneeExtension_TorqueMuscle.calcJointTorque(q[4], qdot[4], a[ControlRightKneeExtensionRotY]);
  tau_expanded[ControlRightKneeFlexionRotY] 	= RightKneeFlexion_TorqueMuscle.calcJointTorque(q[4], qdot[4], a[ControlRightKneeFlexionRotY]);    
  tau_expanded[ControlRightAnkleExtensionRotY] 	= RightAnkleExtension_TorqueMuscle.calcJointTorque(q[5], qdot[5], a[ControlRightAnkleExtensionRotY]);
  tau_expanded[ControlRightAnkleFlexionRotY] 	= RightAnkleFlexion_TorqueMuscle.calcJointTorque(q[5], qdot[5], a[ControlRightAnkleFlexionRotY]);  
  tau_expanded[ControlLeftHipExtensionRotY] 	= LeftHipExtension_TorqueMuscle.calcJointTorque(q[6], qdot[6], a[ControlLeftHipExtensionRotY]);
  tau_expanded[ControlLeftHipFlexionRotY] 	= LeftHipFlexion_TorqueMuscle.calcJointTorque(q[6], qdot[6], a[ControlLeftHipFlexionRotY]);  
  tau_expanded[ControlLeftKneeExtensionRotY] 	= LeftKneeExtension_TorqueMuscle.calcJointTorque(q[7], qdot[7], a[ControlLeftKneeExtensionRotY]);
  tau_expanded[ControlLeftKneeFlexionRotY] 	= LeftKneeFlexion_TorqueMuscle.calcJointTorque(q[7], qdot[7], a[ControlLeftKneeFlexionRotY]);    
  tau_expanded[ControlLeftAnkleExtensionRotY] 	= LeftAnkleExtension_TorqueMuscle.calcJointTorque(q[8], qdot[8], a[ControlLeftAnkleExtensionRotY]);
  tau_expanded[ControlLeftAnkleFlexionRotY] 	= LeftAnkleFlexion_TorqueMuscle.calcJointTorque(q[8], qdot[8], a[ControlLeftAnkleFlexionRotY]);  
  tau_expanded[ControlTorsoExtensionRotY] 	= TorsoExtension_TorqueMuscle.calcJointTorque(q[9], qdot[9], a[ControlTorsoExtensionRotY]);
  tau_expanded[ControlTorsoFlexionRotY] 	= TorsoFlexion_TorqueMuscle.calcJointTorque(q[9], qdot[9], a[ControlTorsoFlexionRotY]);
}

void WalkerModel::computeExcitation (const double *u) { 
  // For now control = excitation, so no dynamics are needed
  for (unsigned int ctrl_idx = 0; ctrl_idx < ControlNameLast; ctrl_idx++) {
    e[ctrl_idx] = u[ctrl_idx];
  }
}

// From current controls and activation, compute first derivative as per Thelen 2013 model
void WalkerModel::computeActivationDynamics () {
  // Activation-Deactivation dynamics and parameters from Winter et al. 
  double tc_a = 0.011, tc_d = 0.068;
  for (unsigned int ctrl_idx = 0; ctrl_idx < ControlNameLast; ctrl_idx++) {
    if (e[ctrl_idx] >= a[ctrl_idx]) {
      adot[ctrl_idx] = (e[ctrl_idx] - a[ctrl_idx])*(e[ctrl_idx]/tc_a + (1-e[ctrl_idx])/tc_d);
    }
    else {
      adot[ctrl_idx] = (e[ctrl_idx] - a[ctrl_idx])/tc_d;
    }
  }
}

void WalkerModel::updateKinematics () {
	UpdateKinematics (model, q, qdot, qddot);
	kinematicsUpdated = true;
}

void WalkerModel::calcForwardDynamicsRhs (double *res) {
	ForwardDynamicsContactsKokkevis (model, q, qdot, tau_condensed, constraintSets[activeConstraintSet], qddot);
	
	dynamicsComputed = true;
	unsigned int k = 0;
	for (unsigned int i = 0; i < PosNameLast; i++) {
	  res[i] = qdot[i];
	  res[i + PosNameLast] = qddot[i];
	  k = i + PosNameLast + 1;
	}
	for (unsigned int j = 0; j < ControlNameLast; j++) {
	  res[k] = adot[j];
	  k = k + 1;
	}
}

Vector3d WalkerModel::getPointPosition (const PointName &point_name) {
	if (!kinematicsUpdated) {
		updateKinematics();
	}

	unsigned int body_id = pointInfos[point_name].body_id;
	Vector3d point_local = pointInfos[point_name].point_local;

	return CalcBodyToBaseCoordinates (model, q, body_id, point_local, false);
}

Vector3d WalkerModel::getPointVelocity (const PointName &point_name) {
	if (!kinematicsUpdated) {
		updateKinematics();
	}

	unsigned int body_id = pointInfos[point_name].body_id;
	Vector3d point_local = pointInfos[point_name].point_local;

	return CalcPointVelocity (model, q, qdot, body_id, point_local, false);
}

Vector3d WalkerModel::getPointForce (const PointName &point_name) {
	if (!dynamicsComputed) {
		ForwardDynamicsContactsKokkevis (model, q, qdot, tau_condensed, constraintSets[activeConstraintSet], qddot);
		dynamicsComputed = true;
	}

	Vector3d result (0., 0., 0.);
	bool found = false;

	const RigidBodyDynamics::ConstraintSet &active_constraint_set = constraintSets[activeConstraintSet];
	const ConstraintSetInfo constraint_set_info = constraintSetInfos[activeConstraintSet];

	std::vector<ConstraintInfo>::const_iterator constraint_iter = constraintSetInfos[activeConstraintSet].constraints.begin();

	for (unsigned int ci = 0; ci < constraint_set_info.constraints.size(); ci++) {
		const ConstraintInfo& constraint_info = constraint_set_info.constraints[ci];

		if (constraint_info.point_id != point_name)
			continue;

		found = true;
		assert (constraint_info.normal == active_constraint_set.normal[ci]);

		result += active_constraint_set.force[ci] * active_constraint_set.normal[ci];
	}

	if (!found) {
		cerr << "Error (" << __func__ << "): Point '" << pointInfos[point_name].name << "' is not constrained in constraint set '" << constraintSetInfos[activeConstraintSet].name << "'!" << endl;
		abort();
	}

	return result;
}

bool WalkerModel::loadFromFile (const char* filename, string datfilename, bool verbose) {
	if (!RigidBodyDynamics::Addons::LuaModelReadFromFile (filename, &model, verbose)) {
	  cerr << "Error loading LuaModel: " << filename << endl;
	  abort();
	}

	nDof = model.dof_count;
	nActuatedDof = nDof - 3;
	
	assert (nActuatedDof >= 1 && nActuatedDof <= nDof);

	if (nDof != PosNameLast) {
	  cerr << "Error: Number of model degrees of freedom (" << nDof << ") does not match number of positional variables (" << PosNameLast << ")!" << endl;
	  abort();
	}

	q 		= VectorNd::Zero (nDof);
	qdot 		= VectorNd::Zero (nDof);
	qddot 		= VectorNd::Zero (nDof);
	tau_condensed 	= VectorNd::Zero (nDof);
	tau_expanded 	= VectorNd::Zero (ControlNameLast);
	a 		= VectorNd::Zero (ControlNameLast);
	adot 		= VectorNd::Zero (ControlNameLast);
	e 		= VectorNd::Zero (ControlNameLast);
	q_max 		= VectorNd::Zero (nDof);
	q_min 		= VectorNd::Zero (nDof);
	qdot_max 	= VectorNd::Zero (nDof);
	
	VectorNd sd_max = datfile_get_vector(datfilename.c_str(), "sd_max");
	VectorNd sd_min = datfile_get_vector(datfilename.c_str(), "sd_min");
	
	for (unsigned int i = 0; i < PosNameLast; i++) {
	  q_max[i] = sd_max[i];
	  q_min[i] = sd_min[i];
	  qdot_max[i] = sd_max[i + PosNameLast];
	}
	
	setupTorqueMuscles(datfilename, verbose);
	
	return true;
}

void WalkerModel::setupTorqueMuscles (std::string datfilename, bool verbose) {
	
	VectorNd jointStrength = datfile_get_vector(datfilename.c_str(), "p");
	
	int gender = 0; //male
	int age = 0; //young
	int hip = 0;
	int knee = 1;
	int ankle = 2;
	int ext = 0;
	int flex = 1;
	double subjectHeight = 1.25;
	double subjectMass = 24.7;
	double jointAngleOffset = 0;
	
	double angleSign = -1.0; double torqueSign = 1.0;
	RightHipExtension_TorqueMuscle 	 = Anderson2007TorqueMuscle (hip, ext, gender, age, subjectHeight, subjectMass, jointAngleOffset, angleSign, torqueSign, "RhipExt");
	angleSign  = -1.0; torqueSign = -1.0;
	RightHipFlexion_TorqueMuscle 	 = Anderson2007TorqueMuscle (hip, flex, gender, age, subjectHeight, subjectMass, jointAngleOffset, angleSign, torqueSign, "RhipFlex");
	angleSign  =  1.0; torqueSign = -1.0;
	RightKneeExtension_TorqueMuscle  = Anderson2007TorqueMuscle (knee, ext, gender, age, subjectHeight, subjectMass, jointAngleOffset, angleSign, torqueSign, "RkneeExt");
	angleSign  =  1.0; torqueSign =  1.0;
	RightKneeFlexion_TorqueMuscle 	 = Anderson2007TorqueMuscle (knee, flex, gender, age, subjectHeight, subjectMass, jointAngleOffset, angleSign, torqueSign, "RkneeFlex");	
	angleSign  = -1.0; torqueSign =  1.0;
	RightAnkleExtension_TorqueMuscle = Anderson2007TorqueMuscle (ankle, ext, gender, age, subjectHeight, subjectMass, jointAngleOffset, angleSign, torqueSign, "RankleExt");
	angleSign  = -1.0; torqueSign = -1.0;
	RightAnkleFlexion_TorqueMuscle 	 = Anderson2007TorqueMuscle (ankle, flex, gender, age, subjectHeight, subjectMass, jointAngleOffset, angleSign, torqueSign, "RankleFlex");	
	angleSign = -1.0; torqueSign = 1.0;
	LeftHipExtension_TorqueMuscle 	 = Anderson2007TorqueMuscle (hip, ext, gender, age, subjectHeight, subjectMass, jointAngleOffset, angleSign, torqueSign, "LhipExt");
	angleSign  = -1.0; torqueSign = -1.0;
	LeftHipFlexion_TorqueMuscle 	 = Anderson2007TorqueMuscle (hip, flex, gender, age, subjectHeight, subjectMass, jointAngleOffset, angleSign, torqueSign, "LhipFlex");	
	angleSign  =  1.0; torqueSign = -1.0;
	LeftKneeExtension_TorqueMuscle 	 = Anderson2007TorqueMuscle (knee, ext, gender, age, subjectHeight, subjectMass, jointAngleOffset, angleSign, torqueSign, "LkneeExt");
	angleSign  =  1.0; torqueSign =  1.0;
	LeftKneeFlexion_TorqueMuscle 	 = Anderson2007TorqueMuscle (knee, flex, gender, age, subjectHeight, subjectMass, jointAngleOffset, angleSign, torqueSign, "LkneeFlex");	
	angleSign  = -1.0; torqueSign =  1.0;
	LeftAnkleExtension_TorqueMuscle  = Anderson2007TorqueMuscle (ankle, ext, gender, age, subjectHeight, subjectMass, jointAngleOffset, angleSign, torqueSign, "LankleExt");
	angleSign  = -1.0; torqueSign = -1.0;
	LeftAnkleFlexion_TorqueMuscle 	 = Anderson2007TorqueMuscle (ankle, flex, gender, age, subjectHeight, subjectMass, jointAngleOffset, angleSign, torqueSign, "LankleFlex");	
	angleSign = -1.0; torqueSign = 1.0;
	TorsoExtension_TorqueMuscle	 = Anderson2007TorqueMuscle (hip, ext, gender, age, subjectHeight, subjectMass, jointAngleOffset, angleSign, torqueSign, "TorsoExt");
	angleSign  = -1.0; torqueSign = -1.0;
	TorsoFlexion_TorqueMuscle 	 = Anderson2007TorqueMuscle (hip, flex, gender, age, subjectHeight, subjectMass, jointAngleOffset, angleSign, torqueSign, "TorsoFlex");
	
	double passiveTorqueScale = 0.0;
	RightHipExtension_TorqueMuscle.setMaximumActiveIsometricTorque(jointStrength[ParamRightHip_Ext]);
	RightHipExtension_TorqueMuscle.setPassiveTorqueScale(passiveTorqueScale);
	RightHipFlexion_TorqueMuscle.setMaximumActiveIsometricTorque(jointStrength[ParamRightHip_Flex]);
	RightHipFlexion_TorqueMuscle.setPassiveTorqueScale(passiveTorqueScale);
	RightKneeExtension_TorqueMuscle.setMaximumActiveIsometricTorque(jointStrength[ParamRightKnee_Ext]);
	RightKneeExtension_TorqueMuscle.setPassiveTorqueScale(passiveTorqueScale);
	RightKneeFlexion_TorqueMuscle.setMaximumActiveIsometricTorque(jointStrength[ParamRightKnee_Flex]);
	RightKneeFlexion_TorqueMuscle.setPassiveTorqueScale(passiveTorqueScale);
	RightAnkleExtension_TorqueMuscle.setMaximumActiveIsometricTorque(jointStrength[ParamRightAnkle_Ext]);
	RightAnkleExtension_TorqueMuscle.setPassiveTorqueScale(passiveTorqueScale);
	RightAnkleFlexion_TorqueMuscle.setMaximumActiveIsometricTorque(jointStrength[ParamRightAnkle_Flex]);
	RightAnkleFlexion_TorqueMuscle.setPassiveTorqueScale(passiveTorqueScale);
	LeftHipExtension_TorqueMuscle.setMaximumActiveIsometricTorque(jointStrength[ParamLeftHip_Ext]);
	LeftHipExtension_TorqueMuscle.setPassiveTorqueScale(passiveTorqueScale);
	LeftHipFlexion_TorqueMuscle.setMaximumActiveIsometricTorque(jointStrength[ParamLeftHip_Flex]);
	LeftHipFlexion_TorqueMuscle.setPassiveTorqueScale(passiveTorqueScale);
	LeftKneeExtension_TorqueMuscle.setMaximumActiveIsometricTorque(jointStrength[ParamLeftKnee_Ext]);
	LeftKneeExtension_TorqueMuscle.setPassiveTorqueScale(passiveTorqueScale);
	LeftKneeFlexion_TorqueMuscle.setMaximumActiveIsometricTorque(jointStrength[ParamLeftKnee_Flex]);
	LeftKneeFlexion_TorqueMuscle.setPassiveTorqueScale(passiveTorqueScale);
	LeftAnkleExtension_TorqueMuscle.setMaximumActiveIsometricTorque(jointStrength[ParamLeftAnkle_Ext]);
	LeftAnkleExtension_TorqueMuscle.setPassiveTorqueScale(passiveTorqueScale);
	LeftAnkleFlexion_TorqueMuscle.setMaximumActiveIsometricTorque(jointStrength[ParamLeftAnkle_Flex]);
	LeftAnkleFlexion_TorqueMuscle.setPassiveTorqueScale(passiveTorqueScale);
	TorsoExtension_TorqueMuscle.setMaximumActiveIsometricTorque(jointStrength[ParamTorso_Ext]);
	TorsoExtension_TorqueMuscle.setPassiveTorqueScale(passiveTorqueScale);
	TorsoFlexion_TorqueMuscle.setMaximumActiveIsometricTorque(jointStrength[ParamTorso_Flex]);
	TorsoFlexion_TorqueMuscle.setPassiveTorqueScale(passiveTorqueScale);
}

bool OrthosisInfo::loadParamsFromFile (const char* filename) {
  LuaTable lua_table = LuaTable::fromFile (filename);
  
  angleCenterOffset 			= lua_table[orthosis_name.c_str()]["angleCenterOffset"];
  angleHardStopInPlantarFlexion 	= lua_table[orthosis_name.c_str()]["angleHardStopInPlantarFlexion"];
  angleHardStopInDorsiFlexion 		= lua_table[orthosis_name.c_str()]["angleHardStopInDorsiFlexion"];
  stiffnessEndRangeInDorsiFlexion 	= lua_table[orthosis_name.c_str()]["stiffnessEndRangeInDorsiFlexion"];
  stiffnessEndRangeInPlantarFlexion 	= lua_table[orthosis_name.c_str()]["stiffnessEndRangeInPlantarFlexion"];
  torquePreload 			= lua_table[orthosis_name.c_str()]["torquePreload"];
  anglePreloadWindow 			= lua_table[orthosis_name.c_str()]["anglePreloadWindow"];
  stiffnessInDorsiFlexion 		= lua_table[orthosis_name.c_str()]["stiffnessInDorsiFlexion"];
  stiffnessInPlantarFlexion		= lua_table[orthosis_name.c_str()]["stiffnessInPlantarFlexion"];
  jointOffsetAngle			= lua_table[orthosis_name.c_str()]["jointOffsetAngle"];
  jointAngleSign			= lua_table[orthosis_name.c_str()]["jointAngleSign"];
  jointTorqueSign			= lua_table[orthosis_name.c_str()]["jointTorqueSign"];
}

bool WalkerModel::loadOrthosisParams (const char* filename, bool verbose) {
  
  orthosisParams_r.orthosis_name = "orthosis_r";
  orthosisParams_r.loadParamsFromFile(filename);
  orthosisParams_l.orthosis_name = "orthosis_l";
  orthosisParams_l.loadParamsFromFile(filename);
     
  orthosis_r = RigidBodyDynamics::Addons::CustomForces::NeuroSwingTorsionSpring (
                              orthosisParams_r.angleHardStopInPlantarFlexion, 
                              orthosisParams_r.angleCenterOffset,             
			      orthosisParams_r.angleHardStopInDorsiFlexion,   
                              orthosisParams_r.stiffnessEndRangeInPlantarFlexion,     
                              orthosisParams_r.stiffnessInPlantarFlexion,     
                              orthosisParams_r.stiffnessInDorsiFlexion,       
			      orthosisParams_r.stiffnessEndRangeInDorsiFlexion,
                              orthosisParams_r.anglePreloadWindow,            
                              orthosisParams_r.torquePreload,                 
                              orthosisParams_r.jointOffsetAngle,              
                              orthosisParams_r.jointAngleSign,                
                              orthosisParams_r.jointTorqueSign,
                              orthosisParams_r.orthosis_name);
  
  orthosis_l = RigidBodyDynamics::Addons::CustomForces::NeuroSwingTorsionSpring (
                              orthosisParams_l.angleHardStopInPlantarFlexion, 
                              orthosisParams_l.angleCenterOffset,             
			      orthosisParams_l.angleHardStopInDorsiFlexion,   
                              orthosisParams_l.stiffnessEndRangeInPlantarFlexion,     
                              orthosisParams_l.stiffnessInPlantarFlexion,     
                              orthosisParams_l.stiffnessInDorsiFlexion,       
			      orthosisParams_l.stiffnessEndRangeInDorsiFlexion,
                              orthosisParams_l.anglePreloadWindow,            
                              orthosisParams_l.torquePreload,                 
                              orthosisParams_l.jointOffsetAngle,              
                              orthosisParams_l.jointAngleSign,                
                              orthosisParams_l.jointTorqueSign,
                              orthosisParams_l.orthosis_name);
  return true;
}

bool WalkerModel::loadPoints (const char* filename, bool verbose) {
	LuaTable lua_table = LuaTable::fromFile (filename);

	int point_count = lua_table["points"].length();

	for (int pi = 1; pi <= point_count; pi++) {
		Point point = lua_table["points"][pi];
		
		PointName point_name = getPointNameFromString (point.name.c_str());

		if (point_name == PointNameLast) {
		  continue;
		}

		point.body_id = model.GetBodyId (point.body_name.c_str());
		pointInfos[point_name] = point;

		if (verbose) {
		  cout << "Point '" << point.name << "' (PointName = " << point_name << ")" << endl;
		  cout << "  body        = " << point.body_name << " (id = " << point.body_id << ")" << endl;
		  cout << "  point_local = '" << point.point_local.transpose() << endl;
		}
	}

	// check whether we missed some points
	Point default_point;
	for (unsigned int i = 0; i < PointNameLast - 1; i++) {
		if (pointInfos[i].name == default_point.name) {
		  cerr << "Error: could not find point info for point '" << PointMap[i].name_str << "' in file " << filename << "." << endl;
		  abort();
		}
	}

	return true;
}

bool WalkerModel::loadConstraintSets (const char* filename, bool verbose) {
	// initialize Constraint Sets
	LuaTable lua_table = LuaTable::fromFile (filename);

	vector<LuaKey> constraint_set_keys = lua_table["constraint_sets"].keys();
	vector<string> constraint_set_names;
	for (size_t i = 0; i < constraint_set_keys.size(); i++) {
		if (constraint_set_keys[i].type == LuaKey::String)
			constraint_set_names.push_back (constraint_set_keys[i].string_value);
		else {
			cerr << "Found invalid constraint set name, string expected!" <<  endl;
			abort();
		}
	}

	constraintSets = vector<RigidBodyDynamics::ConstraintSet> (constraint_set_names.size());

	for (size_t si = 0; si < constraint_set_names.size(); si++) {
		string set_name_str = constraint_set_names[si];
		ConstraintSetName set_name = getConstraintSetNameFromString (set_name_str.c_str());

		if (set_name == ConstraintSetNameLast) {
			continue;
		}

		if (verbose) {
			cout << "ConstraintSet '" << set_name_str << "' (ConstraintSetName = " << set_name << ")" << endl;
		}

		unsigned int constraint_count = lua_table["constraint_sets"][constraint_set_names[si].c_str()].length();
		constraintSetInfos[set_name].name = set_name_str;
		constraintSetInfos[set_name].constraints.resize (constraint_count);

		for (int ci = 0; ci < constraint_count; ci++) {
			ConstraintInfo constraint_info = lua_table["constraint_sets"][set_name_str.c_str()][ci + 1];
			PointName point_name = getPointNameFromString (constraint_info.point_name.c_str());
			constraint_info.point_id = point_name;

			if (verbose) {
				cout << "  Adding Constraint point '" << pointInfos[point_name].name << "' normal = " << constraint_info.normal.transpose() << " body id = " << pointInfos[point_name].body_id << endl;
			}
			constraintSets[set_name].AddConstraint (pointInfos[point_name].body_id, pointInfos[point_name].point_local, constraint_info.normal);
			constraintSetInfos[set_name].constraints[ci] = constraint_info;
		}

		constraintSets[set_name].Bind (model);
	}

	// check whether we missed some sets
	ConstraintSetInfo default_constraint_set;
	for (unsigned int i = 0; i < ConstraintSetNameLast; i++) {
		if (constraintSetInfos[i].name != ConstraintSetMap[i].name_str) {
			cerr << "Error: could not find ConstraintSet info for set '" << ConstraintSetMap[i].name_str << "' in file " << filename << "." << endl;
			abort();
		}
	}

	return true;
}
