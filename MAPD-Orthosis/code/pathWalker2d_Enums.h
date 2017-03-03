#ifndef _PATHWALKER2D_ENUMS_H
#define _PATHWALKER2D_ENUMS_H

enum StageName {
	StageRightFlat = 0,
	StageRightHalx,
	StageRightHalxLeftHeel,
	StageRightHalxLeftFlat,
	StageLeftFlat,
	StageLeftHalx,
	StageLeftHalxRightHeel,
	StageLeftHalxRightFlat,
	StageNameLast
};

enum PosName {
	PosPelvisPosX = 0,
	PosPelvisPosZ,
	PosPelvisRotY,
	PosRightHipRotY,
	PosRightKneeRotY,
	PosRightAnkleRotY,
	PosLeftHipRotY,
	PosLeftKneeRotY,
	PosLeftAnkleRotY,
	PosTorsoRotY,
	PosNameLast
};

enum StateName {
	StatePelvisPosX = 0,
	StatePelvisPosZ,
	StatePelvisRotY,
	StateRightHipRotY,
	StateRightKneeRotY,
	StateRightAnkleRotY,
	StateLeftHipRotY,
	StateLeftKneeRotY,
	StateLeftAnkleRotY,
	StateTorsoRotY,
	StatePelvisVelX,
	StatePelvisVelZ,
	StatePelvisRotVelY,
	StateRightHipRotVelY,
	StateRightKneeRotVelY,
	StateRightAnkleRotVelY,
	StateLeftHipRotVelY,
	StateLeftKneeRotVelY,
	StateLeftAnkleRotVelY,
	StateTorsoRotVelY,	
	StateRightHipExtensionRotActY,
	StateRightHipFlexionRotActY,
	StateRightKneeExtensionRotActY,
	StateRightKneeFlexionRotActY,
	StateRightAnkleExtensionRotActY,
	StateRightAnkleFlexionRotActY,
	StateLeftHipExtensionRotActY,
	StateLeftHipFlexionRotActY,
	StateLeftKneeExtensionRotActY,
	StateLeftKneeFlexionRotActY,
	StateLeftAnkleExtensionRotActY,
	StateLeftAnkleFlexionRotActY,
	StateTorsoExtensionRotActY,
	StateTorsoFlexionRotActY,
	MayerStateActivationSquared,
	StateNameLast
};

enum ControlName {
	ControlRightHipExtensionRotY = 0,  
	ControlRightHipFlexionRotY,
	ControlRightKneeExtensionRotY,
	ControlRightKneeFlexionRotY,
	ControlRightAnkleExtensionRotY,
	ControlRightAnkleFlexionRotY,
	ControlLeftHipExtensionRotY,
	ControlLeftHipFlexionRotY,
	ControlLeftKneeExtensionRotY,
	ControlLeftKneeFlexionRotY,
	ControlLeftAnkleExtensionRotY,
	ControlLeftAnkleFlexionRotY,
	ControlTorsoExtensionRotY,
	ControlTorsoFlexionRotY,
	ControlNameLast
};

enum ParamName {
	ParamRightHip_Ext = 0,
	ParamRightHip_Flex,
	ParamRightKnee_Ext,
	ParamRightKnee_Flex,
	ParamRightAnkle_Ext,
	ParamRightAnkle_Flex,
	ParamLeftHip_Ext,
	ParamLeftHip_Flex,
	ParamLeftKnee_Ext,
	ParamLeftKnee_Flex,
	ParamLeftAnkle_Ext,
	ParamLeftAnkle_Flex,
	ParamTorso_Ext,
	ParamTorso_Flex,
	ParamRightOrthosisStiffnessInDorsiFlexion,  
	ParamRightOrthosisStiffnessInPlantarFlexion,
	ParamLeftOrthosisStiffnessInDorsiFlexion,  
	ParamLeftOrthosisStiffnessInPlantarFlexion,
	ParamNameLast
};

enum ConstraintSetName {
	CSRightFlat = 0,
	CSRightHalx,
	CSRightHalxLeftHeel,
	CSRightHalxLeftFlat,
	CSLeftFlat,
	CSLeftHalx,
	CSLeftHalxRightHeel,
	CSLeftHalxRightFlat,
	ConstraintSetNameLast
};

enum PointName {
	PointRightHeel,
	PointRightHalx,
	PointLeftHeel,
	PointLeftHalx,
	PointNameLast
};

#endif
