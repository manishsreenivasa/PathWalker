#include <cmath>
#include "def_usrmod.hpp"
#include <iostream>
#include <iomanip>
#include <limits>
#include <vector>
#include <fstream>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include "pathWalker2d_Model.h"
#include "SplineInterpolator.h"

using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;
using namespace std;

// datfileutils.h needs VectorNd and therefore RigidBodyDynamics::Math
#include "datfileutils.h"

int nmos   = -1;  /* Number of phases (MOdel Stages) */
int np     = -1;  /* Number of parameters */
int nrc    =  0;  /* Number of coupled constraints */
int nrce   =  0;  /* Number of coupled equality constraints */
int nxd    = -1;  /* Number of differential states */
int nxa    =  0;  /* Number of algebraic states */
int nu     = -1;  /* Number of controls */
int npr    =  0;  /* Number of local parameters */
int nlsq   =  0;  /* Number of least squares fits */

const unsigned int NoOfControls = ControlNameLast;
const unsigned int NoOfStates   = StateNameLast;
const unsigned int NoOfPos      = PosNameLast;
const unsigned int NoOfParams   = ParamNameLast;

WalkerModel walker;
const char* const model_file = "barefeet.lua";
char problem_name[255];
string datfile_name;
bool verbose = false;

// Static friction coefficient based on concrete and rubber contact
double coeff_friction = 0.8;

/* Model */
void LoadModelAndConstraints () {
	if (verbose)
	  cout << "[LoadModelAndConstraints] Loading Model: '" << model_file << "'" << endl;
	if (!walker.loadFromFile (model_file, datfile_name, verbose))
	  abort();

	if (verbose)
	  cout << "[LoadModelAndConstraints] Loading Points" <<  endl;
	if (!walker.loadPoints (model_file, verbose))
	  abort();
	
	if (verbose)
	  cout << "[LoadModelAndConstraints] Loading Constraint Sets" <<  endl;
	if (!walker.loadConstraintSets (model_file, verbose))
	  abort();

	if (verbose)
	  cout << "[LoadModelAndConstraints] Model and Constraint loading successful!" << endl;
}

// Interpolator - read from file - For motion capture data
SplineInterpolator<VectorNd> interpolator;

static void initialize_from_data (
  long   *imos,      /* index of model stage (I) */
  long   *imsn,      /* index of m.s. node on current model stage (I) */
  double *sd,        /* differential states at m.s. node (I/O) */
  double *sa,        /* algebraic states at m.s. node (I/O) */
  double *u,         /* controls at m.s. node (I/O) */
  double *udot,      /* control slopes at m.s. node (I/O) */
  double *ue,        /* controls at end of m.s. interval (I/O) */
  double *uedot,     /* control slopes at end of m.s. interval (I/O) */
  double *p,         /* global model parameters (I/O) */
  double *h,         /* model stage durations (I/O) */
  double *pr         /* local i.p.c. parameters (I/O) */
  )
{ 
  // compute the current time based on the information in the datfile
  VectorNd h_datfile = datfile_get_vector(datfile_name.c_str(), "h");
  VectorNd nshoot_datfile = datfile_get_vector(datfile_name.c_str(), "nshoot");

  double t = 0.;

  for (unsigned int i = 0; i < *imos; i++)
    t += h_datfile[i];

  t += static_cast<double>(*imsn) * h_datfile[*imos] / nshoot_datfile[*imos];
  
  if (sd) {
    VectorNd q = interpolator.getValues(t);
  
    for (size_t i = 0; i < walker.nDof; i++)
      sd[i] = q[i];
  }

  if (sd && u && p)
    walker.updateState (sd, u, p, CSRightFlat);

  VectorNd target_pos = interpolator.getValues (t);
}

// Objective function (Lagrangian type)
void lfcn_reg(double *t, double *xd, double *xa, double *u, double *p, double *lval, double *rwh, long *iwh, InfoPtr *info) {
  unsigned int i;
  *lval = 0.;

  double regFactor = 1e-3;
  for (i = 0; i < NoOfControls; i++) {
    *lval = u[i]*u[i];
  }
  *lval*=regFactor;
}

// Objective function (Least-Squares)
void lsqfcn_trackMotion(double *t, double *sd, double *sa, double *u,
			double *p, double *pr, double *res, long *dpnd, InfoPtr *info) {
  if (*dpnd) {
    *dpnd = RFCN_DPND(0, *sd, 0, 0, 0, 0);
    return;
  }
  VectorNd path_values = interpolator.getValues (*t);
  for (unsigned int i = 0; i < NoOfPos; i++) {
    res[i] = (path_values[i] - sd[i]);
  }
}

/* Right Hand Sides */
void ffcn_right_flat (double *t, double *xd, double *xa, double *u, double *p, double *rhs, double *rwh, long *iwh, InfoPtr *info) {
	walker.updateState (xd, u, p, CSRightFlat);
	walker.calcForwardDynamicsRhs (rhs);
}
void ffcn_right_halx (double *t, double *xd, double *xa, double *u, double *p, double *rhs, double *rwh, long *iwh, InfoPtr *info) {
	walker.updateState (xd, u, p, CSRightHalx);
	walker.calcForwardDynamicsRhs (rhs);
}
void ffcn_right_halx_left_heel (double *t, double *xd, double *xa, double *u, double *p, double *rhs, double *rwh, long *iwh, InfoPtr *info) {
	walker.updateState (xd, u, p, CSRightHalxLeftHeel);
	walker.calcForwardDynamicsRhs (rhs);
}
void ffcn_right_halx_left_flat (double *t, double *xd, double *xa, double *u, double *p, double *rhs, double *rwh, long *iwh, InfoPtr *info) {
	walker.updateState (xd, u, p, CSRightHalxLeftFlat);
	walker.calcForwardDynamicsRhs (rhs);
}
void ffcn_left_flat (double *t, double *xd, double *xa, double *u, double *p, double *rhs, double *rwh, long *iwh, InfoPtr *info) {
	walker.updateState (xd, u, p, CSLeftFlat);
	walker.calcForwardDynamicsRhs (rhs);
}
void ffcn_left_halx (double *t, double *xd, double *xa, double *u, double *p, double *rhs, double *rwh, long *iwh, InfoPtr *info) {
	walker.updateState (xd, u, p, CSLeftHalx);
	walker.calcForwardDynamicsRhs (rhs);
}
void ffcn_left_halx_right_heel (double *t, double *xd, double *xa, double *u, double *p, double *rhs, double *rwh, long *iwh, InfoPtr *info) {
	walker.updateState (xd, u, p, CSLeftHalxRightHeel);
	walker.calcForwardDynamicsRhs (rhs);
}
void ffcn_left_halx_right_flat (double *t, double *xd, double *xa, double *u, double *p, double *rhs, double *rwh, long *iwh, InfoPtr *info) {
	walker.updateState (xd, u, p, CSLeftHalxRightFlat);
	walker.calcForwardDynamicsRhs (rhs);
}

/* Decoupled Constraints */
static int rdfcn_right_flat_s_n = 10, rdfcn_right_flat_s_ne = 6;
void rdfcn_right_flat_s(double *ts, double *sd, double *sa, double *u, double *p, double *pr, double *res, long *dpnd, InfoPtr *info) {
  
  if (*dpnd) {
    *dpnd = RFCN_DPND(NULL, *sd, 0, *u, *p, 0);
    return;
  }

  walker.updateState (sd, u, p, CSRightFlat);

  Vector3d RightHeelPos = walker.getPointPosition (PointRightHeel);
  Vector3d RightHalxPos = walker.getPointPosition (PointRightHalx);
  Vector3d RightHeelVel = walker.getPointVelocity (PointRightHeel);
  Vector3d RightHalxVel = walker.getPointVelocity (PointRightHalx);
  Vector3d LeftHalxPos = walker.getPointPosition (PointLeftHalx);
  Vector3d RightHeelForce = walker.getPointForce (PointRightHeel);
  Vector3d RightHalxForce = walker.getPointForce (PointRightHalx);
  double limitingFriction = coeff_friction*(RightHeelForce[2]+RightHalxForce[2]);

  // Maintain Right Heel and Halx vertical positions on the ground
  res[0] = RightHeelPos[2];
  res[1] = RightHalxPos[2];
  
  // Maintain Right Heel velocities to zero
  res[2] = RightHeelVel[0];
  res[3] = RightHeelVel[2];
  
  // Maintain Right Halx anterior-posterior velocity to zero
  res[4] = RightHalxVel[2];  
  
  // Maintain Left Halx vertical position on the ground
  res[5] = LeftHalxPos[2];  
  
  // Maintain positive vertical force at Right Heel and Halx
  res[6] = RightHeelForce[2];
  res[7] = RightHalxForce[2];
  
  // Maintain anterior-posterior forces at Right Heel within friction cone
  res[8] = (limitingFriction-RightHeelForce[0]);
  res[9] = (limitingFriction+RightHeelForce[0]);
}

static int rdfcn_right_flat_i_n = 6, rdfcn_right_flat_i_ne = 0;
void rdfcn_right_flat_i(double *ts, double *sd, double *sa, double *u, double *p, double *pr, double *res, long *dpnd, InfoPtr *info) {
  if (*dpnd) {
    *dpnd = RFCN_DPND(0, *sd, 0, *u, *p, 0);
    return;
  }

  walker.updateState (sd, u, p, CSRightFlat);

  Vector3d LeftHalxPos = walker.getPointPosition (PointLeftHalx);
  Vector3d LeftHeelPos = walker.getPointPosition (PointLeftHeel);
  Vector3d RightHeelForce = walker.getPointForce (PointRightHeel);
  Vector3d RightHalxForce = walker.getPointForce (PointRightHalx);
  double limitingFriction = coeff_friction*(RightHeelForce[2]+RightHalxForce[2]);

  // Maintain Left Heel and Halx vertical positions off the ground
  res[0] = LeftHeelPos[2];
  res[1] = LeftHalxPos[2];
  
  // Maintain positive vertical force at Right Heel and Halx
  res[2] = RightHeelForce[2];
  res[3] = RightHalxForce[2];
  
  // Maintain anterior-posterior forces at Right Heel within friction cone
  res[4] = (limitingFriction-RightHeelForce[0]);
  res[5] = (limitingFriction+RightHeelForce[0]);  
}

static int rdfcn_right_halx_s_n = 7, rdfcn_right_halx_s_ne = 2;
void rdfcn_right_halx_s(double *ts, double *sd, double *sa, double *u, double *p, double *pr, double *res, long *dpnd, InfoPtr *info) {
  if (*dpnd) {
    *dpnd = RFCN_DPND(0, *sd, 0, *u, *p, 0);
    return;
  }

  walker.updateState (sd, u, p, CSRightFlat);

  Vector3d LeftHalxPos = walker.getPointPosition (PointLeftHalx);
  Vector3d LeftHeelPos = walker.getPointPosition (PointLeftHeel);
  Vector3d RightHeelPos = walker.getPointPosition (PointRightHeel);
  Vector3d RightHalxForce = walker.getPointForce (PointRightHalx);
  Vector3d RightHeelForce = walker.getPointForce (PointRightHeel);
  double limitingFriction = coeff_friction*(RightHeelForce[2]+RightHalxForce[2]);

  // Set zero force at Right Heel while on ground
  res[0] = RightHeelForce[2];
  res[1] = RightHeelPos[2];
  
  // Maintain Left Heel and Halx vertical positions off the ground
  res[2] = LeftHeelPos[2];
  res[3] = LeftHalxPos[2];
  
  // Maintain positive vertical force at Right Halx
  res[4] = RightHalxForce[2];
  
  // Maintain anterior-posterior forces at Right Halx within friction cone
  res[5] = (limitingFriction-RightHalxForce[0]);
  res[6] = (limitingFriction+RightHalxForce[0]);  
}

static int rdfcn_right_halx_i_n = 6, rdfcn_right_halx_i_ne = 0;
void rdfcn_right_halx_i(double *ts, double *sd, double *sa, double *u, double *p, double *pr, double *res, long *dpnd, InfoPtr *info) {
  if (*dpnd) {
    *dpnd = RFCN_DPND(0, *sd, 0, *u, *p, 0);
    return;
  }

  walker.updateState (sd, u, p, CSRightHalx);

  Vector3d LeftHalxPos = walker.getPointPosition (PointLeftHalx);
  Vector3d LeftHeelPos = walker.getPointPosition (PointLeftHeel);
  Vector3d RightHeelPos = walker.getPointPosition (PointRightHeel);
  Vector3d RightHalxForce = walker.getPointForce (PointRightHalx);
  double limitingFriction = coeff_friction*RightHalxForce[2];

  // Maintain Left Heel, Left Halx and Right Heel vertical positions off the ground
  res[0] = LeftHeelPos[2];
  res[1] = LeftHalxPos[2];
  res[2] = RightHeelPos[2];
  
  // Maintain positive vertical force at Right Halx
  res[3] = RightHalxForce[2];
  
  // Maintain anterior-posterior forces at Right Halx within friction cone
  res[4] = (limitingFriction-RightHalxForce[0]);
  res[5] = (limitingFriction+RightHalxForce[0]);  
}

static int rdfcn_right_halx_left_heel_s_n = 9, rdfcn_right_halx_left_heel_s_ne = 4;
void rdfcn_right_halx_left_heel_s(double *ts, double *sd, double *sa, double *u, double *p, double *pr, double *res, long *dpnd, InfoPtr *info) {
  if (*dpnd) {
    *dpnd = RFCN_DPND(0, *sd, 0, *u, *p, 0);
    return;
  }

  walker.updateState (sd, u, p, CSRightHalxLeftHeel);

  Vector3d LeftHalxPos = walker.getPointPosition (PointLeftHalx);
  Vector3d LeftHeelPos = walker.getPointPosition (PointLeftHeel);
  Vector3d LeftHeelVel = walker.getPointVelocity (PointLeftHeel);
  Vector3d RightHeelPos = walker.getPointPosition (PointRightHeel);
  Vector3d RightHalxForce = walker.getPointForce (PointRightHalx);
  Vector3d LeftHeelForce = walker.getPointForce (PointLeftHeel);
  double limitingFriction = coeff_friction*RightHalxForce[2];

  // Maintain Left Heel on the ground at zero velocity and zero forces
  res[0] = LeftHeelPos[2];
  res[1] = LeftHeelVel[2];  
  res[2] = LeftHeelForce[2];  
  res[3] = LeftHeelForce[0];
    
  // Maintain Left Halx and Right Heel off the ground
  res[4] = LeftHalxPos[2];
  res[5] = RightHeelPos[2];
  
  // Maintain positive vertical force at Right Halx
  res[6] = RightHalxForce[2]; 
  
  // Maintain anterior-posterior forces at Right Halx within friction cone
  res[7] = (limitingFriction-RightHalxForce[0]);
  res[8] = (limitingFriction+RightHalxForce[0]);
}

static int rdfcn_right_halx_left_heel_i_n = 8, rdfcn_right_halx_left_heel_i_ne = 0;
void rdfcn_right_halx_left_heel_i(double *ts, double *sd, double *sa, double *u, double *p, double *pr, double *res, long *dpnd, InfoPtr *info) {
  if (*dpnd) {
    *dpnd = RFCN_DPND(0, *sd, 0, *u, *p, 0);
    return;
  }

  walker.updateState (sd, u, p, CSRightHalxLeftHeel);
  Vector3d LeftHalxPos = walker.getPointPosition (PointLeftHalx);
  Vector3d RightHeelPos = walker.getPointPosition (PointRightHeel);
  Vector3d RightHalxForce = walker.getPointForce (PointRightHalx);
  Vector3d LeftHeelForce = walker.getPointForce (PointLeftHeel);
  double limitingFrictionRightHalx = coeff_friction*RightHalxForce[2];
  double limitingFrictionLeftHeel = coeff_friction*LeftHeelForce[2];

  // Maintain Left Halx and Right Heel off the ground
  res[0] = LeftHalxPos[2];
  res[1] = RightHeelPos[2];
  
  // Maintain positive vertical forces at Right Halx and Left Heel
  res[2] = RightHalxForce[2];
  res[3] = LeftHeelForce[2];
  
  // Maintain anterior-posterior forces at Right Halx within friction cone
  res[4] = (limitingFrictionRightHalx-RightHalxForce[0]);  
  res[5] = (limitingFrictionRightHalx+RightHalxForce[0]);

  // Maintain anterior-posterior forces at Left Heel within friction cone
  res[6] = (limitingFrictionLeftHeel-LeftHeelForce[0]);  
  res[7] = (limitingFrictionLeftHeel+LeftHeelForce[0]);    
}

static int rdfcn_right_halx_left_flat_s_n = 10, rdfcn_right_halx_left_flat_s_ne = 3;
void rdfcn_right_halx_left_flat_s(double *ts, double *sd, double *sa, double *u, double *p, double *pr, double *res, long *dpnd, InfoPtr *info) {
  if (*dpnd) {
    *dpnd = RFCN_DPND(0, *sd, 0, *u, *p, 0);
    return;
  }

  walker.updateState (sd, u, p, CSRightHalxLeftFlat);

  Vector3d LeftHalxPos = walker.getPointPosition (PointLeftHalx);
  Vector3d LeftHalxVel = walker.getPointVelocity (PointLeftHalx);
  Vector3d RightHeelPos = walker.getPointPosition (PointRightHeel);
  Vector3d RightHalxForce = walker.getPointForce (PointRightHalx);
  Vector3d LeftHeelForce = walker.getPointForce (PointLeftHeel);
  Vector3d LeftHalxForce = walker.getPointForce (PointLeftHalx);
  double limitingFrictionRightHalx = coeff_friction*RightHalxForce[2];
  double limitingFrictionLeft = coeff_friction*(LeftHeelForce[2]+LeftHalxForce[2]);  

  // Maintain Left Halx on ground, with zero vertical velocity, and zero vertical force 
  res[0] = LeftHalxPos[2];
  res[1] = LeftHalxVel[2];
  res[2] = LeftHalxForce[2];
  
  // Maintain Right Heel off ground
  res[3] = RightHeelPos[2];
  
  // Maintain positive vertical forces at Right Halx
  res[4] = RightHalxForce[2];
  
  // Maintain positive vertical forces at Left Heel
  res[5] = LeftHeelForce[2];
  
  // Maintain anterior-posterior forces at Right Halx within friction cone
  res[6] = (limitingFrictionRightHalx-RightHalxForce[0]);  
  res[7] = (limitingFrictionRightHalx+RightHalxForce[0]);
  
  // Maintain anterior-posterior forces at Left Heel within friction cone
  res[8] = (limitingFrictionLeft-LeftHeelForce[0]);  
  res[9] = (limitingFrictionLeft+LeftHeelForce[0]);  
}

static int rdfcn_right_halx_left_flat_i_n = 8, rdfcn_right_halx_left_flat_i_ne = 0;
void rdfcn_right_halx_left_flat_i(double *ts, double *sd, double *sa, double *u, double *p, double *pr, double *res, long *dpnd, InfoPtr *info) {
  if (*dpnd) {
    *dpnd = RFCN_DPND(0, *sd, 0, *u, *p, 0);
    return;
  }

  walker.updateState (sd, u, p, CSRightHalxLeftFlat);

  Vector3d RightHeelPos   = walker.getPointPosition (PointRightHeel);
  Vector3d RightHalxForce = walker.getPointForce (PointRightHalx);
  Vector3d LeftHeelForce  = walker.getPointForce (PointLeftHeel);
  Vector3d LeftHalxForce  = walker.getPointForce (PointLeftHalx);
  double limitingFrictionRightHalx = coeff_friction*RightHalxForce[2];
  double limitingFrictionLeft = coeff_friction*(LeftHeelForce[2]+LeftHalxForce[2]);  

  // Maintain Right Heel off the ground
  res[0] = RightHeelPos[2];
  
  // Maintain positive vertical forces at Right Halx
  res[1] = RightHalxForce[2];
  
  // Maintain positive vertical forces at Left Heel and Left Halx
  res[2] = LeftHeelForce[2];
  res[3] = LeftHalxForce[2];
  
  // Maintain anterior-posterior forces at Right Halx within friction cone
  res[4] = (limitingFrictionRightHalx-RightHalxForce[0]);  
  res[5] = (limitingFrictionRightHalx+RightHalxForce[0]);
  
  // Maintain anterior-posterior forces at Left Heel within friction cone
  res[6] = (limitingFrictionLeft-LeftHeelForce[0]);  
  res[7] = (limitingFrictionLeft+LeftHeelForce[0]);    
}

static int rdfcn_left_flat_s_n = 12, rdfcn_left_flat_s_ne = 8;
void rdfcn_left_flat_s(double *ts, double *sd, double *sa, double *u, double *p, double *pr, double *res, long *dpnd, InfoPtr *info) {
  
  if (*dpnd) {
    *dpnd = RFCN_DPND(NULL, *sd, 0, *u, *p, 0);
    return;
  }

  walker.updateState (sd, u, p, CSRightHalxLeftFlat);

  Vector3d LeftHeelPos = walker.getPointPosition (PointLeftHeel);
  Vector3d LeftHalxPos = walker.getPointPosition (PointLeftHalx);
  Vector3d LeftHeelVel = walker.getPointVelocity (PointLeftHeel);
  Vector3d LeftHalxVel = walker.getPointVelocity (PointLeftHalx);
  Vector3d RightHalxPos = walker.getPointPosition (PointRightHalx);
  Vector3d LeftHeelForce = walker.getPointForce (PointLeftHeel);
  Vector3d LeftHalxForce = walker.getPointForce (PointLeftHalx);
  Vector3d RightHalxForce = walker.getPointForce (PointRightHalx);
  double limitingFrictionLeft = coeff_friction*(LeftHeelForce[2]+LeftHalxForce[2]);  

  // Maintain Left Heel and Halx on the ground
  res[0] = LeftHeelPos[2];
  res[1] = LeftHalxPos[2];  
  
  // Maintain Left Heel and Halx velocities to zero
  res[2] = LeftHeelVel[0];
  res[3] = LeftHeelVel[2];  
  res[4] = LeftHalxVel[2];  
  
  // Maintain Right Halx on the ground with zero forces
  res[5] = RightHalxPos[2];  
  res[6] = RightHalxForce[2];
  res[7] = RightHalxForce[0];
  
  // Maintain positive vertical forces at Left Heel and Left Halx
  res[8] = LeftHeelForce[2];
  res[9] = LeftHalxForce[2];
  
  // Maintain anterior-posterior forces at Left Heel within friction cone
  res[10] = (limitingFrictionLeft-LeftHeelForce[0]);  
  res[11] = (limitingFrictionLeft+LeftHeelForce[0]);    
}

static int rdfcn_left_flat_i_n = 6, rdfcn_left_flat_i_ne = 0;
void rdfcn_left_flat_i(double *ts, double *sd, double *sa, double *u, double *p, double *pr, double *res, long *dpnd, InfoPtr *info) {
  if (*dpnd) {
    *dpnd = RFCN_DPND(0, *sd, 0, *u, *p, 0);
    return;
  }

  walker.updateState (sd, u, p, CSLeftFlat);

  Vector3d RightHalxPos  = walker.getPointPosition (PointRightHalx);
  Vector3d RightHeelPos  = walker.getPointPosition (PointRightHeel);
  Vector3d LeftHeelForce = walker.getPointForce (PointLeftHeel);
  Vector3d LeftHalxForce = walker.getPointForce (PointLeftHalx);
  double limitingFrictionLeft = coeff_friction*(LeftHeelForce[2]+LeftHalxForce[2]);  
    
  // Maintain Right Heel and Halx off the ground
  res[0] = RightHeelPos[2];
  res[1] = RightHalxPos[2];
  
  // Maintain positive vertical forces at Left Heel and Left Halx
  res[2] = LeftHeelForce[2];
  res[3] = LeftHalxForce[2];
  
  // Maintain anterior-posterior forces at Left Heel within friction cone
  res[4] = (limitingFrictionLeft-LeftHeelForce[0]);  
  res[5] = (limitingFrictionLeft+LeftHeelForce[0]);  
}

static int rdfcn_left_halx_s_n = 7, rdfcn_left_halx_s_ne = 2;
void rdfcn_left_halx_s(double *ts, double *sd, double *sa, double *u, double *p, double *pr, double *res, long *dpnd, InfoPtr *info) {
  if (*dpnd) {
    *dpnd = RFCN_DPND(0, *sd, 0, *u, *p, 0);
    return;
  }

  walker.updateState (sd, u, p, CSLeftFlat);

  Vector3d RightHalxPos  = walker.getPointPosition (PointRightHalx);
  Vector3d RightHeelPos  = walker.getPointPosition (PointRightHeel);
  Vector3d LeftHeelPos   = walker.getPointPosition (PointLeftHeel);
  Vector3d LeftHalxForce = walker.getPointForce (PointLeftHalx);
  Vector3d LeftHeelForce = walker.getPointForce (PointLeftHeel);
  double limitingFrictionLeft = coeff_friction*(LeftHeelForce[2]+LeftHalxForce[2]);  

  // Maintain Left Heel on the ground with zero vertical force
  res[0] = LeftHeelForce[2];
  res[1] = LeftHeelPos[2];
  
  // Maintain Right Heel and Halx off the ground
  res[2] = RightHeelPos[2];
  res[3] = RightHalxPos[2];
  
  // Maintain positive vertical forces at Left Halx
  res[4] = LeftHalxForce[2];
  
  // Maintain anterior-posterior forces at Left Halx within friction cone
  res[5] = (limitingFrictionLeft-LeftHalxForce[0]);  
  res[6] = (limitingFrictionLeft+LeftHalxForce[0]);
}

static int rdfcn_left_halx_i_n = 6, rdfcn_left_halx_i_ne = 0;
void rdfcn_left_halx_i(double *ts, double *sd, double *sa, double *u, double *p, double *pr, double *res, long *dpnd, InfoPtr *info) {
  if (*dpnd) {
    *dpnd = RFCN_DPND(0, *sd, 0, *u, *p, 0);
    return;
  }

  walker.updateState (sd, u, p, CSLeftHalx);

  Vector3d RightHalxPos  = walker.getPointPosition (PointRightHalx);
  Vector3d RightHeelPos  = walker.getPointPosition (PointRightHeel);
  Vector3d LeftHeelPos   = walker.getPointPosition (PointLeftHeel);
  Vector3d LeftHalxForce = walker.getPointForce (PointLeftHalx);
  double limitingFrictionLeft = coeff_friction*(LeftHalxForce[2]);  

  // Maintain Right Heel, Right Halx and Left Heel off the ground
  res[0] = RightHeelPos[2];
  res[1] = RightHalxPos[2];
  res[2] = LeftHeelPos[2];
  
  // Maintain positive vertical forces at Left Halx
  res[3] = LeftHalxForce[2];
  
  // Maintain anterior-posterior forces at Left Halx within friction cone
  res[4] = (limitingFrictionLeft-LeftHalxForce[0]);  
  res[5] = (limitingFrictionLeft+LeftHalxForce[0]);
}

static int rdfcn_left_halx_right_heel_s_n = 9, rdfcn_left_halx_right_heel_s_ne = 4;
void rdfcn_left_halx_right_heel_s (double *ts, double *sd, double *sa, double *u, double *p, double *pr, double *res, long *dpnd, InfoPtr *info) {
  if (*dpnd) {
    *dpnd = RFCN_DPND(0, *sd, 0, *u, *p, 0);
    return;
  }

  walker.updateState (sd, u, p, CSLeftHalxRightHeel);  

  Vector3d RightHalxPos  = walker.getPointPosition (PointRightHalx);
  Vector3d RightHeelPos  = walker.getPointPosition (PointRightHeel);
  Vector3d RightHeelVel  = walker.getPointVelocity (PointRightHeel);
  Vector3d LeftHeelPos   = walker.getPointPosition (PointLeftHeel);
  Vector3d LeftHalxForce = walker.getPointForce (PointLeftHalx);
  Vector3d RightHeelForce = walker.getPointForce (PointRightHeel);
  double limitingFrictionLeft = coeff_friction*(LeftHalxForce[2]);
  
  // Maintain Right Heel on ground at zero vertical velocity and zero forces
  res[0] = RightHeelPos[2];
  res[1] = RightHeelVel[2];
  res[2] = RightHeelForce[2];
  res[3] = RightHeelForce[0];
  
  // Maintain Right Halx and Left Heel off the ground
  res[4] = RightHalxPos[2];
  res[5] = LeftHeelPos[2];
  
  // Maintain positive vertical forces at Left Halx
  res[6] = LeftHalxForce[2];
  
  // Maintain anterior-posterior forces at Left Halx within friction cone
  res[7] = (limitingFrictionLeft-LeftHalxForce[0]);  
  res[8] = (limitingFrictionLeft+LeftHalxForce[0]);
}

static int rdfcn_left_halx_right_heel_i_n = 8, rdfcn_left_halx_right_heel_i_ne = 0;
void rdfcn_left_halx_right_heel_i(double *ts, double *sd, double *sa, double *u, double *p, double *pr, double *res, long *dpnd, InfoPtr *info) {
  if (*dpnd) {
    *dpnd = RFCN_DPND(0, *sd, 0, *u, *p, 0);
    return;
  }

  walker.updateState (sd, u, p, CSLeftHalxRightHeel);  
  Vector3d RightHalxPos   = walker.getPointPosition (PointRightHalx);
  Vector3d LeftHeelPos    = walker.getPointPosition (PointLeftHeel);
  Vector3d LeftHalxForce  = walker.getPointForce (PointLeftHalx);
  Vector3d RightHeelForce = walker.getPointForce (PointRightHeel);
  double limitingFrictionLeft = coeff_friction*(LeftHalxForce[2]);
  double limitingFrictionRight = coeff_friction*(RightHeelForce[2]);    

  // Maintain Right Halx and Left Heel off the ground
  res[0] = RightHalxPos[2];
  res[1] = LeftHeelPos[2];
  
  // Maintain positive vertical forces at Left Halx and Right Heel
  res[2] = LeftHalxForce[2];
  res[3] = RightHeelForce[2];
  
  // Maintain anterior-posterior forces at Left Halx within friction cone
  res[4] = (limitingFrictionLeft-LeftHalxForce[0]);  
  res[5] = (limitingFrictionLeft+LeftHalxForce[0]);
  
  // Maintain anterior-posterior forces at Right Heel within friction cone
  res[6] = (limitingFrictionRight-RightHeelForce[0]);  
  res[7] = (limitingFrictionRight+RightHeelForce[0]);  
}

static int rdfcn_left_halx_right_flat_s_n = 10, rdfcn_left_halx_right_flat_s_ne = 3;
void rdfcn_left_halx_right_flat_s (double *ts, double *sd, double *sa, double *u, double *p, double *pr, double *res, long *dpnd, InfoPtr *info) {
  if (*dpnd) {
    *dpnd = RFCN_DPND(0, *sd, 0, *u, *p, 0);
    return;
  }

  walker.updateState (sd, u, p, CSLeftHalxRightFlat);

  Vector3d RightHalxPos   = walker.getPointPosition (PointRightHalx);
  Vector3d RightHalxVel   = walker.getPointVelocity (PointRightHalx);
  Vector3d LeftHeelPos    = walker.getPointPosition (PointLeftHeel);
  Vector3d LeftHalxForce  = walker.getPointForce (PointLeftHalx);
  Vector3d RightHeelForce = walker.getPointForce (PointRightHeel);
  Vector3d RightHalxForce = walker.getPointForce (PointRightHalx);
  double limitingFrictionLeft = coeff_friction*(LeftHalxForce[2]);
  double limitingFrictionRight = coeff_friction*(RightHeelForce[2]+RightHalxForce[2]);    

  // Maintain Right Halx on ground at zero velocity and zero vertical force
  res[0] = RightHalxPos[2];
  res[1] = RightHalxVel[2];
  res[2] = RightHalxForce[2];
  
  // Maintain Left Heel off the ground
  res[3] = LeftHeelPos[2];
  
  // Maintain positive vertical forces at Left Halx and Right Heel
  res[4] = LeftHalxForce[2];
  res[5] = RightHeelForce[2];
  
  // Maintain anterior-posterior forces at Left Halx within friction cone
  res[6] = (limitingFrictionLeft-LeftHalxForce[0]);  
  res[7] = (limitingFrictionLeft+LeftHalxForce[0]);
  
  // Maintain anterior-posterior forces at Right Heel within friction cone
  res[8] = (limitingFrictionRight-RightHeelForce[0]);  
  res[9] = (limitingFrictionRight+RightHeelForce[0]);    
}

static int rdfcn_left_halx_right_flat_i_n = 8, rdfcn_left_halx_right_flat_i_ne = 0;
void rdfcn_left_halx_right_flat_i (double *ts, double *sd, double *sa, double *u, double *p, double *pr, double *res, long *dpnd, InfoPtr *info) {
  if (*dpnd) {
    *dpnd = RFCN_DPND(0, *sd, 0, *u, *p, 0);
    return;
  }

  walker.updateState (sd, u, p, CSLeftHalxRightFlat);

  Vector3d LeftHeelPos    = walker.getPointPosition (PointLeftHeel);
  Vector3d LeftHalxForce  = walker.getPointForce (PointLeftHalx);
  Vector3d RightHeelForce = walker.getPointForce (PointRightHeel);
  Vector3d RightHalxForce = walker.getPointForce (PointRightHalx);
  double limitingFrictionLeft = coeff_friction*(LeftHalxForce[2]);
  double limitingFrictionRight = coeff_friction*(RightHeelForce[2]+RightHalxForce[2]);    

  // Maintain Left Heel off the ground
  res[0] = LeftHeelPos[2];
  
  // Maintain positive vertical forces at Left Halx, Right Heel and Right Halx 
  res[1] = LeftHalxForce[2];
  res[2] = RightHeelForce[2];
  res[3] = RightHalxForce[2];
  
  // Maintain anterior-posterior forces at Left Halx within friction cone
  res[4] = (limitingFrictionLeft-LeftHalxForce[0]);  
  res[5] = (limitingFrictionLeft+LeftHalxForce[0]);
  
  // Maintain anterior-posterior forces at Right Heel within friction cone
  res[6] = (limitingFrictionRight-RightHeelForce[0]);  
  res[7] = (limitingFrictionRight+RightHeelForce[0]);    
}

static int rdfcn_left_halx_right_flat_e_n = 6, rdfcn_left_halx_right_flat_e_ne = 1;
void rdfcn_left_halx_right_flat_e(double *ts, double *sd, double *sa, double *u, double *p, double *pr, double *res, long *dpnd, InfoPtr *info) {
  if (*dpnd) {
    *dpnd = RFCN_DPND(*ts, *sd, 0, *u, *p, 0);
    return;
  }

  walker.updateState (sd, u, p, CSRightFlat);

  Vector3d LeftHeelPos = walker.getPointPosition (PointLeftHeel);
  Vector3d LeftHalxPos = walker.getPointPosition (PointLeftHalx);
  Vector3d RightHeelForce = walker.getPointForce (PointRightHeel);
  Vector3d RightHalxForce = walker.getPointForce (PointRightHalx);
  double limitingFrictionRight = coeff_friction*(RightHeelForce[2]+RightHalxForce[2]);    

  // Maintain Left Halx on the ground
  res[0] = LeftHalxPos[2];
  
  // Maintain Left Heel off the ground
  res[1] = LeftHeelPos[2];
  
  // Maintain positive vertical forces at Right Heel and Right Halx 
  res[2] = RightHeelForce[2];
  res[3] = RightHalxForce[2];
  
  // Maintain anterior-posterior forces at Right Heel within friction cone
  res[4] = (limitingFrictionRight-RightHeelForce[0]);  
  res[5] = (limitingFrictionRight+RightHeelForce[0]);    
}

// Hi-Res Data Output
/* *************************************************************************** */
vector<double> t_values;
vector<VectorNd> sd_values;
vector<VectorNd> u_values;
vector<VectorNd> active_torque_values;
vector<VectorNd> passive_torque_values;
vector<VectorNd> feet_pos_values;
vector<VectorNd> feet_force_values;
vector<int> currStage_value;

void hires_data_out(double *t, double *sd, double *sa, double *u, double *p, double *rwh, long *iwh, InfoPtr *info) {

     if (*t == 0.) {
	
	ofstream meshup_csv_stream;
	ofstream meshup_header_stream;
        string meshup_file_name = "RES/pathWalker2d.csv";
	
	string augmented_data_file_name = "RES/pathWalker2d_augmented.txt";
	ofstream augmented_data_stream;
	
	meshup_header_stream.open (meshup_file_name.c_str(), ios_base::trunc);
	augmented_data_stream.open (augmented_data_file_name.c_str(), ios_base::trunc);
	
	if (!meshup_header_stream || !augmented_data_stream) {
	  cerr << "Error opening file " << meshup_file_name << " or " << augmented_data_file_name << endl;
	  abort();
	}
	
	const char* meshup_header = "COLUMNS: \n\
	time,\n\
	pelvis:T:X,\n\
	pelvis:T:Z,\n\
	pelvis:R:Y:rad,\n\
	thigh_r:R:Y:rad,\n\
	shank_r:R:Y:rad,\n\
	foot_r:R:Y:rad,\n\
	thigh_l:R:Y:rad,\n\
	shank_l:R:Y:rad,\n\
	foot_l:R:Y:rad,\n\
	trunk:R:Y:rad\n\
	DATA:";

	meshup_header_stream << meshup_header << endl;
	
	if (t_values.size() > 0) {
	  for (unsigned int i = 0; i < t_values.size(); i ++) {
	    augmented_data_stream << t_values[i] << ", ";
	    meshup_header_stream << t_values[i] << ", ";
	    
	    // Write q values to CSV
	    for (unsigned int j = 0; j < NoOfPos; j++) {
	      meshup_header_stream << sd_values[i][j];
	      if (j < sd_values[i].size() -1 )
		meshup_header_stream << ", ";
	    }
	    
	    // Write rest of states (q_vel and activations)to Augmented data file
	    for (unsigned int j = NoOfPos; j < sd_values[i].size(); j++) {
	      augmented_data_stream << sd_values[i][j];
	      if (j < sd_values[i].size() -1 )
		augmented_data_stream << ", ";
	    }
	    augmented_data_stream << ", ";
	    
	    // Write controls (excitations) to Augmented data file
	    for (unsigned int j = 0; j < NoOfControls; j++) {
	      augmented_data_stream << u_values[i][j];
	      if (j < u_values[i].size() -1 )
		augmented_data_stream << ", ";
	    }
	    augmented_data_stream << ", ";
	    
	    // Write muscle passive and active torques to Augmented data file
	    for (unsigned int j = 0; j < NoOfControls; j++) {
	      augmented_data_stream << active_torque_values[i][j];
	      if (j < active_torque_values[i].size() -1 )
		augmented_data_stream << ", ";
	    }
	    augmented_data_stream << ", ";
	    for (unsigned int j = 0; j < NoOfControls; j++) {
	      augmented_data_stream << passive_torque_values[i][j];
	      if (j < passive_torque_values[i].size() -1 )
		augmented_data_stream << ", ";
	    }
	    augmented_data_stream << ", ";
	    
	    // Write feet heel and halx positions
	    for (unsigned int j = 0; j < 12; j++) {
	      augmented_data_stream << feet_pos_values[i][j];
	      if (j < feet_pos_values[i].size() -1 )
		augmented_data_stream << ", ";
	    }
	    augmented_data_stream << ", ";
	    
	    // Write feet force values
	    for (unsigned int j = 0; j < 12; j++) {
	      augmented_data_stream << feet_force_values[i][j];
	      if (j < feet_force_values[i].size() -1 )
		augmented_data_stream << ", ";
	    }
	    
	    augmented_data_stream << ", " << currStage_value[i];
	    
	    augmented_data_stream << endl;
	    meshup_header_stream << endl;
	  }
	}
	
	t_values.clear();
	sd_values.clear();
	u_values.clear();
	active_torque_values.clear();
	passive_torque_values.clear();
	feet_pos_values.clear();
	feet_force_values.clear();
	currStage_value.clear();
	
	meshup_header_stream.close();
	augmented_data_stream.close();
    }  
	  
    t_values.push_back (*t);

    VectorNd sd_vec (NoOfStates);
    for (unsigned i = 0; i < NoOfStates; i++)
      sd_vec[i] = sd[i];
    sd_values.push_back (sd_vec);
    
    VectorNd u_vec (NoOfControls);
    for (unsigned i = 0; i < NoOfControls; i++)
      u_vec[i] = u[i];
    u_values.push_back (u_vec);
    
    // Recompute and write Muscle Total Forces and Muscle Passive Forces
    VectorNd active_torque_vec (NoOfControls);
    active_torque_vec[ControlRightHipExtensionRotY] 	= walker.RightHipExtension_TorqueMuscle.calcJointTorque   (sd[StateRightHipRotY],   sd[StateRightHipRotVelY],   sd[StateRightHipExtensionRotActY]);
    active_torque_vec[ControlRightHipFlexionRotY]   	= walker.RightHipFlexion_TorqueMuscle.calcJointTorque     (sd[StateRightHipRotY],   sd[StateRightHipRotVelY],   sd[StateRightHipFlexionRotActY]);
    active_torque_vec[ControlRightKneeExtensionRotY] 	= walker.RightKneeExtension_TorqueMuscle.calcJointTorque  (sd[StateRightKneeRotY],  sd[StateRightKneeRotVelY],  sd[StateRightKneeExtensionRotActY]);
    active_torque_vec[ControlRightKneeFlexionRotY] 	= walker.RightKneeFlexion_TorqueMuscle.calcJointTorque    (sd[StateRightKneeRotY],  sd[StateRightKneeRotVelY],  sd[StateRightKneeFlexionRotActY]);
    active_torque_vec[ControlRightAnkleExtensionRotY] 	= walker.RightAnkleExtension_TorqueMuscle.calcJointTorque (sd[StateRightAnkleRotY], sd[StateRightAnkleRotVelY], sd[StateRightAnkleExtensionRotActY]);
    active_torque_vec[ControlRightAnkleFlexionRotY] 	= walker.RightAnkleFlexion_TorqueMuscle.calcJointTorque   (sd[StateRightAnkleRotY], sd[StateRightAnkleRotVelY], sd[StateRightAnkleFlexionRotActY]);
    active_torque_vec[ControlLeftHipExtensionRotY] 	= walker.LeftHipExtension_TorqueMuscle.calcJointTorque    (sd[StateLeftHipRotY],    sd[StateLeftHipRotVelY],    sd[StateLeftHipExtensionRotActY]);
    active_torque_vec[ControlLeftHipFlexionRotY]   	= walker.LeftHipFlexion_TorqueMuscle.calcJointTorque      (sd[StateLeftHipRotY],    sd[StateLeftHipRotVelY],    sd[StateLeftHipFlexionRotActY]);
    active_torque_vec[ControlLeftKneeExtensionRotY] 	= walker.LeftKneeExtension_TorqueMuscle.calcJointTorque   (sd[StateLeftKneeRotY],   sd[StateLeftKneeRotVelY],   sd[StateLeftKneeExtensionRotActY]);
    active_torque_vec[ControlLeftKneeFlexionRotY] 	= walker.LeftKneeFlexion_TorqueMuscle.calcJointTorque     (sd[StateLeftKneeRotY],   sd[StateLeftKneeRotVelY],   sd[StateLeftKneeFlexionRotActY]);
    active_torque_vec[ControlLeftAnkleExtensionRotY] 	= walker.LeftAnkleExtension_TorqueMuscle.calcJointTorque  (sd[StateLeftAnkleRotY],  sd[StateLeftAnkleRotVelY],  sd[StateLeftAnkleExtensionRotActY]);
    active_torque_vec[ControlLeftAnkleFlexionRotY] 	= walker.LeftAnkleFlexion_TorqueMuscle.calcJointTorque    (sd[StateLeftAnkleRotY],  sd[StateLeftAnkleRotVelY],  sd[StateLeftAnkleFlexionRotActY]);
    active_torque_vec[ControlTorsoExtensionRotY] 	= walker.TorsoExtension_TorqueMuscle.calcJointTorque 	  (sd[StateTorsoRotY],      sd[StateTorsoRotVelY], 	sd[StateTorsoExtensionRotActY]);
    active_torque_vec[ControlTorsoFlexionRotY] 		= walker.TorsoFlexion_TorqueMuscle.calcJointTorque   	  (sd[StateTorsoRotY],      sd[StateTorsoRotVelY], 	sd[StateTorsoFlexionRotActY]);
    active_torque_values.push_back(active_torque_vec);
    
    VectorNd passive_torque_vec (NoOfControls);
    passive_torque_vec[ControlRightHipExtensionRotY] 	= walker.RightHipExtension_TorqueMuscle.calcJointTorque   (sd[StateRightHipRotY],   sd[StateRightHipRotVelY],   0.0);
    passive_torque_vec[ControlRightHipFlexionRotY]   	= walker.RightHipFlexion_TorqueMuscle.calcJointTorque     (sd[StateRightHipRotY],   sd[StateRightHipRotVelY],   0.0);
    passive_torque_vec[ControlRightKneeExtensionRotY] 	= walker.RightKneeExtension_TorqueMuscle.calcJointTorque  (sd[StateRightKneeRotY],  sd[StateRightKneeRotVelY],  0.0);
    passive_torque_vec[ControlRightKneeFlexionRotY] 	= walker.RightKneeFlexion_TorqueMuscle.calcJointTorque    (sd[StateRightKneeRotY],  sd[StateRightKneeRotVelY],  0.0);
    passive_torque_vec[ControlRightAnkleExtensionRotY] 	= walker.RightAnkleExtension_TorqueMuscle.calcJointTorque (sd[StateRightAnkleRotY], sd[StateRightAnkleRotVelY], 0.0);
    passive_torque_vec[ControlRightAnkleFlexionRotY] 	= walker.RightAnkleFlexion_TorqueMuscle.calcJointTorque   (sd[StateRightAnkleRotY], sd[StateRightAnkleRotVelY], 0.0);
    passive_torque_vec[ControlLeftHipExtensionRotY] 	= walker.LeftHipExtension_TorqueMuscle.calcJointTorque    (sd[StateLeftHipRotY],    sd[StateLeftHipRotVelY],    0.0);
    passive_torque_vec[ControlLeftHipFlexionRotY]   	= walker.LeftHipFlexion_TorqueMuscle.calcJointTorque      (sd[StateLeftHipRotY],    sd[StateLeftHipRotVelY],    0.0);
    passive_torque_vec[ControlLeftKneeExtensionRotY] 	= walker.LeftKneeExtension_TorqueMuscle.calcJointTorque   (sd[StateLeftKneeRotY],   sd[StateLeftKneeRotVelY],   0.0);
    passive_torque_vec[ControlLeftKneeFlexionRotY] 	= walker.LeftKneeFlexion_TorqueMuscle.calcJointTorque     (sd[StateLeftKneeRotY],   sd[StateLeftKneeRotVelY],   0.0);
    passive_torque_vec[ControlLeftAnkleExtensionRotY] 	= walker.LeftAnkleExtension_TorqueMuscle.calcJointTorque  (sd[StateLeftAnkleRotY],  sd[StateLeftAnkleRotVelY],  0.0);
    passive_torque_vec[ControlLeftAnkleFlexionRotY] 	= walker.LeftAnkleFlexion_TorqueMuscle.calcJointTorque    (sd[StateLeftAnkleRotY],  sd[StateLeftAnkleRotVelY],  0.0);
    passive_torque_vec[ControlTorsoExtensionRotY] 	= walker.TorsoExtension_TorqueMuscle.calcJointTorque 	  (sd[StateTorsoRotY],      sd[StateTorsoRotVelY], 	0.0);
    passive_torque_vec[ControlTorsoFlexionRotY] 	= walker.TorsoFlexion_TorqueMuscle.calcJointTorque   	  (sd[StateTorsoRotY],      sd[StateTorsoRotVelY], 	0.0);
    passive_torque_values.push_back(passive_torque_vec);
    
    // Get foot heel and halx points
    RigidBodyDynamics::Math::VectorNd q_curr;
    q_curr = VectorNd::Zero (walker.model.dof_count);
    for (unsigned int j = 0; j < NoOfPos; j++)
      q_curr[j] = sd_vec[j];

    // 6 = foot_r, 9 = foot_l
    RigidBodyDynamics::Math::Vector3d RightHeelPos = RigidBodyDynamics::CalcBodyToBaseCoordinates (walker.model, q_curr, 6, walker.pointInfos[PointRightHeel].point_local, true);
    RigidBodyDynamics::Math::Vector3d RightHalxPos = RigidBodyDynamics::CalcBodyToBaseCoordinates (walker.model, q_curr, 6, walker.pointInfos[PointRightHalx].point_local, true);
    RigidBodyDynamics::Math::Vector3d LeftHeelPos = RigidBodyDynamics::CalcBodyToBaseCoordinates (walker.model, q_curr, 9, walker.pointInfos[PointLeftHeel].point_local, true);
    RigidBodyDynamics::Math::Vector3d LeftHalxPos = RigidBodyDynamics::CalcBodyToBaseCoordinates (walker.model, q_curr, 9, walker.pointInfos[PointLeftHalx].point_local, true);
	    
    VectorNd feet_pos_vec (12);
    feet_pos_vec[0] = RightHeelPos[0]; feet_pos_vec[1]  = RightHeelPos[1]; feet_pos_vec[2]  = RightHeelPos[2];
    feet_pos_vec[3] = RightHalxPos[0]; feet_pos_vec[4]  = RightHalxPos[1]; feet_pos_vec[5]  = RightHalxPos[2];
    feet_pos_vec[6] = LeftHeelPos[0];  feet_pos_vec[7]  = LeftHeelPos[1];  feet_pos_vec[8]  = LeftHeelPos[2];
    feet_pos_vec[9] = LeftHalxPos[0];  feet_pos_vec[10] = LeftHalxPos[1];  feet_pos_vec[11] = LeftHalxPos[2];
    feet_pos_values.push_back(feet_pos_vec);
    
    const long currStage = info->cimos;
    currStage_value.push_back(currStage);
    VectorNd feet_force_vec (12);        
    Vector3d RightHeelForce(0.0,0.0,0.0);
    Vector3d RightHalxForce(0.0,0.0,0.0);
    Vector3d LeftHeelForce(0.0,0.0,0.0);
    Vector3d LeftHalxForce(0.0,0.0,0.0);

    switch (currStage) {
	  case 0: 		// CSRightFlat
	    walker.updateState (sd, u, p, CSRightFlat);
	    RightHeelForce = walker.getPointForce (PointRightHeel);
	    RightHalxForce = walker.getPointForce (PointRightHalx);
	    break;
	  case 1:	    	// CSRightHalx
	    walker.updateState (sd, u, p, CSRightHalx);
	    RightHalxForce = walker.getPointForce (PointRightHalx);
	    break;
	  case 2:	    	// CSRightHalxLeftHeel
	    walker.updateState (sd, u, p, CSRightHalxLeftHeel);
	    RightHalxForce = walker.getPointForce (PointRightHalx);
	    LeftHeelForce = walker.getPointForce (PointLeftHeel);
	    break;
	   case 3:	    	// CSRightHalxLeftFlat
	    walker.updateState (sd, u, p, CSRightHalxLeftFlat);
	    RightHalxForce = walker.getPointForce (PointRightHalx);
	    LeftHeelForce = walker.getPointForce (PointLeftHeel);
	    LeftHalxForce = walker.getPointForce (PointLeftHalx);
	    break;
	   case 4:	    	// CSLeftFlat
	    walker.updateState (sd, u, p, CSLeftFlat);
	    LeftHeelForce = walker.getPointForce (PointLeftHeel);
	    LeftHalxForce = walker.getPointForce (PointLeftHalx);
	    break;
	   case 5: 	    	// CSLeftHalx
	    walker.updateState (sd, u, p, CSLeftHalx);
	    LeftHalxForce = walker.getPointForce (PointLeftHalx);
	    break;
	   case 6:		// CSLeftHalxRightHeel
	    walker.updateState (sd, u, p, CSLeftHalxRightHeel);
	    LeftHalxForce = walker.getPointForce (PointLeftHalx);
	    RightHeelForce = walker.getPointForce (PointRightHeel);
	    break;
	   case 7:	    	// CSLeftHalxRightFlat
	    walker.updateState (sd, u, p, CSLeftHalxRightFlat);
	    LeftHalxForce = walker.getPointForce (PointLeftHalx);
	    RightHeelForce = walker.getPointForce (PointRightHeel);
	    RightHalxForce = walker.getPointForce (PointRightHalx);
	    break;
	}
	
	feet_force_vec[0] = RightHeelForce[0]; feet_force_vec[1]  = RightHeelForce[1]; feet_force_vec[2]  = RightHeelForce[2];
	feet_force_vec[3] = RightHalxForce[0]; feet_force_vec[4]  = RightHalxForce[1]; feet_force_vec[5]  = RightHalxForce[2];
	feet_force_vec[6] = LeftHeelForce[0];  feet_force_vec[7]  = LeftHeelForce[1];  feet_force_vec[8]  = LeftHeelForce[2];
	feet_force_vec[9] = LeftHalxForce[0];  feet_force_vec[10] = LeftHalxForce[1];  feet_force_vec[11] = LeftHalxForce[2];
	feet_force_values.push_back(feet_force_vec);
}

void write_ocp_output (long *imos, long *imsn, double *ts, double *te, double *sd, double *sa, double *u, double *udot, double *ue, double *uedot, double *p, double *pr, double *ccxd, double *mul_ccxd,
  #if defined(PRSQP) || defined(EXTPRSQP)
    double *ares,
    double *mul_ares,
  #endif
  double *rd, double *mul_rd, double *rc, double *mul_rc, double *obj, double *rwh, long *iwh) {
  InfoPtr info (0, *imos, *imsn);
  hires_data_out( ts, sd, sa, u, p, rwh, iwh, &info);
}

// \brief Entry point for the muscod application
extern "C" void def_model(void);
void def_model(void) {

	datfile_name = string("DAT/") + string("pathWalker2d") + string(".dat");
	if (verbose)
	  cout << "Using Datfile: " << datfile_name << endl;
  
	LoadModelAndConstraints();
	
	// For read-from-file
	string interpolator_filename = datfile_get_string(datfile_name.c_str(), "interpolator_data");
	if (verbose)
	  cout << "Loading interpolator from " << interpolator_filename << endl;
	interpolator.generateFromCSV (interpolator_filename.c_str());
	
	nmos = StageNameLast;
	np = ParamNameLast;
	nrc = 0;
	nrce = 0;
	nxd = StateNameLast;
	nxa = 0;
	nu = NoOfControls;
	nlsq = NoOfPos;	       
		
	// Define OCP dimensions
	def_mdims(nmos, np, nrc, nrce);

	// Right flat
	def_mstage(
			0,
			nxd, nxa, nu,
			NULL, lfcn_reg,
			0, 0, 0, NULL, ffcn_right_flat, NULL,
			NULL, NULL
			);
	// Right Halx
	def_mstage(
			1,
			nxd, nxa, nu,
			NULL, lfcn_reg,
			0, 0, 0, NULL, ffcn_right_halx, NULL,
			NULL, NULL
			);
	
	// Right Halx Left Heel
	def_mstage(
			2,
			nxd, nxa, nu,
			NULL, lfcn_reg,
			0, 0, 0, NULL, ffcn_right_halx_left_heel, NULL,
			NULL, NULL
			);
		
	// Right Halx Left Flat
	def_mstage(
			3,
			nxd, nxa, nu,
			NULL, lfcn_reg,
			0, 0, 0, NULL, ffcn_right_halx_left_flat, NULL,
			NULL, NULL
			);
	// Left Flat
	def_mstage(
			4,
			nxd, nxa, nu,
			NULL, lfcn_reg,
			0, 0, 0, NULL, ffcn_left_flat, NULL,
			NULL, NULL
			);
	// Left Halx
	def_mstage(
			5,
			nxd, nxa, nu,
			NULL, lfcn_reg,
			0, 0, 0, NULL, ffcn_left_halx, NULL,
			NULL, NULL
			);
		
	// Left Halx Right Heel
	def_mstage(
			6,
			nxd, nxa, nu,
			NULL, lfcn_reg,
			0, 0, 0, NULL, ffcn_left_halx_right_heel, NULL,
			NULL, NULL
			);
		
	// Left Halx Right Flat
	def_mstage(
			7,
			nxd, nxa, nu,
			NULL, lfcn_reg,
			0, 0, 0, NULL, ffcn_left_halx_right_flat, NULL,
			NULL, NULL
			);
	
	// Define Constraints for OCP
	def_mpc(0, "s", npr, rdfcn_right_flat_s_n, rdfcn_right_flat_s_ne, rdfcn_right_flat_s, NULL);
	def_mpc(0, "i", npr, rdfcn_right_flat_i_n, rdfcn_right_flat_i_ne, rdfcn_right_flat_i, NULL);

	def_mpc(1, "s", npr, rdfcn_right_halx_s_n, rdfcn_right_halx_s_ne, rdfcn_right_halx_s, NULL);
 	def_mpc(1, "i", npr, rdfcn_right_halx_i_n, rdfcn_right_halx_i_ne, rdfcn_right_halx_i, NULL);
 
 	def_mpc(2, "s", npr, rdfcn_right_halx_left_heel_s_n, rdfcn_right_halx_left_heel_s_ne, rdfcn_right_halx_left_heel_s, NULL);
 	def_mpc(2, "i", npr, rdfcn_right_halx_left_heel_i_n, rdfcn_right_halx_left_heel_i_ne, rdfcn_right_halx_left_heel_i, NULL);
 
 	def_mpc(3, "s", npr, rdfcn_right_halx_left_flat_s_n, rdfcn_right_halx_left_flat_s_ne, rdfcn_right_halx_left_flat_s, NULL); 
 	def_mpc(3, "i", npr, rdfcn_right_halx_left_flat_i_n, rdfcn_right_halx_left_flat_i_ne, rdfcn_right_halx_left_flat_i, NULL);
	
	def_mpc(4, "s", npr, rdfcn_left_flat_s_n, rdfcn_left_flat_s_ne, rdfcn_left_flat_s, NULL);
	def_mpc(4, "i", npr, rdfcn_left_flat_i_n, rdfcn_left_flat_i_ne, rdfcn_left_flat_i, NULL);

	def_mpc(5, "s", npr, rdfcn_left_halx_s_n, rdfcn_left_halx_s_ne, rdfcn_left_halx_s, NULL);
	def_mpc(5, "i", npr, rdfcn_left_halx_i_n, rdfcn_left_halx_i_ne, rdfcn_left_halx_i, NULL);
 
	def_mpc(6, "s", npr, rdfcn_left_halx_right_heel_s_n, rdfcn_left_halx_right_heel_s_ne, rdfcn_left_halx_right_heel_s, NULL); 
	def_mpc(6, "i", npr, rdfcn_left_halx_right_heel_i_n, rdfcn_left_halx_right_heel_i_ne, rdfcn_left_halx_right_heel_i, NULL);
 
	def_mpc(7, "s", npr, rdfcn_left_halx_right_flat_s_n, rdfcn_left_halx_right_flat_s_ne, rdfcn_left_halx_right_flat_s, NULL);
	def_mpc(7, "i", npr, rdfcn_left_halx_right_flat_i_n, rdfcn_left_halx_right_flat_i_ne, rdfcn_left_halx_right_flat_i, NULL);
	def_mpc(7, "e", npr, rdfcn_left_halx_right_flat_e_n, rdfcn_left_halx_right_flat_e_ne, rdfcn_left_halx_right_flat_e, NULL);

	// Define least squares function to be evaluated at stages
	def_lsq(0, "*", 0, nlsq, lsqfcn_trackMotion);
	def_lsq(1, "*", 0, nlsq, lsqfcn_trackMotion);
	def_lsq(2, "*", 0, nlsq, lsqfcn_trackMotion);
	def_lsq(3, "*", 0, nlsq, lsqfcn_trackMotion);
	def_lsq(4, "*", 0, nlsq, lsqfcn_trackMotion);
	def_lsq(5, "*", 0, nlsq, lsqfcn_trackMotion);
	def_lsq(6, "*", 0, nlsq, lsqfcn_trackMotion);
	def_lsq(7, "*", 0, nlsq, lsqfcn_trackMotion);

	// Output plotting and writing results
	def_mio (initialize_from_data, write_ocp_output, hires_data_out);
}