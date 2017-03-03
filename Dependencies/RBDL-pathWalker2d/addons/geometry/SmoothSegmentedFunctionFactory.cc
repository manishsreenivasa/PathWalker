/* -------------------------------------------------------------------------- *
 *        OpenSim:  SmoothSegmentedFunctionFactory.cpp                        *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                                                            *
 * Copyright (c) 2005-2012 Stanford University and the Authors                *
 * Author(s): Matthew Millard                                                 *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */
 /*
  Update:
    This is a port of the original code so that it will work with
    the multibody code RBDL written by Martin Felis.
  
  Author:
    Matthew Millard
  
  Date:
    Nov 2015

*/
    
//=============================================================================
// INCLUDES
//=============================================================================

#include "SmoothSegmentedFunctionFactory.h"
#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;
using namespace RigidBodyDynamics::Addons::Geometry;

//=============================================================================
// STATICS
//=============================================================================
//using namespace std;


static int NUM_SAMPLE_PTS = 100; //The number of knot points to use to sample
                //each Bezier corner section

static double SMOOTHING = 0;   //The amount of smoothing to use when fitting 
                //3rd order splines to the quintic Bezier
                //functions
static bool DEBUG = true;  //When this is set to true, each function's debug
              //routine will be called, which ususally results
              //in a text file of its output being produced

static double UTOL = (double)std::numeric_limits<double>::epsilon()*1e2;

static double INTTOL = (double)std::numeric_limits<double>::epsilon()*1e4;

static int MAXITER = 20;
//=============================================================================
// UTILITY FUNCTIONS
//=============================================================================

//=============================================================================
// MUSCLE CURVE FITTING FUNCTIONS
//=============================================================================
void  SmoothSegmentedFunctionFactory::createFiberActiveForceLengthCurve(
                      double x0,
                      double x1, 
                      double x2, 
                      double x3, 
                      double ylow,  
                      double dydx, 
                      double curviness,
                      const std::string& curveName,
                      SmoothSegmentedFunction& smoothSegmentedFunctionToUpdate)
{
  //Ensure that the inputs are within a valid range
  double rootEPS = sqrt(std::numeric_limits<double>::epsilon());

  if( (!(x0>=0 && x1>x0+rootEPS && x2>x1+rootEPS && x3>x2+rootEPS) ) ){
    cerr  << "SmoothSegmentedFunctionFactory::"
          << "createFiberActiveForceLengthCurve: "
          << curveName
          << ": This must be true: 0 < lce0 < lce1 < lce2 < lce3"
          << endl;
    assert(0);
    abort();     
              
  }


  if( !(ylow >= 0) ){
    cerr << "SmoothSegmentedFunctionFactory::"
         << "createFiberActiveForceLengthCurve:"
         << curveName
         <<": shoulderVal must be greater than, or equal to 0"
         << endl;
    assert(0);
    abort();

  }

  double dydxUpperBound = (1-ylow)/(x2-x1);


  if( !(dydx >= 0 && dydx < dydxUpperBound) ){
    cerr << "SmoothSegmentedFunctionFactory::"
         << "createFiberActiveForceLengthCurve:"
         << curveName
         << ": plateauSlope must be greater than 0 and less than "
         << dydxUpperBound
         << endl;
    assert(0);
    abort();
  }

  if( !(curviness >= 0 && curviness <= 1) ){
    cerr  << "SmoothSegmentedFunctionFactory::"
          << "createFiberActiveForceLengthCurve:"
          << curveName
          << ": curviness must be between 0 and 1"
          << endl;
    assert(0);
    abort();
  }

  std::string name = curveName;
  name.append(".createFiberActiveForceLengthCurve");



  //Translate the users parameters into Bezier curves
    double c = SegmentedQuinticBezierToolkit::scaleCurviness(curviness);

  //The active force length curve is made up of 5 elbow shaped sections. 
  //Compute the locations of the joining point of each elbow section.

  //Calculate the location of the shoulder
     double xDelta = 0.05*x2; //half the width of the sarcomere 0.0259, 
                 //but TM.Winter's data has a wider shoulder than
                 //this

     double xs  = (x2-xDelta);//x1 + 0.75*(x2-x1);
   
   //Calculate the intermediate points located on the ascending limb
     double y0  = 0;   
     double dydx0 = 0;

     double y1  = 1 - dydx*(xs-x1);
     double dydx01= 1.25*(y1-y0)/(x1-x0);//(y1-y0)/(x1-(x0+xDelta));

     double x01   = x0 + 0.5*(x1-x0); //x0 + xDelta + 0.5*(x1-(x0+xDelta));
     double y01   = y0 + 0.5*(y1-y0);
   
   //Calculate the intermediate points of the shallow ascending plateau
     double x1s   = x1 + 0.5*(xs-x1);
     double y1s   = y1 + 0.5*(1-y1);
     double dydx1s= dydx;
   
     //double dydx01c0 = 0.5*(y1s-y01)/(x1s-x01) + 0.5*(y01-y0)/(x01-x0);
     //double dydx01c1 = 2*( (y1-y0)/(x1-x0));
     //double dydx01(1-c)*dydx01c0 + c*dydx01c1; 
     
     //x2 entered
     double y2 = 1;
     double dydx2 = 0;
   
   //Descending limb
     //x3 entered
     double y3 = 0;
     double dydx3 = 0;
     
     double x23 = (x2+xDelta) + 0.5*(x3-(x2+xDelta)); //x2 + 0.5*(x3-x2);
     double y23 = y2 + 0.5*(y3-y2);
       
     //double dydx23c0 = 0.5*((y23-y2)/(x23-x2)) + 0.5*((y3-y23)/(x3-x23));
     //double dydx23c1 = 2*(y3-y2)/(x3-x2);
     double dydx23   = (y3-y2)/((x3-xDelta)-(x2+xDelta)); 
     //(1-c)*dydx23c0 + c*dydx23c1; 
  
  //Compute the locations of the control points
     RigidBodyDynamics::Math::MatrixNd p0 = SegmentedQuinticBezierToolkit::
    calcQuinticBezierCornerControlPoints(x0,ylow,dydx0,x01,y01,dydx01,c);
     RigidBodyDynamics::Math::MatrixNd p1 = SegmentedQuinticBezierToolkit::
    calcQuinticBezierCornerControlPoints(x01,y01,dydx01,x1s,y1s,dydx1s,c);
     RigidBodyDynamics::Math::MatrixNd p2 = SegmentedQuinticBezierToolkit::
    calcQuinticBezierCornerControlPoints(x1s,y1s,dydx1s,x2, y2, dydx2,c);
     RigidBodyDynamics::Math::MatrixNd p3 = SegmentedQuinticBezierToolkit::
    calcQuinticBezierCornerControlPoints(x2, y2, dydx2,x23,y23,dydx23,c);
     RigidBodyDynamics::Math::MatrixNd p4 = SegmentedQuinticBezierToolkit::
    calcQuinticBezierCornerControlPoints(x23,y23,dydx23,x3,ylow,dydx3,c);
                  
    RigidBodyDynamics::Math::MatrixNd mX(6,5), mY(6,5);
    mX.col(0) = p0.col(0);
    mX.col(1) = p1.col(0);
    mX.col(2) = p2.col(0);
    mX.col(3) = p3.col(0);
    mX.col(4) = p4.col(0);

    mY.col(0) = p0.col(1);
    mY.col(1) = p1.col(1);
    mY.col(2) = p2.col(1);
    mY.col(3) = p3.col(1);
    mY.col(4) = p4.col(1);

    smoothSegmentedFunctionToUpdate.updSmoothSegmentedFunction(
                                      mX,mY,x0,x3,ylow,ylow,0,0,curveName);
}

void SmoothSegmentedFunctionFactory::createFiberForceVelocityCurve(
    double fmaxE, 
    double dydxC, 
    double dydxNearC, 
    double dydxIso, 
    double dydxE, 
    double dydxNearE,
    double concCurviness,
    double eccCurviness,
    const std::string& curveName,
    SmoothSegmentedFunction& smoothSegmentedFunctionToUpdate)
{
  //Ensure that the inputs are within a valid range
  
  if( !(fmaxE > 1.0) ){
    cerr << "SmoothSegmentedFunctionFactory::"
        << "createFiberForceVelocityCurve: "
        << curveName
        <<": fmaxE must be greater than 1"
        << endl;
    assert(0);
    abort(); 
  }
  
  if( !(dydxC >= 0.0 && dydxC < 1) ){ 
    cerr  << "SmoothSegmentedFunctionFactory::"
          << "createFiberForceVelocityCurve: "
          << curveName
          << ": dydxC must be greater than or equal to 0 "
          <<" and less than 1"
          << endl;
    assert(0);
    abort();     
  }

  if( !(dydxNearC > dydxC && dydxNearC <= 1) ){
    cerr  << "SmoothSegmentedFunctionFactory::"
          << "createFiberForceVelocityCurve: "
          << curveName
          << ": dydxNearC must be greater than or equal to 0 "
          << "and less than 1"
          << endl;
    assert(0);
    abort();
  }

  if( !(dydxIso > 1) ){ 
    cerr  << "SmoothSegmentedFunctionFactory::"
          << "createFiberForceVelocityCurve: "
          << curveName
          << ": dydxIso must be greater than (fmaxE-1)/1 ("
          << ((fmaxE-1.0)/1.0)
          << ")"
          << endl;
    assert(0);
    abort();
  }

  if( !(dydxE >= 0.0 && dydxE < (fmaxE-1)) ){
    cerr  << "SmoothSegmentedFunctionFactory::"
          << "createFiberForceVelocityCurve: "
          <<  curveName
          <<": dydxE must be greater than or equal to 0 "
          << "and less than fmaxE-1 ("
          << (fmaxE-1) << ")"
          << endl;
    assert(0);
    abort();          
  }

  if(!(dydxNearE >= dydxE && dydxNearE < (fmaxE-1))){
    cerr  << "SmoothSegmentedFunctionFactory::"
          << "createFiberForceVelocityCurve"
          << curveName
          << ": dydxNearE must be greater than or equal to dydxE "
          << "and less than fmaxE-1 (" << (fmaxE-1)
          << ")"
          << endl;
    assert(0);
    abort();          
  }

  if(! (concCurviness <= 1.0 && concCurviness >= 0)){ 
    cerr  << "SmoothSegmentedFunctionFactory::"
          << "createFiberForceVelocityCurve "
          << curveName
          << ": concCurviness must be between 0 and 1"
          << endl;
    assert(0);
    abort();          
  }

  if(! (eccCurviness <= 1.0 && eccCurviness >= 0)){
    cerr  << "SmoothSegmentedFunctionFactory::"
          << "createFiberForceVelocityCurve "
          << curveName    
          << ": eccCurviness must be between 0 and 1"
          << endl;
    assert(0);
    abort();         
  }

  std::string name = curveName;
  name.append(".createFiberForceVelocityCurve");

  //Translate the users parameters into Bezier point locations
  double cC = SegmentedQuinticBezierToolkit::scaleCurviness(concCurviness);
  double cE = SegmentedQuinticBezierToolkit::scaleCurviness(eccCurviness);
  
  //Compute the concentric control point locations
  double xC   = -1;
  double yC   = 0;
  
  double xNearC = -0.9;
  double yNearC = yC + 0.5*dydxNearC*(xNearC-xC) + 0.5*dydxC*(xNearC-xC);

  double xIso = 0;
  double yIso = 1;

  double xE   = 1;
  double yE   = fmaxE;

  double xNearE = 0.9;
  double yNearE = yE + 0.5*dydxNearE*(xNearE-xE) + 0.5*dydxE*(xNearE-xE);


  RigidBodyDynamics::Math::MatrixNd concPts1 = SegmentedQuinticBezierToolkit::
    calcQuinticBezierCornerControlPoints(  xC,   yC,  dydxC, 
                       xNearC, yNearC,dydxNearC,cC);
  RigidBodyDynamics::Math::MatrixNd concPts2 = SegmentedQuinticBezierToolkit::
    calcQuinticBezierCornerControlPoints(xNearC,yNearC,dydxNearC, 
                         xIso,  yIso,  dydxIso,  cC);
  RigidBodyDynamics::Math::MatrixNd eccPts1 = SegmentedQuinticBezierToolkit::
    calcQuinticBezierCornerControlPoints(  xIso,  yIso,  dydxIso, 
                       xNearE,  yNearE,  dydxNearE, cE);
  RigidBodyDynamics::Math::MatrixNd eccPts2 = SegmentedQuinticBezierToolkit::
    calcQuinticBezierCornerControlPoints(xNearE, yNearE, dydxNearE, 
                         xE,   yE,   dydxE, cE);

  RigidBodyDynamics::Math::MatrixNd mX(6,4), mY(6,4);
  mX.col(0) = concPts1.col(0);
  mX.col(1) = concPts2.col(0);
  mX.col(2) = eccPts1.col(0);
  mX.col(3) = eccPts2.col(0);

  mY.col(0) = concPts1.col(1);
  mY.col(1) = concPts2.col(1);
  mY.col(2) = eccPts1.col(1);
  mY.col(3) = eccPts2.col(1);

  smoothSegmentedFunctionToUpdate.updSmoothSegmentedFunction(
                                    mX,mY,xC,xE,yC,yE,dydxC,dydxE,curveName);
}


void SmoothSegmentedFunctionFactory::createFiberForceVelocityInverseCurve(
                  double fmaxE, 
                  double dydxC, 
                  double dydxNearC, 
                  double dydxIso,
                  double dydxE, 
                  double dydxNearE,
                  double concCurviness, 
                  double eccCurviness,
                  const std::string& curveName,
                  SmoothSegmentedFunction& smoothSegmentedFunctionToUpdate)
{
  //Ensure that the inputs are within a valid range  
  if(! (fmaxE > 1.0 )){ 
    cerr  << "SmoothSegmentedFunctionFactory::"
          << "createFiberForceVelocityInverseCurve: "
          << curveName
          << ": fmaxE must be greater than 1"
          << endl;
    assert(0);
    abort();
  }
  
  double SimTKSignificantReal = 
    pow((double)std::numeric_limits<double>::epsilon(), 7.0/8.0);

  if(! (dydxC > SimTKSignificantReal && dydxC < 1 )){ 
    cerr  << "SmoothSegmentedFunctionFactory::"
          << "createFiberForceVelocityInverseCurve "
          << curveName
          << ": dydxC must be greater than 0"
          << "and less than 1"
          << endl;
    assert(0);
    abort();

  }

  if(! (dydxNearC > dydxC && dydxNearC < 1 )){
    std::stringstream errMsg;
    cerr  << "SmoothSegmentedFunctionFactory::"
          << "createFiberForceVelocityInverseCurve "
          << ": dydxNearC must be greater than 0 "
          << curveName
          << " and less than 1"
          << endl;
    assert(0);
    abort();
  }

  if(! (dydxIso > 1)){ 
    cerr  << "SmoothSegmentedFunctionFactory::"
          << "createFiberForceVelocityInverseCurve "
          << curveName
          << ": dydxIso must be greater than or equal to 1"
          << endl;
    assert(0);
    abort();
  }

  //double SimTKSignificantReal =
  //  pow(std::numeric_limits<double>::epsilon(), 7.0/8.0);

  if(! (dydxE > SimTKSignificantReal && dydxE < (fmaxE-1)) ){
    cerr  << "SmoothSegmentedFunctionFactory::"
          << "createFiberForceVelocityInverseCurve "
          << curveName
          << ": dydxE must be greater than or equal to 0"
          << " and less than fmaxE-1 (" << (fmaxE-1) << ")"
          << endl;
    assert(0);
    abort();
  }

  if(! (dydxNearE >= dydxE && dydxNearE < (fmaxE-1)) ){ 
    cerr  << "SmoothSegmentedFunctionFactory::"
          << "createFiberForceVelocityInverseCurve "
          << curveName
          << ": dydxNearE must be greater than or equal to dydxE"
          << "and less than fmaxE-1 ("<< (fmaxE-1) << ")"
          << endl;
    assert(0);
    abort();
  }

  if(! (concCurviness <= 1.0 && concCurviness >= 0) ){ 
    cerr  << "SmoothSegmentedFunctionFactory::"
          << "createFiberForceVelocityInverseCurve "
          << curveName
          << ": concCurviness must be between 0 and 1"
          << endl;
    assert(0);
    abort();
  }

  if(! (eccCurviness <= 1.0 && eccCurviness >= 0) ){ 
    cerr  << "SmoothSegmentedFunctionFactory::"
          << "createFiberForceVelocityInverseCurve "
          << curveName
          << ": eccCurviness must be between 0 and 1"
          << endl;
    assert(0);
    abort();
  }

  std::string name = curveName;
  name.append(".createFiberForceVelocityInverseCurve");

  //Translate the users parameters into Bezier point locations
  double cC = SegmentedQuinticBezierToolkit::scaleCurviness(concCurviness);
  double cE = SegmentedQuinticBezierToolkit::scaleCurviness(eccCurviness);
  
  //Compute the concentric control point locations
  double xC   = -1;
  double yC   = 0;
  
  double xNearC = -0.9;
  double yNearC = yC + 0.5*dydxNearC*(xNearC-xC) + 0.5*dydxC*(xNearC-xC);

  double xIso = 0;
  double yIso = 1;

  double xE   = 1;
  double yE   = fmaxE;

  double xNearE = 0.9;
  double yNearE = yE + 0.5*dydxNearE*(xNearE-xE) + 0.5*dydxE*(xNearE-xE);


  RigidBodyDynamics::Math::MatrixNd concPts1 = SegmentedQuinticBezierToolkit::
    calcQuinticBezierCornerControlPoints(   xC,   yC,  dydxC, 
                      xNearC, yNearC,dydxNearC,cC);
  RigidBodyDynamics::Math::MatrixNd concPts2 = SegmentedQuinticBezierToolkit::
    calcQuinticBezierCornerControlPoints(xNearC,yNearC,dydxNearC, 
                         xIso,  yIso,  dydxIso,  cC);
  RigidBodyDynamics::Math::MatrixNd eccPts1 = SegmentedQuinticBezierToolkit::
    calcQuinticBezierCornerControlPoints(  xIso,  yIso,  dydxIso, 
                       xNearE,  yNearE,  dydxNearE, cE);
  RigidBodyDynamics::Math::MatrixNd eccPts2 = SegmentedQuinticBezierToolkit::
    calcQuinticBezierCornerControlPoints(xNearE, yNearE, dydxNearE, 
                         xE,   yE,   dydxE, cE);

  RigidBodyDynamics::Math::MatrixNd mX(6,4), mY(6,4);
  mX.col(0) = concPts1.col(0);
  mX.col(1) = concPts2.col(0);
  mX.col(2) =  eccPts1.col(0);
  mX.col(3) =  eccPts2.col(0);

  mY.col(0) = concPts1.col(1);
  mY.col(1) = concPts2.col(1);
  mY.col(2) =  eccPts1.col(1);
  mY.col(3) =  eccPts2.col(1);
  
  smoothSegmentedFunctionToUpdate.updSmoothSegmentedFunction(
                                  mY,mX,yC,yE,xC,xE,1/dydxC,1/dydxE, curveName);

}

void SmoothSegmentedFunctionFactory::createFiberCompressiveForcePennationCurve(
                  double phi0, 
                  double k, 
                  double curviness, 
                  const std::string& curveName,
                  SmoothSegmentedFunction& smoothSegmentedFunctionToUpdate)
{
  //Check the input arguments
  if( !(phi0>0 && phi0<(M_PI/2.0)) ){ 
    cerr  << "SmoothSegmentedFunctionFactory::"
          << "createFiberCompressiveForcePennationCurve "
          << curveName
          << ": phi0 must be greater than 0, and less than Pi/2"
          << endl;
    assert(0);
    abort();
  }

  if( !(k > (1.0/(M_PI/2.0-phi0))) ){ 
    cerr  << "SmoothSegmentedFunctionFactory::"
          << "createFiberCompressiveForcePennationCurve " 
          << curveName
          << ": k must be greater than " << (1.0/(M_PI/2.0-phi0))
          << endl;
    assert(0);
    abort();          
  }

  if( !(curviness>=0 && curviness <= 1)  ){
    cerr << "SmoothSegmentedFunctionFactory::"
         << "createFiberCompressiveForcePennationCurve "
         << curveName 
         << ": curviness must be between 0.0 and 1.0"
          << endl;
    assert(0);
    abort();         
  }

  std::string name=curveName;
  name.append(".createFiberCompressiveForcePennationCurve");

  //Translate the user parameters to quintic Bezier points
  double c = SegmentedQuinticBezierToolkit::scaleCurviness(curviness);
  double x0 = phi0;
  double y0 = 0;
  double dydx0 = 0;
  double x1 = M_PI/2.0;
  double y1 = 1;
  double dydx1 = k;

  RigidBodyDynamics::Math::MatrixNd ctrlPts = SegmentedQuinticBezierToolkit::
    calcQuinticBezierCornerControlPoints(x0,y0,dydx0,x1,y1,dydx1,c);
  
  RigidBodyDynamics::Math::MatrixNd mX(6,1), mY(6,1);
  mX.col(0) = ctrlPts.col(0);
  mY.col(0) = ctrlPts.col(1);

  smoothSegmentedFunctionToUpdate.updSmoothSegmentedFunction(
                                    mX,mY,x0,x1,y0,y1,dydx0,dydx1,curveName);
}

void SmoothSegmentedFunctionFactory::
      createFiberCompressiveForceCosPennationCurve(
        double cosPhi0, 
        double k, 
        double curviness, 
        const std::string& curveName,
        SmoothSegmentedFunction& smoothSegmentedFunctionToUpdate)
{
  //Check the input arguments
  if( !(cosPhi0>0 && cosPhi0 < 1) ){ 
    cerr  << "SmoothSegmentedFunctionFactory::"
          << "createFiberCompressiveForceCosPennationCurve " 
          << curveName
          << ": cosPhi0 must be greater than 0, and less than 1"
          << endl;
    assert(0);
    abort();          
  }

  if( !(k < 1/cosPhi0) ){ 
    cerr  << "SmoothSegmentedFunctionFactory::"
          << "createFiberCompressiveForceCosPennationCurve " 
          << curveName
          << ": k must be less than 0"
          << endl;
    assert(0);
    abort();          
  }

  if( !(curviness>=0 && curviness <= 1) ){ 
    cerr  << "SmoothSegmentedFunctionFactory::"
          << "createFiberCompressiveForceCosPennationCurve"
          << curveName
          << ": curviness must be between 0.0 and 1.0"
          << endl;
    assert(0);
    abort();   
  }

  std::string name=curveName;
  name.append(".createFiberCompressiveForceCosPennationCurve");

  //Translate the user parameters to quintic Bezier points
  double c = SegmentedQuinticBezierToolkit::scaleCurviness(curviness);
  double x0 = 0;
  double y0 = 1;
  double dydx0 = k;
  double x1 = cosPhi0;
  double y1 = 0;
  double dydx1 = 0;

  RigidBodyDynamics::Math::MatrixNd ctrlPts = SegmentedQuinticBezierToolkit::
    calcQuinticBezierCornerControlPoints(x0,y0,dydx0,x1,y1,dydx1,c);
  
  RigidBodyDynamics::Math::MatrixNd mX(6,1), mY(6,1);
  mX.col(0) = ctrlPts.col(0);
  mY.col(0) = ctrlPts.col(1);

  smoothSegmentedFunctionToUpdate.updSmoothSegmentedFunction(
                                    mX,mY,x0,x1,y0,y1,dydx0,dydx1,curveName);
}

void SmoothSegmentedFunctionFactory::createFiberCompressiveForceLengthCurve(
                  double lmax, 
                  double k, 
                  double curviness, 
                  const std::string& curveName,
                  SmoothSegmentedFunction& smoothSegmentedFunctionToUpdate)
{

  if( !(lmax>0) ){ 
    cerr  << "SmoothSegmentedFunctionFactory::"
          << "createFiberCompressiveForceLength "
          << curveName
          << ": l0 must be greater than 0"
          << endl;
    assert(0);
    abort();     
  }

  if( !(k < -(1.0/lmax)) ){ 
    cerr  << "SmoothSegmentedFunctionFactory::"
          << "createFiberCompressiveForceLength "
          << curveName
          << ": k must be less than "
          << -(1.0/lmax)
          << endl;
    assert(0);
    abort();
  }

  if( !(curviness>=0 && curviness <= 1) ){ 
    cerr  << "SmoothSegmentedFunctionFactory::"
          << "createFiberCompressiveForceLength "
          << curveName
          << ": curviness must be between 0.0 and 1.0"
          << endl;
    assert(0);
    abort();
  }

  std::string caller = curveName;
  caller.append(".createFiberCompressiveForceLength");

  //Translate the user parameters to quintic Bezier points
  double c = SegmentedQuinticBezierToolkit::scaleCurviness(curviness);
  double x0 = 0.0;
  double y0 = 1;
  double dydx0 = k;
  double x1 = lmax;
  double y1 = 0;
  double dydx1 = 0;

  RigidBodyDynamics::Math::MatrixNd ctrlPts = SegmentedQuinticBezierToolkit::
    calcQuinticBezierCornerControlPoints(x0,y0,dydx0,x1,y1,dydx1,c);

  RigidBodyDynamics::Math::MatrixNd mX(6,1), mY(6,1);
  mX.col(0) = ctrlPts.col(0);
  mY.col(0) = ctrlPts.col(1);

  smoothSegmentedFunctionToUpdate.updSmoothSegmentedFunction(
                                    mX,mY,x0,x1,y0,y1,dydx0,dydx1,curveName);
}


void  SmoothSegmentedFunctionFactory::createFiberForceLengthCurve(
                  double eZero, 
                  double eIso, 
                  double kLow, 
                  double kIso, 
                  double curviness,
                  const std::string& curveName,
                  SmoothSegmentedFunction& smoothSegmentedFunctionToUpdate)
{
  
  if( !(eIso > eZero) ){ 
    cerr  << "SmoothSegmentedFunctionFactory::"
          << "createFiberForceLength "
          << curveName
          << ": The following must hold: eIso  > eZero"
          << endl;
    assert(0);
    abort();
  }

  if( !(kIso > (1.0/(eIso-eZero))) ){ 
    std::stringstream errMsg;
    cerr  << "SmoothSegmentedFunctionFactory::"
          << "createFiberForceLength " 
          << curveName
          << ": kiso must be greater than 1/(eIso-eZero) ("
          << (1.0/(eIso-eZero)) << ")"
          << endl;
    assert(0);
    abort(); 
  }

  if( !(kLow > 0.0 && kLow < 1/(eIso-eZero)) ){
    cerr  << "SmoothSegmentedFunctionFactory::"
          << "createFiberForceLength " 
          << curveName
          << ": kLow must be greater than 0 and less than"
          << 1.0/(eIso-eZero)
          << endl;
    assert(0);
    abort();
  }

  if( !(curviness>=0 && curviness <= 1) ){ 
    cerr  << "SmoothSegmentedFunctionFactory::"
          << "createFiberForceLength "
          << curveName
          << ": curviness must be between 0.0 and 1.0"
          << endl;
    assert(0);
    abort();
  }

  std::string callerName = curveName;
  callerName.append(".createFiberForceLength");


    //Translate the user parameters to quintic Bezier points
  double c = SegmentedQuinticBezierToolkit::scaleCurviness(curviness);
  double xZero = 1+eZero;
  double yZero = 0;
  
  double xIso = 1 + eIso;
  double yIso = 1;
  
  double deltaX = std::min(0.1*(1.0/kIso), 0.1*(xIso-xZero));

  double xLow   = xZero + deltaX;
  double xfoot  = xZero + 0.5*(xLow-xZero);
  double yfoot  = 0;
  double yLow   = yfoot + kLow*(xLow-xfoot);

  //Compute the Quintic Bezier control points
  RigidBodyDynamics::Math::MatrixNd p0 = SegmentedQuinticBezierToolkit::
   calcQuinticBezierCornerControlPoints(xZero, yZero, 0,
                       xLow, yLow,  kLow,c);
  
  RigidBodyDynamics::Math::MatrixNd p1 = SegmentedQuinticBezierToolkit::
   calcQuinticBezierCornerControlPoints(xLow, yLow, kLow,
                      xIso, yIso, kIso, c);
  RigidBodyDynamics::Math::MatrixNd mX(6,2);
  RigidBodyDynamics::Math::MatrixNd mY(6,2);

  mX.col(0) = p0.col(0);
  mY.col(0) = p0.col(1);

  mX.col(1) = p1.col(0);
  mY.col(1) = p1.col(1);
  
  smoothSegmentedFunctionToUpdate.updSmoothSegmentedFunction(
              mX, mY, xZero, xIso, yZero, yIso, 0.0, kIso, curveName);
}




void SmoothSegmentedFunctionFactory::
      createTendonForceLengthCurve( double eIso, double kIso, 
                    double fToe, double curviness,
                    const std::string& curveName,
                    SmoothSegmentedFunction& smoothSegmentedFunctionToUpdate)
{
  
  if( !(eIso>0) ){ 
    cerr  << "SmoothSegmentedFunctionFactory::"
          << "createTendonForceLengthCurve " 
          << curveName
          << ": eIso must be greater than 0, but "
          << eIso << " was entered"
          << endl;
    assert(0);
    abort();

  }   

  if( !(fToe>0 && fToe < 1) ){ 
    cerr  << "SmoothSegmentedFunctionFactory::"
          << "createTendonForceLengthCurve "
          << curveName 
          << ": fToe must be greater than 0 and less than 1, but "
          << fToe
          << " was entered"
          << endl;
    assert(0);
    abort();
  }

  if( !(kIso > (1/eIso)) ){ 
    cerr  << "SmoothSegmentedFunctionFactory::"
          << "createTendonForceLengthCurve "
          << curveName
          << ": kIso must be greater than 1/eIso, ("
          << (1/eIso) << "), but kIso (" 
          << kIso << ") was entered"
          << endl;
    assert(0);
    abort();
  }


  if( !(curviness>=0 && curviness <= 1) ){ 
    cerr  << "SmoothSegmentedFunctionFactory::"
          << "createTendonForceLengthCurve "
          << curveName
          << " : curviness must be between 0.0 and 1.0, but "
          << curviness << " was entered"
          << endl;
    assert(0);
    abort();
  }

  std::string callerName = curveName;
  callerName.append(".createTendonForceLengthCurve");

  //Translate the user parameters to quintic Bezier points
  double c = SegmentedQuinticBezierToolkit::scaleCurviness(curviness);
  double x0 = 1.0;
  double y0 = 0;
  double dydx0 = 0;

  double xIso = 1.0 + eIso;
  double yIso = 1;
  double dydxIso = kIso;

  //Location where the curved section becomes linear
  double yToe = fToe;
  double xToe = (yToe-1)/kIso + xIso;


  //To limit the 2nd derivative of the toe region the line it tends to
  //has to intersect the x axis to the right of the origin
    double xFoot = 1.0+(xToe-1.0)/10.0;
    double yFoot = 0;
  double dydxToe = (yToe-yFoot)/(xToe-xFoot);

  //Compute the location of the corner formed by the average slope of the
  //toe and the slope of the linear section
  double yToeMid = yToe*0.5;
  double xToeMid = (yToeMid-yIso)/kIso + xIso;
  double dydxToeMid = (yToeMid-yFoot)/(xToeMid-xFoot);

  //Compute the location of the control point to the left of the corner
  double xToeCtrl = xFoot + 0.5*(xToeMid-xFoot); 
  double yToeCtrl = yFoot + dydxToeMid*(xToeCtrl-xFoot);



  //Compute the Quintic Bezier control points
  RigidBodyDynamics::Math::MatrixNd p0 = SegmentedQuinticBezierToolkit::
   calcQuinticBezierCornerControlPoints(x0,y0,dydx0,
                    xToeCtrl,yToeCtrl,dydxToeMid,c);
  RigidBodyDynamics::Math::MatrixNd p1 = SegmentedQuinticBezierToolkit::
   calcQuinticBezierCornerControlPoints(xToeCtrl, yToeCtrl, dydxToeMid,
                        xToe,   yToe,  dydxIso, c);
  RigidBodyDynamics::Math::MatrixNd mX(6,2);
  RigidBodyDynamics::Math::MatrixNd mY(6,2);

  mX.col(0) = p0.col(0);
  mY.col(0) = p0.col(1);

  mX.col(1) = p1.col(0);
  mY.col(1) = p1.col(1);

  smoothSegmentedFunctionToUpdate.updSmoothSegmentedFunction(
                    mX, mY, x0, xToe, y0,  yToe, dydx0, dydxIso, curveName);                    
  
}

//=============================================================================
// Anderson 2007 Active Torque Angle Curve
//=============================================================================


void SmoothSegmentedFunctionFactory::
  createAnderson2007ActiveTorqueAngleCurve(
    double c2, 
    double c3,
    const std::string& curveName,
    SmoothSegmentedFunction& smoothSegmentedFunctionToUpdate)
{
  //Check the input arguments
  if( !(c2 > 0) ){ 
    cerr  << "SmoothSegmentedFunctionFactory::"
          << "createAnderson2007ActiveTorqueAngleCurve " 
          << curveName
          << ": c2 must be greater than 0"
          << endl;
    assert(0);
    abort();          
  }


  std::string name=curveName;
  name.append(".createAnderson2007ActiveTorqueAngleCurve");

  //For now these advanced paramters are hidden. They will only be
  //uncovered if absolutely necessary.
  double minValueAtShoulders = 0;
  double minShoulderSlopeMagnitude = 0;

  double curviness = 0.5;
  double c = SegmentedQuinticBezierToolkit::scaleCurviness(curviness);

  //Translate the user parameters to quintic Bezier points
  double x0 = c3 - 1.05*(0.5*(M_PI/c2));
  double x1 = c3 - 0.95*(0.5*(M_PI/c2));
  double x2 = c3;
  double x3 = c3 + 0.95*(0.5*(M_PI/c2));
  double x4 = c3 + 1.05*(0.5*(M_PI/c2));

  double y0 = minValueAtShoulders;
  double y1 = cos(c2*(x1-c3));
  double y2 = cos(c2*(x2-c3));
  double y3 = cos(c2*(x3-c3));
  double y4 = minValueAtShoulders;

  double dydx0 =  minShoulderSlopeMagnitude;
  double dydx1 = -sin(c2*(x1-c3))*c2;
  double dydx2 = -sin(c2*(x2-c3))*c2;
  double dydx3 = -sin(c2*(x3-c3))*c2;
  double dydx4 = -minShoulderSlopeMagnitude;


  //Compute the Quintic Bezier control points
  RigidBodyDynamics::Math::MatrixNd p0 = 
    SegmentedQuinticBezierToolkit::
       calcQuinticBezierCornerControlPoints(x0,y0,dydx0,
                          x1,y1,dydx1,c);

  RigidBodyDynamics::Math::MatrixNd p1 = 
    SegmentedQuinticBezierToolkit::
       calcQuinticBezierCornerControlPoints(x1,y1,dydx1,
                          x2,y2,dydx2,c);

  RigidBodyDynamics::Math::MatrixNd p2 = 
    SegmentedQuinticBezierToolkit::
       calcQuinticBezierCornerControlPoints(x2,y2,dydx2,
                          x3,y3,dydx3,c);

  RigidBodyDynamics::Math::MatrixNd p3 = 
    SegmentedQuinticBezierToolkit::
       calcQuinticBezierCornerControlPoints(x3,y3,dydx3,
                          x4,y4,dydx4,c);                                

  RigidBodyDynamics::Math::MatrixNd mX(6,4);
  RigidBodyDynamics::Math::MatrixNd mY(6,4);

  mX.col(0) = p0.col(0);
  mY.col(0) = p0.col(1);
  mX.col(1) = p1.col(0);
  mY.col(1) = p1.col(1);
  mX.col(2) = p2.col(0);
  mY.col(2) = p2.col(1);
  mX.col(3) = p3.col(0);
  mY.col(3) = p3.col(1);


  smoothSegmentedFunctionToUpdate.updSmoothSegmentedFunction(
              mX, mY, x0, x4, y0, y4, dydx0, dydx4, curveName);
}

//=============================================================================
// ANDERSON 2007 Active Torque Angular Velocity Curve
//=============================================================================
void SmoothSegmentedFunctionFactory::
  createAnderson2007ActiveTorqueVelocityCurve(
    double c4, 
    double c5,
    double c6,
    double minEccentricMultiplier,
    double maxEccentricMultiplier,
    const std::string& curveName,
    SmoothSegmentedFunction& smoothSegmentedFunctionToUpdate)
{
  //Check the input arguments
  if( !(c4 < c5) ){ 
    cerr  << "SmoothSegmentedFunctionFactory::"
          << "createAndersonActiveTorqueVelocityCurve " 
          << curveName
          << ": c4 must be greater than c5"
          << endl;
    assert(0);
    abort();          
  }

  if( !((c4 > 0)) ){ 
    cerr  << "SmoothSegmentedFunctionFactory::"
          << "createAndersonActiveTorqueVelocityCurve " 
          << curveName
          << ": c4 must be greater than 0"
          << endl;
    assert(0);
    abort();          
  }

  if( !(c6 > 0.0) ){ 
    cerr  << "SmoothSegmentedFunctionFactory::"
          << "createAndersonActiveTorqueVelocityCurve " 
          << curveName
          << ": c6 must be greater than 1.0"
          << endl;
    assert(0);
    abort();          
  }

  if( !(minEccentricMultiplier > 1.0) ){ 
    cerr  << "SmoothSegmentedFunctionFactory::"
          << "createAndersonActiveTorqueVelocityCurve " 
          << curveName
          << ": minEccentricMultiplier must be greater than 1.0"
          << endl;
    assert(0);
    abort(); 
  }

  if( !(maxEccentricMultiplier > minEccentricMultiplier) ){ 
    cerr  << "SmoothSegmentedFunctionFactory::"
          << "createAndersonActiveTorqueVelocityCurve " 
          << curveName
          << ": maxEccentricMultiplier must be greater than "
          << " minEccentricMultiplier"
          << endl;
    assert(0);
    abort();
  }

  //Advanced settings that we'll hide for now
  double minShoulderSlopeMagnitude  = 0;
  double curviness = 0.75;
  double c         = SegmentedQuinticBezierToolkit::scaleCurviness(curviness);

  //Go and get the value of the curve that is closest to 
  //the maximum contraction velocity by setting rhs of Eqn. 9
  //to 0 and solving
  double dthMaxConc = fabs( 2.0*c4*c5/(c5-3.0*c4) );

  //Depending on the selection of parameters, the original Anderson 
  //torque-velocity curves won't actually go to 0 on the concentric
  //side. This only show up in the case of the plantar flexors.

  //To get the maximum concentric contraction velocity we'll back 
  //up from the maximum calculated above and linearly extrapolate
  //to zero.
  double x1  =   dthMaxConc*0.80;

  double y1den =  (2*c4*c5 + x1*(2*c5-4*c4));            
  double y1  =   (2*c4*c5 + x1*(c5-3*c4))/y1den;

  double dydx1 =   (c5-3*c4)/(2*c4*c5 + x1*(2*c5-4*c4)) 
          -(2*c4*c5 + x1*(c5-3*c4))*(2*c5-4*c4)
          /( y1den*y1den );

  if( !(dydx1 < 0.0) ){ 
    cerr  << "SmoothSegmentedFunctionFactory::"
          << "createAndersonActiveTorqueVelocityCurve " 
          << curveName
          << ": choice of c3,c5, and c6 is such that "
          << " the curve will never cross the x axis."
          << endl;
    assert(0);
    abort();   
  }

  //Linearly extrapolate to get the maximum velocity, and then go just a
  //little futher - so that curve stays smooth
  double x0 = (x1-y1/dydx1)*1.1;

  double x2 = c5; // 50% iso
  double x3 = c4; // 75% iso
  double x4 = 0; 
  double x5 = -x0;//-60*pi/180; 
  //Anderson et al. 2007: pg 3107, column 2, paragraph 2.
  //they made 1 measurement at 60 degree's sec.

  double y0 = 0;

  double y2 =  (2*c4*c5 + x2*(c5-3*c4))
        /(2*c4*c5 + x2*(2*c5-4*c4));

  double y3 =  (2*c4*c5 + x3*(c5-3*c4))
        /(2*c4*c5 + x3*(2*c5-4*c4));

  double y4 = 1;

  double y5 = minEccentricMultiplier; 
  double x60  = -60*M_PI/180;
  double y5at60 = (2*c4*c5 - x60*(  c5-3*c4))*(1.0-c6*x60)
           /(2*c4*c5 - x60*(2*c5-4*c4));

  //This expression to evaluate y5 given the value of the
  //eccentric curve at an angular velocity of 60/s
  //(where Anderson actually took a measurement) is a
  //hueristic that fits most of Anderson's data that I have
  //observed
  y5 = 1+3*(y5at60-1);

  if(y5 < minEccentricMultiplier){
     y5 = minEccentricMultiplier;
  }
  if(y5 > maxEccentricMultiplier){
     y5 = maxEccentricMultiplier; 
  }



  double dydx0 = -minShoulderSlopeMagnitude;

  double tmpD =   (2*c4*c5 + x2*(2*c5-4*c4));
  double dydx2 =  (c5-3*c4)/(tmpD)
          -(2*c4*c5 + x2*(c5-3*c4))*(2*c5-4*c4)/(tmpD*tmpD);
      
  tmpD     =  (2*c4*c5 + x3*(2*c5-4*c4));      
  double dydx3 =  (c5-3*c4)/(tmpD)
          -(2*c4*c5 + x3*(c5-3*c4))*(2*c5-4*c4)/(tmpD*tmpD);
      
  tmpD     =  (2*c4*c5 + x4*(2*c5-4*c4));     
  double dydx4 =  (c5-3*c4)/(tmpD)
          -(2*c4*c5 + x4*(c5-3*c4))*(2*c5-4*c4)
          /(tmpD*tmpD);

  double dydx5 = -minShoulderSlopeMagnitude;


  RigidBodyDynamics::Math::MatrixNd p0 =
    SegmentedQuinticBezierToolkit::
       calcQuinticBezierCornerControlPoints(x5,y5,dydx5,
                          x4,y4,dydx4,
                          c);
  RigidBodyDynamics::Math::MatrixNd p1 =
    SegmentedQuinticBezierToolkit::
       calcQuinticBezierCornerControlPoints(x4,y4,dydx4,
                          x3,y3,dydx3,
                          c);
  RigidBodyDynamics::Math::MatrixNd p2 =
    SegmentedQuinticBezierToolkit::
       calcQuinticBezierCornerControlPoints(x3,y3,dydx3,
                          x2,y2,dydx2,c);
  RigidBodyDynamics::Math::MatrixNd p3 =
    SegmentedQuinticBezierToolkit::
       calcQuinticBezierCornerControlPoints(x2,y2,dydx2,
                          x1,y1,dydx1,c);
  RigidBodyDynamics::Math::MatrixNd p4 =
    SegmentedQuinticBezierToolkit::
       calcQuinticBezierCornerControlPoints(x1,y1,dydx1,
                          x0,y0,dydx0,c);


  RigidBodyDynamics::Math::MatrixNd mX(6,5);
  RigidBodyDynamics::Math::MatrixNd mY(6,5);

  mX.col(0) = p0.col(0);
  mY.col(0) = p0.col(1);
  mX.col(1) = p1.col(0);
  mY.col(1) = p1.col(1);
  mX.col(2) = p2.col(0);
  mY.col(2) = p2.col(1);
  mX.col(3) = p3.col(0);
  mY.col(3) = p3.col(1);
  mX.col(4) = p4.col(0);
  mY.col(4) = p4.col(1);

  smoothSegmentedFunctionToUpdate.updSmoothSegmentedFunction(
              mX, mY, x5, x0, y5, y0, dydx5, dydx0, curveName);
}
//=============================================================================
// ANDERSON 2007 Passive Torque Angle Curve
//=============================================================================
void SmoothSegmentedFunctionFactory:: 
  createAnderson2007PassiveTorqueAngleCurve(
    double scale,
    double c1,
    double b1,
    double k1,
    double b2,
    double k2,                    
    const std::string& curveName,
    SmoothSegmentedFunction& smoothSegmentedFunctionToUpdate)
{

  if( !(scale > 0) ){ 
    cerr  << "SmoothSegmentedFunctionFactory::"
          << "createAnderson2007PassiveTorqueAngleCurve " 
          << curveName
          << ": scale must be greater than 0"
          << endl;
    assert(0);
    abort();         
  }

  if( !(c1 > 0) ) { 
    cerr  << "SmoothSegmentedFunctionFactory::"
          << "createAnderson2007PassiveTorqueAngleCurve " 
          << curveName
          << ": c1 must be greater than 0"
          << endl;
    assert(0);
    abort();          
  }

  //Advanced settings that we'll hide for now
  double curviness = 0.75;
  double c = SegmentedQuinticBezierToolkit::scaleCurviness(curviness);
  double minShoulderSlopeMagnitude = 0;
  
  //Zero out the coefficients associated with a 
  //the side of the curve that goes negative.
  bool flag_oneSided = true;

  if(flag_oneSided){
    if(fabs(b1) > 0){
      if(b1 > 0){
        b2 = 0;
        k2 = 0;
      }else{
        b1 = 0; 
        k1 = 0;
      }

    }else if(fabs(b2) > 0){
      if(b2 > 0){
        b1 = 0; 
        k1 = 0;   
      }else{
        b2 = 0;
        k2 = 0;
      }
    }   
  }
  //Divide up the curve into a left half
  //and a right half, rather than 1 and 2. 
  //Why? These two different halves require different
  //   Bezier curves.

  double c1Scale = c1*scale;
  double thL    = 0.; //left
  double thR    = 0.; //right
  double DtauDthL = 0.;
  double DtauDthR = 0.;
  double bL     = 0.;
  double kL     = 0.;
  double bR     = 0.;
  double kR     = 0.;

  int curveType   = 0; //flat curve
  int flag_thL  = 0;
  int flag_thR  = 0;

  if(fabs(k1)>0 && fabs(b1)>0){
    //The value of theta where the passive force generated by the
    //muscle is equal to 1 maximum isometric contraction.
    thL     = (1/k1)*log(fabs( c1Scale/b1 ));
    DtauDthL  = b1*k1*exp(thL*k1);      
    bL      = b1;
    kL      = k1;
    flag_thL  = 1;
  }
  
  if(fabs(k2)>0 && fabs(b2)>0){
    //The value of theta where the passive force generated by the
    //muscle is equal to 1 maximum isometric contraction.
    thR     = (1/k2)*log(fabs( c1Scale/b2 ));
    DtauDthR  = b2*k2*exp(thR*k2);  
    bR      = b2;
    kR      = k2;   
    flag_thR  = 1;
  }

  //A 'left' curve has a negative slope,
  //A 'right' curve has a positive slope.
  if(DtauDthL > DtauDthR){
    double tmpD = thL;
    thL   = thR;
    thR   = tmpD;

    tmpD  =  bR;
    bR    = bL;
    bL    = tmpD;

    tmpD  = kR;
    kR    = kL;
    kL    = tmpD;

    tmpD    = DtauDthL;
    DtauDthL  = DtauDthR;
    DtauDthR  = tmpD;

    int tmpI = flag_thL;
    flag_thL = flag_thR;
    flag_thR = tmpI;    
  }


  if(flag_thL){
    curveType   = curveType + 1;
  }
  if(flag_thR){
    curveType   = curveType + 2;
  }

  RigidBodyDynamics::Math::MatrixNd mX(6,1);
  RigidBodyDynamics::Math::MatrixNd mY(6,1);

  double xStart  = 0;
  double xEnd    = 0;
  double yStart  = 0;
  double yEnd    = 0;
  double dydxStart = 0;
  double dydxEnd   = 0;


  switch(curveType){

    case 0:
      {//No curve - it's a flat line
        RigidBodyDynamics::Math::MatrixNd p0 = 
          SegmentedQuinticBezierToolkit::
             calcQuinticBezierCornerControlPoints(0.,0.,0.,
                                1.,0.,0.,c); 

        mX.col(0) = p0.col(0);
        mY.col(0) = p0.col(1);

      }break;
    case 1:
      {        
        //Get a point on the curve that is close to 0.
        double x1  = (1/kL)*log(fabs(0.01*c1Scale/bL) );
        double y1  = bL*exp(kL*x1);
        double dydx1 = bL*kL*exp(kL*x1);

        //Get a point that is at 1 maximum isometric torque
        double x3  = thL;
        double y3  = bL*exp(kL*x3);
        double dydx3 = bL*kL*exp(kL*x3);
           
        //Get a mid-point         
        double x2  = 0.5*(x1+x3);
        double y2  = bL*exp(kL*x2);
        double dydx2 = bL*kL*exp(kL*x2);
        
        //Past the crossing point of the linear extrapolation
        double x0   = x1 - 2*y1/dydx1;
        double y0   = 0;
        double dydx0  = minShoulderSlopeMagnitude*copysign(1.0,dydx1);
        
        xStart    = x3;
        xEnd    = x0;
        yStart    = y3;
        yEnd    = y0;
        dydxStart   = dydx3;
        dydxEnd   = dydx0;

        RigidBodyDynamics::Math::MatrixNd p0 = 
          SegmentedQuinticBezierToolkit::
             calcQuinticBezierCornerControlPoints(x3,y3,dydx3,
                                x2,y2,dydx2,c); 
        RigidBodyDynamics::Math::MatrixNd p1 = 
          SegmentedQuinticBezierToolkit::
             calcQuinticBezierCornerControlPoints(x2,y2,dydx2,
                                x1,y1,dydx1,c); 
        RigidBodyDynamics::Math::MatrixNd p2 = 
          SegmentedQuinticBezierToolkit::
             calcQuinticBezierCornerControlPoints(x1,y1,dydx1,
                                x0,y0,dydx0,c);

        mX.resize(6,3);
        mY.resize(6,3);

        mX.col(0) = p0.col(0);
        mY.col(0) = p0.col(1);     
        mX.col(1) = p1.col(0);
        mY.col(1) = p1.col(1); 
        mX.col(2) = p2.col(0);
        mY.col(2) = p2.col(1); 

      }break;
    case 2:
      {
        //Get a point on the curve that is close to 0.
        double x1  = (1/kR)*log(fabs(0.01*c1Scale/bR) );
        double y1  = bR*exp(kR*x1);
        double dydx1 = bR*kR*exp(kR*x1);

        //Go just past the crossing point of the linear extrapolation
        double x0   = x1 - 2*y1/dydx1;
        double y0   = 0;
        double dydx0 = minShoulderSlopeMagnitude*copysign(1.0,dydx1);

        //Get a point close to 1 maximum isometric torque
        double x3  = thR;
        double y3  = bR*exp(kR*x3);
        double dydx3 = bR*kR*exp(kR*x3);

        //Get a mid point.
        double x2   = 0.5*(x1+x3);
        double y2   = bR*exp(kR*x2);
        double dydx2  = bR*kR*exp(kR*x2);    

        xStart    = x0;
        xEnd    = x3;
        yStart    = y0;
        yEnd    = y3;
        dydxStart   = dydx0;
        dydxEnd   = dydx3;

        RigidBodyDynamics::Math::MatrixNd p0 = 
          SegmentedQuinticBezierToolkit::
             calcQuinticBezierCornerControlPoints(x0,y0,dydx0,
                                x1,y1,dydx1,c); 
        RigidBodyDynamics::Math::MatrixNd p1 = 
          SegmentedQuinticBezierToolkit::
             calcQuinticBezierCornerControlPoints(x1,y1,dydx1,
                                x2,y2,dydx2,c); 
        RigidBodyDynamics::Math::MatrixNd p2 = 
          SegmentedQuinticBezierToolkit::
             calcQuinticBezierCornerControlPoints(x2,y2,dydx2,
                                x3,y3,dydx3,c);
        mX.resize(6,3);
        mY.resize(6,3);

        mX.col(0) = p0.col(0);
        mY.col(0) = p0.col(1);     
        mX.col(1) = p1.col(0);
        mY.col(1) = p1.col(1); 
        mX.col(2) = p2.col(0);
        mY.col(2) = p2.col(1); 

      }break;  
    case 3:
      {
        //Only when flag_oneSided = false;
        double x0 = thL;
        double x4 = thR;

        double x2 = 0.5*(x0 + x4);
        double x1 = 0.5*(x0 + x2);
        double x3 = 0.5*(x2 + x4);

        double y0 = b1*exp(k1*x0)
              + b2*exp(k2*x0);
        double y1 = b1*exp(k1*x1)
              + b2*exp(k2*x1);        
        double y2 = b1*exp(k1*x2)
              + b2*exp(k2*x2);
        double y3 = b1*exp(k1*x3)
              + b2*exp(k2*x3);
        double y4 = b1*exp(k1*x4)
              + b2*exp(k2*x4);    

        double dydx0 =   b1*k1*exp(k1*x0)
                 + b2*k2*exp(k2*x0);
        double dydx1 =   b1*k1*exp(k1*x1)
                 + b2*k2*exp(k2*x1);
        double dydx2 =   b1*k1*exp(k1*x2)
                 + b2*k2*exp(k2*x2);
        double dydx3 =   b1*k1*exp(k1*x3)
                 + b2*k2*exp(k2*x3);
        double dydx4 =   b1*k1*exp(k1*x4)
                 + b2*k2*exp(k2*x4);  

        xStart    = x0;
        xEnd    = x4;
        yStart    = y0;
        yEnd    = y4;
        dydxStart   = dydx0;
        dydxEnd   = dydx4;                   

        RigidBodyDynamics::Math::MatrixNd p0 = 
          SegmentedQuinticBezierToolkit::
             calcQuinticBezierCornerControlPoints(x0,y0,dydx0,
                                x1,y1,dydx1,c); 
        RigidBodyDynamics::Math::MatrixNd p1 = 
          SegmentedQuinticBezierToolkit::
             calcQuinticBezierCornerControlPoints(x1,y1,dydx1,
                                x2,y2,dydx2,c); 
        RigidBodyDynamics::Math::MatrixNd p2 = 
          SegmentedQuinticBezierToolkit::
             calcQuinticBezierCornerControlPoints(x2,y2,dydx2,
                                x3,y3,dydx3,c);
        RigidBodyDynamics::Math::MatrixNd p3 = 
          SegmentedQuinticBezierToolkit::
             calcQuinticBezierCornerControlPoints(x3,y3,dydx3,
                                x4,y4,dydx4,c);

        mX.resize(6,4);
        mY.resize(6,4);

        mX.col(0) = p0.col(0);
        mY.col(0) = p0.col(1);     
        mX.col(1) = p1.col(0);
        mY.col(1) = p1.col(1); 
        mX.col(2) = p2.col(0);
        mY.col(2) = p2.col(1);    
        mX.col(3) = p3.col(0);
        mY.col(3) = p3.col(1); 

      }break;
    default:
    {
      cerr  << "SmoothSegmentedFunctionFactory::"
            << "createAnderson2007PassiveTorqueAngleCurve " 
            << curveName
            << ": undefined curveType"
            << endl;
    assert(0);
    abort();   
    }

  };

  //Normalize the y values.
  mY = mY*(1/c1Scale);

  smoothSegmentedFunctionToUpdate.updSmoothSegmentedFunction(
    mX,                mY, 
    xStart,            xEnd, 
    yStart/c1Scale,    yEnd/c1Scale,
    dydxStart/c1Scale, dydxEnd/c1Scale,
    curveName);
}

void SmoothSegmentedFunctionFactory::
  createNeuroSwingSpringProfile(
    RigidBodyDynamics::Math::VectorNd& controlPointAngles,
    RigidBodyDynamics::Math::VectorNd& controlPointTorques,
    RigidBodyDynamics::Math::VectorNd& controlPointSlopes,
    const std::string& curveName,
    SmoothSegmentedFunction& smoothSegmentedFunctionToUpdate)
{

  if(   !(  (controlPointAngles.rows()  == 5)
      &&  (controlPointTorques.rows() == 5)
      &&  (controlPointSlopes.rows()  == 5))  ) {
    cerr  << "SmoothSegmentedFunctionFactory::"
          << "createNeuroSwingSpringProfile "
          << curveName
          << ": the vectors controlPointAngles, controlPointTorques, "
          << " and controlPointSlopes must have 5 elements."
          << endl;
    assert(0);
    abort();
  }


  double slopeSign = 1.0;
  double maxSlope = 0.0;
  for(int i=0;i<5;++i){
    if(fabs(maxSlope) < fabs(controlPointSlopes(i))){
      maxSlope = controlPointSlopes(i);
    }
  }
  if(maxSlope >= 0){
    slopeSign = 1.0;
  }else{
    slopeSign = -1;
  }


  for(int i=0; i < 5; ++i ){
      if( (controlPointSlopes(i)*slopeSign < 0) ) {
        cerr  << "SmoothSegmentedFunctionFactory::"
              << "createNeuroSwingSpringProfile "
              << curveName
              << ": all slopes must have the same sign"
              << endl;
        assert(0);
        abort();
      }
  }

  for(int i=1; i < 5; ++i ){
      if( !(controlPointAngles(i) >= controlPointAngles(i-1)) ){
        cerr  << "SmoothSegmentedFunctionFactory::"
              << "createNeuroSwingSpringProfile "
              << curveName
              << ": controlPointAngle elements must be ordered from "
              << " smallest to largest."
              << endl;
        assert(0);
        abort();
      }
  }

  for(int i=1; i < 5; ++i ){
      if( !(controlPointTorques(i)*slopeSign >= controlPointTorques(i-1)*slopeSign) ){
        cerr  << "SmoothSegmentedFunctionFactory::"
              << "createNeuroSwingSpringProfile "
              << curveName
              << ": controlPointTorque elements cannot be decreasing."
              << endl;
        assert(0);
        abort();
      }
  }

  double curviness = 1.0;
  double c = SegmentedQuinticBezierToolkit::scaleCurviness(curviness);

  RigidBodyDynamics::Math::MatrixNd mX(6,4);
  RigidBodyDynamics::Math::MatrixNd mY(6,4);

  for(int i=1; i<5; ++i){
    RigidBodyDynamics::Math::MatrixNd pTmp =
      SegmentedQuinticBezierToolkit::
         calcQuinticBezierCornerControlPoints(
          controlPointAngles(i-1),
          controlPointTorques(i-1),
          controlPointSlopes(i-1),
            controlPointAngles(i),
            controlPointTorques(i),
            controlPointSlopes(i),
            c);
    mX.col(i-1) = pTmp.col(0);
    mY.col(i-1) = pTmp.col(1);
  }

   smoothSegmentedFunctionToUpdate.updSmoothSegmentedFunction(  
                    mX,                     mY,
                    controlPointAngles(0),  controlPointAngles(4),
                    controlPointTorques(0), controlPointTorques(4),
                    controlPointSlopes(0),  controlPointSlopes(4),
                    curveName);
}
