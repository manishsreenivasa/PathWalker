//==============================================================================
/* 
 * RBDL - Rigid Body Dynamics Library: Addon : forceElements
 * Copyright (c) 2016 Matthew Millard <matthew.millard@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#include "NeuroSwingTorsionSpring.h"
#include "../geometry/SegmentedQuinticBezierToolkit.h"
#include "csvtools.h"

#include <limits>

#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>
 
static double EPSILON      = std::numeric_limits<double>::epsilon();
static double SQRTEPSILON  = sqrt(std::numeric_limits<double>::epsilon());
static double ROOT4EPSILON = sqrt(SQRTEPSILON);

using namespace RigidBodyDynamics::Math;
using namespace RigidBodyDynamics::Addons::CustomForces;
using namespace RigidBodyDynamics::Addons::Geometry;

using namespace std;

NeuroSwingTorsionSpring::
    NeuroSwingTorsionSpring( )
        :_jointAngleOffset(1.0),
        _jointAngleSign(1.0),
        _jointTorqueSign(1.0),
        _springName("empty")
{
    
}

NeuroSwingTorsionSpring::NeuroSwingTorsionSpring( 
             double angleHardStopInPlantarFlexion, 
             double angleCenterOffset,
             double angleHardStopInDorsiFlexion,
             double stiffnessEndRangeInPlantarFlexion,
             double stiffnessInPlantarFlexion,
             double stiffnessInDorsiFlexion,       
             double stiffnessEndRangeInDorsiFlexion,
             double anglePreloadWindow,
             double torquePreload,
             double jointAngleOffset,
             double jointAngleSign,                 
             double jointTorqueSign,               
            const std::string& springName)
    :_angleHardStopInPlantarFlexion(angleHardStopInPlantarFlexion),
     _angleCenterOffset(angleCenterOffset),
     _angleHardStopInDorsiFlexion(angleHardStopInDorsiFlexion),
     _stiffnessEndRangeInPlantarFlexion(stiffnessEndRangeInPlantarFlexion),
     _stiffnessInPlantarFlexion(stiffnessInPlantarFlexion),              
     _stiffnessInDorsiFlexion(stiffnessInDorsiFlexion),
     _stiffnessEndRangeInDorsiFlexion(stiffnessEndRangeInDorsiFlexion),
     _anglePreloadWindow(anglePreloadWindow),
     _torquePreload(torquePreload),
     _jointAngleOffset(jointAngleOffset),
     _jointAngleSign(jointAngleSign),
     _jointTorqueSign(jointTorqueSign),
     _springName(springName)
{
    _neuroSwingCurve = SmoothSegmentedFunction();
    _neuroSwingCurve.setName(_springName);
                 
    _parametersAreDirty = true;

    buildNeuroSwingTorsionSpringCurve();
}

void NeuroSwingTorsionSpring::buildNeuroSwingTorsionSpringCurve() 
{

  //Check the new input arguments
  if(!( (_angleHardStopInDorsiFlexion < _angleCenterOffset) 
     && (_angleCenterOffset < _angleHardStopInPlantarFlexion) )){

     cerr << "NeuroSwingTorsionSpring::buildNeuroSwingTorsionSpringCurve "
          << _springName.c_str() 
          << ": Incompatible angle settings. Note that "
          << "angleHardStopInDorsiFlexion < angleCenterOffset "
          << "< angleHardStopInPlantarFlexion" 
          << endl;
      printCurveParametersUsingCerr();
    assert(0);
    abort();
  }

  if(!(     (_stiffnessEndRangeInDorsiFlexion >= 0) 
        &&  (_stiffnessInDorsiFlexion >= 0) 
        &&  (_stiffnessInPlantarFlexion >= 0)
        &&  (_stiffnessEndRangeInPlantarFlexion >= 0) )){

     cerr <<"NeuroSwingTorsionSpring::buildNeuroSwingTorsionSpringCurve "
          << _springName.c_str() 
          <<":Invalid stiffness settings. All stiffness parameters must be >= 0"
          << endl;
     printCurveParametersUsingCerr();
    assert(0);
    abort();
  }

  if(!(     (_torquePreload >= 0) 
        &&  (_anglePreloadWindow >= ROOT4EPSILON) )){

     cerr << "NeuroSwingTorsionSpring::buildNeuroSwingTorsionSpringCurve " 
          << _springName.c_str() 
          << ": Preload setting error. Check that anglePreloadWindow > epsilion^(1/4)"
          <<" which for a double means anglePreloadWindow > "
          << ROOT4EPSILON
          << endl;
    printCurveParametersUsingCerr();
    assert(0);
    abort();
  }


  if(!(_anglePreloadWindow <= 0.25*(_angleHardStopInPlantarFlexion
                                  -_angleHardStopInDorsiFlexion)))
  {
    cerr << "NeuroSwingTorsionSpring::buildNeuroSwingTorsionSpringCurve "
         << _springName.c_str()
         << ": Preload setting error. Check that anglePreloadWindow < "
         << "0.25*(_angleHardStopInPlantarFlexion-_angleHardStopInDorsiFlexion)"
         << "which in this case means anglePreloadWindow < "
         << 0.25*(_angleHardStopInPlantarFlexion -_angleHardStopInDorsiFlexion)
         << endl;
   printCurveParametersUsingCerr();
   assert(0);
   abort();
  }

  if( abs(abs(_jointAngleSign)-1) >  EPSILON){
    std::cerr   << "NeuroSwingTorsionSpring::"
           << "NeuroSwingTorsionSpring:"
           << _springName
           << ": jointAngleSign must be -1 or 1  not "
           << _jointAngleSign;
    printCurveParametersUsingCerr();
    assert(0);
    abort();
  }


  if( abs(abs(_jointTorqueSign)-1) >  EPSILON){
      std::cerr   << "NeuroSwingTorsionSpring::"
             << "NeuroSwingTorsionSpring:"
             << _springName
             << ": jointTorqueSign must be -1 or 1 not "
             << _jointTorqueSign;
      printCurveParametersUsingCerr();
      assert(0);
      abort();
  }

                                
  //Go and constuct the curve
  double dorsiFlexionWindow = fabs(_angleHardStopInDorsiFlexion
                            - (_angleCenterOffset - 0.5*_anglePreloadWindow));


  double plantarFlexionWindow = fabs(_angleHardStopInPlantarFlexion
                             - (_angleCenterOffset + 0.5*_anglePreloadWindow));



  RigidBodyDynamics::Math::VectorNd xPts(7);
  RigidBodyDynamics::Math::VectorNd yPts(7);
  RigidBodyDynamics::Math::VectorNd dydxPts(7);

  //double elbowWidth = max( _anglePreloadWindow*0.5/dorsiFlexionWindow,
  //                         _anglePreloadWindow*0.5/plantarFlexionWindow);

 //if(elbowWidth < 0.05){
  //  elbowWidth = 0.05;
  //}
  double elbowWidth = 0.1;

  double deltaY0, deltaY1 = 0;
  double deltaX0 = dorsiFlexionWindow*elbowWidth;
  double deltaX1 = plantarFlexionWindow*elbowWidth;

  if(   _stiffnessEndRangeInDorsiFlexion > 0){
    deltaY0 = _torquePreload
            + dorsiFlexionWindow*_stiffnessInDorsiFlexion;
    deltaY0 = deltaY0*elbowWidth;
    if(deltaY0 < SQRTEPSILON){
      deltaY0 = SQRTEPSILON;
    }

    deltaX0 = deltaY0/_stiffnessEndRangeInDorsiFlexion;

  }

  if(   _stiffnessEndRangeInPlantarFlexion > 0){
    deltaY1 = _torquePreload
            + plantarFlexionWindow*_stiffnessInPlantarFlexion;
    deltaY1 = deltaY1*elbowWidth;
    if(deltaY1 < SQRTEPSILON){
      deltaY1 = SQRTEPSILON;
    }
    deltaX1 = deltaY1/_stiffnessEndRangeInPlantarFlexion;
  }



  xPts(0)    = _angleHardStopInDorsiFlexion
               - deltaX0;
  xPts(1)    = _angleHardStopInDorsiFlexion               
               + dorsiFlexionWindow*elbowWidth;//_stiffnessInDorsiFlexion;
  xPts(2)    = _angleCenterOffset - _anglePreloadWindow;
  xPts(3)    = _angleCenterOffset;
  xPts(4)    = _angleCenterOffset + _anglePreloadWindow;
  xPts(5)    = _angleHardStopInPlantarFlexion
              - plantarFlexionWindow*elbowWidth;///_stiffnessInPlantarFlexion;
  xPts(6)    = _angleHardStopInPlantarFlexion
               +deltaX1;
               //+ plantarFlexionWindow*elbowWidth/_stiffnessEndRangeInPlantarFlexion;

  dydxPts(0) = _stiffnessEndRangeInDorsiFlexion;
  dydxPts(1) = _stiffnessInDorsiFlexion;
  dydxPts(2) = _stiffnessInDorsiFlexion;
  dydxPts(3) = 2.0*_torquePreload/(0.5*_anglePreloadWindow);
  dydxPts(4) = _stiffnessInPlantarFlexion;
  dydxPts(5) = _stiffnessInPlantarFlexion;
  dydxPts(6) = _stiffnessEndRangeInPlantarFlexion;

 
  yPts(0) = -_torquePreload
            - dorsiFlexionWindow*_stiffnessInDorsiFlexion
            - deltaY0;

  yPts(1) = -_torquePreload
            - dorsiFlexionWindow*(1-elbowWidth)*_stiffnessInDorsiFlexion;

  yPts(2) = -_torquePreload - 0.5*_anglePreloadWindow*_stiffnessInDorsiFlexion;
  yPts(3) = 0;
  yPts(4) = _torquePreload + 0.5*_anglePreloadWindow*_stiffnessInPlantarFlexion;

  yPts(5) = _torquePreload
            + plantarFlexionWindow*(1-elbowWidth)*_stiffnessInPlantarFlexion;

  yPts(6) =  _torquePreload
            + _stiffnessInPlantarFlexion*plantarFlexionWindow
            + deltaY1;
  
  RigidBodyDynamics::Math::MatrixNd mX(6,6);
  RigidBodyDynamics::Math::MatrixNd mY(6,6);
  RigidBodyDynamics::Math::MatrixNd xyPts(6,2);

  double curviness        = 0.75;
  double scaledCurviness  = SegmentedQuinticBezierToolkit::
                              scaleCurviness(curviness);


  for(int i=1; i < 7; ++i){
    xyPts = SegmentedQuinticBezierToolkit::
              calcQuinticBezierCornerControlPoints(
                  xPts(i-1),yPts(i-1),dydxPts(i-1),
                  xPts(i)  ,yPts(i)  ,dydxPts(i)  , 
                  scaledCurviness);

    mX.col(i-1) = xyPts.col(0);
    mY.col(i-1) = xyPts.col(1);

  }

  _neuroSwingCurve.updSmoothSegmentedFunction(mX,     mY,
                                              xPts(0),xPts(6),
                                              yPts(0),yPts(6),
                                              dydxPts(0),dydxPts(6),
                                              _springName);

  _parametersAreDirty = false;

}

const RigidBodyDynamics::Addons::Geometry::SmoothSegmentedFunction&
NeuroSwingTorsionSpring::getSmoothSegmentedFunction() const
{
  return _neuroSwingCurve;
}

void NeuroSwingTorsionSpring::printCurveParametersUsingCerr()
{
  std::cerr
  << "_angleHardStopInDorsiFlexion      :"
  << _angleHardStopInDorsiFlexion      << endl
  << "_angleCenterOffset                :"
  << _angleCenterOffset                << endl
  << "_anglePreloadWindow               :"
  << _anglePreloadWindow               << endl
  << "_angleHardStopInPlantarFlexion    :"
  << _angleHardStopInPlantarFlexion    << endl
  << "_stiffnessEndRangeInDorsiFlexion  :"
  << _stiffnessEndRangeInDorsiFlexion  << endl
  << "_stiffnessInDorsiFlexion          :"
  << _stiffnessInDorsiFlexion          << endl
  << "_torquePreload                    :"
  << _torquePreload                    << endl
  << "_stiffnessInPlantarFlexion        :"
  << _stiffnessInPlantarFlexion        << endl
  << "_stiffnessEndRangeInPlantarFlexion:"
  << _stiffnessEndRangeInPlantarFlexion<< endl
  << "_jointAngleOffset                 :"
  << _jointAngleOffset                 << endl
  << "_jointAngleSign                   :"
  << _jointAngleSign                   << endl
  << "_jointTorqueSign                  :"
  << _jointTorqueSign                  << endl;
}

double NeuroSwingTorsionSpring::calcNeuroSwingAngle(double jointAngle) const{
  return _jointAngleSign*(jointAngle - _jointAngleOffset);
}

double NeuroSwingTorsionSpring::calcJointAngle(double neuroSwingAngle) const{
  return neuroSwingAngle*_jointAngleSign + _jointAngleOffset;
}

double NeuroSwingTorsionSpring::calcJointTorque(                        
        double jointAngle) const
{

  //Here we are doing a lazy update of the curve because
  //if you are changing several parameters you can have
  //intermediate settings for the curve that are illeagle.
  //To maintain the const status of this function we're using
  //a const_cast
  if(_parametersAreDirty == true){
    NeuroSwingTorsionSpring* mutableThis =
        const_cast<NeuroSwingTorsionSpring*>(this);
    mutableThis->buildNeuroSwingTorsionSpringCurve();
  }

  double neuroAngle =  calcNeuroSwingAngle(jointAngle);
  return _neuroSwingCurve.calcValue(neuroAngle)*_jointTorqueSign;

}

double NeuroSwingTorsionSpring::calcJointStiffness(
        double jointAngle) const
{
  //See comment on calcJointTorque regarding the const_cast
  if(_parametersAreDirty == true){
    NeuroSwingTorsionSpring* mutableThis =
        const_cast<NeuroSwingTorsionSpring*>(this);
    mutableThis->buildNeuroSwingTorsionSpringCurve();
  }

  double neuroAngle     =  calcNeuroSwingAngle(jointAngle);
  double jointStiffness =  _neuroSwingCurve.calcDerivative(neuroAngle,1) 
                            * (_jointTorqueSign*_jointAngleSign);
  return jointStiffness;                            
}


void NeuroSwingTorsionSpring::printJointTorqueProfileToFile(
        const std::string& path,
        const std::string& fileNameWithoutExtension) const
{
  //See comment on calcJointTorque regarding the const_cast
  if(_parametersAreDirty == true){
    NeuroSwingTorsionSpring* mutableThis =
        const_cast<NeuroSwingTorsionSpring*>(this);
    mutableThis->buildNeuroSwingTorsionSpringCurve();
  }

  RigidBodyDynamics::Math::VectorNd xRange(2);
  xRange          = _neuroSwingCurve.getCurveDomain();
  double range    = (xRange(1)-xRange(0));
  double angle0   = (xRange(0)-0.1*range);
  double angle1   = (xRange(1)+0.1*range);



  RigidBodyDynamics::Math::MatrixNd curveSample
      = _neuroSwingCurve.calcSampledCurve(1,angle0,angle1);


  std::vector< std::vector < double > > matrix;
  std::vector < double > row(6);
  std::string header("angle,"
                     "torque,"
                     "stiffness,"
                     "jointAngle,"
                     "jointTorque,"
                     "jointStiffness");
  double angle;
  double torque;
  double stiffness;
  double jointAngle;
  double jointTorque;
  double jointStiffness;


  for(int i=0; i<curveSample.rows(); ++i){
    angle       = curveSample(i,0);
    torque      = curveSample(i,1);
    stiffness   = curveSample(i,2);

    jointAngle      = calcJointAngle(angle);
    jointTorque     = calcJointTorque(jointAngle);
    jointStiffness  = calcJointStiffness(jointAngle);

    row.at(0) = angle;
    row.at(1) = torque;
    row.at(2) = stiffness;
    row.at(3) = jointAngle;
    row.at(4) = jointTorque;
    row.at(5) = jointStiffness;
    matrix.push_back(row);
  }

  std::string fullFilePath = path;
  if(!path.empty()){
      fullFilePath.append("/");
  }
  fullFilePath.append(fileNameWithoutExtension);
  fullFilePath.append(".csv");
  printMatrixToFile(matrix,header,fullFilePath);

  matrix.clear();

  //Print the control points.

  std::vector< double > row2(2);
  RigidBodyDynamics::Math::MatrixNd ctrlX(1,1);
  RigidBodyDynamics::Math::MatrixNd ctrlY(1,1);

  _neuroSwingCurve.getXControlPoints(ctrlX);
  _neuroSwingCurve.getYControlPoints(ctrlY);

  for(int i=0; i<ctrlX.rows();++i){
    for(int j=0; j<ctrlX.cols();++j){
      row2.at(0) = ctrlX(i,j);
      row2.at(1) = ctrlY(i,j);
      matrix.push_back(row2);
    }
  }

  fullFilePath = path;
  if(!path.empty()){
      fullFilePath.append("/");
  }
  fullFilePath.append(fileNameWithoutExtension);
  fullFilePath.append("ControlPoints.csv");
  header = "x,y";
  printMatrixToFile(matrix,header,fullFilePath);
}

std::string NeuroSwingTorsionSpring::getName() const
{
  return _springName;
}
void NeuroSwingTorsionSpring::setName(std::string& name)
{
  _neuroSwingCurve.setName(name);
  _springName = name;

}

double NeuroSwingTorsionSpring::getAngleHardStopInDorsiFlexion() const
{
  return _angleHardStopInDorsiFlexion;
}

double NeuroSwingTorsionSpring::getAngleCenterOffset()  const
{
  return _angleCenterOffset;
}               
double NeuroSwingTorsionSpring::getAnglePreloadWindow() const
{
  return _anglePreloadWindow;
}           

double NeuroSwingTorsionSpring::getAngleHardStopInPlantarFlexion() const
{
  return _angleHardStopInPlantarFlexion;
} 

double NeuroSwingTorsionSpring::getStiffnessEndRangeInDorsiFlexion() const
{
  return _stiffnessEndRangeInDorsiFlexion;
}

double NeuroSwingTorsionSpring::getStiffnessInDorsiFlexion() const
{
  return _stiffnessInDorsiFlexion;
}  

double NeuroSwingTorsionSpring::getTorquePreload() const
{
  return _torquePreload;
}                  
double NeuroSwingTorsionSpring::getStiffnessInPlantarFlexion() const
{
  return _stiffnessInPlantarFlexion;
}      
double NeuroSwingTorsionSpring::getStiffnessEndRangeInPlantarFlexion() const
{
  return _stiffnessEndRangeInPlantarFlexion;
}    
double NeuroSwingTorsionSpring::getJointAngleOffset() const
{
  return _jointAngleOffset;
}              

double NeuroSwingTorsionSpring::getJointAngleSign() const
{
  return _jointAngleSign;
}                 


double NeuroSwingTorsionSpring::getJointTorqueSign() const
{
  return _jointTorqueSign;
}

void NeuroSwingTorsionSpring::setAngleHardStopInDorsiFlexion(double val)
{
  _angleHardStopInDorsiFlexion = val;
  _parametersAreDirty = true;

}

void NeuroSwingTorsionSpring::setAngleCenterOffset(double val)
{
   _angleCenterOffset = val;
   _parametersAreDirty = true;
}

void NeuroSwingTorsionSpring::setAnglePreloadWindow(double val)
{
  _anglePreloadWindow = val;
  _parametersAreDirty = true;
}

void NeuroSwingTorsionSpring::setAngleHardStopInPlantarFlexion(double val)
{
  _angleHardStopInPlantarFlexion = val;
  _parametersAreDirty = true;
}

void NeuroSwingTorsionSpring::setStiffnessEndRangeInDorsiFlexion(double val)
{
  _stiffnessEndRangeInDorsiFlexion = val;
  _parametersAreDirty = true;
}

void NeuroSwingTorsionSpring::setStiffnessInDorsiFlexion(double val)
{
  _stiffnessInDorsiFlexion = val;
  _parametersAreDirty = true;
}

void NeuroSwingTorsionSpring::setTorquePreload(double val)
{
  _torquePreload = val;
  _parametersAreDirty = true;
}

void NeuroSwingTorsionSpring::setStiffnessInPlantarFlexion(double val)
{
  _stiffnessInPlantarFlexion = val;
  _parametersAreDirty = true;
}

void NeuroSwingTorsionSpring::setStiffnessEndRangeInPlantarFlexion(double val)
{
  _stiffnessEndRangeInPlantarFlexion = val;
  _parametersAreDirty = true;
}

void NeuroSwingTorsionSpring::setJointAngleOffset(double val)
{
  _jointAngleOffset = val;
  _parametersAreDirty = true;
}

void NeuroSwingTorsionSpring::setJointAngleSign(double val)
{
  _jointAngleSign = val;
  _parametersAreDirty = true;
}

void NeuroSwingTorsionSpring::setJointTorqueSign(double val)
{
  _jointTorqueSign = val;
  _parametersAreDirty = true;
}
