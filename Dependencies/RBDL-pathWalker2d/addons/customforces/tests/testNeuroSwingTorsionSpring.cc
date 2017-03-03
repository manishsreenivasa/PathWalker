/*                                                                             *
 * testNeuroSwingSpring
 * Copyright (c) 2016 Matthew Millard <matthew.millard@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */


//==============================================================================
// INCLUDES
//==============================================================================



#include "../NeuroSwingTorsionSpring.h"
#include "../csvtools.h"

#include <UnitTest++.h>
#include <rbdl/rbdl_math.h>
#include <ctime>
#include <string>
#include <ostream>
#include <sstream>
#include <stdio.h>
#include <exception>
#include <cassert>
#include <vector>

static double EPSILON     = std::numeric_limits<double>::epsilon();
static double SQRTEPSILON = sqrt(std::numeric_limits<double>::epsilon());
static double TOL         = std::numeric_limits<double>::epsilon()*1e7;

using namespace RigidBodyDynamics::Addons::CustomForces;
using namespace RigidBodyDynamics::Addons::Geometry;
using namespace std;



TEST(NeuroSwingCalcJointTorqueCorrectnessTests){

  double angleHardStopInPlantarFlexion  =  9.599311e-02;
  double angleCenterOffset              = -1.745329e-02;
  double angleHardStopInDorsiFlexion    = -1.308997e-01;

  double stiffnessEndRangeInPlantarFlexion = 100;
  double stiffnessInPlantarFlexion      =  10;//15;
  double stiffnessInDorsiFlexion        =  55;//35;
  double stiffnessEndRangeInDorsiFlexion=  100;

  double anglePreloadWindow             =  1e-2;
  double torquePreload                  =  1.0;


  string springName("test");

  NeuroSwingTorsionSpring leftAFO =
      NeuroSwingTorsionSpring(angleHardStopInPlantarFlexion,
                              angleCenterOffset,
                              angleHardStopInDorsiFlexion,
                              stiffnessEndRangeInPlantarFlexion,
                              stiffnessInPlantarFlexion,
                              stiffnessInDorsiFlexion,
                              stiffnessEndRangeInDorsiFlexion,
                              anglePreloadWindow,
                              torquePreload,
                              0.0,
                              1.0,
                              1.0,
                              springName);

  leftAFO.printJointTorqueProfileToFile("","NeuroSwingProfile");

  const SmoothSegmentedFunction& ssf0 = leftAFO.getSmoothSegmentedFunction();
  RigidBodyDynamics::Math::MatrixNd ssfXpts(1,1);
  ssf0.getXControlPoints(ssfXpts);

  //Go and calculate all of the key point locations and
  //then test to make sure tMatrhey have the proper values

  RigidBodyDynamics::Math::VectorNd x(7);

  x(0) = ssfXpts(0,0);
  x(1) = ssfXpts(1,0);
  x(2) = ssfXpts(2,0);
  x(3) = ssfXpts(3,0);
  x(4) = ssfXpts(4,0);
  x(5) = ssfXpts(5,0);
  x(6) = ssfXpts(5,5);

  double err = 0;

  //==================================================================
  // Check for calculation correctness
  //==================================================================
  //torquePreload
  double test_torquePreload =
      leftAFO.calcJointTorque(x(2));
  err = test_torquePreload
      + stiffnessInDorsiFlexion*0.5*anglePreloadWindow
      + torquePreload;
  CHECK(fabs(test_torquePreload
             + stiffnessInDorsiFlexion*0.5*anglePreloadWindow
             + torquePreload)    < TOL);

  test_torquePreload =
      leftAFO.calcJointTorque(x(4));
  err = test_torquePreload
      - stiffnessInPlantarFlexion*0.5*anglePreloadWindow
      - torquePreload;
  CHECK(fabs(test_torquePreload
             - stiffnessInPlantarFlexion*0.5*anglePreloadWindow
             - torquePreload)    < TOL);

  //Now test that the curve has the correct updated numerical values
  //stiffnessEndRangeInDorsiFlexion
  double test_stiffnessEndRangeInDorsiFlexion = leftAFO.calcJointStiffness(x(0));
  err = test_stiffnessEndRangeInDorsiFlexion -stiffnessEndRangeInDorsiFlexion;
  CHECK(fabs(test_stiffnessEndRangeInDorsiFlexion
             -stiffnessEndRangeInDorsiFlexion)    < TOL);

  //stiffnessInDorsiFlexion
  double test_stiffnessInDorsiFlexion = leftAFO.calcJointStiffness(x(1));
  err = test_stiffnessInDorsiFlexion - stiffnessInDorsiFlexion;
  CHECK(fabs(test_stiffnessInDorsiFlexion
             -stiffnessInDorsiFlexion)    < TOL);



  //stiffnessInPlantarFlexion
  double test_stiffnessInPlantarFlexion = leftAFO.calcJointStiffness(x(4));
  err = test_stiffnessInPlantarFlexion -stiffnessInPlantarFlexion;
  CHECK(fabs(test_stiffnessInPlantarFlexion
             -stiffnessInPlantarFlexion)    < TOL);

  //stiffnessEndRangeInPlantarFlexion
  double test_stiffnessEndRangeInPlantarFlexion = leftAFO.calcJointStiffness(x(6));
  err = test_stiffnessEndRangeInPlantarFlexion
       -stiffnessEndRangeInPlantarFlexion;
  CHECK(fabs(test_stiffnessEndRangeInPlantarFlexion
             -stiffnessEndRangeInPlantarFlexion)    < TOL);


  double tqCenter = leftAFO.calcJointTorque(x(3));
  CHECK(tqCenter < TOL);

  err = test_stiffnessEndRangeInDorsiFlexion
      -stiffnessEndRangeInDorsiFlexion;
  CHECK(fabs(test_stiffnessEndRangeInDorsiFlexion
                 -stiffnessEndRangeInDorsiFlexion)    < TOL);

  CHECK(fabs(test_stiffnessInDorsiFlexion
                 -stiffnessInDorsiFlexion)            < TOL);
  CHECK(fabs(test_stiffnessInPlantarFlexion
                 -stiffnessInPlantarFlexion)          < TOL);

  err = test_stiffnessEndRangeInPlantarFlexion
      -stiffnessEndRangeInPlantarFlexion;
  CHECK(fabs(test_stiffnessEndRangeInPlantarFlexion
                 -stiffnessEndRangeInPlantarFlexion)  < TOL);

  CHECK(fabs(test_stiffnessInDorsiFlexion
             -leftAFO.calcJointStiffness(x(2)))       < TOL);
  CHECK(fabs(test_stiffnessInPlantarFlexion
             -leftAFO.calcJointStiffness(x(4)))       < TOL);


  //==================================================================
  // Check that the get and set functions work
  //==================================================================
  angleHardStopInDorsiFlexion    = leftAFO.getAngleHardStopInDorsiFlexion()*0.5;
  leftAFO.setAngleHardStopInDorsiFlexion(angleHardStopInDorsiFlexion);
  CHECK(fabs(leftAFO.getAngleHardStopInDorsiFlexion()
             - angleHardStopInDorsiFlexion)    < TOL);

  angleCenterOffset    = leftAFO.getAngleCenterOffset()*0.5;
  leftAFO.setAngleCenterOffset(angleCenterOffset);
  CHECK(fabs(leftAFO.getAngleCenterOffset()
             - angleCenterOffset)    < TOL);

  anglePreloadWindow    = leftAFO.getAnglePreloadWindow()*0.5;
  leftAFO.setAnglePreloadWindow(anglePreloadWindow);
  CHECK(fabs(leftAFO.getAnglePreloadWindow()
             - anglePreloadWindow)    < TOL);

  angleHardStopInPlantarFlexion    =
      leftAFO.getAngleHardStopInPlantarFlexion()*0.5;
  leftAFO.setAngleHardStopInPlantarFlexion(angleHardStopInPlantarFlexion);
  CHECK(fabs(leftAFO.getAngleHardStopInPlantarFlexion()
             - angleHardStopInPlantarFlexion)    < TOL);

  stiffnessEndRangeInDorsiFlexion =
      leftAFO.getStiffnessEndRangeInDorsiFlexion()*0.5;
  leftAFO.setStiffnessEndRangeInDorsiFlexion(
        stiffnessEndRangeInDorsiFlexion);
  CHECK(fabs(leftAFO.getStiffnessEndRangeInDorsiFlexion()
             - stiffnessEndRangeInDorsiFlexion)    < TOL);

  stiffnessInDorsiFlexion =
      leftAFO.getStiffnessInDorsiFlexion()*0.5;
  leftAFO.setStiffnessInDorsiFlexion(
        stiffnessInDorsiFlexion);
  CHECK(fabs(leftAFO.getStiffnessInDorsiFlexion()
             - stiffnessInDorsiFlexion)    < TOL);

  stiffnessInPlantarFlexion =
      leftAFO.getStiffnessInPlantarFlexion()*0.5;
  leftAFO.setStiffnessInPlantarFlexion(
        stiffnessInPlantarFlexion);
  CHECK(fabs(leftAFO.getStiffnessInPlantarFlexion()
             - stiffnessInPlantarFlexion)    < TOL);

  stiffnessEndRangeInPlantarFlexion =
      leftAFO.getStiffnessEndRangeInPlantarFlexion()*0.5;
  leftAFO.setStiffnessEndRangeInPlantarFlexion(
        stiffnessEndRangeInPlantarFlexion);
  CHECK(fabs(leftAFO.getStiffnessEndRangeInPlantarFlexion()
             - stiffnessEndRangeInPlantarFlexion)    < TOL);


  torquePreload = leftAFO.getTorquePreload()*0.5;
  leftAFO.setTorquePreload(torquePreload);
  CHECK(fabs(leftAFO.getTorquePreload() -torquePreload)    < TOL);


  //Now make sure that a legal curve can be created for a wide
  //range of legal parameter values

  double stiffnessEndRangeInPlantarFlexion_min = 0+SQRTEPSILON;
  double stiffnessInPlantarFlexion_min = 0+SQRTEPSILON;
  double stiffnessInDorsiFlexion_min = 0+SQRTEPSILON;
  double stiffnessEndRangeInDorsiFlexion_min = 0+SQRTEPSILON;
  double anglePreloadWindow_min = sqrt(SQRTEPSILON);
  double torquePreload_min = 0+SQRTEPSILON;

  double stiffnessEndRangeInPlantarFlexion_max = 100-SQRTEPSILON;
  double stiffnessInPlantarFlexion_max = 50-SQRTEPSILON;
  double stiffnessInDorsiFlexion_max = 50-SQRTEPSILON;
  double stiffnessEndRangeInDorsiFlexion_max = 100-SQRTEPSILON;
  double anglePreloadWindow_max = 0.25*(angleHardStopInPlantarFlexion
                                   -angleHardStopInDorsiFlexion)-SQRTEPSILON;
  double torquePreload_max = 10-SQRTEPSILON;

  int steps = 3;


  for(int i0=0;i0<steps;++i0){
      stiffnessEndRangeInPlantarFlexion =
          stiffnessEndRangeInPlantarFlexion_min
          + i0*(stiffnessEndRangeInPlantarFlexion_max
                -stiffnessEndRangeInPlantarFlexion_min)/((double)steps-1.0);

    for(int i1=0;i1<steps;++i1){
      stiffnessInPlantarFlexion =
          stiffnessInPlantarFlexion_min
          + i1*(stiffnessInPlantarFlexion_max
                -stiffnessInPlantarFlexion_min)/((double)steps-1.0);

      for(int i2=0;i2<steps;++i2){
        stiffnessEndRangeInDorsiFlexion =
            stiffnessEndRangeInDorsiFlexion_min
            + i2*(stiffnessEndRangeInDorsiFlexion_max
                  -stiffnessEndRangeInDorsiFlexion_min)/((double)steps-1.0);

        for(int i3=0;i3<steps;++i3){
          stiffnessInDorsiFlexion =
              stiffnessInDorsiFlexion_min
              + i3*(stiffnessInDorsiFlexion_max
                    -stiffnessInDorsiFlexion_min)/((double)steps-1.0);

          for(int i4=0;i4<steps;++i4){
            anglePreloadWindow = anglePreloadWindow_min
                + i4*(anglePreloadWindow_max-anglePreloadWindow_min)/((double)steps-1.0);

            for(int i5=0;i5<steps;++i5){
              torquePreload = torquePreload_min
                  + i5*(torquePreload_max-torquePreload_min)/((double)steps-1.0);

              NeuroSwingTorsionSpring tempAFO =
                  NeuroSwingTorsionSpring(angleHardStopInPlantarFlexion,
                                          angleCenterOffset,
                                          angleHardStopInDorsiFlexion,
                                          stiffnessEndRangeInPlantarFlexion,
                                          stiffnessInPlantarFlexion,
                                          stiffnessInDorsiFlexion,
                                          stiffnessEndRangeInDorsiFlexion,
                                          anglePreloadWindow,
                                          torquePreload,
                                          0.0,
                                          1.0,
                                          1.0,
                                          springName);
            }

          }

        }

      }

    }

  }




}






