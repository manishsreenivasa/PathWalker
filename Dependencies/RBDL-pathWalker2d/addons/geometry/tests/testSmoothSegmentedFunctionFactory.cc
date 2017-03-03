/* -------------------------------------------------------------------------- *
 *        OpenSim:  testSmoothSegmentedFunctionFactory.cpp                    *
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
    
     This also includes additional curves (the Anderson2007 curves)
     which are not presently in OpenSim.

    Author:
     Matthew Millard
    
    Date:
     Nov 2015

*/
/* 
Below is a basic bench mark simulation for the SmoothSegmentedFunctionFactory
class, a class that enables the easy generation of C2 continuous curves 
that define the various characteristic curves required in a muscle model
 */

// Author:  Matthew Millard

//==============================================================================
// INCLUDES
//==============================================================================


//#include <OpenSim/Common/Exception.h>
//#include <OpenSim/Common/SegmentedQuinticBezierToolkit.h>
//#include <OpenSim/Common/SmoothSegmentedFunctionFactory.h>
//#include "SegmentedQuinticBezierToolkit.h"
//#include "SmoothSegmentedFunctionFactory.h"

#include "geometry.h"


//#include <SimTKsimbody.h>
#include <UnitTest++.h>
#include <rbdl/rbdl_math.h>
#include <ctime>
#include <string>
#include <stdio.h>
#include <exception>
#include <cassert>

using namespace RigidBodyDynamics::Addons::Geometry;

using namespace std;
//using namespace OpenSim;
//using namespace SimTK;

static double EPSILON = numeric_limits<double>::epsilon();

static bool FLAG_PLOT_CURVES    = false;
static string FILE_PATH      = "";
static double TOL_DX      = 5e-3;
static double TOL_DX_BIG     = 1e-2;
static double TOL_BIG     = 1e-6;
static double TOL_SMALL      = 1e-12;

/**
This function will print cvs file of the column vector col0 and the matrix 
 data

@params col0: A vector that must have the same number of rows as the data matrix
       This column vector is printed as the first column
@params data: A matrix of data
@params filename: The name of the file to print
*/
void printMatrixToFile(
    const RigidBodyDynamics::Math::VectorNd& col0, 
    const RigidBodyDynamics::Math::MatrixNd& data, 
    string& filename)
{
    
    ofstream datafile;
    datafile.open(filename.c_str());

    for(int i = 0; i < data.rows(); i++){
     datafile << col0(i) << ",";
     for(int j = 0; j < data.cols(); j++){
      if(j<data.cols()-1)
       datafile << data(i,j) << ",";
      else
       datafile << data(i,j) << "\n";
     }   
    }
    datafile.close();
} 


/**
This function will print cvs file of the matrix 
 data

@params data: A matrix of data
@params filename: The name of the file to print
*/
void printMatrixToFile( 
    const RigidBodyDynamics::Math::MatrixNd& data, 
    string& filename)
{
    ofstream datafile;
    datafile.open(filename.c_str());

    for(int i = 0; i < data.rows(); i++){
     for(int j = 0; j < data.cols(); j++){
      if(j<data.cols()-1)
       datafile << data(i,j) << ",";
      else
       datafile << data(i,j) << "\n";
     }   
    }
    datafile.close();
}

/**
    This function computes a standard central difference dy/dx. If 
    extrap_endpoints is set to 1, then the derivative at the end points is 
    estimated by linearly extrapolating the dy/dx values beside the end points

 @param x domain vector
 @param y range vector
 @param extrap_endpoints: (false)   Endpoints of the returned vector will be 
               zero, because a central difference
               is undefined at these endpoints
             (true)  Endpoints are computed by linearly 
               extrapolating using a first difference from 
               the neighboring 2 points
 @returns dy/dx computed using central differences
*/
RigidBodyDynamics::Math::VectorNd 
    calcCentralDifference(  RigidBodyDynamics::Math::VectorNd& x, 
             RigidBodyDynamics::Math::VectorNd& y, 
             bool extrap_endpoints){
 
    RigidBodyDynamics::Math::VectorNd dy(x.size());
    double dx1,dx2;
    double dy1,dy2;
    int size = x.size();
    for(int i=1; i<size-1; i++){
     dx1 = x(i)-x(i-1);
     dx2 = x(i+1)-x(i);
     dy1 = y(i)-y(i-1);
     dy2 = y(i+1)-y(i);
     dy(i)= 0.5*dy1/dx1 + 0.5*dy2/dx2;
    }

    if(extrap_endpoints == true){
     dy1 = dy(2)-dy(1);
     dx1 = x(2)-x(1);
     dy(0) = dy(1) + (dy1/dx1)*(x(0)-x(1));

     dy2 = dy(size-2)-dy(size-3);
     dx2 = x(size-2)-x(size-3);
     dy(size-1) = dy(size-2) + (dy2/dx2)*(x(size-1)-x(size-2));
    }
    return dy;
}

/**
    This function computes a standard central difference dy/dx at each point in
    a vector x, for a SmoothSegmentedFunction mcf, to a desired tolerance. This 
    function will take the best step size at each point to minimize the 
    error caused by taking a numerical derivative, and the error caused by
    numerical rounding error:

    For a step size of h/2 to the left and to the right of the point of 
    interest the error is
    error = 1/4*h^2*c3 + r*f(x)/h,         (1)
      
    Where c3 is the coefficient of the 3rd order Taylor series expansion
    about point x. Thus c3 can be computed if the order + 2 derivative is
    known
     
     c3 = (d^3f(x)/dx^3)/(6)            (2)
     
    And r*f(x)/h is the rounding error that occurs due to the central 
    difference.

    Taking a first derivative of 1 and solving for h yields

    h = (r*f(x)*2/c3)^(1/3)

    Where r is EPSILON

  @param x domain vector
  @param mcf the SmoothSegmentedFunction of interest
  @param order the order of the numerical derivative
  @param tolerance desired tolerance on the returned result
  @returns dy/dx computed using central differences
*/
RigidBodyDynamics::Math::VectorNd 
calcCentralDifference(  RigidBodyDynamics::Math::VectorNd& x, 
            SmoothSegmentedFunction& mcf,
            double tol, int order){
 

    RigidBodyDynamics::Math::VectorNd dyV(x.size());
    RigidBodyDynamics::Math::VectorNd yV(x.size());

    double y = 0;
    double dy = 0;
    double dyNUM = 0;
    double err= 0;
    double h = 0;
    double xL = 0;
    double xR = 0;

    double c3 = 0;
    double fL = 0;
    double fR = 0;
    double rootEPS = sqrt(EPSILON);

    double y_C3min = 1e-10;
    double y_C3max = 1e1;


    for(int i=0; i<x.size(); i++){
     yV(i) = mcf.calcDerivative(x(i),order-1);
    }
   

    for(int i=0; i< x.size(); i++){
     
     c3 = abs(mcf.calcDerivative(x(i),order+2));
     
     //singularity prevention
     if(abs(c3) < y_C3min)
      c3 = y_C3min;
     //Compute h
     y  = abs(mcf.calcDerivative(x(i), order-1));
     //preventing 0 from being assigned to y
     if(y < y_C3min)
      y = y_C3min;

     //Dumb check
     if(y/c3 < y_C3min){
      c3 = 1;
      y = y_C3min;
     }
     if(y/c3 > y_C3max){
      c3 = 1;
      y = y_C3max;
     }

     h  = pow( ( (EPSILON*y*2.0)/(c3) ) , 1.0/3.0);
    
     //Now check that h to the left and right are at least similar
     //If not, take the smallest one.
     xL = x(i)-h/2;
     xR = x(i)+h/2;

     fL = mcf.calcDerivative(xL, order-1);
     fR = mcf.calcDerivative(xR, order-1);

     //Just for convenience checking ...
     dyNUM = (fR-fL)/h;
     dy    = mcf.calcDerivative(x(i),order);
     err   = abs(dy-dyNUM);

     /*if(err > tol && abs(dy) > rootEPS && order <= 2){
      err = err/abs(dy);
      if(err > tol)
       cout << "rel tol exceeded" << endl;     
     }*/

     dyV(i) = dyNUM;

    }


    return dyV;
}

/**
    This function tests numerically for continuity of a curve. The test is 
    performed by taking a point on the curve, and then two points (called the 
    shoulder points) to the left and right of the point in question. The value
    of the functions derivative is evaluated at each of the shoulder points and
    used to linearly extrapolate from the shoulder points back to the original 
    point. If the original point and the linear extrapolations of each of the 
    shoulder points agree within tol, then the curve is assumed to be 
    continuous.


    @param x     Values to test for continuity
    @param yx    The SmoothSegmentedFunction to test
    @param order    The order of the curve of SmoothSegmentedFunction to test
        for continuity
    @param minTol   The minimum error allowed - this prevents the second order
        error term from going to zero
    @param taylorErrorMult  This scales the error tolerance. The default error
             tolerance is the the 2nd order Taylor series
             term.
*/
bool isFunctionContinuous(  RigidBodyDynamics::Math::VectorNd& xV, 
             SmoothSegmentedFunction& yV, 
             int order, 
             double minTol,
             double taylorErrorMult)
{
    bool flag_continuous = true;

    double xL = 0;   // left shoulder point
    double xR = 0;   // right shoulder point
    double yL = 0;   // left shoulder point function value
    double yR = 0;   // right shoulder point function value
    double dydxL = 0;   // left shoulder point derivative value
    double dydxR = 0;   // right shoulder point derivative value

    double xVal = 0;    //x value to test
    double yVal = 0;    //Y(x) value to test

    double yValEL = 0;  //Extrapolation to yVal from the left
    double yValER = 0;  //Extrapolation to yVal from the right

    double errL = 0;
    double errR = 0;

    double errLMX = 0;
    double errRMX = 0;


    for(int i =1; i < xV.size()-1; i++){
     xVal = xV(i);
     yVal = yV.calcDerivative(xVal, order);

     xL = 0.5*(xV(i)+xV(i-1));
     xR = 0.5*(xV(i)+xV(i+1));

     yL = yV.calcDerivative(xL,order);
     yR = yV.calcDerivative(xR,order);

     dydxL = yV.calcDerivative(xL,order+1);
     dydxR = yV.calcDerivative(xR,order+1);

     
     yValEL = yL + dydxL*(xVal-xL);
     yValER = yR - dydxR*(xR-xVal);

     errL = abs(yValEL-yVal);
     errR = abs(yValER-yVal);

     errLMX = abs(yV.calcDerivative(xL,order+2)*0.5*(xVal-xL)*(xVal-xL));
     errRMX = abs(yV.calcDerivative(xR,order+2)*0.5*(xR-xVal)*(xR-xVal));

     errLMX*=taylorErrorMult;
     errRMX*=taylorErrorMult;

     if(errLMX < minTol)
      errLMX = minTol;

     if(errRMX < minTol)
      errRMX = minTol; // to accomodate numerical
              //error in errL

     if(errL > errLMX || errR > errRMX){      
      flag_continuous = false;
     }
    }

    return flag_continuous;
}


/**
This function will scan through a vector and determine if it is monotonic or
not

@param y the vector of interest
@param multEPS The tolerance on the monotoncity check, expressed as a scaling of
       EPSILON
@return true if the vector is monotonic, false if it is not
*/
bool isVectorMonotonic( RigidBodyDynamics::Math::VectorNd& y, 
            int multEPS)
{
    double dir = y(y.size()-1)-y(0);
    bool isMonotonic = true;

    if(dir < 0){
     for(int i =1; i <y.size(); i++){
      if(y(i) > y(i-1)+EPSILON*multEPS){
       isMonotonic = false;
      //printf("Monotonicity broken at idx %i, since %fe-16 > %fe-16\n",
       //     i,y(i)*1e16,y(i-1)*1e16);
       printf("Monotonicity broken at idx %i, since "
        "y(i)-y(i-1) < tol, (%f*EPSILON < EPSILON*%i) \n",
            i,((y(i)-y(i-1))/EPSILON), multEPS);
      }
     }
    }
    if(dir > 0){
     for(int i =1; i <y.size(); i++){
      if(y(i) < y(i-1)-EPSILON*multEPS){
       isMonotonic = false;
       printf("Monotonicity broken at idx %i, since "
        "y(i)-y(i-1) < -tol, (%f*EPSILON < -EPSILON*%i) \n",
            i,((y(i)-y(i-1))/EPSILON), multEPS);
      }
     }
    }
    if(dir == 0){
     isMonotonic = false;
    }

    return isMonotonic;
}

/**
This function will compute the numerical integral of y(x) using the trapezoidal
method

@param x the domain vector
@param y the range vector, of y(x), evaluated at x
@param flag_TrueIntForward_FalseIntBackward 
    When this flag is set to true, the integral of y(x) will be evaluated from
    left to right, starting with int(y(0)) = 0. When this flag is false, then
    y(x) will be evaluated from right to left with int(y(n)) = 0, where n is 
    the maximum number of elements.                    
@return the integral of y(x)
*/

RigidBodyDynamics::Math::VectorNd calcTrapzIntegral(
            RigidBodyDynamics::Math::VectorNd& x, 
            RigidBodyDynamics::Math::VectorNd& y, 
            bool flag_TrueIntForward_FalseIntBackward)
{
    RigidBodyDynamics::Math::VectorNd inty 
     = RigidBodyDynamics::Math::VectorNd::Zero(y.size());
    //inty = 0;


    int startIdx = 1;
    int endIdx = y.size()-1;

    if(flag_TrueIntForward_FalseIntBackward == true){
      
     double width = 0;
     for(int i = 1; i <= endIdx; i=i+1){
      width = abs(x(i)-x(i-1));
      inty(i) = inty(i-1) +  width*(0.5)*(y(i)+y(i-1));
     }

    }else{
     
     double width = 0;      
     for(int i = endIdx-1; i >= 0; i=i-1){
      width = abs(x(i)-x(i+1));
      inty(i) = inty(i+1) +  width*(0.5)*(y(i)+y(i+1));
     }
    }

  
    return inty;
}

/**
   @param a The first vector
   @param b The second vector
   @return Returns the maximum absolute difference between vectors a and b
*/
double calcMaximumVectorError(RigidBodyDynamics::Math::VectorNd& a, 
            RigidBodyDynamics::Math::VectorNd& b)
{
    double error = 0;
    double cerror=0;
    for(int i = 0; i< a.size(); i++)
    {
     cerror = abs(a(i)-b(i));
     if(cerror > error){
      error = cerror;
     }     
    }
    return error;
}


/*
void testQuinticBezier_Exceptions(){
    //cout <<"**************************************************"<<endl;
    //cout << "   TEST: Bezier Curve Exceptions" << endl;
    string name  = "testQuinticBezier_Exceptions()";

    //Generate a Bezier curve
    RigidBodyDynamics::Math::VectorNd xPts(6);
    RigidBodyDynamics::Math::VectorNd yPts(6);

    RigidBodyDynamics::Math::MatrixNd xMPts(6,1);
    RigidBodyDynamics::Math::MatrixNd yMPts(6,1);

    xPts(0) = 0;
    xPts(1) = 0.5;
    xPts(2) = 0.5;
    xPts(3) = 0.75;
    xPts(4) = 0.75;
    xPts(5) = 1;
    xMPts.col(0) = xPts;

    yPts(0) = 0;
    yPts(1) = 0.125;
    yPts(2) = 0.125;
    yPts(3) = 0.5;
    yPts(4) = 0.5;
    yPts(5) = 1;
    yMPts.col(0) = yPts;

    //RigidBodyDynamics::Math::VectorNd u(100);
    //RigidBodyDynamics::Math::VectorNd x(100);
    //SimTK::Array_< SimTK::Spline > aSplineUX(1);
    //for(int i=0; i<100; i++){
    //   u(i) = ((double)i)/((double)99);
    //    x(i)=SegmentedQuinticBezierToolkit::
    //         calcQuinticBezierCurveVal(u(i), xPts);
    //}

    //aSplineUX[0] = SimTK::SplineFitter<double>::
    //     fitForSmoothingParameter(3,x,u,0).getSpline();
    //Now we have a curve.

    //=========================================================================
    //Test exceptions for calcQuinticBezierCornerControlPoints
    //=========================================================================
    double x0 = 0;
    double x1 = 1;
    double y0 = 0;
    double y1 = 1;
    double dydx0 = 0;
    double dydx1 = 1;
    double curviness = 0;

    double curvinessEX1 = 1.01; //illegal value 
    double curvinessEX2 = -0.01; //illegal value 


    double dydx0EX1 = 0; //illeagle pair
    double dydx1EX1 = 0.1;

    bool exceptionThrown = false;
    //try{
    CHECK_ASSERT(
     RigidBodyDynamics::Math::MatrixNd test1 =
     SegmentedQuinticBezierToolkit::
     calcQuinticBezierCornerControlPoints(   x0, y0, dydx0, 
                     x1, y1, dydx1, 
                      curvinessEX1);
    );
    //}catch(...){
    //    exceptionThrown = true;
    //}
    CHECK(exceptionThrown);

    exceptionThrown = false;
    try{
     RigidBodyDynamics::Math::MatrixNd test2 = 
     SegmentedQuinticBezierToolkit::
     calcQuinticBezierCornerControlPoints(x0, y0, dydx0, 
                  x1, y1,dydx1, 
                  curvinessEX2);
    }catch(...){
     exceptionThrown = true;
    }
    CHECK(exceptionThrown);


    exceptionThrown = false;
    try{
     RigidBodyDynamics::Math::MatrixNd test2 = 
     SegmentedQuinticBezierToolkit::
     calcQuinticBezierCornerControlPoints(x0, y0, dydx0EX1, 
                  x1, y1, dydx1EX1, 
                  curviness);
    }catch(...){
     exceptionThrown = true;
    }
    CHECK(exceptionThrown);

    //=========================================================================
    //Test exceptions for calcIndex
    //=========================================================================

    double xEX1 = xPts(0)-0.01; //This is not in the set.
    double xEX2 = xPts(5)+0.01; //This is not in the set.

    exceptionThrown = false;
    try{
     int t = SegmentedQuinticBezierToolkit::calcIndex(xEX1, xMPts);
    }catch(...){
     exceptionThrown = true;
    }
    CHECK(exceptionThrown);

    exceptionThrown = false;
    try{    
     int t = SegmentedQuinticBezierToolkit::calcIndex(xEX2, xMPts);
    }catch(...){
     exceptionThrown = true;
    }
    CHECK(exceptionThrown);

    //=========================================================================
    //Test exceptions for calcU
    //=========================================================================
    
    //xEX1 is not within the curve, and so the Newton iteration will not
    //converge
    exceptionThrown = false;
    try{
     double uPt = 
     SegmentedQuinticBezierToolkit::calcU(xEX1, xPts, 1e-8, 10);
    }catch(...){
     exceptionThrown = true;
    }
    CHECK(exceptionThrown);

    //=========================================================================
    //Test exceptions for calcQuinticBezierCurveVal
    //=========================================================================

    double uEX1 = -0.01; //illeagle
    double uEX2 = 1.01;  //illeagle
    RigidBodyDynamics::Math::VectorNd xPtsEX 
     = RigidBodyDynamics::Math::VectorNd::Zero(5);
    //xPtsEX = 0;

    exceptionThrown = false;
    try{double tst = SegmentedQuinticBezierToolkit::
         calcQuinticBezierCurveVal(uEX1,xPts);
    }catch(...){
     exceptionThrown = true;
    }
    CHECK(exceptionThrown);

    exceptionThrown = false;
    try{double tst = SegmentedQuinticBezierToolkit::
        calcQuinticBezierCurveVal(uEX2,xPts);
    }catch(...){
     exceptionThrown = true;
    }
    CHECK(exceptionThrown);     

    exceptionThrown = false;
    try{
     double tst = SegmentedQuinticBezierToolkit::
        calcQuinticBezierCurveVal(0.5,xPtsEX);
    }catch(...){
     exceptionThrown = true;
    }
    CHECK(exceptionThrown);

    //=========================================================================
    //Test exceptions for calcQuinticBezierCurveDerivU
    //=========================================================================
    try{
     double tst = SegmentedQuinticBezierToolkit::
     calcQuinticBezierCurveDerivU(uEX1, xPts, (int)1);
    }catch(...){
     exceptionThrown = true;
    }
    CHECK(exceptionThrown);     

    exceptionThrown = false;
    try{    
     double tst = SegmentedQuinticBezierToolkit::
     calcQuinticBezierCurveDerivU(uEX2, xPts, (int)1);
    }catch(...){
     exceptionThrown = true;
    }
    CHECK(exceptionThrown);     

    exceptionThrown = false;
    try{
     double tst = SegmentedQuinticBezierToolkit::
     calcQuinticBezierCurveDerivU(0.5, xPtsEX, (int)1);
    }catch(...){
     exceptionThrown = true;
    }
    CHECK(exceptionThrown);     

    exceptionThrown = false;
    try{     
     double tst = SegmentedQuinticBezierToolkit::
     calcQuinticBezierCurveDerivU(0.5, xPts, 0);
    }catch(...){
     exceptionThrown = true;
    }
    CHECK(exceptionThrown);     

    //=========================================================================
    //Test exceptions for calcQuinticBezierCurveDerivDYDX
    //=========================================================================
    
    try{double test= SegmentedQuinticBezierToolkit::
     calcQuinticBezierCurveDerivDYDX(uEX1,xPts,yPts,1);
    }catch(...){
     exceptionThrown = true;
    }
    CHECK(exceptionThrown);     

    exceptionThrown = false;
    try{double test= SegmentedQuinticBezierToolkit::
     calcQuinticBezierCurveDerivDYDX(uEX2,xPts,yPts,1);
    }catch(...){
     exceptionThrown = true;
    }
    CHECK(exceptionThrown);     

    exceptionThrown = false;
    try{double test= SegmentedQuinticBezierToolkit::
     calcQuinticBezierCurveDerivDYDX(0.5,xPts,yPts,0);
    }catch(...){
     exceptionThrown = true;
    }
    CHECK(exceptionThrown);     

    exceptionThrown = false;
    try{double test= SegmentedQuinticBezierToolkit::
     calcQuinticBezierCurveDerivDYDX(0.5,xPts,yPts,7);
    }catch(...){
     exceptionThrown = true;
    }
    CHECK(exceptionThrown);     

    exceptionThrown = false;
    try{double test= SegmentedQuinticBezierToolkit::
     calcQuinticBezierCurveDerivDYDX(0.5,xPtsEX,yPts,3);
    }catch(...){
     exceptionThrown = true;
    }
    CHECK(exceptionThrown);     

    exceptionThrown = false;
    try{double test= SegmentedQuinticBezierToolkit::
     calcQuinticBezierCurveDerivDYDX(0.5,xPts,xPtsEX,3);
    }catch(...){
     exceptionThrown = true;
    }
    CHECK(exceptionThrown);     
    //=========================================================================
    //Test exceptions for calcNumIntBezierYfcnX
    //=========================================================================
    //There are none.
    //cout << "    passed. All exceptions are functioning as intended" << endl;
    //cout <<"**************************************************"<<endl;
}
*/


/**
    This function will create a quintic Bezier curve y(x) and sample it, its 
    first derivative w.r.t. U (dx(u)/du and dy(u)/du), and its first derivative
    w.r.t. to X and print it to the screen.
*/
void testQuinticBezier_DU_DYDX(int maximumNumberOfToleranceViolations)
{
    //cout <<"**************************************************"<<endl;
    //cout << "   TEST: Bezier Curve Derivative DU" << endl;
    string name  = "testQuinticBezier_DU_DYDX()";
     RigidBodyDynamics::Math::VectorNd xPts(6);
     RigidBodyDynamics::Math::VectorNd yPts(6);
     xPts(0) = 0;
     xPts(1) = 0.5;
     xPts(2) = 0.5;
     xPts(3) = 0.75;
     xPts(4) = 0.75;
     xPts(5) = 1;

     yPts(0) = 0;
     yPts(1) = 0.125;
     yPts(2) = 0.125;
     yPts(3) = 0.5;
     yPts(4) = 0.5;
     yPts(5) = 1;

     double val = 0;
     double d1 = 0;
     double d2 = 0;
     double d3 = 0;
     double d4 = 0;
     double d5 = 0;
     double d6 = 0;

     double u = 0;

     int steps = 100;

     RigidBodyDynamics::Math::MatrixNd analyticDerXU(steps,8);
     RigidBodyDynamics::Math::MatrixNd analyticDerYU(steps,8);
     RigidBodyDynamics::Math::VectorNd uV(steps);
     for(int i = 0; i<steps; i++){
      //int i = 10;
      u = (double)i/(steps-1);
      uV(i) = u;

      val= SegmentedQuinticBezierToolkit::
       calcQuinticBezierCurveVal(u,xPts);
      d1 = SegmentedQuinticBezierToolkit::
       calcQuinticBezierCurveDerivU(u,xPts,1);
      d2 = SegmentedQuinticBezierToolkit::
       calcQuinticBezierCurveDerivU(u,xPts,2);
      d3 = SegmentedQuinticBezierToolkit::
       calcQuinticBezierCurveDerivU(u,xPts,3);
      d4 = SegmentedQuinticBezierToolkit::
       calcQuinticBezierCurveDerivU(u,xPts,4);
      d5 = SegmentedQuinticBezierToolkit::
       calcQuinticBezierCurveDerivU(u,xPts,5);
      d6 = SegmentedQuinticBezierToolkit::
       calcQuinticBezierCurveDerivU(u,xPts,6);

      analyticDerXU(i,0) = u;
      analyticDerXU(i,1) = val;
      analyticDerXU(i,2) = d1;
      analyticDerXU(i,3) = d2;
      analyticDerXU(i,4) = d3;
      analyticDerXU(i,5) = d4;
      analyticDerXU(i,6) = d5;
      analyticDerXU(i,7) = d6;

      val= SegmentedQuinticBezierToolkit::
       calcQuinticBezierCurveVal(u,yPts);
      d1 = SegmentedQuinticBezierToolkit::
       calcQuinticBezierCurveDerivU(u,yPts,1);
      d2 = SegmentedQuinticBezierToolkit::
       calcQuinticBezierCurveDerivU(u,yPts,2);
      d3 = SegmentedQuinticBezierToolkit::
       calcQuinticBezierCurveDerivU(u,yPts,3);
      d4 = SegmentedQuinticBezierToolkit::
       calcQuinticBezierCurveDerivU(u,yPts,4);
      d5 = SegmentedQuinticBezierToolkit::
       calcQuinticBezierCurveDerivU(u,yPts,5);
      d6 = SegmentedQuinticBezierToolkit::
       calcQuinticBezierCurveDerivU(u,yPts,6);

      analyticDerYU(i,0) = u;
      analyticDerYU(i,1) = val;
      analyticDerYU(i,2) = d1;
      analyticDerYU(i,3) = d2;
      analyticDerYU(i,4) = d3;
      analyticDerYU(i,5) = d4;
      analyticDerYU(i,6) = d5;
      analyticDerYU(i,7) = d6;

     }

     int mxDU = 6-1;
     RigidBodyDynamics::Math::MatrixNd numericDer(analyticDerXU.rows(), mxDU);
     RigidBodyDynamics::Math::MatrixNd errorDer(analyticDerXU.rows(), mxDU);

     double tol = (double)(1.0/steps);
     tol = tol*tol*50; 
     //Numerical error in a central difference increases with the 
     //square of h. 
     //http://en.wikipedia.org/wiki/Finite_difference

     RigidBodyDynamics::Math::VectorNd domainX = 
      RigidBodyDynamics::Math::VectorNd::Zero(analyticDerXU.rows());
     RigidBodyDynamics::Math::VectorNd rangeY = 
      RigidBodyDynamics::Math::VectorNd::Zero(analyticDerXU.rows());

     RigidBodyDynamics::Math::VectorNd analyticDerYX = 
      RigidBodyDynamics::Math::VectorNd::Zero(analyticDerXU.rows());
     for(int j=0; j<analyticDerXU.rows(); j++){
      domainX(j) = analyticDerXU(j,0);
     }


     for(int i=0;i<mxDU;i++){

      for(int j=0; j<analyticDerXU.rows(); j++){
       rangeY(j)     = analyticDerXU(j,i+1);
       analyticDerYX(j)    = analyticDerXU(j,i+2);
      }


      numericDer.col(i) = calcCentralDifference(domainX,
                    rangeY,true); 
      errorDer.col(i) = analyticDerYX-numericDer.col(i);


      //errorDer(i)= abs( errorDer(i).elementwiseDivide(numericDer(i)) );
      //The end points can't be tested because a central difference
      //cannot be accurately calculated at these locations
      for(int j=1; j<analyticDerXU.rows()-1; j++){
       assert( abs(errorDer(j,i))<tol );
       //if(errorDer(j,i)>tol)
       //printf("Error > Tol: (%i,%i): %f > %f\n",j,i,errorDer(j,i),tol);
      }
     }     
     errorDer.cwiseAbs();
     //cout << errorDer << endl;

    //printf("...absolute tolerance of %f met\n", tol);

     //cout << "   TEST: Bezier Curve Derivative DYDX to d6y/dx6" << endl;
     RigidBodyDynamics::Math::MatrixNd numericDerXY(analyticDerXU.rows(), 6);
     RigidBodyDynamics::Math::MatrixNd analyticDerXY(analyticDerXU.rows(),6);

     for(int i=0; i< analyticDerXU.rows(); i++)
     {
      analyticDerXY(i,0) = SegmentedQuinticBezierToolkit::
       calcQuinticBezierCurveDerivDYDX(uV(i),xPts,yPts,1);
      analyticDerXY(i,1) = SegmentedQuinticBezierToolkit::
       calcQuinticBezierCurveDerivDYDX(uV(i),xPts,yPts,2);
      analyticDerXY(i,2) = SegmentedQuinticBezierToolkit::
       calcQuinticBezierCurveDerivDYDX(uV(i),xPts,yPts,3);
      analyticDerXY(i,3) = SegmentedQuinticBezierToolkit::
       calcQuinticBezierCurveDerivDYDX(uV(i),xPts,yPts,4);
      analyticDerXY(i,4) = SegmentedQuinticBezierToolkit::
       calcQuinticBezierCurveDerivDYDX(uV(i),xPts,yPts,5);
      analyticDerXY(i,5) = SegmentedQuinticBezierToolkit::
       calcQuinticBezierCurveDerivDYDX(uV(i),xPts,yPts,6);
     }


     for(int j=0; j<numericDerXY.cols();j++){

      for(int k=0; k<numericDerXY.rows(); k++){
       domainX(k) = analyticDerXU(k,1);
       if(j == 0){
        rangeY(k)  = analyticDerYU(k,1);
       }else{
        rangeY(k)  = analyticDerXY(k,j-1);
       }
      }
      numericDerXY.col(j) = calcCentralDifference(domainX,
                       rangeY,true);
     
     }


     //Generate numerical derivative curves for the first 3 derivatives
     /*
     numericDerXY.col(0) = calcCentralDifference(analyticDerXU.col(1),
                     analyticDerYU.col(1),true);
     numericDerXY.col(1) = calcCentralDifference(analyticDerXU.col(1),
                     analyticDerXY.col(0),true);
     numericDerXY.col(2) = calcCentralDifference(analyticDerXU.col(1),
                     analyticDerXY.col(1),true);
     numericDerXY.col(3) = calcCentralDifference(analyticDerXU.col(1),
                     analyticDerXY.col(2),true);
     numericDerXY.col(4) = calcCentralDifference(analyticDerXU.col(1),
                     analyticDerXY.col(3),true);
     numericDerXY.col(5) = calcCentralDifference(analyticDerXU.col(1),
                     analyticDerXY.col(4),true);
     */

     //Create the matrix of errors
     RigidBodyDynamics::Math::MatrixNd errorDerXYNum(analyticDerXU.rows(), 6);
     RigidBodyDynamics::Math::MatrixNd errorDerXYDen(analyticDerXU.rows(), 6);
     RigidBodyDynamics::Math::MatrixNd errorDerXY(analyticDerXU.rows(), 6);
     errorDerXYNum = analyticDerXY-numericDerXY;
     errorDerXYNum.cwiseAbs();
     errorDerXYDen = analyticDerXY+numericDerXY;
     errorDerXYDen.cwiseAbs();
     errorDerXY    = errorDerXYNum.cwiseQuotient(errorDerXYDen);

     double relTol = 5e-2;

     int relTolExceeded = 0;

     for(int j=0;j<6;j++){
      //can't test the first and last entries because a central diff.
      //cannot calculate these values accurately.
      for(int i=1;i<analyticDerXU.rows()-1;i++){
       if(errorDerXY(i,j)>relTol){
        //printf("Error > Tol: (%i,%i): %f > %f\n",i,j,
        //          errorDerXY(i,j),relTol);
        relTolExceeded++;
       }
      }
     }
     //cout << relTolExceeded << endl;

     //The relative tolerance gets exceeded occasionally in locations of
     //rapid change in the curve. Provided there are only a few locations
     //where the relative tolerance of 5% is broken, the curves should be
     //regarded as being good. Ten errors out of a possible 100*6 data points
     //seems relatively small.
     CHECK(relTolExceeded < maximumNumberOfToleranceViolations);


     //std::string fname = "analyticDerXY.csv";
     //printMatrixToFile(analyticDerXY,fname);
     //fname = "numericDerXY.csv";
     //printMatrixToFile(numericDerXY,fname);
     //fname = "errorDerXY.csv";
     //printMatrixToFile(errorDerXY,fname);
     //printf("   ...relative tolerance of %f not exceeded more than %i times\n"
     //    "   across all 6 derivatives, with 100 samples each\n",
     //        relTol, 10);
     //cout <<"**************************************************"<<endl;


}

/**
    This function will generate a series of quintic Bezier control points 
    that form a C2 corner and print the results to the command window
*/
void sampleBezierCornerGeneration()
{
    string name = "testBezierCornerGeneration()";
    double x0 = 0;
    double y0 = 0;
    double dydx0 = 0.01;
    double x1 = 1;
    double y1 = 1;
    double dydx1 = 2;
    double curviness = 0.5;

    
    
    
    RigidBodyDynamics::Math::MatrixNd xyPts = 
     SegmentedQuinticBezierToolkit::
     calcQuinticBezierCornerControlPoints(x0,y0,dydx0, 
                  x1,y1,dydx1,
                  curviness);

    cout << "XY Corner Control Points" << endl;
    cout << xyPts << endl;


}
/**
    This function will sample and print a Function to file
*/
void samplePrintExtrapolatedFunction(const Function_<double>& fcn, 
                double x0, double x1, 
                int num, string& name)
{
    RigidBodyDynamics::Math::VectorNd tcX(num);
    RigidBodyDynamics::Math::MatrixNd tcY(num,3);
    std::vector<int> dx(1), ddx(2);
    dx[0]=0;
    ddx[0]=0;
    ddx[1]=0;

    double tcD = (x1-x0)/(num-1);
    RigidBodyDynamics::Math::VectorNd x(1);

    for(int i=0; i<num; i++){
     tcX(i) = x0 + i*tcD;
     x(0)   = tcX(i);
     tcY(i,0)  = fcn.calcValue(x);
     tcY(i,1)  = fcn.calcDerivative(dx,x);
     tcY(i,2)  = fcn.calcDerivative(ddx,x);
    }
    printMatrixToFile(tcX,tcY,name);
}



/**
 1. The SmoothSegmentedFunction's derivatives will be compared against 
    numerically calculated derivatives to ensure that the errors between the 
    2 are small

*/
bool areCurveDerivativesCloseToNumericDerivatives(
     SmoothSegmentedFunction& mcf,
     RigidBodyDynamics::Math::MatrixNd& mcfSample,
     double tol)
{
    //cout << "   TEST: Derivative correctness " << endl;
    int maxDer = 4;//mcf.getMaxDerivativeOrder() - 2;
   
    RigidBodyDynamics::Math::MatrixNd numSample(mcfSample.rows(),maxDer);  
    RigidBodyDynamics::Math::MatrixNd relError(mcfSample.rows(),maxDer);
    

    RigidBodyDynamics::Math::VectorNd domainX = 
     RigidBodyDynamics::Math::VectorNd::Zero(mcfSample.rows());

    for(int j=0; j<mcfSample.rows(); j++)
     domainX(j) = mcfSample(j,0);

    for(int i=0; i < maxDer; i++){
     //Compute the relative error
     numSample.col(i)=calcCentralDifference(domainX,mcf,tol,i+1);
     relError.col(i)= mcfSample.col(i+2)-numSample.col(i);  

     //compute a relative error where possible
     for(int j=0; j < relError.rows(); j++){
      if(abs(mcfSample(j,i+2)) > tol){
       relError(j,i) = relError(j,i)/mcfSample(j,i+2);
      }
     }

    }

    RigidBodyDynamics::Math::VectorNd errRelMax =
     RigidBodyDynamics::Math::VectorNd::Zero(6);
    RigidBodyDynamics::Math::VectorNd errAbsMax = 
     RigidBodyDynamics::Math::VectorNd::Zero(6);

    double absTol = 5*tol;

    bool flagError12=false;
    RigidBodyDynamics::Math::VectorNd tolExceeded12V =
     RigidBodyDynamics::Math::VectorNd::Zero(mcfSample.rows());
    
    int tolExceeded12 = 0;
    int tolExceeded34 = 0;
    for(int j=0;j<maxDer;j++){
     
     for(int i=0; i<mcfSample.rows(); i++){
      if(relError(i,j) > tol && mcfSample(i,j+2) > tol){
       if(j <= 1){
        tolExceeded12++;
        tolExceeded12V(i)=1;
        flagError12=true;
       }
       if(j>=2)
        tolExceeded34++;       
      }
      if(mcfSample(i,j+2) > tol)
      if(errRelMax(j) < abs(relError(i,j)))
        errRelMax(j) = abs(relError(i,j));

      //This is a harder test: here we're comparing absolute error
      //so the tolerance margin is a little higher
      if(relError(i,j) > absTol && mcfSample(i,j+2) <= tol){
       if(j <= 1){
        tolExceeded12++;
        tolExceeded12V(i)=1;
        flagError12=true;
       }
       if(j>=2)
        tolExceeded34++;            
      }

      if(mcfSample(i,j+2) < tol)
      if(errAbsMax(j) < abs(relError(i,j)))
        errAbsMax(j) = abs(relError(i,j));
      
     
     }

     /*
     if(flagError12 == true){
      printf("Derivative %i Rel Error Exceeded:\n",j);
      printf("x     dx_relErr dx_calcVal dx_sample"
        " dx2_relErr dx2_calcVal dx2_sample\n");
      for(int i=0; i<mcfSample.rows(); i++){
       if(tolExceeded12V(i) == 1){
          printf("%f %f %f  %f %f   %f    %f",
        mcfSample(i,0),relError(i,0),mcfSample(i,2),numSample(i,0),
                relError(i,1),mcfSample(i,3),numSample(i,1));
        
       }
      }
     }
     flagError12=false;*/
     //tolExceeded12V = 
     //RigidBodyDynamics::Math::VectorNd::Zero(mcfSample.rows());
    }
    




    //CHECK(tolExceeded12 == 0);

  //   printMatrixToFile(mcfSample,"analyticDerivatives.csv");
  //   printMatrixToFile(numSample,"numericDerivatives.csv");
  //   printMatrixToFile(numError,"numAnalyticError.csv");
  //   cout << "Matricies Printed" << endl;
    /*
   printf("passed: A tolerance of %fe-3 reached with a maximum relative error\n"
       "     of %fe-3 and %fe-3 for the first two derivatives\n"
       "     at points where it is possible to compute the relative\n"
       "     error.                     \n\n"
       "     At points where the relative error couldn't be computed,\n"
       "     due to divide by 0, the first two derivatives met an \n"
       "     absolute tolerance of %fe-3, with a maximum value of \n"
       "     %fe-3 and %fe-3. \n\n"
       "     Derivatives of order 3 through %i exceeded the \n"
       "     relative or absolute tolerance (as appropriate) %i times.\n",
       tol*1e3, errRelMax(0)*1e3, errRelMax(1)*1e3, absTol*1e3,
       errAbsMax(0)*1e3, errAbsMax(1)*1e3, maxDer, tolExceeded34);
  */
  //cout << endl;

   return (tolExceeded12 == 0);
    
}
 /*
 2. The SmoothSegmentedFunction's integral, when calculated, will be compared 
    against a numerically computed integral (using the trapezoidal method)
 */

/*
void testMuscleCurveIntegral(SmoothSegmentedFunction& mcf,
              RigidBodyDynamics::Math::MatrixNd& mcfSample)
{


    //2. Integral test
    if(mcf.isIntegralAvailable()){
     cout << "   TEST: Integral correctness " << endl;
     RigidBodyDynamics::Math::VectorNd intyTRAPZ 
            = calcTrapzIntegral(mcfSample(0), mcfSample(1), 
                mcf.isIntegralComputedLeftToRight());
     //cout << intyTRAPZ << endl;
     //The error of the trapezoidal integration method is
     //no more than (width^2/2)*sup(f'(x)), which we approximate by
     //http://en.wikipedia.org/wiki/Numerical_integration#Conservative_.28a
     //_priori.29_error_estimation

     //Compute the width between each of the samples, because the width
     //of this sample dictates how accurate the numerical derivative is
     RigidBodyDynamics::Math::VectorNd xWidth(mcfSample.rows());

     double xWidthMax = 0;
     for(int i=0; i<mcfSample.rows()-1; i++){
      xWidth(i) = mcfSample(i+1,0)-mcfSample(i,0);
      if(xWidth(i) > xWidthMax){
       xWidthMax = xWidth(i);
      }

     }
     xWidth(mcfSample.rows()-1) = xWidth(mcfSample.rows()-2);

     //cout << xWidth << endl;
     //cout <<endl;

     RigidBodyDynamics::Math::VectorNd intyCumError(xWidth.size());
     double supdf = 0;

     if(mcf.isIntegralComputedLeftToRight()){
      //Get the vector of accumulated errors left to rigth
      intyCumError(0) = xWidth(0)*xWidth(0)*0.5*mcfSample(0,1);
      for(int i=1; i< intyCumError.size(); i++){
       supdf = max(mcfSample(i,1),mcfSample(i-1,1));
       intyCumError(i) = intyCumError(i-1)
                + xWidth(i)*xWidth(i)*0.5*supdf;
      }
      //margin of error
      intyCumError = intyCumError*2.0;
      intyCumError += intyCumError(intyCumError.size()-1)*1.05;
      
     }
     else{
      //Get the vector of accumulated errors right to left
      int eidx = intyCumError.size()-1;
      intyCumError(eidx) = xWidth(eidx)*xWidth(eidx)*0.5*mcfSample(eidx,1);
      for(int i=eidx-1; i>=0; i--){
       supdf = max(mcfSample(i,1),mcfSample(i+1,1));
       intyCumError(i) = intyCumError(i+1)
                + xWidth(i)*xWidth(i)*0.5*supdf;
      }
      //margin of error
      intyCumError = intyCumError*2.0;
      intyCumError += intyCumError(0)*1.05;
      
     }
     

     //cout << intyCumError << endl;
     //cout << endl;

     double maxError = 0;
     double mxErr = 0;
     double mnErr = 1e10;

     RigidBodyDynamics::Math::VectorNd error(mcfSample.rows());
     //double relErrorTol = 0.01;

     int intyCol = 2+6;//mcf.getMaxDerivativeOrder();
     int errCtr = 0;
     for(int i=0; i< intyCumError.size(); i++){
      error(i)=abs( intyTRAPZ(i)-mcfSample(i,intyCol) );
      maxError = intyCumError(i);

      if(error(i)>mxErr)
       mxErr=error(i);
      if(mnErr>error(i))
       mnErr=error(i);


      
      if(error(i)>maxError){
       printf("Tol exceeded (%i), %f e-6> %f e-6\n",i,error(i)*1e6,
         maxError*1e6);
       errCtr++;
      }
     }
     SimTK_TEST(errCtr==0);
     //cout << error << endl;
     printf("   passed: integral agrees with trapz within a factor of 2.05 \n"
      "   of the theoretical accuracy of trapz, with a maximum \n"
      "   error of %fe-6\n",mxErr*1e6);
     cout << endl;
    }

}
*/

/*
 3. The MuscleCurveFunctions function value, first and second derivative curves
    will be numerically tested for continuity.
*/ 
bool isCurveC2Continuous(SmoothSegmentedFunction& mcf,
              RigidBodyDynamics::Math::MatrixNd& mcfSample)
{
    //cout << "   TEST: C2 Continuity " << endl;

    int multC0 = 5;
    int multC1 = 50;
    int multC2 = 100;

    RigidBodyDynamics::Math::VectorNd fcnSample = 
     RigidBodyDynamics::Math::VectorNd::Zero(mcfSample.rows());

    for(int i=0; i < mcfSample.rows(); i++){
     fcnSample(i) = mcfSample(i,0);
    }


    bool c0 = isFunctionContinuous(fcnSample, mcf, 0, 1e-6, multC0);
    bool c1 = isFunctionContinuous(fcnSample, mcf, 1, 1e-6, multC1);
    bool c2 = isFunctionContinuous(fcnSample, mcf, 2, 1e-6, multC2);



    return (c0 && c1 && c2);
    //printf( "   passed: C2 continuity established to a multiple\n"
    //     "     of the next Taylor series error term.\n "
    //     "     C0,C1, and C2 multiples: %i,%i and %i\n",
    //        multC0,multC1,multC2);
    //cout << endl;
}

/*
 4. The MuscleCurveFunctions which are supposed to be monotonic will be
    tested for monotonicity.
*/
bool isCurveMontonic(RigidBodyDynamics::Math::MatrixNd mcfSample)
{
    //cout << "   TEST: Monotonicity " << endl;
    int multEps = 10;

    RigidBodyDynamics::Math::VectorNd fcnSample = 
     RigidBodyDynamics::Math::VectorNd::Zero(mcfSample.rows());

    for(int i=0; i < mcfSample.rows(); i++){
     fcnSample(i) = mcfSample(i,1);
    }

    bool monotonic = isVectorMonotonic(fcnSample,10);
    return monotonic;
    //printf("   passed: curve is monotonic to %i*EPSILON",multEps);
    //cout << endl;
}


   
TEST(QuinticBezierToolKitDerivatives)
{
    int maximumNumberOfToleranceViolations = 10;
    testQuinticBezier_DU_DYDX(10);
}


TEST(tendonCurve)
{
    //cout <<"**************************************************"<<endl;
    //cout <<"TENDON CURVE TESTING            "<<endl;
    double e0   = 0.04;
    double kiso = 1.5/e0;
    double c    = 0.5;//0.75;    
    double ftoe = 1.0/3.0;

    SmoothSegmentedFunction tendonCurve = SmoothSegmentedFunction();
    SmoothSegmentedFunctionFactory::createTendonForceLengthCurve(
         e0,kiso,ftoe,c, "test_tendonCurve", tendonCurve);
    

    RigidBodyDynamics::Math::MatrixNd tendonCurveSample
     =tendonCurve.calcSampledCurve(6,1.0,1+e0);
    //tendonCurve.printMuscleCurveToCSVFile(FILE_PATH);

//0. Test that each curve fulfills its contract at the end points.
    //cout << "   Keypoint Testing" << endl;
    RigidBodyDynamics::Math::VectorNd tendonCurveDomain = 
     tendonCurve.getCurveDomain();
    CHECK(abs(tendonCurve.calcValue(tendonCurveDomain(0)))<TOL_SMALL);      
    CHECK(abs(tendonCurve.calcValue(tendonCurveDomain(1))-ftoe)<TOL_SMALL);

    CHECK(abs(tendonCurve.calcValue(1.0)    )< TOL_SMALL);  
    CHECK(abs(tendonCurve.calcDerivative(1.0,1))< TOL_BIG);
    CHECK(abs(tendonCurve.calcDerivative(1.0,2))< TOL_BIG);

    CHECK(abs(tendonCurve.calcValue(1+e0)    -1.0 ) < TOL_SMALL);  
    CHECK(abs(tendonCurve.calcDerivative(1+e0,1)-kiso) < TOL_BIG);
    CHECK(abs(tendonCurve.calcDerivative(1+e0,2)-0   ) < TOL_BIG);
    //cout << "   passed" << endl;
    //cout << endl;
//1. Test each derivative sample for correctness against a numerically
//   computed version
    bool areCurveDerivativesGood = 
      areCurveDerivativesCloseToNumericDerivatives(
       tendonCurve,
       tendonCurveSample,TOL_DX_BIG);

    CHECK(areCurveDerivativesGood);
//2. Test each integral, where computed for correctness.
    //testMuscleCurveIntegral(tendonCurve, tendonCurveSample);

//3. Test numerically to see if the curve is C2 continuous
    bool curveIsContinuous = isCurveC2Continuous(tendonCurve,
                      tendonCurveSample);

    CHECK(curveIsContinuous);
//4. Test for montonicity where appropriate
    bool curveIsMonotonic = isCurveMontonic(tendonCurveSample);
    CHECK(curveIsMonotonic);

    if(FLAG_PLOT_CURVES){
     tendonCurve.printCurveToCSVFile(FILE_PATH,
                 "tendonCurve",
                 1.0-(e0/10),
                 1+e0);
    }
    //cout << "    passed" << endl;

}

TEST(activeForceLengthCurve)
{
    //cout << endl;
    //cout << endl;
    //cout <<"**************************************************"<<endl;
    //cout <<"FIBER ACTIVE FORCE LENGTH CURVE TESTING     "<<endl;
    double lce0 = 0.4;
    double lce1 = 0.75;
    double lce2 = 1;
    double lce3 = 1.6;
    double shoulderVal  = 0.05;
    double plateauSlope = 0.75;//0.75;
    double curviness    = 0.75;
    SmoothSegmentedFunction fiberfalCurve = SmoothSegmentedFunction();

    SmoothSegmentedFunctionFactory::
     createFiberActiveForceLengthCurve(lce0, lce1, lce2, lce3, 
          shoulderVal, plateauSlope, curviness,
          "test_fiberActiveForceLengthCurve", fiberfalCurve);
    //fiberfalCurve.printMuscleCurveToCSVFile(FILE_PATH);

    RigidBodyDynamics::Math::MatrixNd fiberfalCurveSample 
        = fiberfalCurve.calcSampledCurve(6,0,lce3);

    //0. Test that each curve fulfills its contract.
    //cout << "   Keypoint Testing" << endl;

    CHECK(abs(fiberfalCurve.calcValue(lce0) - shoulderVal) < TOL_SMALL);  
    CHECK(abs(fiberfalCurve.calcDerivative(lce0,1))     < TOL_BIG);
    CHECK(abs(fiberfalCurve.calcDerivative(lce0,2))     < TOL_BIG);
   
    //lce2 isn't the location of the end of a quintic Bezier curve
    //so I can't actually do any testing on this point.
    //SimTK_TEST_EQ_TOL(fiberfalCurve.calcValue(lce0),shoulderVal,TOL_SMALL);  
    //SimTK_TEST_EQ_TOL(
    //fiberfalCurve.calcDerivative(lce2,1),plateauSlope,TOL_BIG);
    //SimTK_TEST_EQ_TOL(fiberfalCurve.calcDerivative(lce2,2),0.0,TOL_BIG);

    CHECK(abs(fiberfalCurve.calcValue(lce2) - 1.0)    <  TOL_SMALL);  
    CHECK(abs(fiberfalCurve.calcDerivative(lce2,1))   <  TOL_BIG);
    CHECK(abs(fiberfalCurve.calcDerivative(lce2,2))   <  TOL_BIG);
    CHECK(abs(fiberfalCurve.calcValue(lce3)-shoulderVal) <  TOL_SMALL);
    CHECK(abs(fiberfalCurve.calcDerivative(lce3,1))   <  TOL_BIG);
    CHECK(abs(fiberfalCurve.calcDerivative(lce3,2))   <  TOL_BIG);

    //cout << "   passed" << endl;
    //cout << endl;
//1. Test each derivative sample for correctness against a numerically
//   computed version
    bool areCurveDerivativesGood = 
     areCurveDerivativesCloseToNumericDerivatives(
      fiberfalCurve,
      fiberfalCurveSample,
      TOL_DX);
    CHECK(areCurveDerivativesGood);
//2. Test each integral, where computed for correctness.
    //testMuscleCurveIntegral(fiberfalCurve,fiberfalCurveSample);

//3. Test numerically to see if the curve is C2 continuous
    bool curveIsContinuous = isCurveC2Continuous(fiberfalCurve,
                     fiberfalCurveSample);
    CHECK(curveIsContinuous);

    //fiberfalCurve.MuscleCurveToCSVFile("C:/mjhmilla/Stanford/dev");

    if(FLAG_PLOT_CURVES){
     fiberfalCurve.printCurveToCSVFile(FILE_PATH,
                  "fiberFalCurve",
                  0.3,
                  1.8);
    }




}

TEST(ForceVelocityCurve)
{
    //cout <<"**************************************************"<<endl;
    //cout <<"FIBER FORCE VELOCITY CURVE TESTING       "<<endl;

    double fmaxE     = 1.8;
    double dydxC     = 0.1;
    double dydxNearC    = 0.15;
    double dydxE     = 0.1;
    double dydxNearE    = 0.1+0.0001;
    double dydxIso   = 5;
    double concCurviness= 0.1;
    double eccCurviness = 0.75;

    SmoothSegmentedFunction fiberFVCurve = SmoothSegmentedFunction();
    SmoothSegmentedFunctionFactory::
     createFiberForceVelocityCurve(fmaxE, 
            dydxC,
            dydxNearC,
            dydxIso,
            dydxE,
            dydxNearE,
            concCurviness,
            eccCurviness,
            "test_fiberForceVelocityCurve",
            fiberFVCurve);
    //fiberFVCurve.printMuscleCurveToCSVFile(FILE_PATH);

    RigidBodyDynamics::Math::MatrixNd fiberFVCurveSample 
        = fiberFVCurve.calcSampledCurve(6,-1.0,1.0);

//0. Test that each curve fulfills its contract.
    //cout << "   Keypoint Testing" << endl;

    assert(abs(fiberFVCurve.calcValue(-1)      ) < TOL_SMALL);  
    assert(abs(fiberFVCurve.calcDerivative(-1,1)-dydxC  ) < TOL_BIG  );
    assert(abs(fiberFVCurve.calcDerivative(-1,2)     ) < TOL_BIG  );
    assert(abs(fiberFVCurve.calcValue(0)    -1.0  ) < TOL_SMALL);  
    assert(abs(fiberFVCurve.calcDerivative(0,1)-dydxIso ) < TOL_BIG  );
    assert(abs(fiberFVCurve.calcDerivative(0,2)      ) < TOL_BIG  );
    assert(abs(fiberFVCurve.calcValue(1)    -fmaxE   ) < TOL_SMALL);  
    assert(abs(fiberFVCurve.calcDerivative(1,1)-dydxE   ) < TOL_BIG  );
    assert(abs(fiberFVCurve.calcDerivative(1,2)      ) < TOL_BIG  );

    //cout << "   passed" << endl;
    //cout << endl;
//1. Test each derivative sample for correctness against a numerically
//   computed version
    bool areCurveDerivativesGood = 
     areCurveDerivativesCloseToNumericDerivatives(
      fiberFVCurve,
      fiberFVCurveSample,
      TOL_DX);
    CHECK(areCurveDerivativesGood);

//2. Test each integral, where computed for correctness.
    //testMuscleCurveIntegral(fiberFVCurve,fiberFVCurveSample);

//3. Test numerically to see if the curve is C2 continuous
    bool curveIsContinuous = isCurveC2Continuous(
              fiberFVCurve,
              fiberFVCurveSample);
    CHECK(curveIsContinuous);
//4. Test for montonicity where appropriate

    isCurveMontonic(fiberFVCurveSample);
    CHECK(curveIsContinuous);

    if(FLAG_PLOT_CURVES){
     fiberFVCurve.printCurveToCSVFile(FILE_PATH,
                 "fiberFvCurve",
                 -1.1,
                  1.1);
    }

    //cout << "    passed" << endl;

}

TEST(ForceVelocityInverseCurve)
{
    //cout <<"**************************************************"<<endl;
    //cout <<"FIBER FORCE VELOCITY INVERSE CURVE TESTING     "<<endl;
    double fmaxE     = 1.8;
    double dydxC     = 0.1;
    double dydxNearC    = 0.15;
    double dydxE     = 0.1;
    double dydxNearE    = 0.1+0.0001;
    double dydxIso   = 5;
    double concCurviness= 0.1;
    double eccCurviness = 0.75;

    SmoothSegmentedFunction fiberFVInvCurve = SmoothSegmentedFunction();    
    SmoothSegmentedFunctionFactory::
      createFiberForceVelocityInverseCurve( 
        fmaxE, 
        dydxC, 
        dydxNearC, 
        dydxIso, 
        dydxE, 
        dydxNearE,
        concCurviness,  
        eccCurviness,
        "test_fiberForceVelocityInverseCurve",
        fiberFVInvCurve);

    SmoothSegmentedFunction fiberFVCurve = SmoothSegmentedFunction();
    SmoothSegmentedFunctionFactory::
     createFiberForceVelocityCurve(
      fmaxE, 
      dydxC,
      dydxNearC,
      dydxIso,
      dydxE,
      dydxNearE,
      concCurviness,
      eccCurviness,
      "test_fiberForceVelocityCurve",
      fiberFVCurve);    
    //fiberFVInvCurve.printMuscleCurveToFile(FILE_PATH);
    

    RigidBodyDynamics::Math::MatrixNd fiberFVCurveSample
     = fiberFVCurve.calcSampledCurve(6, -1.0, 1.0);

    RigidBodyDynamics::Math::MatrixNd fiberFVInvCurveSample 
        = fiberFVInvCurve.calcSampledCurve(6,0,fmaxE);

//0. Test that each curve fulfills its contract.
    //cout << "   Keypoint Testing" << endl;

    CHECK(abs( fiberFVInvCurve.calcValue(0)    +1    ) < TOL_SMALL);
    CHECK(abs( fiberFVInvCurve.calcDerivative(0,1)-1/dydxC    ) < TOL_BIG);
    CHECK(abs( fiberFVInvCurve.calcDerivative(0,2)      ) < TOL_BIG);               
    CHECK(abs( fiberFVInvCurve.calcValue(1)          ) < TOL_SMALL);  
    CHECK(abs( fiberFVInvCurve.calcDerivative(1,1)-1/dydxIso  ) < TOL_BIG);
    CHECK(abs( fiberFVInvCurve.calcDerivative(1,2)      ) < TOL_BIG);                
    CHECK(abs( fiberFVInvCurve.calcValue(fmaxE)   -1    ) < TOL_SMALL);  
    CHECK(abs( fiberFVInvCurve.calcDerivative(fmaxE,1)-1/dydxE) < TOL_BIG);
    CHECK(abs( fiberFVInvCurve.calcDerivative(fmaxE,2)     ) <  TOL_BIG);

    //cout << "   passed" << endl;
    //cout << endl;
//1. Test each derivative sample for correctness against a numerically
//   computed version
    bool areCurveDerivativesGood = 
     areCurveDerivativesCloseToNumericDerivatives(
      fiberFVInvCurve,
      fiberFVInvCurveSample,
      TOL_DX);
    CHECK(areCurveDerivativesGood);

//2. Test each integral, where computed for correctness.
    //testMuscleCurveIntegral(fiberFVInvCurve,fiberFVInvCurveSample);

//3. Test numerically to see if the curve is C2 continuous
    bool curveIsContinuous = isCurveC2Continuous(
              fiberFVInvCurve,
              fiberFVInvCurveSample);
    CHECK(curveIsContinuous);
//4. Test for montonicity where appropriate

    bool curveIsMonotonic = isCurveMontonic(fiberFVInvCurveSample);
    CHECK(curveIsMonotonic);

//5. Testing the inverse of the curve - is it really an inverse?
    //cout << endl;
    //cout << "   TEST: Inverse correctness:fv(fvinv(fv)) = fv" << endl;
    



    double fv = 0;
    double dlce = 0;
    double fvCalc = 0;
    double fvErr = 0;
    double fvErrMax = 0;
    for(int i = 0; i < fiberFVInvCurveSample.rows(); i++){
     fv   = fiberFVCurveSample(i,0);
     dlce    = fiberFVInvCurve.calcValue(fv);
     fvCalc  = fiberFVCurve.calcValue(dlce);
     fvErr   = abs(fv-fvCalc);
     if(fvErrMax < fvErr)
      fvErrMax = fvErr;

     CHECK( fvErr < TOL_BIG);
    }

    if(FLAG_PLOT_CURVES){
     fiberFVInvCurve.printCurveToCSVFile(FILE_PATH,
                  "fiberFvInvCurve",
                  -0.1,
                  fmaxE+0.1);
    }

    //printf("   passed with a maximum error of %fe-12",fvErrMax*1e12);

}

TEST(passiveForceLengthCurve)
{
    double e0f   = 0.6;
    double kisof    = 8.389863790885878;
    double cf    = 0.65;
    double klow  = 0.5*(1.0/e0f);

    SmoothSegmentedFunction fiberFLCurve = SmoothSegmentedFunction();
    SmoothSegmentedFunctionFactory::
             createFiberForceLengthCurve(
              0.0, e0f,klow,kisof,cf,
             "test_fiberForceLength",
             fiberFLCurve);

    RigidBodyDynamics::Math::MatrixNd fiberFLCurveSample 
        = fiberFLCurve.calcSampledCurve(6,1.0,1.0+e0f);

//0. Test that each curve fulfills its contract.
    //cout << "   Keypoint Testing" << endl;

    RigidBodyDynamics::Math::VectorNd fiberFLCurveDomain 
    = fiberFLCurve.getCurveDomain();

    CHECK(abs(fiberFLCurveDomain(0) -1)     < TOL_SMALL);
    CHECK(abs(fiberFLCurveDomain(1) - (1+e0f)) < TOL_SMALL);
    CHECK(abs(fiberFLCurve.calcValue(1.0)  )   < TOL_SMALL);  
    CHECK(abs(fiberFLCurve.calcDerivative(1.0,1)) < TOL_BIG);
    CHECK(abs(fiberFLCurve.calcDerivative(1.0,2)) < TOL_BIG);
    CHECK(abs(fiberFLCurve.calcValue(1+e0f) -1.0)      < TOL_SMALL);  
    CHECK(abs(fiberFLCurve.calcDerivative(1+e0f,1)-kisof) < TOL_BIG);
    CHECK(abs(fiberFLCurve.calcDerivative(1+e0f,2))    < TOL_BIG);
    //cout << "   passed" << endl;
    //cout << endl;
//1. Test each derivative sample for correctness against a numerically
//   computed version
    bool areCurveDerivativesGood = 
     areCurveDerivativesCloseToNumericDerivatives(
      fiberFLCurve,
      fiberFLCurveSample,
      TOL_DX);
    CHECK(areCurveDerivativesGood);
//2. Test each integral, where computed for correctness.
    //testMuscleCurveIntegral(fiberFLCurve,fiberFLCurveSample);

//3. Test numerically to see if the curve is C2 continuous
    bool curveIsContinuous = isCurveC2Continuous(fiberFLCurve,
                      fiberFLCurveSample);
    CHECK(curveIsContinuous);
//4. Test for montonicity where appropriate

    bool curveIsMonotonic = isCurveMontonic(fiberFLCurveSample);
    CHECK(curveIsMonotonic);

    if(FLAG_PLOT_CURVES){
     fiberFLCurve.printCurveToCSVFile(FILE_PATH,
                 "fiberFLCurve",
                 1.0-0.1,
                 1.0+e0f+0.1);
    }
}

TEST(compressiveForceLengthCurve)
{
///////////////////////////////////////
//FIBER COMPRESSIVE FORCE LENGTH
///////////////////////////////////////
    //cout <<"**************************************************"<<endl;
    //cout <<"FIBER COMPRESSIVE FORCE LENGTH CURVE TESTING   "<<endl;


    double lmax = 0.6;
    double kce  = -8.389863790885878;
    double cce  = 0.5;//0.0;

    SmoothSegmentedFunction fiberCECurve = SmoothSegmentedFunction();
    SmoothSegmentedFunctionFactory::
    createFiberCompressiveForceLengthCurve(
      lmax,
      kce,
      cce,
      "test_fiberCompressiveForceLengthCurve",
      fiberCECurve);

    //fiberCECurve.printMuscleCurveToFile("C:/mjhmilla/Stanford/dev"
    //    "/OpenSim_LOCALPROJECTS/MuscleLibrary_Bench_20120210/build");
    RigidBodyDynamics::Math::MatrixNd fiberCECurveSample 
        = fiberCECurve.calcSampledCurve(6,0,lmax);

//0. Test that each curve fulfills its contract.
    //cout << "   Keypoint Testing" << endl;

    RigidBodyDynamics::Math::VectorNd fiberCECurveDomain 
     = fiberCECurve.getCurveDomain();
    CHECK(abs(fiberCECurveDomain(0))          < TOL_SMALL);
    CHECK(abs(fiberCECurveDomain(1)- lmax)       < TOL_SMALL);
    CHECK(abs(fiberCECurve.calcValue(lmax))      < TOL_SMALL);  
    CHECK(abs(fiberCECurve.calcDerivative(lmax,1))  < TOL_BIG);
    CHECK(abs(fiberCECurve.calcDerivative(lmax,2))  < TOL_BIG);
    CHECK(abs(fiberCECurve.calcValue(0) - 1.0)      < TOL_SMALL);  
    CHECK(abs(fiberCECurve.calcDerivative(0,1)-kce)    < TOL_BIG);
    CHECK(abs(fiberCECurve.calcDerivative(0,2))     < TOL_BIG);
    //cout << "   passed" << endl;
    //cout << endl;
//1. Test each derivative sample for correctness against a numerically
//   computed version
    bool areCurveDerivativesGood = 
     areCurveDerivativesCloseToNumericDerivatives(
      fiberCECurve,
      fiberCECurveSample,
      TOL_DX);
    CHECK(areCurveDerivativesGood);

//2. Test each integral, where computed for correctness.
    //testMuscleCurveIntegral(fiberCECurve,fiberCECurveSample);

//3. Test numerically to see if the curve is C2 continuous
    bool curveIsContinuous = isCurveC2Continuous(fiberCECurve,
                      fiberCECurveSample);
    CHECK(curveIsContinuous);
//4. Test for montonicity where appropriate

    bool curveIsMonotonic = isCurveMontonic(fiberCECurveSample);
    CHECK(curveIsMonotonic);

//5. Testing Exceptions
    //cout << endl;
    //cout << "   Exception Testing" << endl;

    if(FLAG_PLOT_CURVES){
     fiberCECurve.printCurveToCSVFile(FILE_PATH,
                 "fiberCECurve",
                 -0.1,
                 lmax + 0.1);
    }

    //cout << "    passed" << endl;

}

TEST(compressivePhiCurve)
{
    //cout <<"**************************************************"<<endl;
    //cout <<"FIBER COMPRESSIVE FORCE PHI CURVE TESTING   "<<endl;

    double phi0 = (M_PI/2)*(1.0/2.0);
    double phi1 = M_PI/2;
    double kphi  = 8.389863790885878;
    double cphi  = 0.0;  

    SmoothSegmentedFunction fiberCEPhiCurve = SmoothSegmentedFunction();
    SmoothSegmentedFunctionFactory::
     createFiberCompressiveForcePennationCurve(phi0,kphi,cphi,
         "test_fiberCompressiveForcePennationCurve", fiberCEPhiCurve);    


    RigidBodyDynamics::Math::MatrixNd fiberCEPhiCurveSample 
        = fiberCEPhiCurve.calcSampledCurve(6,phi0,phi1);

//0. Test that each curve fulfills its contract.
    //cout << "   Keypoint Testing" << endl;

    CHECK(abs(fiberCEPhiCurve.calcValue(phi0))     < TOL_SMALL);  
    CHECK(abs(fiberCEPhiCurve.calcDerivative(phi0,1)) < TOL_BIG);
    CHECK(abs(fiberCEPhiCurve.calcDerivative(phi0,2)) < TOL_BIG);
    CHECK(abs(fiberCEPhiCurve.calcValue(phi1)) -1     < TOL_SMALL);
    CHECK(abs(fiberCEPhiCurve.calcDerivative(phi1,1)-kphi)  < TOL_BIG);
    CHECK(abs(fiberCEPhiCurve.calcDerivative(phi1,2))    < TOL_BIG);
    //cout << "   passed" << endl;
    //cout << endl;
//1. Test each derivative sample for correctness against a numerically
//   computed version
    bool areCurveDerivativesGood = 
     areCurveDerivativesCloseToNumericDerivatives(
      fiberCEPhiCurve,
      fiberCEPhiCurveSample,
      TOL_DX);
    CHECK(areCurveDerivativesGood);
//2. Test each integral, where computed for correctness.
    //testMuscleCurveIntegral(fiberCEPhiCurve,fiberCEPhiCurveSample);

//3. Test numerically to see if the curve is C2 continuous
    //cout << "Attention: Removed test for the Cos Phi Compressive Curve"<<endl;
    bool curveIsContinuous = isCurveC2Continuous(fiberCEPhiCurve,
                      fiberCEPhiCurveSample);
    CHECK(curveIsContinuous);
//4. Test for montonicity where appropriate
    bool curveIsMonotonic = isCurveMontonic(fiberCEPhiCurveSample);
    CHECK(curveIsMonotonic);


    if(FLAG_PLOT_CURVES){
     fiberCEPhiCurve.printCurveToCSVFile(FILE_PATH,
                 "fiberCEPhiCurve",
                 phi0-0.1,
                 phi1+0.1);
    }

    //cout << "    passed" << endl;

}

TEST(compressiveCosPhiCurve)
{
    //cout <<"**************************************************"<<endl;
    //cout <<"FIBER COMPRESSIVE FORCE COSPHI CURVE TESTING   "<<endl;

    double cosPhi0 = cos( (80.0/90.0)*M_PI/2);
    double kcosPhi  = -1.2/(cosPhi0);
    double ccosPhi  = 0.5;
    SmoothSegmentedFunction fiberCECosPhiCurve = SmoothSegmentedFunction();

    SmoothSegmentedFunctionFactory::
        createFiberCompressiveForceCosPennationCurve(
           cosPhi0,kcosPhi,ccosPhi,
           "test_fiberCompressiveForceCosPennationCurve",
           fiberCECosPhiCurve);

    RigidBodyDynamics::Math::MatrixNd fiberCECosPhiCurveSample 
        = fiberCECosPhiCurve.calcSampledCurve(6,0,cosPhi0);

//0. Test that each curve fulfills its contract.
    //cout << "   Keypoint Testing" << endl;

    CHECK(abs(fiberCECosPhiCurve.calcValue(cosPhi0)    )  < TOL_SMALL);  
    CHECK(abs(fiberCECosPhiCurve.calcDerivative(cosPhi0,1)   )  < TOL_BIG);
    CHECK(abs(fiberCECosPhiCurve.calcDerivative(cosPhi0,2)   )  < TOL_BIG);
    CHECK(abs(fiberCECosPhiCurve.calcValue(0)    - 1.0    )  < TOL_SMALL);  
    CHECK(abs(fiberCECosPhiCurve.calcDerivative(0,1) -kcosPhi)  < TOL_BIG);
    CHECK(abs(fiberCECosPhiCurve.calcDerivative(0,2)      )  < TOL_BIG);
    //cout << "   passed" << endl;
    //cout << endl;
//1. Test each derivative sample for correctness against a numerically
//   computed version
    bool areCurveDerivativesGood = 
     areCurveDerivativesCloseToNumericDerivatives(
      fiberCECosPhiCurve,
      fiberCECosPhiCurveSample,
      TOL_DX);
    CHECK(areCurveDerivativesGood);
//2. Test each integral, where computed for correctness.
    //testMuscleCurveIntegral(fiberCECosPhiCurve,fiberCECosPhiCurveSample);

//3. Test numerically to see if the curve is C2 continuous

    bool curveIsContinuous = isCurveC2Continuous(
              fiberCECosPhiCurve,
              fiberCECosPhiCurveSample);
    CHECK(curveIsContinuous);
//4. Test for montonicity where appropriate

    bool curveIsMonotonic = isCurveMontonic(fiberCECosPhiCurveSample);
    CHECK(curveIsMonotonic);


    //cout << "    passed" << endl;
    if(FLAG_PLOT_CURVES){
     fiberCECosPhiCurve.printCurveToCSVFile(FILE_PATH,
                 "fiberCECosPhiCurve",
                 -0.1,
                 cosPhi0+0.1);
    }

}

TEST(Anderson2007ActiveTorqueAngleCurve)
{
    double subjectWeight    = 75.0*9.81;
    double subjectHeight    = 1.75;
    double scale      = subjectHeight*subjectWeight; 

    //These parameters are taken from table 3 for hip extension for
    //men between the ages of 18-25
    double c1      = 0.161; //normalized maximum hip joint torque
    double c2      = 0.958; // pi/(theta_max - theta_min)
    double c3      = 0.932; //theta_max_iso_torque
    double c4      = 1.578; //omega_1: angular velocity at 75% tq_iso_max
    double c5      = 3.190; //omega_2: angular velocity at 50% tq_iso_max
    double c6      = 0.242; //E, where eccentric slope = (1+E)*concentricSlope
                //Passive torque angle curve parameters
    double b1      =-1.210; // torque_passive = b1*exp(k1*theta) 
    double k1      =-6.351; //        +b2*exp(k2*theta)
    double b2      = 0.476;            
    double k2      = 5.910;



    //cout <<endl;      
    //cout <<endl;
    //cout <<"**************************************************"<<endl;
    //cout <<"ANDERSON 2007 ACTIVE TORQUE ANGLE CURVE TESTING   "<<endl;

    SmoothSegmentedFunction andersonTaCurve = SmoothSegmentedFunction();
    SmoothSegmentedFunctionFactory::
     createAnderson2007ActiveTorqueAngleCurve(c2,c3,
          "test_Anderson2007TorqueAngleCurve",
          andersonTaCurve);      

    double angleRange    = (M_PI/c4); 
    double angleActiveMin   = c3 - angleRange*0.75;
    double angleActiveMax   = c3 + angleRange*0.75;

    RigidBodyDynamics::Math::MatrixNd andersonTaCurveSample 
        = andersonTaCurve.calcSampledCurve( 6,
                       angleActiveMin,
                       angleActiveMax);    

    //cout << "   Keypoint Testing" << endl;

    CHECK(abs(andersonTaCurve.calcValue(c3) - 1.0) < TOL_SMALL);  
    CHECK(abs(andersonTaCurve.calcDerivative(c3,1))     < TOL_BIG);
    CHECK(abs(andersonTaCurve.calcDerivative(c3,2))     < TOL_BIG);

    RigidBodyDynamics::Math::VectorNd curveDomain 
     = andersonTaCurve.getCurveDomain();

    CHECK(abs(andersonTaCurve.calcValue(curveDomain(0)))     < TOL_SMALL);
    CHECK(abs(andersonTaCurve.calcDerivative(curveDomain(0),1))    < TOL_BIG);
    CHECK(abs(andersonTaCurve.calcValue(curveDomain(1)))     < TOL_SMALL);
    CHECK(abs(andersonTaCurve.calcDerivative(curveDomain(1),1))    < TOL_BIG);
    //cout << "   passed " << endl;

    //cout << "    Continuity and Smoothness Testing" << endl;
    bool areCurveDerivativesGood = 
     areCurveDerivativesCloseToNumericDerivatives(
      andersonTaCurve,  
      andersonTaCurveSample, 
      TOL_DX);
    CHECK(areCurveDerivativesGood);

    bool curveIsContinuous = isCurveC2Continuous(andersonTaCurve, 
                      andersonTaCurveSample);
    CHECK(curveIsContinuous);

    if(FLAG_PLOT_CURVES){
     andersonTaCurve.printCurveToCSVFile(
        FILE_PATH,
        "anderson2007ActiveTorqueAngleCurve",
        angleActiveMin,
        angleActiveMax);
    }


}

TEST(Anderson2007PassiveTorqueAngleCurve)
{
    double subjectWeight    = 75.0*9.81;
    double subjectHeight    = 1.75;
    double scale      = subjectHeight*subjectWeight; 

    //These parameters are taken from table 3 for hip extension for
    //men between the ages of 18-25
    double c1      = 0.161; //normalized maximum hip joint torque
    double c2      = 0.958; // pi/(theta_max - theta_min)
    double c3      = 0.932; //theta_max_iso_torque
    double c4      = 1.578; //omega_1: angular velocity at 75% tq_iso_max
    double c5      = 3.190; //omega_2: angular velocity at 50% tq_iso_max
    double c6      = 0.242; //E, where eccentric slope = (1+E)*concentricSlope
                //Passive torque angle curve parameters
    double b1      =-1.210; // torque_passive = b1*exp(k1*theta) 
    double k1      =-6.351; //        +b2*exp(k2*theta)
    double b2      = 0.476;            
    double k2      = 5.910;

    
    //cout <<endl;
    //cout <<endl;
    //cout <<"**************************************************"<<endl;
    //cout <<"ANDERSON 2007 PASSIVE TORQUE ANGLE CURVE TESTING   "<<endl;

    double curveSign = 1.0;



    for(int z = 0; z<2; ++z){

     if(z == 0){
      curveSign = 1.0;
      //cout <<"    TESTING SIDE 1"<<endl;
     }else{
      curveSign = -1.0;
      //cout <<"    TESTING SIDE 2"<<endl;

     }
     SmoothSegmentedFunction andersonTpCurve = SmoothSegmentedFunction();
     SmoothSegmentedFunctionFactory::
       createAnderson2007PassiveTorqueAngleCurve(
            scale,
            c1,
            b1*curveSign,
            k1,
            b2*curveSign,
            k2,
            "test_passiveTorqueAngleCurve",
            andersonTpCurve);

     RigidBodyDynamics::Math::VectorNd curveDomain 
      = andersonTpCurve.getCurveDomain();

     double angleMin = curveDomain(0);
     double angleMax = curveDomain(1);

     RigidBodyDynamics::Math::MatrixNd andersonTpCurveSample
            = andersonTpCurve.calcSampledCurve( 6,
                        angleMin-0.1,
                        angleMax+0.1);

     //cout << "   Keypoint Testing" << endl;

     double tauMin = andersonTpCurve.calcValue(angleMin);
     double tauMax = andersonTpCurve.calcValue(angleMax);
     double tauMinAngle = angleMin;

     if(tauMin > tauMax){
      double tmp  = tauMin;
      tauMin   = tauMax;
      tauMax   = tmp;
      tauMinAngle = angleMax;
     }

     CHECK( abs(tauMin) < TOL_SMALL);
     CHECK( abs(tauMax - 1.0) < TOL_SMALL);
     CHECK( abs(andersonTpCurve.calcDerivative(tauMinAngle,1)) < TOL_SMALL);

     //cout << "   passed " << endl;

     //cout << "   Continuity and Smoothness Testing " << endl;
     bool areCurveDerivativesGood = 
      areCurveDerivativesCloseToNumericDerivatives( 
       andersonTpCurve,  
       andersonTpCurveSample, 
       TOL_DX);
     CHECK(areCurveDerivativesGood);

     bool curveIsContinuous = isCurveC2Continuous(andersonTpCurve,  
                       andersonTpCurveSample);
     CHECK(curveIsContinuous);

     bool curveIsMonotonic = isCurveMontonic(andersonTpCurveSample);
     CHECK(curveIsMonotonic);
     //cout << "   passed " << endl;


    }

    SmoothSegmentedFunction andersonTpCurve = SmoothSegmentedFunction();
    SmoothSegmentedFunctionFactory::
      createAnderson2007PassiveTorqueAngleCurve(
        scale,
        c1,
        b1,
        k1,
        b2,
        k2,
        "test_passiveTorqueAngleCurve",
        andersonTpCurve);

    if(FLAG_PLOT_CURVES){
     andersonTpCurve.printCurveToCSVFile(
        FILE_PATH,
        "anderson2007PassiveTorqueAngleCurve",
        andersonTpCurve.getCurveDomain()[0]-0.1,
        andersonTpCurve.getCurveDomain()[1]+0.1);
    }

}

TEST(Anderson2007ActiveTorqueVelocityCurve)
{
    double subjectWeight    = 75.0*9.81;
    double subjectHeight    = 1.75;
    double scale      = subjectHeight*subjectWeight; 

    //These parameters are taken from table 3 for hip extension for
    //men between the ages of 18-25
    double c1      = 0.161; //normalized maximum hip joint torque
    double c2      = 0.958; // pi/(theta_max - theta_min)
    double c3      = 0.932; //theta_max_iso_torque
    double c4      = 1.578; //omega_1: angular velocity at 75% tq_iso_max
    double c5      = 3.190; //omega_2: angular velocity at 50% tq_iso_max
    double c6      = 0.242; //E, where eccentric slope = (1+E)*concentricSlope
                //Passive torque angle curve parameters
    double b1      =-1.210; // torque_passive = b1*exp(k1*theta) 
    double k1      =-6.351; //        +b2*exp(k2*theta)
    double b2      = 0.476;            
    double k2      = 5.910;

    //cout <<endl;
    //cout <<endl;
    //cout <<"**************************************************"<<endl;
    //cout <<"ANDERSON 2007 ACTIVE TORQUE VELOCITY CURVE TESTING"<<endl;

    double minEccentricMultiplier = 1.1;
    double maxEccentricMultiplier = 1.4;

    SmoothSegmentedFunction andersonTvCurve = SmoothSegmentedFunction();
    SmoothSegmentedFunctionFactory::
    createAnderson2007ActiveTorqueVelocityCurve(
     c4,c5,c6,minEccentricMultiplier,maxEccentricMultiplier,
     "test_anderson2007ActiveTorqueVelocityCurve",
     andersonTvCurve);

    RigidBodyDynamics::Math::VectorNd curveDomain 
     = andersonTvCurve.getCurveDomain();

    double angularVelocityMin = curveDomain(0);
    double angularVelocityMax = curveDomain(1);


    RigidBodyDynamics::Math::MatrixNd andersonTvCurveSample
        = andersonTvCurve.calcSampledCurve( 6,
                       angularVelocityMin-0.1,
                       angularVelocityMax+0.1);

    //cout << "   Keypoint Testing" << endl;

    CHECK(abs(andersonTvCurve.calcValue(0) - 1.0)        < TOL_SMALL);
    CHECK(abs(andersonTvCurve.calcValue(c4) - 0.75)    < TOL_BIG);
    CHECK(abs(andersonTvCurve.calcValue(c5) - 0.5)     < TOL_BIG);
    CHECK(abs(andersonTvCurve.calcValue(angularVelocityMax))   < TOL_BIG);

    double maxTv = andersonTvCurve.calcValue(angularVelocityMin);
    CHECK(maxTv >= minEccentricMultiplier);
    CHECK(maxTv <= maxEccentricMultiplier);

    CHECK(abs(andersonTvCurve.calcDerivative(angularVelocityMin,1))<TOL_SMALL);
    CHECK(abs(andersonTvCurve.calcDerivative(angularVelocityMax,1))<TOL_SMALL);

    //cout << "   passed " << endl;
    //cout << "    Continuity and Smoothness Testing" << endl;

    bool areCurveDerivativesGood = 
     areCurveDerivativesCloseToNumericDerivatives(
      andersonTvCurve,  
      andersonTvCurveSample, 
      TOL_DX);

    CHECK(areCurveDerivativesGood);

    bool curveIsContinuous = isCurveC2Continuous(andersonTvCurve,
                      andersonTvCurveSample);
    CHECK(curveIsContinuous);

    bool curveIsMonotonic = isCurveMontonic(andersonTvCurveSample);
    CHECK(curveIsMonotonic);


    if(FLAG_PLOT_CURVES){
     andersonTvCurve.printCurveToCSVFile(
        FILE_PATH,
        "anderson2007ActiveTorqueVelocityCurve",
        angularVelocityMin,
        angularVelocityMax);
    }

}


TEST(NeuroSwingTorqueAngleCurve)
{

    //cout <<"**************************************************"<<endl;
   //cout <<"NeuroSwing TORQUE ANGLE CURVE TESTING   "<<endl;

    //Coefficients that Manish extracted out of a manual
    //testing trial completed by Julia at Schlearbach.
    RigidBodyDynamics::Math::VectorNd controlPointAngles(5);
    RigidBodyDynamics::Math::VectorNd controlPointAnglesEX(5);
    RigidBodyDynamics::Math::VectorNd controlPointTorques(5);
    RigidBodyDynamics::Math::VectorNd controlPointTorquesEX(5);
    RigidBodyDynamics::Math::VectorNd controlPointSlopes(5);
    RigidBodyDynamics::Math::VectorNd controlPointSlopesEX(5);

    RigidBodyDynamics::Math::VectorNd xpoints;
    xpoints = RigidBodyDynamics::Math::VectorNd::Zero (6);

    double centerOffsetInRadians = -1.745329e-02;
    double hardStopPlantarFlexion = 9.599311e-02;
    double hardStopDorsiFlexion = -1.308997e-01;
    double orthosisStiffnessDorsiFlexionEndRange = 100;
    double orthosisStiffnessPlantarFlexionEndRange = 100;
    double torquePreloadInNm = 1;
    double torquePreloadRangeInRadians = 1.000000e-02;
    double springStiffnessForefoot = 35;
    double springStiffnessHeel = 15;

    xpoints[0] = hardStopDorsiFlexion-0.05;
    xpoints[1] = hardStopDorsiFlexion;
    xpoints[2] = -torquePreloadRangeInRadians + centerOffsetInRadians;
    xpoints[3] = torquePreloadRangeInRadians + centerOffsetInRadians;
    xpoints[4] = hardStopPlantarFlexion;
    xpoints[5] = hardStopPlantarFlexion+0.05;

    RigidBodyDynamics::Math::VectorNd ypoints;
    ypoints = RigidBodyDynamics::Math::VectorNd::Zero (6);
    ypoints[0] = torquePreloadInNm + springStiffnessForefoot*(fabs(xpoints[2] - xpoints[1])) + orthosisStiffnessDorsiFlexionEndRange*0.05;
    ypoints[1] = torquePreloadInNm + springStiffnessForefoot*(fabs(xpoints[2] - xpoints[1]));
    ypoints[2] = torquePreloadInNm;
    ypoints[3] = -torquePreloadInNm;
    ypoints[4] = -torquePreloadInNm - springStiffnessHeel*(fabs(xpoints[4]-xpoints[3]));
    ypoints[5] = -torquePreloadInNm - springStiffnessHeel*(fabs(xpoints[4]-xpoints[3])) - orthosisStiffnessPlantarFlexionEndRange*0.05;

    controlPointAngles[0] = (xpoints[0] + xpoints[1])/2;
    controlPointAngles[1] = (xpoints[1] + xpoints[2])/2;
    controlPointAngles[2] = centerOffsetInRadians;
    controlPointAngles[3] = (xpoints[3] + xpoints[4])/2;
    controlPointAngles[4] = (xpoints[4] + xpoints[5])/2;

    controlPointTorques[0] = (ypoints[0] + ypoints[1])/2;
    controlPointTorques[1] = (ypoints[1] + ypoints[2])/2;
    controlPointTorques[2] = 0.0;
    controlPointTorques[3] = (ypoints[3] + ypoints[4])/2;
    controlPointTorques[4] = (ypoints[4] + ypoints[5])/2;

    controlPointSlopes[0] = (ypoints[1] - ypoints[0])/(xpoints[1] - xpoints[0]);
    controlPointSlopes[1] = (ypoints[2] - ypoints[1])/(xpoints[2] - xpoints[1]);
    controlPointSlopes[2] = (ypoints[3] - ypoints[2])/(xpoints[3] - xpoints[2]);
    controlPointSlopes[3] = (ypoints[4] - ypoints[3])/(xpoints[4] - xpoints[3]);
    controlPointSlopes[4] = (ypoints[5] - ypoints[4])/(xpoints[5] - xpoints[4]);

    controlPointAnglesEX = controlPointAngles;
    controlPointTorquesEX = controlPointTorques;
    controlPointSlopesEX = controlPointSlopes;

    SmoothSegmentedFunction neuroSwingCurve = SmoothSegmentedFunction();
    SmoothSegmentedFunctionFactory::
      createNeuroSwingSpringProfile(  controlPointAngles,
                    controlPointTorques,
                    controlPointSlopes,
                    "neuroSwingSpring",
                    neuroSwingCurve);
    
    RigidBodyDynamics::Math::VectorNd neuroSwingDomain =
      neuroSwingCurve.getCurveDomain();

    RigidBodyDynamics::Math::MatrixNd neuroSwingCurveSample
        = neuroSwingCurve.calcSampledCurve( 
            6,
            neuroSwingDomain(0) - 0.1,
            neuroSwingDomain(1) + 0.1);

    //cout << "   Keypoint Testing" << endl;
    double angle = 0.;

    for(int i=0; i<5;++i){
     angle = controlPointAngles(i);
     CHECK(abs(  neuroSwingCurve.calcValue(angle)
          - controlPointTorques(i))    < TOL_SMALL);
     CHECK(abs(  neuroSwingCurve.calcDerivative(angle,1)
          - controlPointSlopes(i))    < TOL_BIG);
    }

    //cout << "   passed " << endl;

    //cout << "   Continuity and Smoothness Testing " << endl;
    bool areCurveDerivativesGood = 
     areCurveDerivativesCloseToNumericDerivatives( 
      neuroSwingCurve,  
      neuroSwingCurveSample, 
      TOL_DX_BIG);
    CHECK(areCurveDerivativesGood);

    bool curveIsContinuous = isCurveC2Continuous(neuroSwingCurve,  
                      neuroSwingCurveSample);
    CHECK(curveIsContinuous);

    //if(FLAG_PLOT_CURVES){
     neuroSwingCurve.printCurveToCSVFile(FILE_PATH,
                    "neuroSwingTorqueAngleCurve",
                    neuroSwingDomain(0)-0.1,
                    neuroSwingDomain(1)+0.1);
    //}
}


int main (int argc, char *argv[])
{
    return UnitTest::RunAllTests ();
}
