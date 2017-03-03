#ifndef NEUROSWINGTORSIONSPRING_H_
#define NEUROSWINGTORSIONSPRING_H_

/* 
 * RBDL - Rigid Body Dynamics Library: Addon : forceElements
 * Copyright (c) 2016 Matthew Millard <matthew.millard@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#include <vector>
#include <rbdl/rbdl_math.h>
#include "../geometry/SmoothSegmentedFunction.h"
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>

namespace RigidBodyDynamics {
    namespace Addons {
      namespace CustomForces{


      class NeuroSwingTorsionSpring {

        public:
        /**
            Creates an empty NeuroSwingTorsionSpring. 
        */
        NeuroSwingTorsionSpring();


        /**
        This class implements torsion spring profile that has the characteristic
        shape of an ankle-foot-orthoses (AFO) that has a Neuroswing ankle joint.
        The profile curve is represented using \f$5^{th}\f$ order 2D Bezier 
        splines which are continuous to the \f$2^{nd}\f$ derivative. 

        \image html fig_NeuroSwingTorsionSpring.png "NeuroSwingTorsionSpring profile parameters" 

        <b>Coordinate Mapping</b>

        The variables jointAngleOffset, jointAngleSign, jointTorqueSign are used
        so that the coordinate system used by your model's ankle joint (\f$q\f$) 
        can be mapped to the coordinate system displayed in the figure. 
        The angles are mapped using this linear equation

        \f[
            neuroSwingAngle =  jointAngleSign*(q-jointAngleOffset)
        \f]

        while Neuroswing torques are converted to joint torques (\f$\tau\f$) 
        simply using

        \f[
            \tau =  jointTorqueSign * neuroSwingTorque
        \f]


        @param angleHardStopInPlantarFlexion
                The maximum dorsiflexion angle permitted by the Neuroswing
                (radians)

        @param angleCenterOffset            
                The angle when the torque provided by the Neuroswing goes to 
                zero (radians). 

        @param angleHardStopInDorsiFlexion
                The dorsiflextion angle at which the NeuroSwing joint hits its
                hardstop and the stiffness transitions from the stiffness
                of the NeuroSwing joint to that of the orthoses.

        @param stiffnessEndRangeInPlantarFlexion
                The stiffness of the ankle-foot-orthoses when the Neuroswing
                joint has hit the plantarflextion hard-stop. Note that this
                stiffness is due to the material used to make the AFO.

        @param stiffnessInPlantarFlexion
                The stiffness of the Neuroswing spring when the ankle angle
                is between angleCenterOffset+0.5*anglePreloadWindow and
                angleHardStopInPlantarFlexion  (Nm).


        @param stiffnessInDorsiFlexion
                The stiffness of the Neuroswing spring when the ankle angle
                is between the angleHardStopInDorsiFlexion and
                angleCenterOffset - 0.5*anglePreloadWindow (Nm).


        @param stiffnessEndRangeInDorsiFlexion
                The stiffness of the ankle-foot-orthoses when the Neuroswing
                joint has hit the dorsiflexion end stop. Note that this 
                stiffness is due to the material used to make the AFO (Nm).


        @param anglePreloadWindow
                The angular width of the window in which the spring's torque
                transitions from a value of negative preload to positive
                preload (radians).


        @param torquePreload
                The torque preload of the Neuroswing joint when it is at an
                angle of (angleCenterOffset +/- 0.5*anglePreloadWindow).


        @param jointAngleOffset
                Offset angle between your model's ankle joint and the
                curve shown in the figure. 

        @param jointAngleSign
                The sign convention that converts your model's ankle joint
                angles to the convention shown in the figure. This variable 
                must be either -1 or 1.

        @param jointTorqueSign
                The sign that maps torque convention illustrated in the figure
                to the correct convention for your model's ankle joint. This
                variable must be either -1 or 1.

        @param springName
                The name of the spring. This is needed to do useful
                things like provide error messages that are human
                readable.

        */
        NeuroSwingTorsionSpring( 
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
             const std::string& springName);



        double getAngleHardStopInDorsiFlexion() const;    
        double getAngleCenterOffset() const;               
        double getAnglePreloadWindow() const;            
        double getAngleHardStopInPlantarFlexion() const; 
        double getStiffnessEndRangeInDorsiFlexion() const;
        double getStiffnessInDorsiFlexion() const;       
        double getTorquePreload() const;                  
        double getStiffnessInPlantarFlexion() const;      
        double getStiffnessEndRangeInPlantarFlexion() const;
        double getJointAngleOffset() const;              
        double getJointAngleSign() const;                 
        double getJointTorqueSign() const;

        void setAngleHardStopInDorsiFlexion(double val);   
        void setAngleCenterOffset(double val);             
        void setAnglePreloadWindow(double val);            
        void setAngleHardStopInPlantarFlexion(double val); 
        void setStiffnessEndRangeInDorsiFlexion(double val);
        void setStiffnessInDorsiFlexion(double val);       
        void setTorquePreload(double val);                 
        void setStiffnessInPlantarFlexion(double val);     
        void setStiffnessEndRangeInPlantarFlexion(double val);
        void setJointAngleOffset(double val);              
        void setJointAngleSign(double val);                
        void setJointTorqueSign(double val);

        std::string getName() const;
        void setName(std::string& name);


        /**
            @param jointAngle the angle of the joint that has the NeuroSwing
                   spring on it (radians)
            @return the joint torque developed by the NeuroSwing spring (Nm)
        */
        double calcJointTorque(                        
                double jointAngle) const;

        /**
            @param jointAngle the angle of the joint that has the NeuroSwing
                   spring on it (radians)
            @return the joint stiffness of the NeuroSwing spring (Nm/rad)
        */
        double calcJointStiffness(
                double jointAngle) const;

        /**
            This function will densely sample the NeuroSwing profile and will
            write the results to a csv file. The file has 6 columns: angle,
            torque, stiffness, jointAngle, jointTorque, jointStiffness. Here
            angle, torque, and stiffness are reported using the convention
            angle, torque and stiffness convention shown in the image that is 
            in the class description. The 'joint' quantities are resolved using
            the values of jointAngleOffset, jointAngleSign and jointTorqueSign.


            @param path : the path to the folder you want the the file stored on
            @return fileNameWithoutExtension : the name of the file but without
                    any extension.
        */
        void printJointTorqueProfileToFile(
                const std::string& path,
                const std::string& fileNameWithoutExtension) const;

        /**
           @return A const reference to the underlying SmoothSegmentedFunction
        */
        const RigidBodyDynamics::Addons::Geometry::SmoothSegmentedFunction&
        getSmoothSegmentedFunction() const;


        private:
                    
        /**Rebuilds the entire curve*/
        void buildNeuroSwingTorsionSpringCurve();
        /**Prints all of the member variable names and values using Cerr*/
        void printCurveParametersUsingCerr();

        /**Maps from the jointAngle to the NeuroSwing angle*/
        double calcNeuroSwingAngle(double jointAngle) const;
        /**Maps from the NeuroSwing angle to the jointAngle*/
        double calcJointAngle(double neuroSwingAngle) const;


        /**The quintic Bezier spline that represents this curve*/
        RigidBodyDynamics::Addons::Geometry::SmoothSegmentedFunction 
                _neuroSwingCurve;

        double _jointAngleOffset;
        double _jointAngleSign;
        double _jointTorqueSign;  

        double _angleHardStopInDorsiFlexion;                                 
        double _angleCenterOffset;                                            
        double _anglePreloadWindow;                                          
        double _angleHardStopInPlantarFlexion;                              
        double _stiffnessEndRangeInDorsiFlexion;                             
        double _stiffnessInDorsiFlexion;                                    
        double _torquePreload;                                               
        double _stiffnessInPlantarFlexion;                                   
        double _stiffnessEndRangeInPlantarFlexion;                             
        
        /**If the parameters of the curve have been changed this is set
        to true. Once the command buildNeuroSwingTorsionSpringCurve has
        been called this parameter is set to false*/
        bool _parametersAreDirty;

        std::string _springName;

    };



}
}
}

#endif
