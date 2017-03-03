#include <vector>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include <btkAcquisitionFileReader.h>	
#include <btkAcquisition.h>
#include <btkForcePlatformsExtractor.h>
#include <btkGroundReactionWrenchFilter.h>
#include "btkProcessObject.h"
#include "btkForcePlatformCollection.h"
#include "btkForcePlatformTypes.h"
#include <btkMergeAcquisitionFilter.h>
#include <btkAcquisitionUnitConverter.h>
#include <btkConfigure.h>
#include "btkMetaData.h"
#include "btkConvert.h"
#include <btkMacro.h>

#include <rbdl/rbdl.h>
#include "rbdl/Logging.h"
#include "addons/luamodel/luamodel.h"
#include "addons/luamodel/luatables.h"
#include <algorithm>   

using namespace std;
using namespace boost;
using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;

class ReadAnimationCSV {                  
	vector<VectorNd> Q;
	const char* FileName;
	int DegreesOfFreedom;
	double StartTime;    
	double EndTime;
	int TimeSteps;

	// Three variables to localize the header. The headers for Meshup typically end with the word "DATA"
	// If the word "DATA" is present, its position is detected and the text only readed after this position
	int LineWith_DATA;
	bool LineWith_DATAFOUND;
	string Line;	

	public:
	  ReadAnimationCSV(const char* NameOfFile);
	  void Read();
	  int GetTimeSteps();
	  int GetDegreesOfFreedom();
	  double GetStartTime();
	  double GetEndTime();
	  vector<VectorNd> GetQ();   

	private:	
	  void Initiate();        
	  template<class string> void tokenizeV(const std::string &s,std::vector<string> &o);
};

ReadAnimationCSV::ReadAnimationCSV(const char* NameOfFile) { FileName = NameOfFile;}
		
void ReadAnimationCSV::Read() {
	Initiate(); 
	ifstream Animationfile(FileName);

	if(LineWith_DATAFOUND) {
		for(int i=0;i<=LineWith_DATA;i++)
			getline(Animationfile, Line);
	} 

	getline(Animationfile, Line);
	std::vector<double> VectorFirstRow;
	tokenizeV(Line, VectorFirstRow);
	StartTime = VectorFirstRow[0];
	DegreesOfFreedom = VectorFirstRow.size()-1;
	Q.push_back(VectorNd::Zero(DegreesOfFreedom));
	for(int i=0;i<DegreesOfFreedom;i++) 
	  Q[0][i]=VectorFirstRow[i+1]; 
	TimeSteps=1; 

	while(getline(Animationfile, Line)) {   
	  std::vector<double> Row;
	  tokenizeV(Line, Row);
	  EndTime=Row[0];
	  Q.push_back(VectorNd::Zero(DegreesOfFreedom));
	  for(int i=0;i<DegreesOfFreedom;i++) 
	    Q[TimeSteps][i]=Row[i+1];
	  TimeSteps++;
	}
}

void ReadAnimationCSV::Initiate() {
	ifstream Animationfile(FileName);
	LineWith_DATA = 0;
	LineWith_DATAFOUND = false;
	while(getline(Animationfile, Line)) {  
		if(!strstr(&Line[0],"DATA")) {
			LineWith_DATA++;}
		else {
		  LineWith_DATAFOUND=true;
		  break;
		}
	}
}

template<class string> void ReadAnimationCSV::tokenizeV(const std::string &s,std::vector<string> &o) { 
	typedef boost::tokenizer<boost::escaped_list_separator<char> >  tok_t;
	tok_t tok(s);
	for(tok_t::iterator j (tok.begin());j != tok.end();++j) {
		std::string f(*j);
		boost::trim(f);
		o.push_back(boost::lexical_cast<string>(f));
	}
}

// Functions to extract the variables that have been read
int ReadAnimationCSV::GetTimeSteps() { return TimeSteps; }
int ReadAnimationCSV::GetDegreesOfFreedom() { return DegreesOfFreedom; }
double ReadAnimationCSV::GetStartTime() { return StartTime; }
double ReadAnimationCSV::GetEndTime() { return EndTime; }
vector<VectorNd> ReadAnimationCSV::GetQ() { return Q; }

class Dynamics{
friend class ForcePlatformsExtractor;
	vector<VectorNd> Q; 
	vector<VectorNd> QDot;			
	vector<VectorNd> QDDot;		
	vector<VectorNd> Tau;   	
	double T0;		
	double Tf;
	double delta_t;	
	int NumT;       
 	Vector3d Footstartposition_r;
        Vector3d Footstartposition_l;
	int NumQ;       
	int FootPrintLeftInterval;	
	int FootPrintRightInterval;		
	int InitialLeftStep;		
	int InitialRightStep;	
	int EndRightStep;
	int EndLeftStep;
	unsigned int Id_Left;
	unsigned int Id_Right;                   
	vector<VectorNd> LForcePlate;	  
	vector<VectorNd> RForcePlate;	         
	vector<VectorNd> LInverseForceMomentum;
	vector<VectorNd> RInverseForceMomentum;	
	Model* model;	                     

	public:
		Dynamics(const char *FilenameLuamodel);
		void ReadJointAnglesFromFile(const char *filename);
		void ReadForcesFromFile(const string& filename); 
		void ComputeJointTorques();

		void WriteQValuesToFile(const char *filenameOfOutputFile);
		void WriteQDotValuesToFile(const char *filenameOfOutputFile);
		void WriteQDDotValuesToFile(const char *filenameOfOutputFile);
		void WriteTauValuesToFile(const char *filenameOfOutputFile);

		void PrintLabelsFromc3dFile(const string& filename);   
   		~Dynamics();	

	private:
		void GetDerivativeswithFiniteDifferencing(); 
		void PrintObjectofTypevectorVectorNDtoFile(vector<VectorNd>* PointerToValues, const char *filename);
		
		// Some BTK functions
		btk::Acquisition::Pointer readAcquisition(const string& filename);  
		btk::WrenchCollection::Pointer readbtkForces(const string& filename);  
		std::string GetPointTypeAsString(btk::Point::Type t);		  

		int ReadStepInterval(btk::Acquisition::Pointer Acq); 
		
		vector<SpatialVector>  SetupExternalForces(VectorNd ForcesMomentumsLeft, VectorNd ForcesMomentumsRight);
};

Dynamics::Dynamics(const char *FilenameLuamodel) {
	model = new Model();
	Addons::LuaModelReadFromFile(FilenameLuamodel, model, false);	
	NumQ = model->dof_count;
	Id_Left  = model->GetBodyId("foot_l");
	Id_Right = model->GetBodyId("foot_r");
}

void Dynamics::ReadJointAnglesFromFile(const char *animationFile) { 
	ReadAnimationCSV Animation(animationFile); 
	Animation.Read();
	T0 = Animation.GetStartTime();
	Tf = Animation.GetEndTime();
	NumT = Animation.GetTimeSteps();
	Q = Animation.GetQ();	
	
	if (NumQ !=Animation.GetDegreesOfFreedom()) {   
		cerr << "		*** WARNING: Animation Degrees Of Freedom: " << Animation.GetDegreesOfFreedom()<< endl;
		cerr << "		*** Number of Q in the model: " << NumQ << endl;
	}
	delta_t=(Tf-T0)/(NumT-1);

	// Initialize force vectors	
	for (int i=0;i<NumT;i++) {
		QDot.push_back(VectorNd::Zero(NumQ));
		QDDot.push_back(VectorNd::Zero(NumQ));
		LForcePlate.push_back(VectorNd::Zero(6));
		RForcePlate.push_back(VectorNd::Zero(6));
		LInverseForceMomentum.push_back(VectorNd::Zero(6));
		RInverseForceMomentum.push_back(VectorNd::Zero(6));
	}
	GetDerivativeswithFiniteDifferencing();
}                   

void Dynamics::WriteQValuesToFile(const char *filenameOfOutputFile) {
	PrintObjectofTypevectorVectorNDtoFile(&Q, filenameOfOutputFile);
}
void Dynamics::WriteQDotValuesToFile(const char *filenameOfOutputFile) {
	PrintObjectofTypevectorVectorNDtoFile(&QDot, filenameOfOutputFile);
}
void Dynamics::WriteQDDotValuesToFile(const char *filenameOfOutputFile) {
	PrintObjectofTypevectorVectorNDtoFile(&QDDot, filenameOfOutputFile);
}
void Dynamics::WriteTauValuesToFile(const char *filenameOfOutputFile) {
	PrintObjectofTypevectorVectorNDtoFile(&Tau, filenameOfOutputFile);
}

void Dynamics::ComputeJointTorques() {
	for(int i = 0; i < NumT; i++) {
		Tau.push_back(VectorNd::Zero(Q[0].size()));
		vector<Math::SpatialVector> f = SetupExternalForces(LInverseForceMomentum[i],RInverseForceMomentum[i]);   
		InverseDynamics(*model, Q[i], QDot[i], QDDot[i], Tau[i], &f);
	}
}   

vector<Math::SpatialVector> Dynamics::SetupExternalForces(VectorNd ForcesMomentumsLeft,VectorNd ForcesMomentumsRight) {   
	vector<SpatialVector> fext;
	
	for(int i=0; i<model->mBodies.size(); i++) {
		fext.push_back(SpatialVector::Zero(6));
	}	
	for (int i=0;i<6;i++) {
		fext[Id_Left][i]=ForcesMomentumsLeft[i];                   
		fext[Id_Right][i]=ForcesMomentumsRight[i];
	}
	return fext;
}

Dynamics::~Dynamics(){}

void Dynamics::ReadForcesFromFile(const string& filename) { 

	btk::Acquisition::Pointer acq = readAcquisition(filename);
	ReadStepInterval(acq);

	int Nframes=(acq->GetPoint("SACR"))->GetFrameNumber();
	if (Nframes != NumT)
	   cerr << "		*** WARNING: The number of frames in the c3d file does not match the animation" << endl;

	// Read Ground Reaction Forces if available, else set to Zero
	Eigen::Matrix<double,Eigen::Dynamic,3> LForce;
	Eigen::Matrix<double,Eigen::Dynamic,3> RForce;

	btk::WrenchCollection::Pointer fpws = readbtkForces(filename);
	btk::WrenchCollection::ConstIterator itWrench = fpws->Begin();
	int numberOfForcePlates = fpws->GetItemNumber();
	int numberOfForcePlateFrames = (*itWrench)->GetPosition()->GetFrameNumber();
	int SamplePerFrame = static_cast<double>(acq->GetNumberAnalogSamplePerFrame());
	if (numberOfForcePlateFrames != Nframes*SamplePerFrame)
	  cerr << "		*** WARNING: The number of frames in the c3d file does not match the the resampled forces" << endl;
	
	Eigen::MatrixXd force_array[numberOfForcePlates];
	Eigen::MatrixXd moment_array[numberOfForcePlates];
	Eigen::MatrixXd position_array[numberOfForcePlates];
	Eigen::MatrixXd resized_force_array[numberOfForcePlates];
	Eigen::MatrixXd resized_moment_array[numberOfForcePlates];
	Eigen::MatrixXd resized_position_array[numberOfForcePlates];
	Eigen::MatrixXd resized_globalposition_array[numberOfForcePlates];
	Eigen::MatrixXd Rstep_PlateValues[numberOfForcePlates];
	Eigen::MatrixXd Lstep_PlateValues[numberOfForcePlates];
	int R_plate_id = 0, L_plate_id = 0;
	double Rstep_Norm[numberOfForcePlates];
	double Lstep_Norm[numberOfForcePlates];

	// Initialize all arrays
	for (int i = 0 ; i < numberOfForcePlates ; ++i) {
		force_array[i]  = Eigen::MatrixXd::Zero(numberOfForcePlateFrames,3);
		moment_array[i] = Eigen::MatrixXd::Zero(numberOfForcePlateFrames,3);
		position_array[i]  = Eigen::MatrixXd::Zero(numberOfForcePlateFrames,3);
		resized_force_array[i]  = Eigen::MatrixXd::Zero(Nframes,3);
		resized_moment_array[i] = Eigen::MatrixXd::Zero(Nframes,3);
		resized_position_array[i] = Eigen::MatrixXd::Zero(Nframes,3);
		resized_globalposition_array[i]      = Eigen::MatrixXd::Zero(Nframes,3);
		Rstep_PlateValues[i]   = Eigen::MatrixXd::Zero(FootPrintRightInterval,3);
		Lstep_PlateValues[i]   = Eigen::MatrixXd::Zero(FootPrintLeftInterval,3);
	}
	
	// Read raw data from C3D
	for (int i = 0 ; i < numberOfForcePlates ; ++i) {
		if (itWrench->get() != 0) {
			force_array[i]  = (*itWrench)->GetForce()->GetValues();
			moment_array[i] = (*itWrench)->GetMoment()->GetValues();
			position_array[i]  = (*itWrench)->GetPosition()->GetValues();
		}
		else
		  cerr << "		*** WARNING: Force platform wrench #"+btk::ToString(i)+" is empty." << endl;
		++itWrench;
	}
	
	// Resize data to match motion capture frames
	for (int i = 0 ; i < numberOfForcePlates ; ++i) {
		for (int k = 0; k < 3; ++k) {
			for (int j = 1; j < (Nframes-1); ++j) {
				for (int n = -(SamplePerFrame/2); n < (SamplePerFrame-SamplePerFrame/2); ++n) {
					resized_force_array[i](j,k)    = resized_force_array[i](j,k)    + (force_array[i](SamplePerFrame*j+n,k))/SamplePerFrame;
					resized_moment_array[i](j,k)   = resized_moment_array[i](j,k)   + (moment_array[i](SamplePerFrame*j+n,k))/SamplePerFrame;
					resized_position_array[i](j,k) = resized_position_array[i](j,k) + (position_array[i](SamplePerFrame*j+n,k))/SamplePerFrame;
				}
			}
			for (int n = 0; n < SamplePerFrame/2; ++n) {
				resized_force_array[i](0,k)            = resized_force_array[i](0,k)            + (force_array[i](n,k))/(SamplePerFrame/2);
				resized_moment_array[i](0,k)           = resized_moment_array[i](0,k)           + (moment_array[i](n,k))/(SamplePerFrame/2);
				resized_position_array[i](0,k)         = resized_position_array[i](0,k)         + (position_array[i](n,k))/(SamplePerFrame/2);
				resized_force_array[i](Nframes-1,k)    = resized_force_array[i](Nframes-1,k)    + (force_array[i](Nframes-1-n,k))/(SamplePerFrame/2);
				resized_moment_array[i](Nframes-1,k)   = resized_moment_array[i](Nframes-1,k)   + (moment_array[i](Nframes-1-n,k))/(SamplePerFrame/2);
				resized_position_array[i](Nframes-1,k) = resized_position_array[i](Nframes-1,k) + (position_array[i](Nframes-1-n,k))/(SamplePerFrame/2); 
			}
		}
		
		// Compute global position of the impact from resized force and moment values
		for (int m = 0; m < Nframes; ++m) {
			Vector3d FF = resized_force_array[i].row(m);
			Vector3d MM = resized_moment_array[i].row(m);
			Vector3d OO = resized_position_array[i].row(m);
			resized_globalposition_array[i].row(m) = FF.cross(MM)/FF.squaredNorm() + OO;
			resized_globalposition_array[i](m,2)=0;
		}
	}
	
	// Check if foot fell on multiple plates
	for (int i = 0 ; i < numberOfForcePlates ; ++i) {
		for (int k = 0; k < 3; ++k) {
		  for (int j = 0; j < FootPrintRightInterval; ++j)
		    Rstep_PlateValues[i](j,k) = resized_force_array[i](j+InitialRightStep,k);
		  for (int j = 0; j < FootPrintLeftInterval; ++j)
		    Lstep_PlateValues[i](j,k) = resized_force_array[i](j+InitialLeftStep,k);
		}
		Rstep_Norm[i] = Rstep_PlateValues[i].squaredNorm();
		Lstep_Norm[i] = Lstep_PlateValues[i].squaredNorm();
		
		if (Rstep_Norm[i] > Rstep_Norm[R_plate_id])
		  R_plate_id=i;
		if (Lstep_Norm[i] > Lstep_Norm[L_plate_id])
		  L_plate_id=i;
	}
	if (R_plate_id == L_plate_id)
	  cerr << "		*** WARNING: Force platform wrench #"+ btk::ToString(R_plate_id) +  " identified as left and right at the same time" << endl;

	
	
	// The forces for these files need to be reversed, as the marker positions (in c3d file) were manually switched around
	// This requires a more elegant solution   
	std::vector<string> reverseFilenames;
	reverseFilenames.push_back(filename);
	reverseFilenames.push_back("../Data/c3d_barefeet/3053176.c3d");
	reverseFilenames.push_back("../Data/c3d_barefeet/3053178.c3d");
	reverseFilenames.push_back("../Data/c3d_barefeet/3053183.c3d");
	
	int reverseSign = -1;
	for (int i = 1; i < reverseFilenames.size(); i++) {
	  if (strcmp(reverseFilenames[0].c_str(), reverseFilenames[i].c_str()) == 0)     
	     reverseSign = 1;
	}
	
	// Shift forces and moments to local frame
	for(int i=0; i < Nframes; i++) { 
		Vector3d FFR = resized_force_array[R_plate_id].row(i);
		FFR[0] = reverseSign*FFR[1]; FFR[1] = 0.0;
		Vector3d FFL = resized_force_array[L_plate_id].row(i);
		FFL[0] = reverseSign*FFL[1]; FFL[1] = 0.0;
		Vector3d PPR = resized_globalposition_array[R_plate_id].row(i);
		PPR[0] = reverseSign*PPR[1]; PPR[1] = 0.0;
		Vector3d PPL = resized_globalposition_array[L_plate_id].row(i);
		PPL[0] = reverseSign*PPL[1]; PPL[1] = 0.0;
		Vector3d MMR = PPR.cross(FFR)/1000;
		Vector3d MML = PPL.cross(FFL)/1000;
		for (int j = 0; j < 3; j++) {
			RInverseForceMomentum[i](j) = MMR(j);
			LInverseForceMomentum[i](j) = MML(j);
			RInverseForceMomentum[i](j+3) = FFR(j);
			LInverseForceMomentum[i](j+3) = FFL(j);
		}
	}
}

int Dynamics::ReadStepInterval(btk::Acquisition::Pointer Acq) { 
  
	InitialRightStep = -1;
	InitialLeftStep  = -1;
	EndRightStep = -1;
	EndLeftStep  = -1;
	
	for (int i = 0; i < Acq->GetEventNumber(); i++) {
		const std::string & label = (Acq->GetEvent(i))->GetLabel();
		const std::string & side = (Acq->GetEvent(i))->GetContext();
		int Frame = (Acq->GetEvent(i))->GetFrame();
		if(label == "Foot Strike") {
			if(side == "Right") {
				if((InitialRightStep == -1) or (Frame < InitialRightStep))
				  InitialRightStep = Frame;
			}
			else if(side == "Left") {
				if((InitialLeftStep  == -1) or (Frame < InitialLeftStep))
				  InitialLeftStep = Frame;
			}
		}
	}
	for(int i = 0; i < Acq->GetEventNumber(); i++) {
		const std::string & label = (Acq->GetEvent(i))->GetLabel();
		const std::string & side = (Acq->GetEvent(i))->GetContext();
		int Frame = (Acq->GetEvent(i))->GetFrame();
		if(label == "Foot Off") {
			if(side == "Right") {
				if (((EndRightStep == -1) or (Frame < EndRightStep)) and (Frame > InitialRightStep))
				  EndRightStep = Frame;
			}
			else if(side=="Left") {
				if (((EndLeftStep  == -1) or (Frame < EndLeftStep)) and (Frame > InitialLeftStep))
				  EndLeftStep = Frame;
			}
		}
	}
	
	int FirstFrame         = Acq->GetFirstFrame();
	EndRightStep           = EndRightStep     - FirstFrame;
	EndLeftStep            = EndLeftStep      - FirstFrame;
	InitialRightStep       = InitialRightStep - FirstFrame;
	InitialLeftStep        = InitialLeftStep  - FirstFrame;
	FootPrintRightInterval = EndRightStep     - InitialRightStep;
	FootPrintLeftInterval  = EndLeftStep      - InitialLeftStep;
}

void Dynamics::GetDerivativeswithFiniteDifferencing() {
       	for(int i = 0; i < NumQ; i++) {
		for (int j = 0; j < NumT-1; j++)
			(QDot[j])[i] = ((Q[j+1])[i]-(Q[j])[i])/(delta_t);
		QDot[NumT-1][i] = QDot[NumT-2][i];
	} 
  	for(int i = 0; i < NumQ; i++) {
		for (int j = 0; j < NumT-1; j++) 
		  (QDDot[j])[i] = ((QDot[j+1])[i]-(QDot[j])[i])/(delta_t);
		QDDot[NumT-1][i]=QDDot[NumT-2][i];
	}
}

void Dynamics::PrintObjectofTypevectorVectorNDtoFile(vector<VectorNd>* PointerToValues, const char * filename) {
	ofstream output_file (filename, ios::out);
	vector<VectorNd> Values=*PointerToValues;
	if (!output_file) {cerr << "Error: could not open file " << filename << "." << endl;abort();}
	for (int i = 0; i < Values.size(); i++) {    
		output_file << i << ", ";
		for(int j = 0; j<Values[i].size(); j++){
			output_file << (Values[i])[j];
			if (j != Values[i].size() - 1){output_file << ", ";}}
		output_file << endl;}
	output_file.close();}

// BTK functions 
btk::Acquisition::Pointer Dynamics::readAcquisition(const std::string& filename) {
  	btk::AcquisitionFileReader::Pointer reader = btk::AcquisitionFileReader::New();
	reader->SetFilename(filename);
	reader->Update();
	return reader->GetOutput();
}
	
btk::WrenchCollection::Pointer Dynamics::readbtkForces(const std::string& filename) {
	btk::Acquisition::Pointer reader = readAcquisition(filename);
	btk::ForcePlatformsExtractor::Pointer fpExtractor = btk::ForcePlatformsExtractor::New();
	fpExtractor->SetInput(reader);
	btk::ForcePlatformWrenchFilter::Pointer fpwFilter = btk::ForcePlatformWrenchFilter::New();
	fpwFilter->SetInput(fpExtractor->GetOutput());
	btk::WrenchCollection::Pointer fpwrs = fpwFilter->GetOutput();	
	fpwFilter->SetTransformToGlobalFrame(true);
	fpwrs->Update();
	return fpwrs;
}

std::string Dynamics::GetPointTypeAsString(btk::Point::Type t) {
	if (t == btk::Point::Marker) return "Marker";
	else if (t == btk::Point::Angle) return "Angle";
	else if (t == btk::Point::Force) return "Force";
	else if (t == btk::Point::Moment) return "Moment";
	else if (t == btk::Point::Power) return "Power";
	else if (t == btk::Point::Scalar) return "Scalar";
}

void Dynamics::PrintLabelsFromc3dFile(const string& filename) {  
	btk::Acquisition::Pointer acq = readAcquisition(filename);
	for (btk::Acquisition::PointConstIterator it = acq->BeginPoint() ; it != acq->EndPoint() ; ++it) 
	  std::cout << (*it)->GetLabel() << " (" << (*it)->GetDescription() << "): " << GetPointTypeAsString((*it)->GetType())<<endl;
}

void app_analyze_options(int argc, char** argv, char** filenames) {

	if(argc < 3) {
	  cerr << "There have to be exactly 3 arguments with the formats .lua, .csv, .c3d" << endl;
	  exit(2);
	}    
	
	for(int c=1; c<argc; c++) {
	  if(strstr(argv[c], ".lua")>0)
	    filenames[1] = argv[c];
	  else if(strstr(argv[c], ".csv")>0)
	    filenames[2] = argv[c];
	  else if(strstr(argv[c], ".c3d")>0)
	    filenames[3] = argv[c];
	  else {
	    cerr << "*** WARNING: undefined option: " << argv[c] << endl;
	    exit(2);
	  }
	}
}

int main(int argc, char* argv[]) {
	
	char *filenames[argc];                                                 
	app_analyze_options(argc, argv, filenames);    
	
	Dynamics ID_OBJ(filenames[1]);
	ID_OBJ.ReadJointAnglesFromFile(filenames[2]);
	ID_OBJ.ReadForcesFromFile(filenames[3]);
	ID_OBJ.ComputeJointTorques();
	ID_OBJ.WriteTauValuesToFile("./id_res.txt");
	ID_OBJ.WriteQValuesToFile("./id_q.txt");
	ID_OBJ.WriteQDotValuesToFile("./id_qdot.txt");
	return 0;
}
