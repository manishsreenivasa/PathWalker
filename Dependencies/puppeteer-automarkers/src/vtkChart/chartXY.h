#ifndef CHARTXY_H
#define CHARTXY_H

#include "SimpleMath/SimpleMath.h"
#include "simpleInterpolation/SplineInterpolator.h"

#include <vector>

#include <vtkVersion.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderWindow.h>
#include <vtkSmartPointer.h>
#include <vtkChartXY.h>
#include <vtkAxis.h>
#include <vtkTable.h>
#include <vtkPlot.h>
#include <vtkPlotPoints.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkContextView.h>
#include <vtkContextScene.h>
#include <vtkPen.h>

#include "QVTKWidget.h"

  
struct ChartColor {
  unsigned char red;
  unsigned char green;
  unsigned char blue;
  unsigned char alpha;

ChartColor() : red(0),
    green(0),
    blue(0),
    alpha(0){
	   }
	 ChartColor(const unsigned char _red, const unsigned char _green, const unsigned char _blue, const unsigned char _alpha) : red(_red),
	   green(_green),
	   blue(_blue),
	   alpha(_alpha) {
	   }
	   };
	 class ChartContainer {
	 protected:
	   std::string chartTitle;
	   double timePtr;
	   std::vector<SplineInterpolator<VectorNd> > dataContainerVec;
	   vtkSmartPointer<vtkContextView> viewPtr;
	   vtkSmartPointer<vtkChartXY> chartPtr;
    
   
	 public:
	   ChartContainer(const std::string _chartTitle, const std::string _xAxisTitle, const std::string _yAxisTitle, const bool legend = false);
	   void registerParent(QVTKWidget* _parent);
	   void setXaxisTitle(const std::string _title);
	   void setYaxisTitle(const std::string _title);
	   void setTitle(const std::string _chartTitle);
	   void reset();
	   void pushData(const std::string _title, const VectorNd _Tdata, const VectorNd _Xdata, const double _width, const ChartColor _color, const int _lineType = 1);
	   void setTimePtr(const double _timePtr);
	   void update();

	 };

#endif
