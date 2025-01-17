/*==============================================================================

  Program: 3D Slicer

  Portions (c) Copyright Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

==============================================================================*/

#ifndef __qSlicerHBVDSciVisPModuleWidget_h
#define __qSlicerHBVDSciVisPModuleWidget_h

// Slicer includes
#include "qSlicerAbstractModuleWidget.h"
#include "qSlicerHBVDSciVisPModuleExport.h"
#include "vtkSlicerApplicationLogic.h"
#include <vtkSlicerVolumeRenderingLogic.h>
#include <qSlicerLayoutManager.h>


// MRML includes
#include <vtkMRMLScalarVolumeNode.h>
#include "vtkMRMLCrosshairNode.h"//12.July.2023

//qMRML includes
#include "qMRMLWidget.h"
#include <qMRMLThreeDWidget.h>
#include <qMRMLThreeDView.h>
#include <qMRMLSliceControllerWidget.h>
#include "qMRMLSliceWidget.h"


// vtk includes
#include <vtkNew.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkActor.h>
#include <vtkActor2D.h>
#include <vtkNamedColors.h>
#include <vtkPolyDataMapper.h>
#include <vtkCollisionDetectionFilter.h> //5.July.2023
#include <vtkTransform.h>//5.July.2023
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkTransform.h>//26.July.2023

class qSlicerHBVDSciVisPModuleWidgetPrivate;
class vtkMRMLNode;
class vtkTimerCallbackHaptic;

/// \ingroup Slicer_QtModules_ExtensionTemplate
class Q_SLICER_QTMODULES_HBVDSCIVISP_EXPORT qSlicerHBVDSciVisPModuleWidget :
  public qSlicerAbstractModuleWidget
{
  Q_OBJECT
  QVTK_OBJECT

public:
    
  bool bDRD_HBVDSV_glbl = false;
  typedef qSlicerAbstractModuleWidget Superclass;
  qSlicerHBVDSciVisPModuleWidget(QWidget *parent=0);
  virtual ~qSlicerHBVDSciVisPModuleWidget();
  vtkMRMLScalarVolumeNode* loadVolume(const char* volume, vtkMRMLScene* scene);
  double dHapticWorldSpaceRadius = 0;
  double dXYZPtr[3];
  void getXYZBasedonHapticWorkSpaceRadius(double dHapticWorldSpaceRadiusINm,double* dXYZPtr, char UpDownLeftRightInOutinUDLRIO);

  vtkSmartPointer<vtkMRMLCrosshairNode> CrosshairNode;

  

public slots:
    void onbtnStartHapticLoopCycle();
    void onbtnStopHapticLoopCycle();
    void onGetNodeData();
    void onInitializeHaptics();
    void setMRMLScene(vtkMRMLScene* varMRMLScene = nullptr);
    void onHapticThreadCalled();
    void onbtnLeftMoveHaptic();
    void onbtnRightMoveHaptic();
    void onbtnUpMoveHaptic();
    void onbtnDownMoveHaptic();
    void onbtnInMoveHaptic();
    void onbtnOutMoveHaptic();
    void onbtnForcesOnOffClicked(bool bClickedState);
    void onbtnTestHapticDeviceMovementClicked(bool bClickedState);
    void onbtnTestVTKPrimitiveSurfaceClicked(bool bClickedState);
    void CallHapticThreadandUIUpdateTimer(float xHapt, float yHapt, float zHapt);
    static void onHapticPositionUpdateCursor(vtkObject* caller, unsigned long eid, void* clientData, void* vtkNotUsed(callData));



protected:
  QScopedPointer<qSlicerHBVDSciVisPModuleWidgetPrivate> d_ptr;

  void setup() override;

  vtkNew<vtkNamedColors> colorHapticCursor;
  vtkNew<vtkSphereSource> vtkSphereHaptCursor; //VTK Cursor
  vtkNew<vtkSphereSource> vtkSphereHaptGoalCursor; //VTK Cursor
  vtkNew<vtkPolyDataMapper> vtkSphereHapticMapper;
  vtkNew<vtkPolyDataMapper> vtkSphereHapticGoalMapper;
  vtkNew<vtkActor> vtkActorHaptCursor;// Gives error for 2 --3D Views as this actor gets deleted -
                                      // 3D View Id -- Accesssible Name :   --
                                      // 3D View Id : vtkMRMLViewNode2      -- Generic Warning: In C:\D\S5D\VTK\Common\Core\vtkDebugLeaks.cxx, line 245
                                      // Deleting unknown object : vtkHapticCursorActorClass
                                         

                                      // 3D View Id : vtkMRMLViewNode1
                                      //3D View Actor Coordinate System : World--
                                      // Generic Warning : In C : \D\S5D\VTK\Common\Core\vtkDebugLeaks.cxx, line 245
                                      // Deleting unknown object : vtkHapticCursorActorClass
  vtkNew<vtkActor> vtkActorHaptGoalCursor;

  static double getHapticDevicePositionInDouble(std::string hapticDevicePositionStr3, char cXYZ);

  qMRMLThreeDWidget* threeDWidget = nullptr;
  vtkMRMLViewLogic* threeDWidgetViewLogic = nullptr;
  qMRMLThreeDView* threeDView = nullptr;
  vtkRenderWindow* renderWindow = nullptr;
  vtkRenderer* renderer = nullptr;
  vtkRenderWindowInteractor* interactor = nullptr;
  qMRMLSliceWidget* sliceWidget = nullptr;
  qMRMLSliceControllerWidget* sliceControlWidget = nullptr;
  vtkSlicerApplicationLogic* appLogicLocal = nullptr;
  vtkSlicerVolumeRenderingLogic* appVolRendMod = nullptr;
  vtkMRMLSelectionNode* appSelNode = nullptr;
  qSlicerLayoutManager* layoutManager = nullptr;
  vtkMRMLScalarVolumeNode* scalarNode = nullptr;

  vtkRenderWindowInteractor* interactorSliceRed = nullptr; // 11.Aug.2023


  


private:
  Q_DECLARE_PRIVATE(qSlicerHBVDSciVisPModuleWidget);
  Q_DISABLE_COPY(qSlicerHBVDSciVisPModuleWidget);
};

#endif

