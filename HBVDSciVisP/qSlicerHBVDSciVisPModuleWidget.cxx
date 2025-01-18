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

// Qt includes
#include <QDebug>
#include <QTimer> // DRD To Update Position of Haptic after Thread is called https://discourse.slicer.org/t/place-fiducial-every-n-seconds/11764/2
#include <qlabel.h>//21.July.2023
#include <QPointer>

#include <itksys/SystemTools.hxx>

//chai3d includes
#include <chai3d.h>

// qMRML includes
#include "qMRMLNodeFactory.h"
#include "qMRMLSceneFactoryWidget.h"
#include "qMRMLSliceView.h"
#include "qMRMLUtils.h"
#include "qMRMLLayoutManager.h"

// Slicer includes
#include "qSlicerHBVDSciVisPModuleWidget.h"
#include "ui_qSlicerHBVDSciVisPModuleWidget.h"
#include <qSlicerApplication.h>

// MRML includes
#include <vtkMRMLScalarVolumeDisplayNode.h>
#include "vtkMRMLVolumeArchetypeStorageNode.h"
#include <vtkMRMLVolumeRenderingDisplayNode.h>
#include <vtkMRMLThreeDViewInteractorStyle.h>
#include <vtkMRMLSelectionNode.h>
#include <vtkMRMLScene.h>
#include <vtkMRMLViewNode.h>
#include <vtkMRMLTransformNode.h> //3.July.2023
#include <vtkMRMLSubjectHierarchyNode.h>
#include <vtkMRMLCameraNode.h> //17.July.2023
#include <vtkMRMLApplicationLogic.h>
#include <vtkMRMLColorTableNode.h>
#include <vtkMRMLModelNode.h>//8.Aug.2023
#include <vtkMRMLMarkupsCurveNode.h>//9.Aug.2023

#include <vtkMRMLAbstractSliceViewDisplayableManager.h>//11.Aug.2023
#include <vtkMRMLCrosshairDisplayableManager.h>//11.Aug.2023
#include <vtkMRMLMarkupsFiducialNode.h>//16.August.2023

// VTK includes
#include <vtkCollection.h>
#include <vtkRendererCollection.h> // for renderWindow->GetRenderers()
#include <vtkCallbackCommand.h>
#include <vtkProperty.h>
#include <vtkMatrix4x4.h>
#include <vtkNumberToString.h> //3.July.2023
#include <vtkGeneralTransform.h> //3.July.2023
#include <vtkCodedEntry.h>//13.July.2023
#include <vtkDoubleArray.h>//24.July.2023
#include <vtkCellData.h>//18.July.2023
#include <vtkPointData.h>//14.July.2023
#include <vtkGeneralTransform.h> //12.July.2023
#include <vtkOrientedImageData.h> // To get imageData
#include <vtkGeometryFilter.h>//4.Aug.2023

#include <vtkTriangleFilter.h>//7.Aug.2023
#include <vtkConeSource.h>//8.Aug.2023
#include <vtkOrientedGlyphContourRepresentation.h>//11.Aug.2023
#include <vtkContourWidget.h>//11.Aug.2023
#include <vtkBezierContourLineInterpolator.h>//11.Aug.2023
#include <vtkPolygonalSurfacePointPlacer.h>//11.Aug.2023
#include <vtkPolygonalSurfaceContourLineInterpolator.h>//11.Aug.2023
#include <vtkInformation.h>//15.August.2023
#include <vtkIdTypeArray.h>//15.Aug.2023
#include <vtkImageData.h>
//#include <vtkMRMLTransformableNode.h>//21.Aug.2023
#include <vtkScalarsToColors.h>//22.Aug.2023
#include <vtkMRMLModelDisplayNode.h>//23.Aug.2023
#include <vtkMRMLMarkupsNode.h>//24.Aug.2023
#include <vtkSlicerMarkupsLogic.h>//24.Aug.2023
#include <vtkBinaryLabelmapToClosedSurfaceConversionRule.h>//24.Aug.2023
#include <vtkPiecewiseFunction.h>//8.Sep.2023 -- For Making the volume translucent
#include <vtkMRMLVolumePropertyNode.h>
#include <vtkDiscreteFlyingEdges3D.h>
#include <vtkCellPicker.h>//10.Nov.2023
#include <vtkPicker.h>//10.Nov.2023
#include <vtkPointPlacer.h>//10.Nov.2023
#include <vtkCamera.h>//10.Nov.2023
#include <vtkImplicitPolyDataDistance.h>//13.Nov.2023
#include <vtkImageDataGeometryFilter.h>//13.Nov.2023
//#if Slicer_MAIN_PROJECT_VERSION_MINOR <= 3
#include <vtkFlyingEdges3D.h>//17.Nov.2023
//#else
//#include <vtkSurfaceNets3D.h>//23.Nov.2023
//#endif
#include <vtkPolyDataNormals.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkGlyph3D.h>
#include <vtkMaskPoints.h>
//#include <vtkMaskPoints.h>
//#include <vtkMaskPoints.h>
#include <vtkFloatArray.h>

// CTK includes -- DRD Added
#include <ctkMessageBox.h>

// Generic includes
#include <exception>

// Includes for csv reading
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>

// Includes for timestamps
#include <chrono>
#include <ctime>


//These variable are used for generating debug information.
//----------------------------------------------------------------------------
static const std::string SERIALIZED_GEOMETRY_SEPARATOR = ";";
static const std::string SERIALIZATION_SEPARATOR = "&";
static const std::string SERIALIZATION_SEPARATOR_INNER = "|";

//------------------------------------------------------------------------------
using namespace chai3d;



/// <summary>
/// This is the class storing pointer for actor.  It is helpful in getting the widget reference.
/// </summary>
class vtkHapticHBVDSciVisPEventCallbackCommand : public vtkCallbackCommand
{
public:
	static vtkHapticHBVDSciVisPEventCallbackCommand* New()
	{
		return new vtkHapticHBVDSciVisPEventCallbackCommand;
	}
	/// Segment editor widget observing the event
	QPointer<qSlicerHBVDSciVisPModuleWidget> EditorWidget;
	/// Slice widget or 3D widget
	QPointer<qMRMLWidget> ViewWidget;
	cVector3d cVect3dhapticDevicePositionClassPVar;
	vtkSmartPointer<vtkActor> vtkCursorActor;
};

/// <summary>
/// Structure for connecting various class variables.
///</summary>
//-----------------------------------------------------------------------------
struct HapticEditorEventObservation
{
	vtkSmartPointer<vtkHapticHBVDSciVisPEventCallbackCommand> CallbackCommand;
	vtkWeakPointer<vtkObject> ObservedObject;
	QVector<int> ObservationTags;
};


/// <summary>
/// Global variables required for integration between 3D Slicer  and Chai3d
/// </summary>
int contactMode = 0;
vtkNew<vtkSphereSource> sphere1;
vtkNew<vtkConeSource> objVTKCone;
vtkNew<vtkMatrix4x4> matrix1;
vtkNew<vtkTransform> transform0;
vtkSmartPointer<vtkActor> ptrNRRDActor = vtkSmartPointer<vtkActor>::New();
vtkSmartPointer<vtkPolyData> ptrNRRDPolyData;
vtkSmartPointer<vtkPolyDataMapper> ptrNRRDPolyDataMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
vtkSmartPointer<vtkGeometryFilter> ptrNRRDGeometryFilter;
//vtkSmartPointer<vtkGeometryFilter> oGeomFltr = vtkSmartPointer<vtkGeometryFilter>::New();//13.Nov.2023
vtkSmartPointer<vtkImageDataGeometryFilter> oGeomFltr = vtkSmartPointer<vtkImageDataGeometryFilter>::New();//13.Nov.2023

vtkSmartPointer<vtkImplicitPolyDataDistance> implicitPolyDataDistance = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();//13.Nov.2023
vtkSmartPointer<vtkFloatArray> oFArrNormal = vtkSmartPointer<vtkFloatArray>::New();
vtkSmartPointer<vtkFloatArray> oFArrPDPtNormal = vtkSmartPointer<vtkFloatArray>::New();
vtkSmartPointer<vtkPolyData> oPDVol = vtkSmartPointer<vtkPolyData>::New();//Moved to Global 29.Nov.2023

vtkNew<vtkCollisionDetectionFilter> collide;
vtkNew<vtkPolyDataMapper> mapper1;
vtkNew<vtkActor> actor1;
vtkNew<vtkPolyDataMapper> mapper2;
vtkNew<vtkActor> actor2;
//vtkNew<vtkCellPicker> oCellPicker;


vtkNew<vtkImageData> imageData;
vtkSmartPointer<vtkMatrix4x4> volumeRAStoIJKMatrixFromLoadedVolume = vtkSmartPointer<vtkMatrix4x4>::New();
bool bDRD_DEBUG_MODE_DRDVAR = false; // Prints or Generates information for debug
bool bDRD_Collision_Detect_Polygon = false; // Adds polygon to detect collision
bool bDRD_GiveonlySurfaceForce = false; // Ignores volume voxel scalar values
bool bDRD_Collision_Detect_Cone_Construct = false; // If false takes test volume for polygon construction -- If true uses generic VTK Cone as actor into consideration
bool bDRD_Compute_Normal_Glyph = true; // If false takes test volume for polygon construction -- If true uses generic VTK Cone as actor into consideration
double drangeTmpGlbl[2]; //DRD 24.July.2023
double dScalarValueCurrentGlbl;
double dScalarValuePastGlbl;
double dScalarValuePushGlbl;
double dDistance = 0.0; // Distance from Surface of volume -- 28.Nov.2023
vtkIdType oprevCellId = 0;
std::string sinfoHDModelName = "";
std::string sinfoHDManufacturerName = "";

std::string slblRates = "";

std::string sForceDataAtPoint = "";


vtkSmartPointer<vtkMatrix4x4> rasToxyz = vtkSmartPointer<vtkMatrix4x4>::New();
vtkSmartPointer<vtkMatrix4x4> ptrDblRotationMatrix = vtkSmartPointer<vtkMatrix4x4>::New();
vtkSmartPointer<vtkMatrix4x4> ptrDblRotationMatrixRev = vtkSmartPointer<vtkMatrix4x4>::New();


//vtkNew<vtkMRMLMarkupsCurveNode> curveNode;
vtkSmartPointer<vtkMRMLMarkupsCurveNode> curveNode = vtkSmartPointer<vtkMRMLMarkupsCurveNode>::New();


int numSteps = 2; //was 100
double dx = 1.0 / static_cast<double>(numSteps) * 2.0;

int iForceCounter = 0;
/// <summary>
/// This is generic namespace for using Chai3d global variables
/// </summary>
namespace
{

	// a haptic device handler
	cHapticDeviceHandler* handler;

	// a pointer to the current haptic device
	cGenericHapticDevicePtr hapticDevice;

	// a frequency counter to measure the simulation graphic rate
	cFrequencyCounter freqCounterGraphics;

	// a frequency counter to measure the simulation haptic rate
	cFrequencyCounter freqCounterHaptics;

	// haptic thread
	cThread* hapticsThread;

	cVector3d positionGlobalVar;

	// a flag for using damping (ON/OFF)
	bool useDamping = false;

	// a flag for using force field (ON/OFF)
	bool useForceField = true;

	// a global variable to store the position [m] of the haptic device
	cVector3d hapticDevicePosition;

	// a flag to indicate if the haptic simulation currently running
	bool simulationRunning = false;

	// a flag to indicate if the haptic simulation has terminated
	bool simulationFinished = true;

	std::string dataloadTimestamp;

	void chai3dRunHapticLoop();

	//bool AccuratePick(double x, double y, double* pickPoint, double* dSelectionPt, vtkRenderer* renderer, double pickNormal[3] = nullptr);

	cVector3d desiredPosition; // Making this also global as unable to pass argument to thread method
	cVector3d oSurfaceClosestPosition;// For setting VTK surfaces closest position

	const char* cptrFileName = "example.nrrd";

	bool hdbutton0, hdbutton1, hdbutton2, hdbutton3;

	bool bSendForces;
	bool btestDeviceMovement;
	bool bVTKPrimitiveSurface;
	bool bSlipOnSurface = true;
	cVector3d force;
	std::ofstream forceLog;
	// Controls whether force magnitude and haptic loop frequency are logged
	bool bForceLog = true;

	int iImageExtArr6[6]; // Avoid VTK error showing up in the loop  when function vtkIdType vtkImageData::GetScalarIndex(int coordinate[3]) is called from file -- C:\D\S5D\VTK\Common\DataModel\vtkImageData.cxx



	/// <summary>
	/// This is the haptic function.  This will be called in a seperate thread.
	/// </summary>
	void chai3dRunHapticLoop()
	{

		// simulation in now running
		simulationRunning = true;
		simulationFinished = false;

		while (simulationRunning)
		{

			/////////////////////////////////////////////////////////////////////
			// READ HAPTIC DEVICE
			/////////////////////////////////////////////////////////////////////

			// read position 
			cVector3d position;

			cVector3d hdPrevPos;

			cVector3d hdPushPos;
			const double dReturnPosition = 0.33333333; // Delta position to which the push to take place

			hapticDevice->getPosition(position);

			hapticDevicePosition = position;
			//if (dScalarValuePastGlbl < dScalarValueCurrentGlbl && dScalarValuePastGlbl < 7)
			if (dScalarValuePastGlbl < dScalarValueCurrentGlbl && dScalarValuePastGlbl < -1) // Modify this line based on noise present in the data
			{
				hdPushPos = hdPrevPos;
			}
			else if (dScalarValuePastGlbl == 0)
			{
				hdPushPos = position;
			}
			else
			{
				hdPushPos = cSub(hdPrevPos, hapticDevicePosition);
				hdPushPos = cMul(dReturnPosition, hdPushPos);
			}
			// read orientation 
			cMatrix3d rotation;
			hapticDevice->getRotation(rotation);

			// read gripper position
			double gripperAngle;
			hapticDevice->getGripperAngleRad(gripperAngle);

			// read linear velocity 
			cVector3d linearVelocity;
			hapticDevice->getLinearVelocity(linearVelocity);

			// read angular velocity
			cVector3d angularVelocity;
			hapticDevice->getAngularVelocity(angularVelocity);

			// read gripper angular velocity
			double gripperAngularVelocity;
			hapticDevice->getGripperAngularVelocity(gripperAngularVelocity);

			// read user-switch status (button 0)
			bool button0, button1, button2, button3;
			button0 = false;
			button1 = false;
			button2 = false;
			button3 = false;

			hapticDevice->getUserSwitch(0, button0);
			hapticDevice->getUserSwitch(1, button1);
			hapticDevice->getUserSwitch(2, button2);
			hapticDevice->getUserSwitch(3, button3);
			if (button0 == true)
			{
				hdbutton0 = true;
			}
			else
			{
				hdbutton0 = false;
			}
			if (button1 == true)
			{
				hdbutton1 = true;
			}
			else
			{
				hdbutton1 = false;
			}
			if (button2 == true)
			{
				hdbutton2 = true;
			}
			else
			{
				hdbutton2 = false;
			}
			if (button3 == true)
			{
				hdbutton3 = true;
			}
			else
			{
				hdbutton3 = false;
			}
			/////////////////////////////////////////////////////////////////////
			// COMPUTE AND APPLY FORCES
			/////////////////////////////////////////////////////////////////////

			// desired orientation
			cMatrix3d desiredRotation;
			desiredRotation.identity();


			// variables for forces
			force = cVector3d(0, 0, 0);
			cVector3d torque(0, 0, 0);
			double gripperForce = 0.0;
			cVector3d forceField(0, 0, 0);
			// apply force field
			if (useForceField)
			{
				// compute linear force
				double Kp = 25; // [N/m]
				if (bDRD_GiveonlySurfaceForce && !btestDeviceMovement && !bVTKPrimitiveSurface)
				{
					Kp = 25; // [N/m]
					forceField = Kp * (hdPushPos - position);
				}
				else if (btestDeviceMovement && !bVTKPrimitiveSurface)
				{
					forceField = Kp * (desiredPosition - position);
				}
				else if (bVTKPrimitiveSurface && !btestDeviceMovement)
				{
					if (dDistance < 0.0)
					{
						double dPosDistance = 0.0;
						dPosDistance = (dDistance * -1); // As distance inside volume is negative
						if (dPosDistance < 10)
						{
							forceField = ((Kp / 10) * dPosDistance) * (oSurfaceClosestPosition - position); // Relative force till 10.0 inside volume
						}
						else
						{
							forceField = ((Kp / 10) * 10) * (oSurfaceClosestPosition - position); // Full force if its more than 10.0 inside volume
						}

					}


				}
				else if (!btestDeviceMovement && !bVTKPrimitiveSurface)
				{
					Kp = (Kp * 3) / drangeTmpGlbl[1]; // DRD 24.July.2023
					//cVector3d forceField = Kp * (desiredPosition - position);
					Kp = Kp * dScalarValueCurrentGlbl;
					if (Kp > 25) { Kp = 25; }
					forceField = Kp * (hdPushPos - position);
					//cVector3d forceField = Kp * (desiredPosition - position);
				}
				force.add(forceField);
				// compute angular torque
				double Kr = 0.05; // [N/m.rad]
				cVector3d axis;
				double angle;
				cMatrix3d deltaRotation = cTranspose(rotation) * desiredRotation;
				deltaRotation.toAxisAngle(axis, angle);
				torque = rotation * ((Kr * angle) * axis);
			}
			if (!btestDeviceMovement && !bVTKPrimitiveSurface)
			{
				double Kv = 1.0 * 20.0;
				if ((sinfoHDModelName == "Touch" && sinfoHDManufacturerName == "Sensable Technologies"))
				{
					Kv = 1.0 * 5.0;
				}
				else
				{
					Kv = 1.0 * 20.0;
				}
				cVector3d forceDamping = -Kv * linearVelocity;
				force.add(forceDamping);
			}
			// apply damping term

			double dForceLength_Magnitude = force.length();
			sForceDataAtPoint = "Force :[" + std::to_string(dForceLength_Magnitude) + "] - ForceField :[" + std::to_string(forceField.length()) + "] - Surface Closest Pos :[" + oSurfaceClosestPosition.str(3) + "]";

			if (useDamping)
			{
				cHapticDeviceInfo info = hapticDevice->getSpecifications();  // Put this as global -- to avoid repeated call -- Generally not used. When using keep in mind and do appropriately.

				// compute linear damping force
				double Kv = 1.0 * info.m_maxLinearDamping;
				cVector3d forceDamping = -Kv * linearVelocity;
				force.add(forceDamping);

				// compute angular damping force
				double Kvr = 1.0 * info.m_maxAngularDamping;
				cVector3d torqueDamping = -Kvr * angularVelocity;
				torque.add(torqueDamping);

				// compute gripper angular damping force
				double Kvg = 1.0 * info.m_maxGripperAngularDamping;
				gripperForce = gripperForce - Kvg * gripperAngularVelocity;
			}

			if (bSendForces)
			{
				hapticDevice->setForceAndTorqueAndGripperForce(force, torque, gripperForce);
			}
			hdPrevPos = hapticDevicePosition;
			freqCounterHaptics.signal(1);
		}

		// exit haptics thread
		simulationFinished = true;
	}


	struct CSVData {
		std::vector<std::unordered_map<std::string, std::string>> data;
		std::vector<std::string> columnNames;
	};

	CSVData loadCSV(const std::string& filename) {
		CSVData csvData;
		std::ifstream file(filename);
		std::string line, cell;

		// Read the first line to get the column names
		if (std::getline(file, line)) {
			std::stringstream lineStream(line);
			while (std::getline(lineStream, cell, ',')) {
				csvData.columnNames.push_back(cell);
			}
		}

		while (std::getline(file, line)) {
			std::unordered_map<std::string, std::string> row;
			std::stringstream lineStream(line);
			int columnIndex = 0;
			while (std::getline(lineStream, cell, ',')) {
				if (columnIndex < csvData.columnNames.size()) {
					// Handle quoted values
					if (!cell.empty() && cell.front() == '"') {
						std::string quotedCell = cell;
						while (std::getline(lineStream, cell, ',')) {
							quotedCell += "," + cell;
							if (!cell.empty() && cell.back() == '"') {
								break;
							}
						}
						cell = quotedCell;
					}
					// Remove surrounding quotes if present
					if (!cell.empty() && cell.front() == '"' && cell.back() == '"') {
						cell = cell.substr(1, cell.size() - 2);
					}
					row[csvData.columnNames[columnIndex]] = cell;
				}
				columnIndex++;
			}
			csvData.data.push_back(row);
		}

		return csvData;
	}

	// Function to retrieve a specific record
	std::unordered_map<std::string, std::string> retrieveRecord(const CSVData& csvData, int uiParticipant, int uiTrial) {
		for (const auto& row : csvData.data) {
			int participant = std::stoi(row.at("Participant"));
			int trial = std::stoi(row.at("Trial"));
			if (participant == uiParticipant && trial == uiTrial) {
				return row; // Return the matched row's variables
			}
		}

		// If no match is found, return an empty map
		return {};
	}

	void saveToCSV(int participant, int trial, const std::string& condition, const std::string& sample,
		const std::string& markerPosition, const std::string& fiducialPosition,
		const std::string& dataloadTimestamp, const std::string& placementTimestamp) {
		std::string filename = std::to_string(participant) + ".csv";
		bool fileExists = std::ifstream(filename).good();

		std::ofstream file(filename, std::ios_base::app); // Open file for appending

		if (!fileExists) {
			// Write headers if the file doesn't exist
			file << "Participant,Trial,Condition,Sample,MarkerPosition,FiducialPosition,DataloadTimestamp,PlacementTimestamp\n";
		}

		file << participant << "," << trial << "," << condition << "," << sample << ","
			<< "\"" << markerPosition << "\"," << "\"" << fiducialPosition << "\","
			<< dataloadTimestamp << "," << placementTimestamp << "\n";

		file.close();
	}
}



////-------------------------------------------------------------------------
//bool AccuratePickHBVDSciVis(double x, double y, double* pickPoint, double* dSelectionPt, vtkRenderer* renderer)
//{
//
//	double dSelectionPtXY[2] = { x,y };
//	double dWorldPtPt[3] = { 0.0 };
//	double dWorldOrPt[9] = { 0.0 };
//	vtkSmartPointer<vtkPointPlacer> oPointPlacer = vtkSmartPointer<vtkPointPlacer>::New();
//
//	oPointPlacer->ComputeWorldPosition(renderer, dSelectionPtXY , dWorldPtPt, dWorldOrPt);
//	//vtkSmartPointer<vtkCellPicker> oCellPicker = vtkSmartPointer<vtkCellPicker>::New();
//	////bool success = oCellPicker->Pick3DPoint(dSelectionPt, renderer);
//	//bool success = oCellPicker->Pick(x, y, 0, renderer);
//	//if (!success)
//	//{
//	//	return false;
//	//}
//	//vtkSmartPointer<vtkPicker> oPicker = vtkSmartPointer<vtkPicker>::New();
//	//vtkPoints* pickPositions = oPicker->GetPickedPositions();
//	//vtkIdType numberOfPickedPositions = pickPositions->GetNumberOfPoints();
//	//if (numberOfPickedPositions < 1)
//	//{
//	//	return false;
//	//}
//	//pickPositions->GetPoint(0, pickPoint);
//
//
//
//
//
//	//bool success = this->AccuratePicker->Pick(x, y, 0, this->Renderer);
//	//if (pickNormal)
//	//{
//	//	this->AccuratePicker->GetPickNormal(pickNormal);
//	//}
//	//if (!success)
//	//{
//	//	return false;
//	//}
//
//	//vtkPoints* pickPositions = this->AccuratePicker->GetPickedPositions();
//	//vtkIdType numberOfPickedPositions = pickPositions->GetNumberOfPoints();
//	//if (numberOfPickedPositions < 1)
//	//{
//	//	return false;
//	//}
//
//	//// There may be multiple picked positions, choose the one closest to the camera
//	double cameraPosition[3] = { 0,0,0 };
//	renderer->GetActiveCamera()->GetPosition(cameraPosition);
//	//pickPositions->GetPoint(0, pickPoint);
//	double minDist2 = vtkMath::Distance2BetweenPoints(dWorldPtPt, cameraPosition);
//	if (minDist2 < 1)
//	{
//		std::cout << "Mindistance is less than one :" << minDist2 << endl;
//	}
//	else
//	{
//		std::cout << "Mindistance :" << minDist2 << endl;
//	}
//	//for (vtkIdType i = 1; i < numberOfPickedPositions; i++)
//	//{
//	//	double currentMinDist2 = vtkMath::Distance2BetweenPoints(pickPositions->GetPoint(i), cameraPosition);
//	//	if (currentMinDist2 < minDist2)
//	//	{
//	//		pickPositions->GetPoint(i, pickPoint);
//	//		minDist2 = currentMinDist2;
//	//	}
//	//}
//	return true;
//}




/// <summary>
/// This is the main haptic loop. Also for any button functionality add the code here.  For modification to coordinates, make them here.  This is can be termed as the graphic loop.
/// </summary>
class vtkTimerCallbackHaptic : public vtkCommand
{
public:
	vtkTimerCallbackHaptic() = default;

	static vtkTimerCallbackHaptic* New()
	{
		vtkTimerCallbackHaptic* cb = new vtkTimerCallbackHaptic;
		cb->TimerCount = 0;
		return cb;
	}

	virtual void Execute(vtkObject* caller, unsigned long eventId,
		void* vtkNotUsed(callData))
	{
		double dLoadedVolumeDiagonal = 1.0; //400.0
		if (vtkCommand::TimerEvent == eventId)
		{
			//	++this->TimerCount;
		}
		if (bDRD_DEBUG_MODE_DRDVAR)
		{
			std::cout << this->TimerCount << std::endl;
		}
		std::string dArrXYZ = hapticDevicePosition.str(3);

		std::string sHapticPos = dArrXYZ;
		std::size_t pos = sHapticPos.find(",");      // position of "live" in str
		std::string sHapticPosX = sHapticPos.substr(0, pos);

		std::string sHapticPosXCut = sHapticPos.substr(pos + 1);
		std::size_t posY = sHapticPosXCut.find(",");

		std::string sHapticPosY = sHapticPos.substr(pos + 1, posY);

		std::string sHapticPosYCut = sHapticPosXCut.substr(posY + 1);
		std::string sHapticPosZ = sHapticPosYCut;

		vtkImageData* ptrScalarNodeImgData = scalarNodeTmp->GetImageData();
		double ijk3Func[4] = { 0.0, 0.0, 0.0, 1.0 };
		if (bDRD_DEBUG_MODE_DRDVAR)
		{
			std::cout << " vtkImageData : [ " << ptrScalarNodeImgData << " ] " << std::endl;
			std::cout << " vtkImageData ptrScalarNodeImgData->GetLength() -- Length of Diagonal of Bounding box: [ " << ptrScalarNodeImgData->GetLength() << " ] " << std::endl;
		}
		try
		{
			dLoadedVolumeDiagonal = ptrScalarNodeImgData->GetLength();
		}
		catch (...)
		{
			std::cout << " Unable to get Length : [" << dLoadedVolumeDiagonal << " ]" << endl;
		}
		//std::cout << " dLoadedVolumeDiagonal : [" << dLoadedVolumeDiagonal << " ]" << endl;
		//dLoadedVolumeDiagonal = (dLoadedVolumeDiagonal * 10.0);
		//dLoadedVolumeDiagonal = (dLoadedVolumeDiagonal * 8.0 * 2.0);// 28.Nov.2023 Had multi factor of 8.0 now doubling it
		dLoadedVolumeDiagonal = (dLoadedVolumeDiagonal * 8.0);// Using Touch Device, so multi factor need not be doubled.
		if (bDRD_DEBUG_MODE_DRDVAR)
		{
			std::cout << " Multiplication factor : dLoadedVolumeDiagonal -- Length of Diagonal of Bounding box: [ " << dLoadedVolumeDiagonal << " ] " << std::endl;
		}
		double dHapticPosCallBackValueX = std::stod(sHapticPosX);
		double dHapticPosCallBackValueY = std::stod(sHapticPosY);
		double dHapticPosCallBackValueZ = std::stod(sHapticPosZ);

		double dCurX = (dHapticPosCallBackValueX * dLoadedVolumeDiagonal);//383.0 //4000.0
		double dCurY = (dHapticPosCallBackValueY * dLoadedVolumeDiagonal);
		double dCurZ = (dHapticPosCallBackValueZ * dLoadedVolumeDiagonal);
		double dHD_XYZ3[4] = { 0.0, 0.0, 0.0, 0.0 }; //Quaternion used in 3D Rotation
		double dHD_XYZRev[4] = { 0.0, 0.0, 0.0, 0.0 }; //Quaternion used in 3D Rotation

		if (bDRD_DEBUG_MODE_DRDVAR)
		{
			std::cout << "vtkTimerCallbackHaptic - setMRMLScene : ["
				<< ptrTmpMRMLScene
				<< "]"
				<< std::endl;
		}

		//double* dCrossHairRAS;

		freqCounterGraphics.signal(1);



		CrosshairNodeHD = vtkMRMLCrosshairNode::SafeDownCast(ptrTmpMRMLScene->GetNthNodeByClass(0, "vtkMRMLCrosshairNode"));

		if (CrosshairNodeHD != nullptr)
		{

			scalarNodeTmp->GetRASToIJKMatrix(vtk4x4VariableCrosshairTmp);
			scalarNodeTmp->GetParentTransformNode();
			if (bDRD_DEBUG_MODE_DRDVAR)
			{
				vtk4x4VariableCrosshairTmp->Print(std::cout);
				double dRASBoundsOfImageVol[6];
				scalarNodeTmp->GetRASBounds(dRASBoundsOfImageVol);
				std::cout << " RAS Bounds of Image Loaded : [xMin " << dRASBoundsOfImageVol[0] << ",  xMax " << dRASBoundsOfImageVol[1] << ",  yMin " << dRASBoundsOfImageVol[2] << ",  yMax " << dRASBoundsOfImageVol[3] << ",  zMin " << dRASBoundsOfImageVol[4] << ",  zMax " << dRASBoundsOfImageVol[5] << " ] " << std::endl;
				std::cout << " RAS Bounds of Image Loaded Diagonals: [xMax - xMin   (" << dRASBoundsOfImageVol[1] - dRASBoundsOfImageVol[0] << ")" << " , yMax - yMin (" << dRASBoundsOfImageVol[3] - dRASBoundsOfImageVol[2] << ") ,   " << "zMax -  zMin (" << dRASBoundsOfImageVol[5] - dRASBoundsOfImageVol[4] << ") ] " << std::endl;
			}
			// Get voxel array
			vtkImageData* imageData = scalarNodeTmp ? scalarNodeTmp->GetImageData() : nullptr;
			vtkPointData* pointData = imageData ? imageData->GetPointData() : nullptr;

			int numberOfCells = imageData->GetNumberOfCells();
			if (bDRD_DEBUG_MODE_DRDVAR)
			{
				std::cout << " numberOfCells : [ " << numberOfCells << " ] " << endl;
			}
			int icellDimsOfImageData[3];
			imageData->GetCellDims(icellDimsOfImageData);
			if (bDRD_DEBUG_MODE_DRDVAR)
			{
				std::cout << " imageData->GetCellDims(icellDimsOfImageData : i,j,k [ " << icellDimsOfImageData[0] << "," << icellDimsOfImageData[1] << "," << icellDimsOfImageData[2] << " ] " << endl;
			}
			vtkDataSetAttributes* ptrCellDataofImage = nullptr;
			ptrCellDataofImage = imageData->GetCellData();
			if (bDRD_DEBUG_MODE_DRDVAR)
			{
				std::cout << "  imageData->GetCellData() [ -------------------------------------------" << endl;
				ptrCellDataofImage->Print(std::cout);
				std::cout << "   ] -------------------------------------------" << endl;
			}

			if (bDRD_DEBUG_MODE_DRDVAR)
			{
				std::cout << " scalarNodeTmp->GetMaxSpacing() : [ " << scalarNodeTmp->GetMaxSpacing() << " ] " << std::endl;
				std::cout << " scalarNodeTmp->GetMinSpacing() : [ " << scalarNodeTmp->GetMinSpacing() << " ] " << std::endl;
			}
			vtkStringArray* ptrVTKStrArr;
			vtkSmartPointer<vtkCodedEntry> cVarVtkCodedEntryTmp;
			cVarVtkCodedEntryTmp = scalarNodeTmp->GetVoxelValueQuantity();
			if (bDRD_DEBUG_MODE_DRDVAR)
			{
				std::cout << " scalarNodeTmp->GetVoxelValueQuantity()->GetAsString().c_str() : [ " << cVarVtkCodedEntryTmp << " ] " << std::endl;
			}
			cVarVtkCodedEntryTmp = scalarNodeTmp->GetVoxelValueUnits();
			if (bDRD_DEBUG_MODE_DRDVAR)
			{
				std::cout << " scalarNodeTmp->GetVoxelValueUnits() : [ " << cVarVtkCodedEntryTmp << " ] " << std::endl;
				std::cout << " scalarNodeTmp->GetVoxelVectorType() : [ " << scalarNodeTmp->GetVoxelVectorType() << " ] " << std::endl;
			}
			const char* ccPtrTmp;
			ccPtrTmp = scalarNodeTmp->GetVoxelVectorTypeAsString(0);
			if (bDRD_DEBUG_MODE_DRDVAR)
			{
				std::cout << " scalarNodeTmp->GetVoxelVectorTypeAsString() : [ " << ccPtrTmp << " ] " << std::endl;
			}
			scalarNodeTmp->GetRASToIJKMatrix(volumeRAStoIJKMatrixFromLoadedVolume);

			double ras34[4] = { dCurX, dCurY, dCurZ, 1.0 };

			ptrDblRotationMatrix->MultiplyPoint(ras34, dHD_XYZ3);

			double ijk3[4] = { 0.0, 0.0, 0.0, 1.0 };
			if (bDRD_DEBUG_MODE_DRDVAR)
			{
				std::cout << " RAS to IJK  matrix: " << std::endl;
			}
			volumeRAStoIJKMatrixFromLoadedVolume->MultiplyPoint(ras34, ijk3); // Required for getting IJK coordinates from RAS coordinates of the Haptic Cursor
			ijk3Func[0] = ijk3[0];
			ijk3Func[1] = ijk3[1];
			ijk3Func[2] = ijk3[2];
			ijk3Func[3] = ijk3[3];
			if (bDRD_DEBUG_MODE_DRDVAR)
			{
				std::cout << "RAS: [" << ras34[0] << ", " << ras34[1] << ", " << ras34[2] << "] -> IJK:[" << ijk3[0] << ", " << ijk3[1] << ", " << ijk3[2] << "]" << std::endl;
			}
			actor->SetPosition(dHD_XYZ3[0], dHD_XYZ3[1], dHD_XYZ3[2]);
			if (bDRD_DEBUG_MODE_DRDVAR)
			{
				std::cout << "RAS: [" << ras34[0] << ", " << ras34[1] << ", " << ras34[2] << "] -> dHD_XYZ3:[" << dHD_XYZ3[0] << ", " << dHD_XYZ3[1] << ", " << dHD_XYZ3[2] << "]" << std::endl;
				std::cout << "Multiply Matrix Point 2 -------------- End" << std::endl;
				std::cout << "scalarNode->GetImageData()->GetScalarTypeMax() : [ " << scalarNodeTmp->GetImageData()->GetScalarTypeMax() << " ] " << std::endl;
			}
			//double dScalarValue = 0.0;


		//	try
		//	{
				////const int* extent = this->Extent;
				//for (int idx = 0; idx < 3; ++idx)
				//{
				//	if (ijk3[idx] < iImageExtArr6[idx * 2] || ijk3[idx] > iImageExtArr6[idx * 2 + 1])
				//	{
				//		std::cout << " Did this remove the VTK Out of memory error?" << endl;
				//	}
				//	else {
					//	dScalarValue = scalarNodeTmp->GetImageData()->GetScalarComponentAsDouble(ijk3[0], ijk3[1], ijk3[2], 0);
				//	}
				//}
			//}
			//catch (...)
			//{
			//	std::cout << " Unable to get scalarNodeTmp->GetImageData()->GetScalarComponentAsDouble : [" << dScalarValue << " ]" << endl;
			//}

			double dScalarValue = scalarNodeTmp->GetImageData()->GetScalarComponentAsDouble(ijk3[0], ijk3[1], ijk3[2], 0);
			if (bDRD_DEBUG_MODE_DRDVAR)
			{
				std::cout << "scalarNode->GetImageData()->GetScalarComponentAsDouble() : [ " << dScalarValue << " ] " << std::endl;
			}
			dScalarValueCurrentGlbl = dScalarValue < 0 ? ((-1) * dScalarValue) : dScalarValue; // We are setting scalar values all positive so as to give positive force always
			if (dScalarValuePastGlbl != dScalarValueCurrentGlbl)
			{
				dScalarValuePastGlbl = dScalarValueCurrentGlbl;
			}
			std::string sDoubleScalarValue = std::to_string(dScalarValueCurrentGlbl); // Made this change on 15.August.2023
			if (bDRD_DEBUG_MODE_DRDVAR)
			{
				std::cout << "scalarNodeTmp->GetImageData()->GetScalarTypeMin() : [ " << scalarNodeTmp->GetImageData()->GetScalarTypeMin() << " ] " << std::endl;
			}
			double drangeTmp[2];
			scalarNodeTmp->GetImageData()->GetScalarRange(drangeTmp);
			ptrQLabelScalarValue->setText("Current Haptic Device Position scalar value : ijk [ " + QString::fromStdString(std::to_string(ijk3[0])) + " ," + QString::fromStdString(std::to_string(ijk3[1])) + "," + QString::fromStdString(std::to_string(ijk3[2])) + " : " + QString::fromStdString(sDoubleScalarValue) + "]");



			vtkSmartPointer <vtkInformation>  ptrVTKTmpInfoDRD = nullptr;
			vtkImageData* imageDataTmpDRD = scalarNodeTmp->GetImageData();

		}
		else
		{
			std::cout << "this->CrosshairNode.GetPointer() null:[" << CrosshairNodeHD << "]" << std::endl;
		}

		if (bVTKPrimitiveSurface)
		{
			double dClosestPoint[3] = { 0.0 };
			double dClosestPointijk[4] = { 0.0 };
			int ijk3FnLocal[3] = { 0 };

			//dHD_XYZ3[0], dHD_XYZ3[1], dHD_XYZ3[2],
			//AccuratePickHBVDSciVis(dHD_XYZ3[0], dHD_XYZ3[1], dpickPoint,dHD_XYZ3, renderer);

			//oFArrNormal->
			//oFArrPDPtNormal->get
			//oPDVol->GetCell(1, 1, 1)->GetPoints()

			dDistance = oIPDD->EvaluateFunctionAndGetClosestPoint(dHD_XYZ3, dClosestPoint);

			//dDistance = oIPDD->EvaluateFunction(dHD_XYZ3);

			std::cout << "Distance to Surface : [ " << dDistance << " ]" << " Haptic Point Loc : [" << dHD_XYZ3[0] << "," << dHD_XYZ3[1] << "," << dHD_XYZ3[2] << "] - Closest Point : [" << dClosestPoint[0] << "," << dClosestPoint[1] << "," << dClosestPoint[2] << "]" << endl;
			//ptrQLabelVTKSurfaceClosePoint->setText("Current ClosestPt Value :  [ " + QString::fromStdString(std::to_string(dClosestPoint[0])) + " ," + QString::fromStdString(std::to_string(dClosestPoint[1])) + "," + QString::fromStdString(std::to_string(dClosestPoint[2])) + " : " + QString::fromStdString(std::to_string(dDistance)) + "]");

			if (dDistance < 0.0)
			{
				double dClosestSurfPos[4] = { dClosestPoint[0], dClosestPoint[1], dClosestPoint[2], 1.0 };

				ptrDblRotationMatrixRev->MultiplyPoint(dClosestSurfPos, dHD_XYZRev);

				double dHaptClosestPosX = (dHD_XYZRev[0] / dLoadedVolumeDiagonal);
				double dHaptClosestPosY = (dHD_XYZRev[1] / dLoadedVolumeDiagonal);
				double dHaptClosestPosZ = (dHD_XYZRev[2] / dLoadedVolumeDiagonal);
				if (iForceCounter < 250)//500
				{
					iForceCounter++; //For increasing force fealt by the user.
				}
				else if (iForceCounter == 250 && bSlipOnSurface)
				{
					oIPDD->EvaluateFunctionAndGetClosestPoint(dHD_XYZ3, dClosestPoint);
					double dClosestSurfPos[4] = { dClosestPoint[0], dClosestPoint[1], dClosestPoint[2], 1.0 };

					ptrDblRotationMatrixRev->MultiplyPoint(dClosestSurfPos, dHD_XYZRev);

					double dHaptClosestPosX = (dHD_XYZRev[0] / dLoadedVolumeDiagonal);
					double dHaptClosestPosY = (dHD_XYZRev[1] / dLoadedVolumeDiagonal);
					double dHaptClosestPosZ = (dHD_XYZRev[2] / dLoadedVolumeDiagonal);
					oSurfaceClosestPosition.set(dHaptClosestPosX, dHaptClosestPosY, dHaptClosestPosZ);
					actorGoalCursor->SetPosition(dClosestPoint[0], dClosestPoint[1], dClosestPoint[2]);
					iForceCounter = 0;
				}
				else {
					oSurfaceClosestPosition.set(dHaptClosestPosX, dHaptClosestPosY, dHaptClosestPosZ);
					actorGoalCursor->SetPosition(dClosestPoint[0], dClosestPoint[1], dClosestPoint[2]);
					iForceCounter = 0;
				}
				ptrQLabelGoalPosition->setText("Current Goal Pos :  [ " + QString::fromStdString(std::to_string(dClosestPoint[0])) + " , " + QString::fromStdString(std::to_string(dClosestPoint[1])) + "," + QString::fromStdString(std::to_string(dClosestPoint[2])) + "]");
				ptrQLabelVTKSurfaceClosePoint->setText("Current ClosestPt Value :  [ " + QString::fromStdString(std::to_string(dHaptClosestPosX)) + " ," + QString::fromStdString(std::to_string(dHaptClosestPosY)) + "," + QString::fromStdString(std::to_string(dHaptClosestPosZ)) + " : " + QString::fromStdString(std::to_string(dDistance)) + "]" + QString::fromStdString(sForceDataAtPoint) + QString::fromStdString("iForceCounter : ") + QString::fromStdString(std::to_string(iForceCounter)));
				volumeRAStoIJKMatrixFromLoadedVolume->MultiplyPoint(dClosestSurfPos, dClosestPointijk);
				ijk3FnLocal[0] = static_cast<int>(dClosestPointijk[0]);
				ijk3FnLocal[1] = static_cast<int>(dClosestPointijk[1]);
				ijk3FnLocal[2] = static_cast<int>(dClosestPointijk[2]);
				std::cout << "XYZ Loc" << ijk3FnLocal[0] << "," << ijk3FnLocal[1] << " , " << ijk3FnLocal[2] << endl;
				try {


					vtkIdType oCellId = ptrScalarNodeImgData->ComputeCellId(ijk3FnLocal);
					std::cout << "Cell Id " << oCellId << endl;
					//ptrQLabelCellIdValue
					ptrQLabelCellIdValue->setText("Cell Id :  [ " + QString::fromStdString(std::to_string(oCellId)) + "] ");

					if (oCellId != oprevCellId && oCellId < 20000)
					{
						oprevCellId = oCellId;
						vtkNew<vtkNamedColors> colorFid1Cursor;
						vtkNew<vtkSphereSource> vtkSphereHaptFiducbt3Tmp; //VTK Cursor
						vtkNew<vtkPolyDataMapper> vtkSphereFid1Mapper;
						vtkNew<vtkActor> vtkActorFidActbt3;

						vtkSphereHaptFiducbt3Tmp->SetRadius(0.9);
						vtkSphereFid1Mapper->SetInputConnection(vtkSphereHaptFiducbt3Tmp->GetOutputPort());
						vtkActorFidActbt3->SetMapper(vtkSphereFid1Mapper);
						vtkActorFidActbt3->GetProperty()->SetColor(colorFid1Cursor->GetColor3d("Red").GetData());
						vtkActorFidActbt3->SetPosition(dHD_XYZ3[0], dHD_XYZ3[1], dHD_XYZ3[2]);

						renderer->AddViewProp(vtkActorFidActbt3);
						renderer->AddActor(vtkActorFidActbt3);
					}
					else if (oCellId != oprevCellId && oCellId >= 20000) {
						oprevCellId = oCellId;
						vtkNew<vtkNamedColors> colorFid1Cursor;
						vtkNew<vtkSphereSource> vtkSphereHaptFiducbt3Tmp; //VTK Cursor
						vtkNew<vtkPolyDataMapper> vtkSphereFid1Mapper;
						vtkNew<vtkActor> vtkActorFidActbt3;

						vtkSphereHaptFiducbt3Tmp->SetRadius(0.9);
						vtkSphereFid1Mapper->SetInputConnection(vtkSphereHaptFiducbt3Tmp->GetOutputPort());
						vtkActorFidActbt3->SetMapper(vtkSphereFid1Mapper);
						vtkActorFidActbt3->GetProperty()->SetColor(colorFid1Cursor->GetColor3d("White").GetData());
						vtkActorFidActbt3->SetPosition(dHD_XYZ3[0], dHD_XYZ3[1], dHD_XYZ3[2]);

						renderer->AddViewProp(vtkActorFidActbt3);
						renderer->AddActor(vtkActorFidActbt3);
					}

					//if (oPDVol->GetCell(oCellId))
					//{
					//	vtkIdList* oIDList = oPDVol->GetCell(oCellId)->GetPointIds(); //Gives ---- vector subscript out of range error.
					//	oIDList->Print(std::cout);
					//	//std::cout << "Point Id Lst" << oIDList[0] << endl;
					//	//if (bDRD_DEBUG_MODE_DRDVAR)
					//	//{
					//	int nc = oIDList->GetNumberOfIds();
					//	for (int i = 0; i < nc; i++)
					//	{
					//		int iIDNo;
					//		iIDNo = oIDList->GetId(i);

					//		std::cout << "Id  ------------ : " << i << " [ " << iIDNo << " ] " << std::endl;
					//	}
					//	//}
					//}

				}
				catch (...)
				{
					std::cout << " Exception Occured ..." << endl;
				}
				vtkNew<vtkNamedColors> colorFidCursor;
				vtkNew<vtkSphereSource> vtkSphereHaptFiducbt2Tmp; //VTK Cursor
				vtkNew<vtkPolyDataMapper> vtkSphereFidMapper;
				vtkNew<vtkActor> vtkActorFidActbt2;

				vtkSphereHaptFiducbt2Tmp->SetRadius(0.6);
				vtkSphereFidMapper->SetInputConnection(vtkSphereHaptFiducbt2Tmp->GetOutputPort());
				vtkActorFidActbt2->SetMapper(vtkSphereFidMapper);
				vtkActorFidActbt2->GetProperty()->SetColor(colorFidCursor->GetColor3d("Green").GetData());
				vtkActorFidActbt2->SetPosition(dClosestPoint[0], dClosestPoint[1], dClosestPoint[2]);

				renderer->AddViewProp(vtkActorFidActbt2);
				renderer->AddActor(vtkActorFidActbt2);


			}
			else if (dDistance > 0.0)
			{
				actorGoalCursor->SetPosition(dClosestPoint[0], dClosestPoint[1], dClosestPoint[2]);
				oSurfaceClosestPosition = hapticDevicePosition;
				ptrQLabelGoalPosition->setText("Current Goal Pos :  [ " + QString::fromStdString(std::to_string(dClosestPoint[0])) + " , " + QString::fromStdString(std::to_string(dClosestPoint[1])) + " , " + QString::fromStdString(std::to_string(dClosestPoint[2])) + "]");
				ptrQLabelVTKSurfaceClosePoint->setText("Current ClosestPt Value :  [ " + QString::fromStdString(hapticDevicePosition.str(3)) + " : " + QString::fromStdString(std::to_string(dDistance)) + "]" + QString::fromStdString(sForceDataAtPoint));
				//oSurfaceClosestPosition.set(dHD_XYZ3[0] , dHD_XYZ3[1] , dHD_XYZ3[2]);

				vtkNew<vtkNamedColors> colorFidCursor;
				vtkNew<vtkSphereSource> vtkSphereHaptFiducbt2Tmp; //VTK Cursor
				vtkNew<vtkPolyDataMapper> vtkSphereFidMapper;
				vtkNew<vtkActor> vtkActorFidActbt2;

				vtkSphereHaptFiducbt2Tmp->SetRadius(0.6);
				vtkSphereFidMapper->SetInputConnection(vtkSphereHaptFiducbt2Tmp->GetOutputPort());
				vtkActorFidActbt2->SetMapper(vtkSphereFidMapper);
				vtkActorFidActbt2->GetProperty()->SetColor(colorFidCursor->GetColor3d("Orange").GetData());
				vtkActorFidActbt2->SetPosition(dClosestPoint[0], dClosestPoint[1], dClosestPoint[2]);

				renderer->AddViewProp(vtkActorFidActbt2);
				renderer->AddActor(vtkActorFidActbt2);

			}
		}

		if (hdbutton0 == true)
		{
			if (bDRD_DEBUG_MODE_DRDVAR)
			{
				std::cout << "Haptic Device Button 0 state : [ " << hdbutton0 << " ]" << endl;
			}

			//transform0->Translate(dx, 0.0, 0.0);

			vtkNew<vtkNamedColors> colorFidCursor;
			vtkNew<vtkSphereSource> vtkSphereHaptFiducbt2Tmp; //VTK Cursor
			vtkNew<vtkPolyDataMapper> vtkSphereFidMapper;
			vtkNew<vtkActor> vtkActorFidActbt2;

			vtkSphereHaptFiducbt2Tmp->SetRadius(5.0);
			vtkSphereFidMapper->SetInputConnection(vtkSphereHaptFiducbt2Tmp->GetOutputPort());
			vtkActorFidActbt2->SetMapper(vtkSphereFidMapper);
			vtkActorFidActbt2->GetProperty()->SetColor(colorFidCursor->GetColor3d("Violet").GetData());
			vtkActorFidActbt2->SetPosition(dHD_XYZ3[0], dHD_XYZ3[1], dHD_XYZ3[2]);

			renderer->AddViewProp(vtkActorFidActbt2);
			renderer->AddActor(vtkActorFidActbt2);


			// -----------------------------------------------------
			// Saving fiducial point to csv
			// -----------------------------------------------------

			// Get the current system time
			std::chrono::system_clock::time_point now = std::chrono::system_clock::now();

			// Convert the time point to a time_t object
			std::time_t time = std::chrono::system_clock::to_time_t(now);

			// Convert the time_t to a string representation
			std::string placementTimestamp = std::ctime(&time);

			int uiParticipant = 1;
			int uiTrial = 2;

			// Formatting the position of the fiducial point as a string
			std::string fiducialPosition = std::to_string(dHD_XYZ3[0]) + "," + std::to_string(dHD_XYZ3[1]) + "," + std::to_string(dHD_XYZ3[2]);

			std::string filename = "trials.csv";
			CSVData csvData = loadCSV(filename);

			std::unordered_map<std::string, std::string> record = retrieveRecord(csvData, uiParticipant, uiTrial);

			if (!record.empty()) {

				std::string condition = record["Condition"];
				std::string sample = record["Sample"];
				std::string markerPosition = record["MarkerPosition"];

				saveToCSV(uiParticipant, uiTrial, condition, sample, markerPosition, fiducialPosition, dataloadTimestamp, placementTimestamp);

			}
			else {
				std::cout << "No matching record found for Participant " << uiParticipant << " and Trial " << uiTrial << std::endl;
			}


		}
		else
		{
			if (bDRD_DEBUG_MODE_DRDVAR)
			{
				std::cout << "Haptic Device Button 0 state : [ " << hdbutton0 << " ]" << endl;
			}
		}

		//25.Aug.2023 Moved this code out of button as Mr Lucas Sir wanted the crosshairs to show as pointer

		// Try to get view group of the 3D view and jump only those slices.
		int viewGroup = -1; // jump all by default
		int viewJumpSliceMode = vtkMRMLSliceNode::OffsetJumpSlice;
		//vtkMRMLSliceNode::OffsetJumpSlice;
		//vtkMRMLSliceNode::DefaultJumpSlice;
		//vtkMRMLSliceNode::CenteredJumpSlice;
		vtkMRMLSliceNode::JumpAllSlices(ptrTmpMRMLScene, dHD_XYZ3[0], dHD_XYZ3[1], dHD_XYZ3[2], viewJumpSliceMode, viewGroup);

		int crosshairMode = 0;
		crosshairMode = CrosshairNodeHD->GetCrosshairMode();
		CrosshairNodeHD->SetCrosshairMode(vtkMRMLCrosshairNode::ShowIntersection);
		/*
		NoCrosshair = 0,
			ShowBasic,
			ShowIntersection,
			ShowHashmarks,
			ShowAll,
			ShowSmallBasic,
			ShowSmallIntersection,
			CrosshairMode_Last
			*/

		CrosshairNodeHD->SetCrosshairRAS(dHD_XYZ3[0], dHD_XYZ3[1], dHD_XYZ3[2]);



		if (hdbutton1 == true)
		{
			std::cout << "Haptic Device Button 1 state : [ " << hdbutton1 << " ]" << endl;


			//This code will help pan camera view.

			vtkMRMLCameraNode* CameraNodeHD = vtkMRMLCameraNode::SafeDownCast(ptrTmpMRMLScene->GetNthNodeByClass(0, "vtkMRMLCameraNode"));
			double dCamFPoint[3];
			CameraNodeHD->GetFocalPoint(dCamFPoint);
			std::cout << "CameraNodeHD->GetPosition(dCamFPoint) : [ " << dCamFPoint[0] << "," << dCamFPoint[1] << "," << dCamFPoint[2] << " ]" << endl;
			CameraNodeHD->SetFocalPoint(dCamFPoint[0] + 5, dCamFPoint[1] + 5, dCamFPoint[2] + 2);
			CameraNodeHD->GetFocalPoint(dCamFPoint);
			std::cout << "CameraNodeHD->GetPosition(dCamFPoint) New : [ " << dCamFPoint[0] << "," << dCamFPoint[1] << "," << dCamFPoint[2] << " ]" << endl;




		}
		else
		{
			if (bDRD_DEBUG_MODE_DRDVAR)
			{
				std::cout << "Haptic Device Button 1 state : [ " << hdbutton1 << " ]" << endl;
			}
		}

		if (hdbutton2 == true)
		{
			if (bDRD_DEBUG_MODE_DRDVAR)
			{
				std::cout << "Haptic Device Button 2 state : [ " << hdbutton2 << " ]" << endl;
			}


		}
		else
		{
			if (bDRD_DEBUG_MODE_DRDVAR)
			{
				std::cout << "Haptic Device Button 2 state : [ " << hdbutton2 << " ]" << endl;
			}
		}

		if (hdbutton3 == true)
		{
			std::cout << "Haptic Device Button 3 state : [ " << hdbutton3 << " ]" << endl;

			vtkMRMLCameraNode* CameraNodeHD = vtkMRMLCameraNode::SafeDownCast(ptrTmpMRMLScene->GetNthNodeByClass(0, "vtkMRMLCameraNode"));
			double dCamFPoint[3];
			CameraNodeHD->GetPosition(dCamFPoint);
			std::cout << "CameraNodeHD->GetPosition(dCamFPoint) : [ " << dCamFPoint[0] << "," << dCamFPoint[1] << "," << dCamFPoint[2] << " ]" << endl;
			CameraNodeHD->SetPosition(dCamFPoint[0] + 5, dCamFPoint[1] + 5, dCamFPoint[2] + 2);
			CameraNodeHD->GetPosition(dCamFPoint);
			std::cout << "CameraNodeHD->GetPosition(dCamFPoint) New : [ " << dCamFPoint[0] << "," << dCamFPoint[1] << "," << dCamFPoint[2] << " ]" << endl;

		}
		else
		{
			if (bDRD_DEBUG_MODE_DRDVAR)
			{
				std::cout << "Haptic Device Button 3 state : [ " << hdbutton3 << " ]" << endl;
			}
		}
		double dActorGoalCurPos[3];
		if (!(actorGoalCursor == nullptr))
		{
			actorGoalCursor->GetPosition(dActorGoalCurPos);
		}

		double dActor0CurPos[3];
		if (!(actor == nullptr))
		{
			actor->GetPosition(dActor0CurPos);
		}
		double dActor1CurPos[3];
		if (!(actor1 == nullptr))
		{
			actor1->GetPosition(dActor1CurPos);
		}
		double dActor2CurPos[3];
		if (!(actor2 == nullptr))
		{
			actor2->GetPosition(dActor2CurPos);
		}
		if (bDRD_DEBUG_MODE_DRDVAR)
		{

			std::cout << " collide->GetCollisionModeAsString()   : [ " << collide->GetCollisionModeAsString()
				<< "] Number of contact cells is  -- collide->GetCollisionModeAsString() : ["
				<< collide->GetNumberOfContacts() << " ] ";
		}

		vtkRenderWindowInteractor* iren =
			dynamic_cast<vtkRenderWindowInteractor*>(caller);
		iren->GetRenderWindow()->Render();
	}

private:
	int TimerCount = 0;
	vtkActor* actorFiducial = nullptr;

public:
	vtkActor* actor = nullptr;
	vtkActor* actorGoalCursor = nullptr;
	vtkActor* actor1 = nullptr;
	vtkActor* actor2 = nullptr;
	vtkActor* actor3 = nullptr;
	vtkSmartPointer<vtkMatrix4x4> vtk4x4VariableCrosshairTmp;
	//vtkSmartPointer<vtkTransform> transformCb;
	vtkSmartPointer<vtkCollisionDetectionFilter> collideCb;
	vtkSmartPointer<vtkMRMLCrosshairNode> CrosshairNodeHD;
	vtkMRMLScalarVolumeNode* scalarNodeTmp = nullptr;
	vtkMRMLScene* ptrTmpMRMLScene = nullptr;
	vtkRenderer* renderer = nullptr;
	vtkSmartPointer<vtkImplicitPolyDataDistance> oIPDD;
	QLabel* ptrQLabelCellIdValue = nullptr;
	QLabel* ptrQLabelScalarValue = nullptr;
	QLabel* ptrQLabelVTKSurfaceClosePoint = nullptr;
	QLabel* ptrQLabelGoalPosition = nullptr;

	int timerId = 0;
	int maxCount = -1;
	std::string sHapticCurPosClsLocal = "0,0,0";

};






//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerHBVDSciVisPModuleWidgetPrivate : public Ui_qSlicerHBVDSciVisPModuleWidget
{
public:
	qSlicerHBVDSciVisPModuleWidgetPrivate();
	QVector<HapticEditorEventObservation> EventObservations;
};

//-----------------------------------------------------------------------------
// qSlicerHBVDSciVisPModuleWidgetPrivate methods

//-----------------------------------------------------------------------------
qSlicerHBVDSciVisPModuleWidgetPrivate::qSlicerHBVDSciVisPModuleWidgetPrivate()
{
}

//-----------------------------------------------------------------------------
// qSlicerHBVDSciVisPModuleWidget methods

//-----------------------------------------------------------------------------
qSlicerHBVDSciVisPModuleWidget::qSlicerHBVDSciVisPModuleWidget(QWidget* _parent)
	: Superclass(_parent)
	, d_ptr(new qSlicerHBVDSciVisPModuleWidgetPrivate)
{
	bDRD_HBVDSV_glbl = true;
}

//-----------------------------------------------------------------------------
qSlicerHBVDSciVisPModuleWidget::~qSlicerHBVDSciVisPModuleWidget()
{
	qSlicerHBVDSciVisPModuleWidget::onbtnStopHapticLoopCycle();
	bDRD_HBVDSV_glbl = false;
}


///<summary>
/// Any initial condition could be set here.
///</summary> 

//-----------------------------------------------------------------------------
void qSlicerHBVDSciVisPModuleWidget::setup()
{
	Q_D(qSlicerHBVDSciVisPModuleWidget);
	d->setupUi(this);
	this->Superclass::setup();

	//--------------------------------------------------------------------------
	// HAPTIC DEVICE
	//--------------------------------------------------------------------------
	// create a haptic device handler
	handler = new cHapticDeviceHandler();

	qSlicerHBVDSciVisPModuleWidget::onInitializeHaptics();


	// create a sphere (cursor) to represent the haptic device
	vtkSphereHaptCursor->SetRadius(0.5); // Changed from 5.0 to 0.5 DRD 7.Aug.2023
	vtkSphereHapticMapper->SetInputConnection(vtkSphereHaptCursor->GetOutputPort());
	vtkActorHaptCursor->SetMapper(vtkSphereHapticMapper);
	vtkSphereHaptGoalCursor->SetRadius(2.5); // Goal Cursor
	vtkSphereHapticGoalMapper->SetInputConnection(vtkSphereHaptGoalCursor->GetOutputPort());
	vtkActorHaptGoalCursor->SetMapper(vtkSphereHapticGoalMapper);
	//vtkActorHaptCursor->GetProperty()->SetColor(colorHapticCursor->GetColor3d("White").GetData());
	cHapticDeviceInfo info = hapticDevice->getSpecifications();

	dHapticWorldSpaceRadius = info.m_workspaceRadius;

	std::cout << "info.m_modelName : [ " << info.m_modelName << " ] " << endl;
	std::cout << "info.m_manufacturerName : [ " << info.m_manufacturerName << " ] " << endl;

	sinfoHDModelName = info.m_modelName;
	sinfoHDManufacturerName = info.m_manufacturerName;



	// if the device has a gripper, enable the gripper to simulate a user switch
	hapticDevice->setEnableGripperUserSwitch(true);

	/////////////////////////////////////////////////////////////////////
	  // READ HAPTIC DEVICE
	  /////////////////////////////////////////////////////////////////////

	  // read position 
	cVector3d position;
	hapticDevice->getPosition(position);
	if (bDRD_DEBUG_MODE_DRDVAR)
	{

		std::cout << "Haptic Device Position : "
			<< " X  , Y , Z ::"
			<< position.str(4)
			<< std::endl;
	}
	d->lblDisplayHapticPosition->setText("Initial Haptic Device Position :  " + QString::fromStdString(position.str(4)));
	hapticDevicePosition = position;




	int uiParticipant = 1; // Replace with the desired participant value
	int uiTrial = 2;       // Replace with the desired trial value

	//retrieveRecord(csvData, uiParticipant, uiTrial);

}

///<summary>
/// Any initial condition could be set here. Also this is called during initialization of the application.  The volume also gets loaded here.
///</summary> 
void qSlicerHBVDSciVisPModuleWidget::setMRMLScene(vtkMRMLScene* varMRMLScene)
{
	Q_D(qSlicerHBVDSciVisPModuleWidget);

	hdbutton0 = false;
	hdbutton1 = false;
	hdbutton2 = false;
	hdbutton3 = false;
	bSendForces = false;
	btestDeviceMovement = false;
	bVTKPrimitiveSurface = false;
	d->btnTestVTKPrimitiveSurface->setStyleSheet("background-color: rgb(0,155,0);");
	d->btnTestVTKPrimitiveSurface->setText("Surface Forces On");
	d->btnTestHapticDeviceMovement->setStyleSheet("background-color: rgb(0,155,0);");
	d->btnTestHapticDeviceMovement->setText("Test Haptic Device Movement On");
	d->btnForcesOnOff->setStyleSheet("background-color: rgb(0,155,0);");
	d->btnForcesOnOff->setText("Switch Forces On");
	d->FooBar->setVisible(false);
	d->btnInitializeHaptics->setVisible(false);
	d->btnGetNodeData->setVisible(false);
	d->btnDownMoveHaptic->setVisible(false);
	d->btnInMoveHaptic->setVisible(false);
	d->btnLeftMoveHaptic->setVisible(false);
	d->btnOutMoveHaptic->setVisible(false);
	d->btnRightMoveHaptic->setVisible(false);
	d->btnUpMoveHaptic->setVisible(false);
	vtkNew<vtkNamedColors> colors;
	double dRotationAngleAroundAxis = 90.0;

	vtkNew<vtkTransform> tr;
	tr->RotateZ(dRotationAngleAroundAxis);
	tr->GetMatrix(ptrDblRotationMatrix);

	double dRtatAngRev = -90.0;

	vtkNew<vtkTransform> trRev;
	trRev->RotateZ(dRtatAngRev);
	trRev->GetMatrix(ptrDblRotationMatrixRev);

	int numberOfVisibleNodes = 0;
	const char* className = "";
	layoutManager = qSlicerApplication::application()->layoutManager();
	scalarNode = loadVolume(cptrFileName, varMRMLScene);  // Volume gets loaded from the file.

	std::cout << "scalarNode->GetVoxelVectorTypeAsString(0) : [ " << scalarNode->GetVoxelVectorTypeAsString(0) << " ] " << endl;
	std::cout << "scalarNode->GetVoxelVectorTypeAsString(0) : [ " << scalarNode->GetVoxelVectorTypeAsString(0) << " ] " << endl;
	std::cout << "scalarNode->GetVoxelVectorTypeAsString(1) : [ " << scalarNode->GetVoxelVectorTypeAsString(1) << " ] " << endl;
	std::cout << "scalarNode->GetVoxelVectorTypeAsString(2) : [ " << scalarNode->GetVoxelVectorTypeAsString(2) << " ] " << endl;
	std::cout << "scalarNode->GetVoxelVectorTypeAsString(3) : [ " << scalarNode->GetVoxelVectorTypeAsString(3) << " ] " << endl;
	double dCenterofScalarNode[3] = { 0.0 };
	double dSpacingofScalarNode[3] = { 0.0 };
	scalarNode->GetOrigin(dCenterofScalarNode);
	scalarNode->GetSpacing(dSpacingofScalarNode);

	if (scalarNode == nullptr)
	{
		std::cerr << "Not a valid volume: " << cptrFileName << std::endl;
		return;
	}
	vtkStdString tmpColorName;

	/// <summary>
	/// For the actors in 3D layout of 3D Slicer the properties are set here.  Here the testing was done only on one 3D view.
	/// </summary>
	/// <param name="varMRMLScene"></param>
	for (int threeDViewId = 0; threeDViewId < layoutManager->threeDViewCount(); ++threeDViewId)
	{

		if (threeDViewId == 0)
		{
			tmpColorName = "Pink";
		}
		else if (threeDViewId == 1)
		{
			tmpColorName = "Aqua";
		}
		else if (threeDViewId == 2)
		{
			tmpColorName = "Gold";
		}
		else
		{
			tmpColorName = "Plum";
		}
		std::cout << "Color for View ID " << threeDViewId << " : [" << tmpColorName << " ]" << std::endl;
		vtkActorHaptCursor->GetProperty()->SetColor(colorHapticCursor->GetColor3d(tmpColorName).GetData());
		vtkActorHaptGoalCursor->GetProperty()->SetColor(colorHapticCursor->GetColor3d("Black").GetData());
		// Create command for 3D view
		threeDWidget = layoutManager->threeDWidget(threeDViewId);
		std::cout << "3D View Id -- Accesssible Name : " << layoutManager->threeDWidget(threeDViewId)->accessibleName().toStdString() << "  -- " << std::endl;
		threeDView = threeDWidget->threeDView();

		renderWindow = threeDWidget->threeDView()->renderWindow();
		renderer = vtkRenderer::SafeDownCast(renderWindow->GetRenderers()->GetItemAsObject(threeDViewId));
		renderer->AddViewProp(vtkActorHaptCursor);
		renderer->AddActor(vtkActorHaptCursor);
		renderer->AddViewProp(vtkActorHaptGoalCursor);
		renderer->AddActor(vtkActorHaptGoalCursor);
		




		//// Connect interactor events
		interactor = threeDView->interactor();
		interactor->Initialize();

		//oGeomFltr = vtkSmartPointer <vtkGeometryFilter>::New();
		//scalarNode->GetImageD

		vtkSmartPointer<vtkImageData> imageDataPoly = scalarNode ? scalarNode->GetImageData() : nullptr;

		//vtkSmartPointer<vtkAlgorithmOutput> oInputVolumeData_AS_AlgoOutpt = scalarNode->GetImageDataConnection();

		vtkSmartPointer<vtkPolyDataMapper> oPDMapper = vtkSmartPointer<vtkPolyDataMapper>::New();

		//oGeomFltr->SetInputConnection(oInputVolumeData_AS_AlgoOutpt);
		//oGeomFltr->Update();
		//oPDMapper->SetInputConnection(oGeomFltr->GetOutputPort());
		//oPDMapper->Update();

		//std::cout << "Main Project Minor version :" << qSlicerCoreApplication::mainApplicationMinorVersion() << endl;

		//#if Slicer_MAIN_PROJECT_VERSION_MINOR <= 3
		vtkSmartPointer<vtkFlyingEdges3D> surface = vtkSmartPointer<vtkFlyingEdges3D>::New();
		//#else
		//		vtkSmartPointer<vtkSurfaceNets3D> surface = vtkSmartPointer<vtkSurfaceNets3D>::New();
		//#endif

		//vtkSmartPointer<vtkSurfaceNets3D> surface = vtkSmartPointer<vtkSurfaceNets3D>::New();


		vtkSmartPointer<vtkPolyDataNormals> normals = vtkSmartPointer<vtkPolyDataNormals>::New();





#if VTK_MAJOR_VERSION <= 5
		surface->SetInput(volume);
#else
		surface->SetInputData(imageDataPoly);
#endif
		if (std::strcmp(cptrFileName, "example.nrrd") == 0) {
			surface->SetValue(0, 0.999); //example.nrrd : 35 -- for flying edges
		}
		else {
			surface->SetValue(0, 0.1);
		}
		surface->Update();

		//vtkSmartPointer<vtkPolyData> oPDVol = vtkSmartPointer<vtkPolyData>::New(); //Moved to Global
		oPDVol = surface->GetOutput();

		oPDMapper->SetInputData(oPDVol);
		if (!bDRD_Compute_Normal_Glyph)
		{
			normals->SetInputData(oPDVol);
			normals->FlipNormalsOn();			
			normals->ComputeCellNormalsOn();
			normals->Update();

			oFArrNormal = vtkFloatArray::SafeDownCast(normals->GetOutput()->GetCellData()->GetNormals());
			//vtkFloatArray* normalsFloatTngt = vtkFloatArray::SafeDownCast(oPDVol->GetCellData()->GetTangents());
			//vtkFloatArray* normalsFloat = vtkFloatArray::SafeDownCast(oPDVol->GetCellData()->GetNormals());
			oFArrPDPtNormal = vtkFloatArray::SafeDownCast(oPDVol->GetPointData()->GetNormals());
			std::cout << " Float Cell Normals 0: " << oFArrNormal[0] << endl;
			std::cout << " Float Point Poly Data Normals 0: " << oFArrPDPtNormal[0] << endl;
			//vtkSmartPointer<vtkTriangleFilter> oTrigFltr = vtkSmartPointer<vtkTriangleFilter>::New(); N-Help
			//std::cout << " Float Poly Data Tngnt 0: " << normalsFloatTngt[0] << endl;
			//std::cout << " Float Poly Data Normals 0: " << normalsFloat[0] << endl;
			//normals->GetOutput()->FindCell()

			if (bDRD_DEBUG_MODE_DRDVAR)
			{
				int nc = oFArrNormal->GetNumberOfTuples();
				for (int i = 0; i < nc; i++)
				{
					double testDouble[3];
					oFArrNormal->GetTuple(i, testDouble);

					std::cout << "Tuple : " << i << " [ " << testDouble[0] << " , " << " " << testDouble[1] << " , "
						<< testDouble[2] << " ] " << std::endl;
				}
			}

			if (bDRD_DEBUG_MODE_DRDVAR)
			{
				int nc = oFArrPDPtNormal->GetNumberOfTuples();
				for (int i = 0; i < nc; i++)
				{
					double testDouble[3];
					oFArrPDPtNormal->GetTuple(i, testDouble);

					std::cout << "Tuple : " << i << " [ " << testDouble[0] << " , " << " " << testDouble[1] << " , "
						<< testDouble[2] << " ] " << std::endl;
				}
			}
		}


		//oTrigFltr->SetInputData(oGeomFltr->GetOutputDataObject(0));
		//oTrigFltr->Update();
		//objVTKCone->SetHeight(90);
		//objVTKCone->SetRadius(50);
		//objVTKCone->Update();
		//
		// //Below code did not help
		//
		//mapper1->SetInputConnection(objVTKCone->GetOutputPort(0));
		//mapper1->ScalarVisibilityOff();
		//actor1->SetMapper(mapper1);

		//actor1->GetProperty()->BackfaceCullingOn();
		//actor1->SetUserTransform(transform0);
		//actor1->GetProperty()->SetDiffuseColor(
		//	colors->GetColor3d("Red").GetData());

		//mapper1->SetInputConnection(


		//mapper1->SetInputConnection(objVTKCone->GetOutputPort());
		//mapper1->ScalarVisibilityOff();
		//actor1->SetMapper(mapper1);
		//actor1->GetProperty()->SetOpacity(0.99);
		//actor1->GetProperty()->SetColor(1, 0, 0);
		//vtkSmartPointer<vtkPolyDataMapper> mapperFlyingEdge = vtkSmartPointer<vtkPolyDataMapper>::New();
		vtkSmartPointer<vtkActor> actorFlyingEdge = vtkSmartPointer<vtkActor>::New();;

		//mapperFlyingEdge->SetInputConnection(oPDMapper->GetOutputPort());
		//mapperFlyingEdge->ScalarVisibilityOff();

		std::cout << "Node Center() Before :[" << static_cast<double>(dCenterofScalarNode[0]) << " , " << static_cast<double>(dCenterofScalarNode[1]) << " , " << static_cast<double>(dCenterofScalarNode[2]) << "]" << std::endl;


		double* dAct1Org;
		double dImgOrg[3] = { 0.0 };
		double dPolOrg[3] = { 0.0 };
		//actor1->SetPosition(dCenterofScalarNode);
		actorFlyingEdge->SetPosition(dCenterofScalarNode); //SetCenter does not help 
		std::cout << "Node Center() After :[" << static_cast<double>(dCenterofScalarNode[0]) << " , " << static_cast<double>(dCenterofScalarNode[1]) << " , " << static_cast<double>(dCenterofScalarNode[2]) << "]" << std::endl;
		dAct1Org = actorFlyingEdge->GetCenter();
		//actorFlyingEdge->SetOrigin(dCenterofScalarNode); //This solves issue with origin but need to test data
		//actorFlyingEdge->SetScale(dSpacingofScalarNode); //This solves issue of scaling but need to test data
		//actorFlyingEdge->RotateX(-90.0);//Changing order of rotation changes the way it gets dispalyed
		//actorFlyingEdge->RotateY(90.0);
		std::cout << "actor1->GetCenter() :[" << static_cast<double>(dAct1Org[0]) << " , " << static_cast<double>(dAct1Org[1]) << " , " << static_cast<double>(dAct1Org[2]) << "]" << std::endl;
		imageDataPoly->GetCenter(dImgOrg);
		//imageDataPoly->Print(std::cout);
		imageDataPoly->GetOrigin(dPolOrg);
		//imageDataPoly->SetOrigin(dCenterofScalarNode); //This solves issue with origin but need to test data
		//imageDataPoly->SetSpacing(dSpacingofScalarNode);//This solves issue of scaling but need to test data
		std::cout << "imageDataPoly->GetCenter() :[" << static_cast<double>(dImgOrg[0]) << " , " << static_cast<double>(dImgOrg[1]) << " , " << static_cast<double>(dImgOrg[2]) << "]" << std::endl;
		std::cout << "imageDataPoly->GetOrigin() :[" << static_cast<double>(dPolOrg[0]) << " , " << static_cast<double>(dPolOrg[1]) << " , " << static_cast<double>(dPolOrg[2]) << "]" << std::endl;
		imageDataPoly->GetOrigin(dPolOrg);
		std::cout << "imageDataPoly->GetOrigin() :[" << static_cast<double>(dPolOrg[0]) << " , " << static_cast<double>(dPolOrg[1]) << " , " << static_cast<double>(dPolOrg[2]) << "]" << std::endl;
		//implicitPolyDataDistance->SetInput(objVTKCone->GetOutput());
		//implicitPolyDataDistance->SetInput(oTrigFltr->GetOutput());
		//implicitPolyDataDistance->SetInput(oPDMapper);



		if (std::strcmp(cptrFileName, "example.nrrd") == 0 || std::strcmp(cptrFileName, "drdModAxialCone_2_mod.nrrd") == 0 || std::strcmp(cptrFileName, "drdModAxialCone_3.nrrd") == 0 || std::strcmp(cptrFileName, "cone_and_sphere.nrrd") == 0)
		{
			actorFlyingEdge->GetProperty()->SetOpacity(0.0001);
		}
		else if (std::strcmp(cptrFileName, "MR-head_LPI_mod.nrrd") == 0)
		{
			actorFlyingEdge->GetProperty()->SetOpacity(0.003);
		}
		else
		{
			actorFlyingEdge->GetProperty()->SetOpacity(0.59);
		}
		actorFlyingEdge->GetProperty()->SetColor(0, 1, 0);
		oPDMapper->Update();




		actorFlyingEdge->SetMapper(oPDMapper);
		//actorFlyingEdge->SetCoordinateSystemToPhysical();
		//actorFlyingEdge->SetCoordinateSystemRenderer(renderer);

		if (bDRD_Compute_Normal_Glyph)
		{
			normals->SetInputData(oPDVol);
			normals->ComputePointNormalsOn();
			normals->ComputeCellNormalsOff();
			//normals->SetFeatureAngle(0.0);
			normals->FlipNormalsOff();/// Flip Normals off is showing cones in outward direction Where as FlipnormalsOn is showing cones facing inwards.


			normals->Update();

			vtkSmartPointer<vtkPolyDataMapper> franMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
			franMapper->SetInputConnection(normals->GetOutputPort());


			vtkSmartPointer<vtkActor> franActor = vtkSmartPointer<vtkActor>::New();
			franActor->SetMapper(franMapper);
			franActor->GetProperty()->SetColor(colors->GetColor3d("Flesh").GetData());
			franActor->GetProperty()->SetOpacity(0.09);
			// We subsample the dataset because we want to glyph just a subset of
			// the points. Otherwise the display is cluttered and cannot be easily
			// read. The RandomModeOn and SetOnRatio combine to random select one out
			// of every 10 points in the dataset.
			//
			vtkSmartPointer<vtkMaskPoints> ptMask = vtkSmartPointer<vtkMaskPoints>::New();
			ptMask->SetInputConnection(normals->GetOutputPort());
			ptMask->SetOnRatio(5);//Was 10 changed to 100
			ptMask->RandomModeOn();
			ptMask->Update();
			// In this case we are using a cone as a glyph. We transform the cone so
			// its base is at 0,0,0. This is the point where glyph rotation occurs.
			vtkSmartPointer<vtkConeSource> cone = vtkSmartPointer<vtkConeSource>::New();
			cone->SetResolution(6);
			cone->SetHeight(5);
			cone->Update();
			vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
			transform->Translate(0.5, 0.0, 0.0);
			transform->Update();
			vtkSmartPointer<vtkTransformPolyDataFilter> transformF = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
			transformF->SetInputConnection(cone->GetOutputPort());
			transformF->SetTransform(transform);
			transformF->Update();
			// vtkGlyph3D takes two inputs: the input point set (SetInputConnection)
			// which can be any vtkDataSet; and the glyph (SetSourceConnection) which
			// must be a vtkPolyData.  We are interested in orienting the glyphs by the
			// surface normals that we previously generated.
			vtkSmartPointer<vtkGlyph3D> glyph = vtkSmartPointer<vtkGlyph3D>::New();
			glyph->SetInputConnection(ptMask->GetOutputPort());
			glyph->SetSourceConnection(transformF->GetOutputPort());
			glyph->SetVectorModeToUseNormal();
			glyph->SetScaleModeToScaleByVector();
			glyph->SetScaleFactor(0.5); //was 0.004 changing to 1
			glyph->Update();
			vtkSmartPointer<vtkPolyDataMapper> spikeMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
			spikeMapper->SetInputConnection(glyph->GetOutputPort());
			spikeMapper->Update();
			vtkSmartPointer<vtkActor> spikeActor = vtkSmartPointer<vtkActor>::New();
			spikeActor->SetMapper(spikeMapper);
			spikeActor->GetProperty()->SetColor(
				colors->GetColor3d("Green").GetData()); // was Emerald_Green .. changed to Green
			renderer->AddViewProp(spikeActor);
			renderer->AddActor(spikeActor);
			renderer->AddViewProp(franActor);
			renderer->AddActor(franActor);
		}


		renderer->AddViewProp(actorFlyingEdge);
		renderer->AddActor(actorFlyingEdge);
		//if (bDRD_Compute_Normal_Glyph)
		//{

		//}
		//implicitPolyDataDistance->se
		implicitPolyDataDistance->SetInput(oPDVol);
		/// <summary>
		/// The logic for construction of collision detection polygon is being done here.  This need to be tested.
		/// </summary>
		/// <param name="varMRMLScene"></param>
		if (bDRD_Collision_Detect_Polygon)
		{
			ptrNRRDGeometryFilter = vtkSmartPointer <vtkGeometryFilter>::New();
			vtkAlgorithmOutput* ptrNRRDAlgorithmOutput = scalarNode->GetImageDataConnection();
			ptrNRRDGeometryFilter->SetPieceInvariant(1);
			ptrNRRDGeometryFilter->SetInputConnection(ptrNRRDAlgorithmOutput);
			objVTKCone->SetHeight(90);
			objVTKCone->SetRadius(50);

			collide->SetInputConnection(0, vtkSphereHaptCursor->GetOutputPort());
			collide->SetTransform(0, transform0);
			if (bDRD_Collision_Detect_Cone_Construct)
			{
				collide->SetInputConnection(1, objVTKCone->GetOutputPort());
			}
			else
			{
				collide->SetInputConnection(1, ptrNRRDGeometryFilter->GetOutputPort());
			}

			collide->SetMatrix(1, matrix1);
			collide->SetBoxTolerance(0.0);
			collide->SetCellTolerance(0.0);
			collide->SetNumberOfCellsPerNode(2);


			if (contactMode == 0)
			{
				collide->SetCollisionModeToAllContacts();
			}
			else if (contactMode == 1)
			{
				collide->SetCollisionModeToFirstContact();
			}
			else
			{
				collide->SetCollisionModeToHalfContacts();
			}
			collide->GenerateScalarsOn();



			mapper1->SetInputConnection(collide->GetOutputPort(0));
			mapper1->ScalarVisibilityOff();
			actor1->SetMapper(mapper1);

			actor1->GetProperty()->BackfaceCullingOn();
			actor1->SetUserTransform(transform0);
			actor1->GetProperty()->SetDiffuseColor(
				colors->GetColor3d("Yellow").GetData());
			actor1->GetProperty()->SetRepresentationToWireframe();

			mapper2->SetInputConnection(collide->GetOutputPort(1));

			actor2->SetMapper(mapper2);
			actor2->GetProperty()->BackfaceCullingOn();
			actor2->GetProperty()->SetDiffuseColor(
				colors->GetColor3d("Green").GetData());

			actor2->SetUserMatrix(matrix1);


			std::cout << "------------------7.Aug.2023----------Start-------" << std::endl;
			vtkImageData* imageData = scalarNode ? scalarNode->GetImageData() : nullptr;

			//scalarNode->GetDisplayNode();
			vtkSmartPointer<vtkMRMLModelDisplayNode> ptrMDLDispNode = vtkSmartPointer<vtkMRMLModelDisplayNode>::New();
			ptrMDLDispNode = vtkMRMLModelDisplayNode::SafeDownCast(scalarNode->GetDisplayNode());



			vtkPolyData* thresholdedSurface = ptrNRRDGeometryFilter->GetOutput();

			vtkSmartPointer<vtkBinaryLabelmapToClosedSurfaceConversionRule> ptrBLMClosedSurface = vtkSmartPointer<vtkBinaryLabelmapToClosedSurfaceConversionRule>::New();
			double dCenterBoundingBoxOfLoadedImage[3] = { 0.0 };
			imageData->Print(std::cout);
			double dOrientation[3];
			double* dCenterTmpDRD;
			int iCenter[3];
			actor1->SetPosition(dCenterofScalarNode);

			std::cout << "actor1->GetCoordinateSystemAsString() :[" << actor1->GetCoordinateSystemAsString() << "]" << std::endl;
			actor1->GetOrientation(dOrientation);
			std::cout << "actor1->GetOrientation() :[" << dOrientation[0] << ", " << dOrientation[1] << ", " << dOrientation[2] << "]" << std::endl;
			dCenterTmpDRD = actor1->GetCenter();
			std::cout << "actor1->GetCenter() :[" << static_cast<double>(dCenterTmpDRD[0]) << " , " << static_cast<double>(dCenterTmpDRD[1]) << " , " << static_cast<double>(dCenterTmpDRD[2]) << "]" << std::endl;

			std::cout << "actor2->GetCoordinateSystemAsString() :[" << actor2->GetCoordinateSystemAsString() << "]" << std::endl;
			actor2->GetOrientation(dOrientation);
			std::cout << "actor2->GetOrientation() :[" << dOrientation[0] << ", " << dOrientation[1] << ", " << dOrientation[2] << "]" << std::endl;
			actor2->SetPosition(dCenterofScalarNode);

			dCenterTmpDRD = actor2->GetCenter();

			actor1->RotateZ(180.0);

			actor2->RotateZ(180.0);

			vtkNew<vtkNamedColors> colorcontour;
			vtkNew<vtkSphereSource> vtkSpherecontour2Tmp; //VTK Cursor
			vtkNew<vtkPolyDataMapper> vtkSpherecontourMapper;
			vtkNew<vtkActor> vtkActorcontourActbt2;

			vtkSpherecontour2Tmp->SetRadius(5.0);
			vtkSpherecontourMapper->SetInputConnection(vtkSpherecontour2Tmp->GetOutputPort());
			vtkActorcontourActbt2->SetMapper(vtkSpherecontourMapper);
			vtkActorcontourActbt2->GetProperty()->SetColor(colorcontour->GetColor3d("Green").GetData());


			std::cout << "------------------7.Aug.2023---------End--------" << std::endl;
			renderer->AddViewProp(actor1);
			renderer->AddActor(actor1);

			renderer->AddViewProp(actor2);
			renderer->AddActor(actor2);
		}
		else
		{

			curveNode->AddToSceneOn();
			std::cout << "---------------curveNode->AddToSceneOn()" << endl;
		}


		this->CrosshairNode = vtkMRMLCrosshairNode::SafeDownCast(varMRMLScene->GetNthNodeByClass(0, "vtkMRMLCrosshairNode"));

		if (this->CrosshairNode != nullptr)
		{
			std::cout << "this->CrosshairNode.GetPointer() :[" << this->CrosshairNode.GetPointer() << "]" << std::endl;

		}
		else
		{
			std::cout << "  this->CrosshairNode  : [ " << this->CrosshairNode << " ] " << std::endl;
		}

		/// <summary>
		/// The GUI loop/callback is initialized here.
		/// </summary>
		/// <param name="varMRMLScene"></param>
		vtkNew<vtkMatrix4x4> vtk4x4VariableCrosshairVar;




		//transform0->Translate(10, 10, -10.8);

		// Sign up to receive TimerEvent
		vtkNew<vtkTimerCallbackHaptic> cb;
		cb->actor = vtkActorHaptCursor;
		cb->actorGoalCursor = vtkActorHaptGoalCursor;
		cb->actor1 = actor1;
		cb->actor2 = actor2;
		//cb->actor3 = actor3;
		//cb->transformCb = transform0;
		cb->collideCb = collide;
		cb->CrosshairNodeHD = this->CrosshairNode;
		cb->vtk4x4VariableCrosshairTmp = vtk4x4VariableCrosshairVar;
		cb->scalarNodeTmp = scalarNode;
		cb->ptrTmpMRMLScene = varMRMLScene;
		cb->renderer = renderer;
		cb->ptrQLabelScalarValue = d->lblDisplayHapticPosAndScalarVal;
		cb->ptrQLabelVTKSurfaceClosePoint = d->lblVTKPrimitiveClosestPosition;
		cb->ptrQLabelGoalPosition = d->lblGoalPosition;
		cb->oIPDD = implicitPolyDataDistance;
		cb->ptrQLabelCellIdValue = d->lblCellId;
		interactor->AddObserver(vtkCommand::TimerEvent, cb); // Haptic Cursor movement

		int timerId = interactor->CreateRepeatingTimer(100); //Works for 30 but not required
		if (bDRD_DEBUG_MODE_DRDVAR)
		{
			std::cout << "timerId: " << timerId << std::endl;
		}
		// Destroy the timer when maxCount is reached.

		cb->maxCount = 30000;
		cb->timerId = timerId;



	}


	/// <summary>
	/// This is code for modifying slice view. The controller could be called here and made to show the slices in the slice views.
	/// </summary>
	/// <param name="varMRMLScene"></param>
	foreach(QString sliceViewName, layoutManager->sliceViewNames())
	{


		sliceWidget = layoutManager->sliceWidget(sliceViewName);
		std::cout << " sliceViewName : [" << sliceViewName.toStdString() << " ] " << std::endl;
		sliceControlWidget = layoutManager->sliceWidget(sliceViewName)->sliceController();

		interactorSliceRed = sliceWidget->sliceView()->interactor();
		vtkCollection* sliceLogicCollection = layoutManager->mrmlSliceLogics();
		if (bDRD_DEBUG_MODE_DRDVAR)
		{
			std::cout << " 0 " << std::endl;
			std::cout << " Logic Class Names " << sliceViewName.toLocal8Bit().constData() << "  ---  " << sliceLogicCollection->GetClassName() << std::endl;
		}

		//25.Aug.2023 -- Commenting Below as Mr Lucas Sir wanted this to be shown on request by user and not by default.
		//sliceControlWidget->setSliceVisible(TRUE); // Show slice views in 3D window --Equivalent to clicking ‘eye’ icon in the slice view controller.  
		// Gives open GL warning nad error here when 2 --3Dviews are open
		if (bDRD_DEBUG_MODE_DRDVAR)
		{
			std::cout << " 1 " << std::endl;
		}
		if (sliceWidget)
		{
			//std::cout << " 2 " << std::endl;
			//sliceWidget->sliceView()->scheduleRender();
		}

		std::cout << " sliceViewName.contains(Red) : [" << sliceViewName.contains("Red") << " ]" << std::endl;


	}



	appLogicLocal = vtkSlicerApplicationLogic::SafeDownCast(this->appLogic());

	appSelNode = appLogicLocal->GetSelectionNode();
	if (appLogicLocal == nullptr)
	{
		std::cerr << "Failed to obtain application logic" << std::endl;
		//return scalarNode.GetPointer();
	}
	appSelNode->SetActiveVolumeID(scalarNode->GetID());
	if (bDRD_DEBUG_MODE_DRDVAR)
	{
		std::cout << " 11 " << std::endl;
	}
	appLogicLocal->PropagateVolumeSelection(); // To show in Slice views // Gives open GL warning nad error here when 2 --3Dviews are open
	if (bDRD_DEBUG_MODE_DRDVAR)
	{
		std::cout << " 12 " << std::endl;
	}
	appLogicLocal->FitSliceToAll();
	if (bDRD_DEBUG_MODE_DRDVAR)
	{
		std::cout << " 13 " << std::endl;
	}


	scalarNode->GetImageData()->GetExtent(iImageExtArr6);

	// Gives open GL warning nad error here when 2 --3Dviews are open -- Even if the above problem causing lines are commented. The error is present after this function.
}


/// <summary>
/// The button calling this code is commented/hidden.  Use if required.
/// </summary>
void qSlicerHBVDSciVisPModuleWidget::onGetNodeData()
{
	Q_D(qSlicerHBVDSciVisPModuleWidget);
	std::cout << "Inside onGetNodeData" << std::endl;
	qSlicerApplication* app = qSlicerApplication::application();
	vtkMRMLScene* varMRMLSceneTmp = qSlicerApplication::application()->mrmlScene();

	QMessageBox msgBox;

	msgBox.setText(" " + QString::fromStdString("Slice Application Version ") + " :  " + qSlicerApplication::application()->applicationVersion());
	//msgBox.exec();

	hapticDevice->getPosition(positionGlobalVar);
	hapticDevicePosition = positionGlobalVar;

	msgBox.setText(" " + QString::fromStdString("Chai3d -- Get Device Position ") + " :  " + QString::fromStdString(hapticDevicePosition.str(4)));
	//msgBox.exec();

	d->lblDisplayHapticPosition->setText("Current Haptic Device Position after GetNodeData:  " + QString::fromStdString(hapticDevicePosition.str(4)));

}


/// <summary>
/// The callback classes are initialized here, make modification if required.
/// </summary>
/// <param name="caller"></param>
/// <param name="eid"></param>
/// <param name="clientData"></param>
/// <param name="vtkNotUsed"></param>
void qSlicerHBVDSciVisPModuleWidget::onHapticPositionUpdateCursor(vtkObject* caller,
	unsigned long eid,
	void* clientData,
	void* vtkNotUsed(callData))
{

	// Get and parse client data
	vtkHapticHBVDSciVisPEventCallbackCommand* callbackCommand = reinterpret_cast<vtkHapticHBVDSciVisPEventCallbackCommand*>(clientData);
	qSlicerHBVDSciVisPModuleWidget* self = callbackCommand->EditorWidget.data();
	qMRMLWidget* viewWidget = callbackCommand->ViewWidget.data();
	qMRMLThreeDWidget* view3DWidget = reinterpret_cast<qMRMLThreeDWidget*>(viewWidget);
	qMRMLThreeDView* threeDView = view3DWidget->threeDView();
	cVector3d  cVect3dhapticDevicePositionGot = callbackCommand->cVect3dhapticDevicePositionClassPVar; ;
	vtkSmartPointer<vtkActor> vtkCursorActor = callbackCommand->vtkCursorActor;
	vtkRenderWindowInteractor* iren =
		static_cast<vtkRenderWindowInteractor*>(caller);

	vtkMRMLThreeDViewInteractorStyle* style = vtkMRMLThreeDViewInteractorStyle::SafeDownCast
	(iren ? iren->GetInteractorStyle() : nullptr);
	if (!style)
	{
		//qCritical() << "qMRMLThreeDView::mouseMoveEvent: no valid interactor style.";
		return;
	}

}

/// <summary>
/// This function calls all subfunctions to start using the haptic device
/// </summary>
void qSlicerHBVDSciVisPModuleWidget::onbtnStartHapticLoopCycle()
{
	//bSendForces = true;
	qSlicerHBVDSciVisPModuleWidget::onbtnStopHapticLoopCycle();
	qSlicerHBVDSciVisPModuleWidget::onInitializeHaptics();
	std::cout << "++++++++++++ Haptic Device Move center is clicked ++++++++++++" << std::endl;

	getXYZBasedonHapticWorkSpaceRadius(dHapticWorldSpaceRadius, dXYZPtr, ' ');
	desiredPosition.set(dXYZPtr[0], dXYZPtr[1], dXYZPtr[2]);
	qSlicerHBVDSciVisPModuleWidget::CallHapticThreadandUIUpdateTimer(dXYZPtr[0], dXYZPtr[1], dXYZPtr[2]); // Go to center

	if (bForceLog)
	{
		forceLog.open("force_log.txt", std::ios_base::app);
	}

	// Move the first object
	//int numSteps = 10; //was 100
	//double dx = 1.0 / static_cast<double>(numSteps) * 2.0;
	//transform0->Translate(-numSteps * dx - .3, 0.0, 0.0); // 8.Aug.2023 DRD

}


/// <summary>
/// Here the GUI labels get updated.
/// </summary>
void qSlicerHBVDSciVisPModuleWidget::onHapticThreadCalled()
{

	//Display position after haptic thread is called
	Q_D(qSlicerHBVDSciVisPModuleWidget);
	hapticDevice->getPosition(positionGlobalVar); //Get Latest position of Haptic device
	hapticDevicePosition = positionGlobalVar;
	d->lblDisplayHapticPosition->setText("Current Haptic Device Position -- QTimer Loop:  " + QString::fromStdString(hapticDevicePosition.str(3)));
	//slblRates = cStr(freqCounterGraphics.getFrequency(), 0) + " Hz / " + cStr(freqCounterHaptics.getFrequency(), 0) + " Hz";
	d->lblRates->setText("Refresh Rates (Graphic[Hz]/Haptics[Hz]) : " + QString::fromStdString(cStr(freqCounterGraphics.getFrequency(), 0)) + " Hz / " + QString::fromStdString(cStr(freqCounterHaptics.getFrequency(), 0)) + " Hz");
	if (bDRD_DEBUG_MODE_DRDVAR)
	{
		std::cout << "Current Haptic Device Position -- QTimer Loop:" << hapticDevicePosition.str(3) << std::endl;
	}

	if (bForceLog)
	{
		double magnitude = force.length(); // compute the magnitude of the force vector

		if (std::isnan(magnitude)) // check if magnitude is NaN
		{
			magnitude = 0; // set magnitude to 0 if it is NaN
		}

		forceLog << magnitude << "," << freqCounterHaptics.getFrequency() << "\n";
	}

}


/// <summary>
/// This is the load volume function used to load the volume data to 3D slicer.  The path of the file is also set here.
/// </summary>
/// <param name="volume"></param>
/// <param name="scene"></param>
/// <returns></returns>
vtkMRMLScalarVolumeNode* qSlicerHBVDSciVisPModuleWidget::loadVolume(const char* volume, vtkMRMLScene* scene)
{
	vtkNew<vtkMRMLScalarVolumeDisplayNode> displayNode;
	vtkNew<vtkMRMLScalarVolumeNode> scalarNode;
	vtkNew<vtkMRMLVolumeArchetypeStorageNode> storageNode;
	vtkSmartPointer<vtkPiecewiseFunction> ptrVtkPiecewiseFuOpacities = vtkSmartPointer<vtkPiecewiseFunction>::New();
	QMessageBox msgBox;
	displayNode->SetAutoWindowLevel(true);//Displays in slice views the visible image based on window and level
	displayNode->SetInterpolate(false);
	std::string sDirPath = itksys::SystemTools::GetCurrentWorkingDirectory();
	sDirPath = itksys::SystemTools::GetParentDirectory(sDirPath);
	std::string sFilePath = sDirPath + "/" + volume;
	// Check that the file exists
	storageNode->SetFileName(sFilePath.c_str());
	//if (vtksys::SystemTools::FileExists(volume) == false)
	if (!itksys::SystemTools::FileExists(sFilePath))
	{
		msgBox.setText(" " + QString::fromStdString("Error! in loadVolume") + " :  " + QString::fromStdString(" Did not find file  ") + cptrFileName + " . \n " + sFilePath.c_str());
		msgBox.exec();
		return nullptr;
	}

	if (storageNode->SupportedFileType(volume) == 0)
	{
		return nullptr;
	}
	std::string sVolumeP1 = volume;
	std::string sVolumeP2 = "HBVDinSciVis - ";
	std::string sVolumeCombine = sVolumeP2 + sVolumeP1;
	scalarNode->SetName(sVolumeCombine.c_str());
	scalarNode->SetScene(scene);
	displayNode->SetScene(scene);
	//scalarNode->SetSpacing(1, 1, 1);//23.Nov.2023
	//scalarNode->SetOrigin(0.0, 0.0, 0.0);//23.Nov.2023
	scene->AddNode(storageNode.GetPointer());
	scene->AddNode(displayNode.GetPointer());
	scalarNode->SetAndObserveStorageNodeID(storageNode->GetID());
	scalarNode->SetAndObserveDisplayNodeID(displayNode->GetID());
	scene->AddNode(scalarNode.GetPointer());

	if (std::strcmp(cptrFileName, "MR-head_LPI_mod.nrrd") == 0)
	{
		// For MR-head_LPI -- Makes it invisible
		ptrVtkPiecewiseFuOpacities->AddPoint(0.0, 0.0);
		ptrVtkPiecewiseFuOpacities->AddPoint(21.21, 0.00);
		ptrVtkPiecewiseFuOpacities->AddPoint(63.64, 0.03);
		ptrVtkPiecewiseFuOpacities->AddPoint(374.12, 0.09);
		ptrVtkPiecewiseFuOpacities->AddPoint(659.53, 0.53);
	}
	else if (std::strcmp(cptrFileName, "example.nrrd") == 0 || std::strcmp(cptrFileName, "drdModAxialCone_2_mod.nrrd") == 0 || std::strcmp(cptrFileName, "drdModAxialCone_3.nrrd") == 0 || std::strcmp(cptrFileName, "cone_and_sphere.nrrd") == 0)
	{
		//For Cone and Sphere or Cone
		ptrVtkPiecewiseFuOpacities->AddPoint(0.0, 0.0);
		ptrVtkPiecewiseFuOpacities->AddPoint(6.07, 0.15);
		ptrVtkPiecewiseFuOpacities->AddPoint(39.84, 0.53);
	}
	else
	{
		//For rest
		ptrVtkPiecewiseFuOpacities->AddPoint(0.0, 0.0);
		ptrVtkPiecewiseFuOpacities->AddPoint(6.07, 0.15);
		//ptrVtkPiecewiseFuOpacities->AddPoint(39.84, 0.53);
	}
	try
	{

		storageNode->ReadData(scalarNode.GetPointer());
	}
	catch (itk::ExceptionObject& excep)
	{


		msgBox.setText(" " + QString::fromStdString("Error!") + " :  " + excep.GetDescription());
		msgBox.exec();
		return nullptr;
	}
	catch (...)
	{
		msgBox.setText(" " + QString::fromStdString("Error! in loadVolume ") + " :  " + QString::fromStdString("loadVolume"));
		msgBox.exec();
		return nullptr;
	}
	std::cout << "Storage node Name :" << storageNode->GetName() << std::endl;
	std::cout << "Storage node - GetCoordinateSystemTypeAsString :" << storageNode->GetCoordinateSystemTypeAsString(0) << std::endl;
	// Default color tables are not present in the scene if there is no vtkMRMLColorLogic,
	// therefore we need to create and set the color node manually.
	vtkNew<vtkMRMLColorTableNode> colorNode;
	colorNode->SetTypeToGrey();
	scene->AddNode(colorNode.GetPointer());
	displayNode->SetAndObserveColorNodeID(colorNode->GetID());
	// Obtain volume rendering logic
	vtkSlicerVolumeRenderingLogic* volumeRenderingLogic = vtkSlicerVolumeRenderingLogic::SafeDownCast(this->moduleLogic("VolumeRendering"));
	// Create volume rendering display node
	vtkMRMLVolumeRenderingDisplayNode* vrDisplayNode = volumeRenderingLogic->CreateDefaultVolumeRenderingNodes(scalarNode.GetPointer()); // Shows in 3d widget
	vtkSmartPointer<vtkMRMLVolumePropertyNode> propertyNode = vrDisplayNode->GetVolumePropertyNode();
	propertyNode->SetScalarOpacity(ptrVtkPiecewiseFuOpacities);

	//vtkImageData* imageData = scalarNode ? scalarNode->GetImageData() : nullptr;

	//vtkNew<vtkDiscreteFlyingEdges3D> marchingCubes;
	//marchingCubes->SetInputData(imageData);
	//marchingCubes->ComputeGradientsOn();//DRD 7.Nov.2023
	////marchingCubes->ComputeNormalsOff(); // While computing normals is faster using the flying edges filter,
	//marchingCubes->ComputeNormalsOn();//DRD 7.Nov.2023
	//marchingCubes->Update();


	return scalarNode.GetPointer();
}










/// <summary>
/// The caliberation and initialization of the haptic device takes place here.
/// </summary>
void qSlicerHBVDSciVisPModuleWidget::onInitializeHaptics()
{
	Q_D(qSlicerHBVDSciVisPModuleWidget);
	// get a handle to the first haptic device
	handler->getDevice(hapticDevice, 0);
	// open a connection to haptic device
	hapticDevice->open();
	// calibrate device (if necessary)
	hapticDevice->calibrate();
	// retrieve information about the current haptic device
	  // read position 
	cVector3d position;
	hapticDevice->getPosition(position);
	std::cout << "Haptic Device Position : "
		<< " X  , Y , Z ::"
		<< position.str()
		<< std::endl;
	d->lblDisplayHapticPosition->setText("Haptic Device Position on Initialization : " + QString::fromStdString(position.str()));
}


/// <summary>
/// Thsis is a convience testing function.  This is used for testing the caliberation of the haptic device.
/// </summary>
/// <param name="dHapticWorldSpaceRadiusINm"></param>
/// <param name="dXYZPtr"></param>
/// <param name="UpDownLeftRightInOutinUDLRIO"></param>
void qSlicerHBVDSciVisPModuleWidget::getXYZBasedonHapticWorkSpaceRadius(double dHapticWorldSpaceRadiusINm, double* dXYZPtr, char UpDownLeftRightInOutinUDLRIO)
{
	//dXYZPtr[3] = {0,0,0};
	switch (UpDownLeftRightInOutinUDLRIO)
	{
	case 'U':
		dXYZPtr[0] = (dHapticWorldSpaceRadiusINm - dHapticWorldSpaceRadiusINm);
		dXYZPtr[1] = (dHapticWorldSpaceRadiusINm - dHapticWorldSpaceRadiusINm);
		dXYZPtr[2] = (dHapticWorldSpaceRadiusINm);
		break;
	case 'D':
		dXYZPtr[0] = (dHapticWorldSpaceRadiusINm - dHapticWorldSpaceRadiusINm);
		dXYZPtr[1] = (dHapticWorldSpaceRadiusINm - dHapticWorldSpaceRadiusINm);
		dXYZPtr[2] = (-1) * (dHapticWorldSpaceRadiusINm);
		break;
	case 'L':
		dXYZPtr[0] = (dHapticWorldSpaceRadiusINm - dHapticWorldSpaceRadiusINm);
		dXYZPtr[1] = (-1) * (dHapticWorldSpaceRadiusINm);
		dXYZPtr[2] = (dHapticWorldSpaceRadiusINm - dHapticWorldSpaceRadiusINm);
		break;
	case 'R':
		dXYZPtr[0] = (dHapticWorldSpaceRadiusINm - dHapticWorldSpaceRadiusINm);
		dXYZPtr[1] = (dHapticWorldSpaceRadiusINm);
		dXYZPtr[2] = (dHapticWorldSpaceRadiusINm - dHapticWorldSpaceRadiusINm);
		break;
	case 'I':
		dXYZPtr[0] = (-1) * (dHapticWorldSpaceRadiusINm);
		dXYZPtr[1] = (dHapticWorldSpaceRadiusINm - dHapticWorldSpaceRadiusINm);
		dXYZPtr[2] = (dHapticWorldSpaceRadiusINm - dHapticWorldSpaceRadiusINm);
		break;
	case 'O':
		dXYZPtr[0] = (dHapticWorldSpaceRadiusINm);
		dXYZPtr[1] = (dHapticWorldSpaceRadiusINm - dHapticWorldSpaceRadiusINm);
		dXYZPtr[2] = (dHapticWorldSpaceRadiusINm - dHapticWorldSpaceRadiusINm);
		break;
	default:
		dXYZPtr[0] = 0.000;
		dXYZPtr[1] = 0.000;
		dXYZPtr[2] = 0.000;
	}
	std::cout << "Calculated X Y Z based on Haptic Work Space : [ " << dXYZPtr[0] << " , " << dXYZPtr[1] << " , " << dXYZPtr[2] << " ] " << std::endl;
	//return dXYZPtr;
}


void qSlicerHBVDSciVisPModuleWidget::onbtnLeftMoveHaptic()
{
	qSlicerHBVDSciVisPModuleWidget::onbtnStopHapticLoopCycle();
	if (!(hapticDevice->open())) qSlicerHBVDSciVisPModuleWidget::onInitializeHaptics();
	getXYZBasedonHapticWorkSpaceRadius(dHapticWorldSpaceRadius, dXYZPtr, 'L');
	//desiredPosition.set(0.001, -0.062, 0.001);
	desiredPosition.set(dXYZPtr[0], dXYZPtr[1], dXYZPtr[2]);
	std::cout << "++++++++++++ Haptic Device Move Left is clicked ++++++++++++ " << std::endl;
	qSlicerHBVDSciVisPModuleWidget::CallHapticThreadandUIUpdateTimer(dXYZPtr[0], dXYZPtr[1], dXYZPtr[2]);// yMax - Left to the user  -- Assuming Right hand user
}


void qSlicerHBVDSciVisPModuleWidget::onbtnRightMoveHaptic()
{
	qSlicerHBVDSciVisPModuleWidget::onbtnStopHapticLoopCycle();
	if (!(hapticDevice->open()))qSlicerHBVDSciVisPModuleWidget::onInitializeHaptics();
	getXYZBasedonHapticWorkSpaceRadius(dHapticWorldSpaceRadius, dXYZPtr, 'R');
	desiredPosition.set(dXYZPtr[0], dXYZPtr[1], dXYZPtr[2]);
	//desiredPosition.set(0.000, 0.065, 0.000);
	std::cout << "++++++++++++ Haptic Device Move Right is clicked ++++++++++++ " << std::endl;
	qSlicerHBVDSciVisPModuleWidget::CallHapticThreadandUIUpdateTimer(dXYZPtr[0], dXYZPtr[1], dXYZPtr[2]); // yMax - Right from the user  -- Assuming Right hand user
}


void qSlicerHBVDSciVisPModuleWidget::onbtnUpMoveHaptic()
{
	qSlicerHBVDSciVisPModuleWidget::onbtnStopHapticLoopCycle();
	if (!(hapticDevice->open()))qSlicerHBVDSciVisPModuleWidget::onInitializeHaptics();
	getXYZBasedonHapticWorkSpaceRadius(dHapticWorldSpaceRadius, dXYZPtr, 'U');
	desiredPosition.set(dXYZPtr[0], dXYZPtr[1], dXYZPtr[2]);
	//desiredPosition.set(0.000, 0.000, 0.067); // zMax - Up from the ground  -- Assuming Right hand user and device is placed in the right orientation
	std::cout << "++++++++++++ Haptic Device Move Up is clicked ++++++++++++ " << std::endl;
	qSlicerHBVDSciVisPModuleWidget::CallHapticThreadandUIUpdateTimer(dXYZPtr[0], dXYZPtr[1], dXYZPtr[2]);// zMax - Up from the ground  -- Assuming Right hand user and device is placed in the right orientation
}


void qSlicerHBVDSciVisPModuleWidget::onbtnDownMoveHaptic()
{
	qSlicerHBVDSciVisPModuleWidget::onbtnStopHapticLoopCycle();
	if (!(hapticDevice->open()))qSlicerHBVDSciVisPModuleWidget::onInitializeHaptics();
	getXYZBasedonHapticWorkSpaceRadius(dHapticWorldSpaceRadius, dXYZPtr, 'D');
	desiredPosition.set(dXYZPtr[0], dXYZPtr[1], dXYZPtr[2]);
	//desiredPosition.set(0.001, 0.004, -0.063);
	std::cout << "++++++++++++ Haptic Device Move Down is clicked ++++++++++++ " << std::endl;
	qSlicerHBVDSciVisPModuleWidget::CallHapticThreadandUIUpdateTimer(dXYZPtr[0], dXYZPtr[1], dXYZPtr[2]);// zMax - Down to the ground  -- Assuming Right hand user and device is placed in the right orientation
}


void qSlicerHBVDSciVisPModuleWidget::onbtnInMoveHaptic()
{
	qSlicerHBVDSciVisPModuleWidget::onbtnStopHapticLoopCycle();
	if (!(hapticDevice->open()))qSlicerHBVDSciVisPModuleWidget::onInitializeHaptics();
	getXYZBasedonHapticWorkSpaceRadius(dHapticWorldSpaceRadius, dXYZPtr, 'I');
	desiredPosition.set(dXYZPtr[0], dXYZPtr[1], dXYZPtr[2]);
	//desiredPosition.set(-0.048, 0.002, 0.002);
	std::cout << "++++++++++++ Haptic Device Move In is clicked ++++++++++++ " << std::endl;
	qSlicerHBVDSciVisPModuleWidget::CallHapticThreadandUIUpdateTimer(dXYZPtr[0], dXYZPtr[1], dXYZPtr[2]); // xMin - In - Away the user towards device -- Assuming Right hand user

}


void qSlicerHBVDSciVisPModuleWidget::onbtnOutMoveHaptic()
{
	qSlicerHBVDSciVisPModuleWidget::onbtnStopHapticLoopCycle();
	if (!(hapticDevice->open()))qSlicerHBVDSciVisPModuleWidget::onInitializeHaptics();
	getXYZBasedonHapticWorkSpaceRadius(dHapticWorldSpaceRadius, dXYZPtr, 'O');
	desiredPosition.set(dXYZPtr[0], dXYZPtr[1], dXYZPtr[2]);
	//desiredPosition.set(0.047, 0.000, 0.000);
	std::cout << "++++++++++++ Haptic Device Move Out is clicked ++++++++++++ " << std::endl;
	qSlicerHBVDSciVisPModuleWidget::CallHapticThreadandUIUpdateTimer(dXYZPtr[0], dXYZPtr[1], dXYZPtr[2]); // xMax - Out Towards the user away from device -- Assuming Right hand user

}


/// <summary>
/// This is the function where the haptic thread and the gui callback are called.
/// </summary>
/// <param name="xHapt"></param>
/// <param name="yHapt"></param>
/// <param name="zHapt"></param>
void qSlicerHBVDSciVisPModuleWidget::CallHapticThreadandUIUpdateTimer(float xHapt, float yHapt, float zHapt)
{
	Q_D(qSlicerHBVDSciVisPModuleWidget);
	//--------------------------------------------------------------------------
	// START SIMULATION
	//--------------------------------------------------------------------------

	// create a thread which starts the main haptics rendering loop
	//hapticsThread->start(chai3dRunHapticLoop, CTHREAD_PRIORITY_HAPTICS, &desiredPositionofHaptic);
	hapticsThread = new cThread();
	hapticsThread->start(chai3dRunHapticLoop, CTHREAD_PRIORITY_HAPTICS);

	QTimer* timer = new QTimer(this);
	connect(timer, &QTimer::timeout, this, &qSlicerHBVDSciVisPModuleWidget::onHapticThreadCalled);
	timer->start(1000);
}


/// <summary>
/// The checking of the haptic loop for running state and it's destruction happens here.
/// </summary>
void qSlicerHBVDSciVisPModuleWidget::onbtnStopHapticLoopCycle()
{
	//bSendForces = false;
	//btestDeviceMovement = false;
	std::cout << "Inside Close Callback function--started now" << std::endl;
	// stop the simulation
	simulationRunning = false;
	// wait for graphics and haptics loops to terminate
	while (!simulationFinished)
	{
		cSleepMs(100);
		std::cout << "Inside Close Callback function--Waiting for graphic simulation to finish" << std::endl;
	}
	if (hapticDevice)
	{
		// close haptic device
		hapticDevice->close();
	}
	std::cout << "Inside Close Callback function--Ended now" << std::endl;

	if (hapticsThread != nullptr)
	{
		// delete resources
		try
		{
			hapticsThread->stop();
			hapticsThread->~cThread();
			//delete hapticsThread;
		}
		catch (itk::ExceptionObject& excep)
		{
			std::cout << "Unable to delete hapticThread : " << excep.GetDescription() << endl;
		}
		catch (...)
		{
			std::cout << "Unable to delete hapticThread unknown exception." << endl;
		}
	}

	if (bForceLog)
	{

		forceLog.close();

	}
}


/// <summary>
/// This is a convience function.  The chai3d stores the coordinate in a different format compared to how its required in 3D slicer, therefore this function is used for simplification.
/// </summary>
/// <param name="hapticDevicePositionStr3"></param>
/// <param name="cXYZ"></param>
/// <returns></returns>
double qSlicerHBVDSciVisPModuleWidget::getHapticDevicePositionInDouble(std::string hapticDevicePositionStr3, char cXYZ)
{
	std::string sHapticPos = hapticDevicePositionStr3;
	std::size_t pos = sHapticPos.find(",");      // position of "live" in str
	std::string sHapticPosX = sHapticPos.substr(0, pos);

	if (cXYZ == 'X')
	{
		return std::stod(sHapticPosX);
	}

	std::string sHapticPosXCut = sHapticPos.substr(pos + 1);
	std::size_t posY = sHapticPosXCut.find(",");

	std::string sHapticPosY = sHapticPos.substr(pos + 1, posY);
	if (cXYZ == 'Y')
	{
		return std::stod(sHapticPosY);
	}
	std::string sHapticPosYCut = sHapticPosXCut.substr(posY + 1);
	std::string sHapticPosZ = sHapticPosYCut;

	if (cXYZ == 'Z')
	{
		return std::stod(sHapticPosZ);
	}

	return 0.0;

}

/// <summary>
/// This gets called on butting click to start or stop haptic forces.
/// </summary>
/// <param name="bClickedState"></param>
void qSlicerHBVDSciVisPModuleWidget::onbtnForcesOnOffClicked(bool bClickedState)
{
	Q_D(qSlicerHBVDSciVisPModuleWidget);
	std::cout << "bClickedState  : [ " << bClickedState << " ] " << endl;
	std::cout << "bSendForces  : [ " << bSendForces << " ] " << endl;
	if (!bSendForces && bClickedState)
	{
		//bSendForces = false;
		d->btnForcesOnOff->setStyleSheet("background-color: rgb(155,0,0);");
		d->btnForcesOnOff->setText("Switch Forces Off");
		bSendForces = true;
	}
	else
	{
		d->btnForcesOnOff->setStyleSheet("background-color: rgb(0,155,0);");
		d->btnForcesOnOff->setText("Switch Forces On");
		bSendForces = false;
	}
}


/// <summary>
/// This code gets called on the button click to test the haptic device caliberation.
/// </summary>
/// <param name="bClickedState"></param>
void qSlicerHBVDSciVisPModuleWidget::onbtnTestHapticDeviceMovementClicked(bool bClickedState)
{
	Q_D(qSlicerHBVDSciVisPModuleWidget);
	if (!btestDeviceMovement && bClickedState)
	{
		d->btnTestHapticDeviceMovement->setStyleSheet("background-color: rgb(155,0,0);");
		d->btnTestHapticDeviceMovement->setText("Test Haptic Device Movement Off");
		btestDeviceMovement = true;
		d->btnDownMoveHaptic->setVisible(true);
		d->btnInMoveHaptic->setVisible(true);
		d->btnLeftMoveHaptic->setVisible(true);
		d->btnOutMoveHaptic->setVisible(true);
		d->btnRightMoveHaptic->setVisible(true);
		d->btnUpMoveHaptic->setVisible(true);
	}
	else
	{
		d->btnTestHapticDeviceMovement->setStyleSheet("background-color: rgb(0,155,0);");
		d->btnTestHapticDeviceMovement->setText("Test Haptic Device Movement On");
		btestDeviceMovement = false;
		d->btnDownMoveHaptic->setVisible(false);
		d->btnInMoveHaptic->setVisible(false);
		d->btnLeftMoveHaptic->setVisible(false);
		d->btnOutMoveHaptic->setVisible(false);
		d->btnRightMoveHaptic->setVisible(false);
		d->btnUpMoveHaptic->setVisible(false);
	}
}




/// <summary>
/// This code gets called on the button click to test the haptic device forces on VTK Primitive surfaces
/// </summary>
/// <param name="bClickedState"></param>
void qSlicerHBVDSciVisPModuleWidget::onbtnTestVTKPrimitiveSurfaceClicked(bool bClickedState)
{
	Q_D(qSlicerHBVDSciVisPModuleWidget);
	if (!bVTKPrimitiveSurface && bClickedState)
	{
		d->btnTestVTKPrimitiveSurface->setStyleSheet("background-color: rgb(155,0,0);");
		d->btnTestVTKPrimitiveSurface->setText("Surface Forces Off");
		bVTKPrimitiveSurface = true;
		//d->btnDownMoveHaptic->setVisible(true);

	}
	else
	{
		d->btnTestVTKPrimitiveSurface->setStyleSheet("background-color: rgb(0,155,0);");
		d->btnTestVTKPrimitiveSurface->setText("Surface Forces On");
		bVTKPrimitiveSurface = false;
		//d->btnDownMoveHaptic->setVisible(false);

	}
}