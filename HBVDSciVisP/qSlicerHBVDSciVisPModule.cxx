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

// HBVDSciVisP Logic includes
#include <vtkSlicerHBVDSciVisPLogic.h>

// HBVDSciVisP includes
#include "qSlicerHBVDSciVisPModule.h"
#include "qSlicerHBVDSciVisPModuleWidget.h"

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerHBVDSciVisPModulePrivate
{
public:
  qSlicerHBVDSciVisPModulePrivate();
};

//-----------------------------------------------------------------------------
// qSlicerHBVDSciVisPModulePrivate methods

//-----------------------------------------------------------------------------
qSlicerHBVDSciVisPModulePrivate::qSlicerHBVDSciVisPModulePrivate()
{
}

//-----------------------------------------------------------------------------
// qSlicerHBVDSciVisPModule methods

//-----------------------------------------------------------------------------
qSlicerHBVDSciVisPModule::qSlicerHBVDSciVisPModule(QObject* _parent)
  : Superclass(_parent)
  , d_ptr(new qSlicerHBVDSciVisPModulePrivate)
{
}

//-----------------------------------------------------------------------------
qSlicerHBVDSciVisPModule::~qSlicerHBVDSciVisPModule()
{
}

//-----------------------------------------------------------------------------
QString qSlicerHBVDSciVisPModule::helpText() const
{
  return "This is a loadable module that can be bundled in an extension.  This helps to connect to a Haptic Device.";
}

//-----------------------------------------------------------------------------
QString qSlicerHBVDSciVisPModule::acknowledgementText() const
{
    QString acknowledgement = QString(
        "<center><table border=\"0\"><tr>"
        "<td><img src=\":Icons/HBVDSciVisP.png\" alt\"HBVDSciVisP\"></td>"
        "<td><img src=\":Icons/ZibLogo.png\" alt\"ZibLogo\"></td>"
        "</tr><tr>"
        "<td><img src=\":Icons/VirtualDissectionIcon.png\" alt\"VirtualDissectionIcon\"></td>"
        "<td><img src=\":Icons/LSRIcon.png\" alt\"LSRIcon\"></td>"
        "</tr></table></center>" +
        tr("This work was partially funded by ZIB grant NXNNXXNNNNNN-NNXN"));
    //std::cout << "Acknowledgement :" << acknowledgement.toStdString() << std::endl;
    return acknowledgement;
}

//-----------------------------------------------------------------------------
QStringList qSlicerHBVDSciVisPModule::contributors() const
{
  QStringList moduleContributors;
  moduleContributors << QString("Dr Lucas Siqueira Rodrigues (Cluster of Excellence Matters of Activity)");
  return moduleContributors;
}

//-----------------------------------------------------------------------------
QIcon qSlicerHBVDSciVisPModule::icon() const
{
  return QIcon(":/Icons/HBVDSciVisP.png");
}

//-----------------------------------------------------------------------------
QStringList qSlicerHBVDSciVisPModule::categories() const
{
  return QStringList() << "Examples";
}

//-----------------------------------------------------------------------------
QStringList qSlicerHBVDSciVisPModule::dependencies() const
{
  return QStringList();
}

//-----------------------------------------------------------------------------
void qSlicerHBVDSciVisPModule::setup()
{
  this->Superclass::setup();
}

//-----------------------------------------------------------------------------
qSlicerAbstractModuleRepresentation* qSlicerHBVDSciVisPModule
::createWidgetRepresentation()
{
  return new qSlicerHBVDSciVisPModuleWidget;
}

//-----------------------------------------------------------------------------
vtkMRMLAbstractLogic* qSlicerHBVDSciVisPModule::createLogic()
{
  return vtkSlicerHBVDSciVisPLogic::New();
}
