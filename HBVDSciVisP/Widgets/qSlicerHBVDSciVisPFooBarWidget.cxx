/*==============================================================================

  Program: 3D Slicer

  Copyright (c) Kitware Inc.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

  This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
  and was partially funded by NIH grant 3P41RR013218-12S1

==============================================================================*/

// FooBar Widgets includes
#include "qSlicerHBVDSciVisPFooBarWidget.h"
#include "ui_qSlicerHBVDSciVisPFooBarWidget.h"

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_HBVDSciVisP
class qSlicerHBVDSciVisPFooBarWidgetPrivate
  : public Ui_qSlicerHBVDSciVisPFooBarWidget
{
  Q_DECLARE_PUBLIC(qSlicerHBVDSciVisPFooBarWidget);
protected:
  qSlicerHBVDSciVisPFooBarWidget* const q_ptr;

public:
  qSlicerHBVDSciVisPFooBarWidgetPrivate(
    qSlicerHBVDSciVisPFooBarWidget& object);
  virtual void setupUi(qSlicerHBVDSciVisPFooBarWidget*);
};

// --------------------------------------------------------------------------
qSlicerHBVDSciVisPFooBarWidgetPrivate
::qSlicerHBVDSciVisPFooBarWidgetPrivate(
  qSlicerHBVDSciVisPFooBarWidget& object)
  : q_ptr(&object)
{
}

// --------------------------------------------------------------------------
void qSlicerHBVDSciVisPFooBarWidgetPrivate
::setupUi(qSlicerHBVDSciVisPFooBarWidget* widget)
{
  this->Ui_qSlicerHBVDSciVisPFooBarWidget::setupUi(widget);
}

//-----------------------------------------------------------------------------
// qSlicerHBVDSciVisPFooBarWidget methods

//-----------------------------------------------------------------------------
qSlicerHBVDSciVisPFooBarWidget
::qSlicerHBVDSciVisPFooBarWidget(QWidget* parentWidget)
  : Superclass( parentWidget )
  , d_ptr( new qSlicerHBVDSciVisPFooBarWidgetPrivate(*this) )
{
  Q_D(qSlicerHBVDSciVisPFooBarWidget);
  d->setupUi(this);
}

//-----------------------------------------------------------------------------
qSlicerHBVDSciVisPFooBarWidget
::~qSlicerHBVDSciVisPFooBarWidget()
{
}
