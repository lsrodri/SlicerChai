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

#ifndef __qSlicerHBVDSciVisPFooBarWidget_h
#define __qSlicerHBVDSciVisPFooBarWidget_h

// Qt includes
#include <QWidget>

// FooBar Widgets includes
#include "qSlicerHBVDSciVisPModuleWidgetsExport.h"

class qSlicerHBVDSciVisPFooBarWidgetPrivate;

/// \ingroup Slicer_QtModules_HBVDSciVisP
class Q_SLICER_MODULE_HBVDSCIVISP_WIDGETS_EXPORT qSlicerHBVDSciVisPFooBarWidget
  : public QWidget
{
  Q_OBJECT
public:
  typedef QWidget Superclass;
  qSlicerHBVDSciVisPFooBarWidget(QWidget *parent=0);
  ~qSlicerHBVDSciVisPFooBarWidget() override;

protected slots:

protected:
  QScopedPointer<qSlicerHBVDSciVisPFooBarWidgetPrivate> d_ptr;

private:
  Q_DECLARE_PRIVATE(qSlicerHBVDSciVisPFooBarWidget);
  Q_DISABLE_COPY(qSlicerHBVDSciVisPFooBarWidget);
};

#endif
