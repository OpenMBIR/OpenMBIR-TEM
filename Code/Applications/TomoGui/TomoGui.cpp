/* ============================================================================
 * Copyright (c) 2010, Michael A. Jackson (BlueQuartz Software)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright notice, this
 * list of conditions and the following disclaimer in the documentation and/or
 * other materials provided with the distribution.
 *
 * Neither the name of Michael A. Jackson nor the names of its contributors may
 * be used to endorse or promote products derived from this software without
 * specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
 * USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "TomoGui.h"

#include <iostream>
#include <sstream>
#include <limits>
#include <fstream>

//-- Qt Includes
#include <QtCore/QPluginLoader>
#include <QtCore/QFileInfo>
#include <QtCore/QFile>
#include <QtCore/QDir>
#include <QtCore/QString>
#include <QtCore/QUrl>
#include <QtCore/QThread>
#include <QtCore/QThreadPool>
#include <QtCore/QFileInfoList>
#include <QtCore/QTimer>

#include <QtGui/QApplication>
#include <QtGui/QFileDialog>
#include <QtGui/QCloseEvent>
#include <QtGui/QMessageBox>
#include <QtGui/QListWidget>
#include <QtGui/QStringListModel>
#include <QtGui/QLineEdit>
#include <QtGui/QDoubleValidator>
#include <QtGui/QImage>

// Our Project wide includes
#include "QtSupport/ApplicationAboutBoxDialog.h"
#include "QtSupport/QRecentFileList.h"
#include "QtSupport/QFileCompleter.h"


//-- TomoEngine Includes
#include "TomoEngine/TomoEngine.h"
#include "TomoEngine/TomoEngineVersion.h"
#include "TomoEngine/Common/EIMMath.h"
#include "TomoEngine/SOC/SOCStructures.h"
#include "TomoEngine/SOC/SOCEngine.h"
#include "TomoEngine/IO/MRCHeader.h"
#include "TomoEngine/IO/MRCReader.h"
#include "TomoEngine/Filters/GainsOffsetsReader.h"
#include "TomoEngine/Filters/TargetGainSigmaXEstimation.h"

#include "License/LicenseFiles.h"

#include "CheckBoxDelegate.h"
#include "LayersDockWidget.h"
#include "GainsOffsetsTableModel.h"
#include "ReconstructionArea.h"


#define GRAY_SCALE 1


#define READ_STRING_SETTING(prefs, var, emptyValue)\
  var->setText( prefs.value(#var).toString() );\
  if (var->text().isEmpty() == true) { var->setText(emptyValue); }


#define READ_SETTING(prefs, var, ok, temp, default, type)\
  ok = false;\
  temp = prefs.value(#var).to##type(&ok);\
  if (false == ok) {temp = default;}\
  var->setValue(temp);

#define READ_VALUE(prefs, var, ok, temp, default, type)\
  ok = false;\
  temp = prefs.value(#var).to##type(&ok);\
  if (false == ok) {temp = default;}\
  var = temp;

#define WRITE_STRING_SETTING(prefs, var)\
  prefs.setValue(#var , this->var->text());

#define WRITE_SETTING(prefs, var)\
  prefs.setValue(#var, this->var->value());

#define READ_BOOL_SETTING(prefs, var, emptyValue)\
  { QString s = prefs.value(#var).toString();\
  if (s.isEmpty() == false) {\
    bool bb = prefs.value(#var).toBool();\
  var->setChecked(bb); } else { var->setChecked(emptyValue); } }

#define WRITE_BOOL_SETTING(prefs, var, b)\
    prefs.setValue(#var, (b) );

#define WRITE_CHECKBOX_SETTING(prefs, var)\
    prefs.setValue(#var, var->isChecked() );

#define WRITE_VALUE(prefs, var)\
    prefs.setValue(#var, var);

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
TomoGui::TomoGui(QWidget *parent) :
QMainWindow(parent),
m_OutputExistsCheck(false),
m_LayersPalette(NULL),
m_StopAnimation(true),
m_CurrentCorner(0),
m_WorkerThread(NULL),
m_MultiResSOC(NULL),
  m_SingleSliceReconstructionActive(false),
  m_FullReconstrucionActive(false),
#if defined(Q_WS_WIN)
m_OpenDialogLastDirectory("C:\\")
#else
m_OpenDialogLastDirectory("~/")
#endif
{
  setupUi(this);
  setupGui();

  m_AnimationTimer = new QTimer(this);
  connect(m_AnimationTimer, SIGNAL(timeout() ), this, SLOT(stepForwardFromTimer() ));

#if defined (Q_OS_MAC)
  QSettings prefs(QSettings::NativeFormat, QSettings::UserScope, QCoreApplication::organizationDomain(), QCoreApplication::applicationName());
#else
  QSettings prefs(QSettings::IniFormat, QSettings::UserScope, QCoreApplication::organizationDomain(), QCoreApplication::applicationName());
#endif
  readSettings(prefs);
  readWindowSettings(prefs);

  QRecentFileList* recentFileList = QRecentFileList::instance();
  connect(recentFileList, SIGNAL (fileListChanged(const QString &)), this, SLOT(updateBaseRecentFileList(const QString &)));
  // Get out initial Recent File List
  this->updateBaseRecentFileList(QString::null);
  qRegisterMetaType<QVector<double> >("QVector<double>");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
TomoGui::~TomoGui()
{

}


// -----------------------------------------------------------------------------
//  Called when the main window is closed.
// -----------------------------------------------------------------------------
void TomoGui::closeEvent(QCloseEvent *event)
{
  qint32 err = checkDirtyDocument();
  if (err < 0)
  {
    event->ignore();
  }
  else
  {
#if defined (Q_OS_MAC)
  QSettings prefs(QSettings::NativeFormat, QSettings::UserScope, QCoreApplication::organizationDomain(), QCoreApplication::applicationName());
#else
  QSettings prefs(QSettings::IniFormat, QSettings::UserScope, QCoreApplication::organizationDomain(), QCoreApplication::applicationName());
#endif
    writeSettings(prefs);
    writeWindowSettings(prefs);
    event->accept();
  }
}


// -----------------------------------------------------------------------------
//  Read the prefs from the local storage file
// -----------------------------------------------------------------------------
void TomoGui::readSettings(QSettings &prefs)
{
  QString val;
  bool ok;
  qint32 i;
  double d;
  prefs.beginGroup("Parameters");

  READ_STRING_SETTING(prefs, inputMRCFilePath, "");
  // This will auto load the MRC File
  on_inputMRCFilePath_textChanged(inputMRCFilePath->text());

  READ_STRING_SETTING(prefs, inputBrightFieldFilePath, "");
  // This will auto load the MRC File
  on_inputBrightFieldFilePath_textChanged(inputBrightFieldFilePath->text());


  READ_STRING_SETTING(prefs, initialReconstructionPath, "");
  READ_STRING_SETTING(prefs, tempDirPath, "");
  READ_STRING_SETTING(prefs, reconstructedVolumeFileName, "");

  READ_STRING_SETTING(prefs, sampleThickness, "");
  READ_STRING_SETTING(prefs, targetGain, "0")
  READ_STRING_SETTING(prefs, sigmaX, ".00000045");

  READ_SETTING(prefs, numResolutions, ok, i, 1, Int);
  READ_SETTING(prefs, finalResolution, ok, i, 1, Int);
  READ_SETTING(prefs, outerIterations, ok, i, 1, Int);
  READ_SETTING(prefs, innerIterations, ok, i, 1, Int);

  READ_BOOL_SETTING(prefs, useDefaultOffset, false);
  READ_STRING_SETTING(prefs, defaultOffset, "0");
  READ_STRING_SETTING(prefs, stopThreshold, "0.001");
  READ_SETTING(prefs, mrf, ok, d, 1.2, Double);
  READ_BOOL_SETTING(prefs, extendObject, false);


//  READ_BOOL_SETTING(prefs, useSubVolume, false);
  READ_STRING_SETTING(prefs, xMin, "0");
  READ_STRING_SETTING(prefs, xMax, "0");
  READ_STRING_SETTING(prefs, yMin, "0");
  READ_STRING_SETTING(prefs, yMax, "0");
//  READ_STRING_SETTING(prefs, zMin, "0");
//  READ_STRING_SETTING(prefs, zMax, "0");



  ok = false;
  i = prefs.value("tiltSelection").toInt(&ok);
  if (false == ok) {i = 0;}
  tiltSelection->setCurrentIndex(i);

}

// -----------------------------------------------------------------------------
//  Write our prefs to file
// -----------------------------------------------------------------------------
void TomoGui::writeSettings(QSettings &prefs)
{
  prefs.beginGroup("Parameters");
  WRITE_STRING_SETTING(prefs, inputMRCFilePath);
  WRITE_STRING_SETTING(prefs, inputBrightFieldFilePath);

  WRITE_STRING_SETTING(prefs, initialReconstructionPath);
  WRITE_STRING_SETTING(prefs, tempDirPath);
  WRITE_STRING_SETTING(prefs, reconstructedVolumeFileName);

  WRITE_STRING_SETTING(prefs, sampleThickness);
  WRITE_STRING_SETTING(prefs, targetGain)
  WRITE_STRING_SETTING(prefs, sigmaX);
  WRITE_SETTING(prefs, numResolutions);
  WRITE_SETTING(prefs, finalResolution);
  WRITE_SETTING(prefs, outerIterations);
  WRITE_SETTING(prefs, innerIterations);

  WRITE_CHECKBOX_SETTING(prefs, useDefaultOffset);
  WRITE_STRING_SETTING(prefs, defaultOffset);
  WRITE_STRING_SETTING(prefs, stopThreshold);

  WRITE_SETTING(prefs, mrf);
  WRITE_CHECKBOX_SETTING(prefs, extendObject);

//  WRITE_BOOL_SETTING(prefs, useSubVolume, useSubVolume->isChecked());
  WRITE_STRING_SETTING(prefs, xMin);
  WRITE_STRING_SETTING(prefs, xMax);
  WRITE_STRING_SETTING(prefs, yMin);
  WRITE_STRING_SETTING(prefs, yMax);
//  WRITE_STRING_SETTING(prefs, zMin);
//  WRITE_STRING_SETTING(prefs, zMax);


  prefs.setValue("tiltSelection", tiltSelection->currentIndex());

  prefs.endGroup();

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::readWindowSettings(QSettings &prefs)
{
  bool ok = false;
  prefs.beginGroup("WindodwSettings");
  if (prefs.contains(QString("Geometry")) )
  {
    QByteArray geo_data = prefs.value(QString("Geometry")).toByteArray();
    ok = restoreGeometry(geo_data);
  }

  if (prefs.contains(QString("Layout")))
  {
    QByteArray layout_data = prefs.value(QString("Layout")).toByteArray();
    restoreState(layout_data);
  }
  prefs.endGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::writeWindowSettings(QSettings &prefs)
{
  prefs.beginGroup("WindodwSettings");
  QByteArray geo_data = saveGeometry();
  QByteArray layout_data = saveState();
  prefs.setValue(QString("Geometry"), geo_data);
  prefs.setValue(QString("Layout"), layout_data);
  prefs.endGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::on_actionSave_Config_File_triggered()
{
  QString proposedFile = m_OpenDialogLastDirectory + QDir::separator() + "Tomo-Config.config";
  QString file = QFileDialog::getSaveFileName(this, tr("Save Tomo Configuration"),
                                              proposedFile,
                                              tr("*.config") );
  if ( true == file.isEmpty() ){ return;  }
  QFileInfo fi(file);
  m_OpenDialogLastDirectory = fi.absolutePath();
  QSettings prefs(file, QSettings::IniFormat, this);
  writeSettings(prefs);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::on_actionLoad_Config_File_triggered()
{
  QString file = QFileDialog::getOpenFileName(this, tr("Select Configuration File"),
                                                 m_OpenDialogLastDirectory,
                                                 tr("Configuration File (*.config)") );
  if ( true == file.isEmpty() ){return;  }
  QFileInfo fi(file);
  m_OpenDialogLastDirectory = fi.absolutePath();
  QSettings prefs(file, QSettings::IniFormat, this);
  readSettings(prefs);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::on_actionParameters_triggered()
{
  parametersDockWidget->show();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::on_actionLayers_Palette_triggered()
{
  m_LayersPalette->show();
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
#define ZOOM_MENU(var, menu, slot)\
  {\
  QAction* action = new QAction(menu);\
  action->setText( #var );\
  QString actionName("action_z" #var "Action");\
  action->setObjectName(actionName);\
  zoomMenu->addAction(action);\
  connect(action, SIGNAL(triggered()), this, SLOT(slot())); \
}

#define ZOOM_MENU_SLOT_DEF(var, index)\
void TomoGui::z##var##_triggered() {\
  zoomButton->setText(#var " % ");\
  m_GraphicsView->setZoomIndex(index);\
}

ZOOM_MENU_SLOT_DEF(10, 0);
ZOOM_MENU_SLOT_DEF(25, 1);
ZOOM_MENU_SLOT_DEF(50, 2);
ZOOM_MENU_SLOT_DEF(100, 3);
ZOOM_MENU_SLOT_DEF(125, 4);
ZOOM_MENU_SLOT_DEF(150, 5);
ZOOM_MENU_SLOT_DEF(200, 6);
ZOOM_MENU_SLOT_DEF(400, 7);
ZOOM_MENU_SLOT_DEF(600, 8);

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::on_fitToWindow_clicked()
{
  m_GraphicsView->setZoomIndex(9);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::setupGui()
{
#ifdef Q_WS_MAC
  // Adjust for the size of the menu bar which is at the top of the screen not in the window
  QSize mySize = size();
  mySize.setHeight(mySize.height() - 30);
  resize(mySize);
#endif

  reconstructedVolumeFileName->setText("");

  QMenu* zoomMenu = new QMenu(this);
  ZOOM_MENU(10, zoomMenu, z10_triggered);
  ZOOM_MENU(25, zoomMenu, z25_triggered);
  ZOOM_MENU(50, zoomMenu, z50_triggered);
  ZOOM_MENU(100, zoomMenu, z100_triggered);
  ZOOM_MENU(125, zoomMenu, z125_triggered);
  ZOOM_MENU(150, zoomMenu, z150_triggered);
  ZOOM_MENU(200, zoomMenu, z200_triggered);
  ZOOM_MENU(400, zoomMenu, z400_triggered);
  ZOOM_MENU(600, zoomMenu, z600_triggered);

  zoomButton->setMenu(zoomMenu);

  m_GraphicsView->setTomoGui(this);
  connect(m_GraphicsView, SIGNAL(fireSingleSliceSelected()),
          this, SLOT(singleSlicePlaneSet()));

  // Just place a really big white image to get our GUI to layout properly
  QImage image(1000, 1000, QImage::Format_ARGB32_Premultiplied);
  image.fill(0);
  m_GraphicsView->loadBaseImageFile(image);

//  if (m_LayersPalette == NULL)
//  {
//    m_LayersPalette = new LayersDockWidget(this);
//    m_LayersPalette->setGraphicsView(m_GraphicsView);
//    m_LayersPalette->setVisible(false);
//    m_LayersPalette->setFeatures(QDockWidget::AllDockWidgetFeatures);
//    m_LayersPalette->setAllowedAreas(Qt::AllDockWidgetAreas);
//    m_LayersPalette->getUseColorTable()->setEnabled(enableUserDefinedAreas->isChecked());
//  }

#if 0
  compositeModeCB->blockSignals(true);

  compositeModeCB->insertItem(0, "Exclusion", QVariant(EmMpm_Constants::Exclusion));
  compositeModeCB->insertItem(1, "Difference", QVariant(EmMpm_Constants::Difference));
  compositeModeCB->insertItem(2, "Alpha Blend", QVariant(EmMpm_Constants::Alpha_Blend));
#endif

#if 0
  compositeModeCB->insertItem(2, "Plus");
  compositeModeCB->insertItem(3, "Multiply");
  compositeModeCB->insertItem(4, "Screen");
  compositeModeCB->insertItem(5, "Darken");
  compositeModeCB->insertItem(6, "Lighten");
  compositeModeCB->insertItem(7, "Color Dodge");
  compositeModeCB->insertItem(8, "Color Burn");
  compositeModeCB->insertItem(9, "Hard Light");
  compositeModeCB->insertItem(10, "Soft Light");

  compositeModeCB->insertItem(12, "Destination");
  compositeModeCB->insertItem(13, "Source Over");
  compositeModeCB->insertItem(14, "Destination Over");
  compositeModeCB->insertItem(15, "Source In");
  compositeModeCB->insertItem(16, "Dest In");

  compositeModeCB->insertItem(17, "Dest Out");
  compositeModeCB->insertItem(18, "Source Atop");
  compositeModeCB->insertItem(19, "Dest Atop");

  compositeModeCB->insertItem(20, "Overlay");
  compositeModeCB->insertItem(21, "Clear");
#endif

#if 0
  compositeModeCB->setCurrentIndex(2);
  compositeModeCB->blockSignals(false);

  compositeModeCB->setEnabled(false);
#endif

//  connect (m_GraphicsView, SIGNAL(fireBaseImageFileLoaded(const QString &)),
//           this, SLOT(baseImageFileLoaded(const QString &)), Qt::QueuedConnection);

  connect(m_GraphicsView, SIGNAL(fireOverlayImageFileLoaded(const QString &)),
          this, SLOT(overlayImageFileLoaded(const QString &)), Qt::QueuedConnection);

  connect(m_GraphicsView, SIGNAL(fireUserInitAreaAdded(ReconstructionArea*)),
          this, SLOT(userInitAreaAdded(ReconstructionArea*)), Qt::QueuedConnection);

  connect(zoomIn, SIGNAL(clicked()), m_GraphicsView, SLOT(zoomIn()), Qt::QueuedConnection);
  connect(zoomOut, SIGNAL(clicked()), m_GraphicsView, SLOT(zoomOut()), Qt::QueuedConnection);

  QFileCompleter* com = new QFileCompleter(this, false);
  inputMRCFilePath->setCompleter(com);
  QObject::connect(com, SIGNAL(activated(const QString &)), this, SLOT(on_inputMRCFilePath_textChanged(const QString &)));

  QFileCompleter* com4 = new QFileCompleter(this, false);
  reconstructedVolumeFileName->setCompleter(com4);
  QObject::connect(com4, SIGNAL(activated(const QString &)), this, SLOT(on_reconstructedVolumeFileName_textChanged(const QString &)));

  // setup the Widget List
  m_WidgetList << inputMRCFilePath << inputMRCFilePathBtn;

  setWidgetListEnabled(true);

  m_ImageWidgets << zoomIn << zoomOut << fitToWindow << layersPalette;
  setImageWidgetsEnabled(false);

  m_GainsOffsetsTableModel = NULL;
#if 1
  // Setup the TableView and Table Models
  QHeaderView* headerView = new QHeaderView(Qt::Horizontal, gainsOffsetsTableView);
  headerView->setResizeMode(QHeaderView::Interactive);
  gainsOffsetsTableView->setHorizontalHeader(headerView);
  headerView->show();

  m_GainsOffsetsTableModel = new GainsOffsetsTableModel;
  m_GainsOffsetsTableModel->setInitialValues();
  gainsOffsetsTableView->setModel(m_GainsOffsetsTableModel);
  QAbstractItemDelegate* idelegate = m_GainsOffsetsTableModel->getItemDelegate();
  gainsOffsetsTableView->setItemDelegateForColumn(GainsOffsetsTableModel::TiltIndex, idelegate);
  gainsOffsetsTableView->setItemDelegateForColumn(GainsOffsetsTableModel::A_Tilt, idelegate);
  gainsOffsetsTableView->setItemDelegateForColumn(GainsOffsetsTableModel::B_Tilt, idelegate);
//  gainsOffsetsTableView->setItemDelegateForColumn(GainsOffsetsTableModel::Gains, idelegate);
 // gainsOffsetsTableView->setItemDelegateForColumn(GainsOffsetsTableModel::Offsets, idelegate);

  QAbstractItemDelegate* cbDelegate = new CheckBoxDelegate;
  gainsOffsetsTableView->setItemDelegateForColumn(GainsOffsetsTableModel::Exclude, cbDelegate);
#endif

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::on_layersPalette_clicked()
{
//  m_LayersPalette->setVisible(true);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
bool TomoGui::sanityCheckOutputDirectory(QLineEdit* le, QString msgTitle)
{

  if (le->text().isEmpty() == true)
  {
    QMessageBox::critical(this, msgTitle,
                          "The output directory has NOT been set. Please set a directory path and try again.",
                          QMessageBox::Ok | QMessageBox::Default);
    return false;
  }

  if (verifyPathExists(le->text(), le) == false)
  {
    QString msg("The Output Directory '");
    msg.append(le->text()).append("'\ndoes not exist. Would you like to create it?");
    int ret = QMessageBox::warning(this, msgTitle,
                                   msg,
                                   QMessageBox::Yes | QMessageBox::Default,
                                   QMessageBox::No);
    if (ret == QMessageBox::No)
    {
      return false;
    }
    else if (ret == QMessageBox::Yes)
    {
      QDir outputDir(le->text());
      if (outputDir.exists() == false)
      {
        bool ok = outputDir.mkpath(".");
        if (ok == false)
        {
          QMessageBox::critical(this,
                                tr("Output Directory Creation"),
                                tr("The output directory could not be created."),
                                QMessageBox::Ok);
          return false;
        }
        else
        {
          return true;
        }
      }
    }
  }
  verifyPathExists(le->text(), le);
  return true;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::on_m_SingleSliceReconstructionBtn_clicked()
{
  // First make sure we are not already running a reconstruction
  if(m_SingleSliceReconstructionBtn->text().compare("Cancel") == 0)
  {
    if(m_MultiResSOC != NULL)
    {
      std::cout << "canceling from GUI...." << std::endl;
      emit cancelPipeline();
    }
    return;
  }

  // Create a Worker Thread that will run the Reconstruction
  if(m_WorkerThread != NULL)
  {
    m_WorkerThread->wait(); // Wait until the thread is complete
    delete m_WorkerThread; // Kill the thread
    m_WorkerThread = NULL;
  }
  m_WorkerThread = new QThread(); // Create a new Thread Resource

  m_MultiResSOC = new QMultiResolutionSOC(NULL);

  // Move the Reconstruction object into the thread that we just created.
  m_MultiResSOC->moveToThread(m_WorkerThread);
  initializeSOCEngine(false);

  /* Connect the signal 'started()' from the QThread to the 'run' slot of the
   * SOCEngine object. Since the SOCEngine object has been moved to another
   * thread of execution and the actual QThread lives in *this* thread then the
   * type of connection will be a Queued connection.
   */
  // When the thread starts its event loop, start the Reconstruction going
  connect(m_WorkerThread, SIGNAL(started()), m_MultiResSOC, SLOT(run()));

  // When the Reconstruction ends then tell the QThread to stop its event loop
  connect(m_MultiResSOC, SIGNAL(finished() ), m_WorkerThread, SLOT(quit()));

  // When the QThread finishes, tell this object that it has finished.
  connect(m_WorkerThread, SIGNAL(finished()), this, SLOT( singleSliceComplete() ));

  // If the use clicks on the "Cancel" button send a message to the Reconstruction object
  // We need a Direct Connection so the
  connect(this, SIGNAL(cancelPipeline() ), m_MultiResSOC, SLOT (on_CancelWorker() ), Qt::DirectConnection);

  // Send Progress from the Reconstruction to this object for display
  connect(m_MultiResSOC, SIGNAL (updateProgress(int)), this, SLOT(pipelineProgress(int) ));

  // Send progress messages from Reconstruction to this object for display
  connect(m_MultiResSOC, SIGNAL (progressMessage(QString)), this, SLOT(addProgressMessage(QString) ));

  // Send progress messages from Reconstruction to this object for display
  connect(m_MultiResSOC, SIGNAL (warningMessage(QString)), this, SLOT(addWarningMessage(QString) ));

  // Send progress messages from Reconstruction to this object for display
  connect(m_MultiResSOC, SIGNAL (errorMessage(QString)), this, SLOT(addErrorMessage(QString) ));

  connect(m_MultiResSOC, SIGNAL(progressImageIsReady(QString)),
          this, SLOT(loadProgressMRCFile(QString) ));

  setWidgetListEnabled(false);
  emit pipelineStarted();
  m_WorkerThread->start();
  m_SingleSliceReconstructionBtn->setText("Cancel");
  m_SingleSliceReconstructionActive = true;
  m_FullReconstrucionActive = false;
}



// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::on_m_GoBtn_clicked()
{
  // First make sure we are not already running a reconstruction
  if(m_GoBtn->text().compare("Cancel") == 0)
  {
    if(m_MultiResSOC != NULL)
    {
      std::cout << "canceling from GUI...." << std::endl;
      emit cancelPipeline();
    }
    return;
  }

  // Make sure we have an output directory setup and created
  if(false == sanityCheckOutputDirectory(tempDirPath, QString("Tomography Reconstruction")))
  {
    return;
  }

  // Now check the input file to make sure it does exist
  QFileInfo fi(inputMRCFilePath->text());
  if(fi.exists() == false)
  {
    QMessageBox::critical(this, tr("Input File File Error"), tr("Input File does not exist. Please check the path."), QMessageBox::Ok);
    return;
  }

  // Make sure we have a name for the output file
  if(reconstructedVolumeFileName->text().isEmpty() == true)
  {
    QMessageBox::critical(this, tr("Output File Error"), tr("Please select a file name for the reconstructed file to be saved as."), QMessageBox::Ok);
    return;
  }
  // We have a name, make sure the user wants to over write the file
  QFile file(reconstructedVolumeFileName->text());
  if(file.exists() == true)
  {
    int ret = QMessageBox::warning(this, tr("EIM Tomo GUI"), tr("The Output File Already Exists\nDo you want to over write the existing file?"), QMessageBox::No
                                       | QMessageBox::Default, QMessageBox::Yes, QMessageBox::Cancel);
    if(ret == QMessageBox::Cancel)
    {
      return;
    }
    else if(ret == QMessageBox::Yes)
    {
      setOutputExistsCheck(true);
    }
    else
    {
      QString outputFile = m_OpenDialogLastDirectory + QDir::separator() + "Untitled.raw";
      outputFile = QFileDialog::getSaveFileName(this, tr("Save Output File As ..."), outputFile, tr("RAW (*.raw)"));
      if(!outputFile.isNull())
      {
        setCurrentProcessedFile("");
        setOutputExistsCheck(true);
      }
      else // The user clicked cancel from the save file dialog

      {
        return;
      }
    }
  }

  // Create a Worker Thread that will run the Reconstruction
  if(m_WorkerThread != NULL)
  {
    m_WorkerThread->wait(); // Wait until the thread is complete
    delete m_WorkerThread; // Kill the thread
    m_WorkerThread = NULL;
  }
  m_WorkerThread = new QThread(); // Create a new Thread Resource

  m_MultiResSOC = new QMultiResolutionSOC(NULL);

  // Move the Reconstruction object into the thread that we just created.
  m_MultiResSOC->moveToThread(m_WorkerThread);
  initializeSOCEngine(true);

  /* Connect the signal 'started()' from the QThread to the 'run' slot of the
   * SOCEngine object. Since the SOCEngine object has been moved to another
   * thread of execution and the actual QThread lives in *this* thread then the
   * type of connection will be a Queued connection.
   */
  // When the thread starts its event loop, start the Reconstruction going
  connect(m_WorkerThread, SIGNAL(started()), m_MultiResSOC, SLOT(run()));

  // When the Reconstruction ends then tell the QThread to stop its event loop
  connect(m_MultiResSOC, SIGNAL(finished() ), m_WorkerThread, SLOT(quit()));

  // When the QThread finishes, tell this object that it has finished.
  connect(m_WorkerThread, SIGNAL(finished()), this, SLOT( pipelineComplete() ));

  // If the use clicks on the "Cancel" button send a message to the Reconstruction object
  // We need a Direct Connection so the
  connect(this, SIGNAL(cancelPipeline() ),
          m_MultiResSOC, SLOT (on_CancelWorker() ), Qt::DirectConnection);

  // Send Progress from the Reconstruction to this object for display
  connect(m_MultiResSOC, SIGNAL (updateProgress(int)),
          this, SLOT(pipelineProgress(int) ));

  // Send progress messages from Reconstruction to this object for display
  connect(m_MultiResSOC, SIGNAL (progressMessage(QString)),
          this, SLOT(addProgressMessage(QString) ));

  // Send progress messages from Reconstruction to this object for display
  connect(m_MultiResSOC, SIGNAL (warningMessage(QString)),
          this, SLOT(addWarningMessage(QString) ));

  // Send progress messages from Reconstruction to this object for display
  connect(m_MultiResSOC, SIGNAL (errorMessage(QString)),
          this, SLOT(addErrorMessage(QString) ));

  connect(m_MultiResSOC, SIGNAL(progressImageIsReady(QString)),
          this, SLOT(loadProgressMRCFile(QString) ));

  setWidgetListEnabled(false);
  emit pipelineStarted();
  m_WorkerThread->start();
  m_GoBtn->setText("Cancel");
  m_SingleSliceReconstructionActive = false;
  m_FullReconstrucionActive = true;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::initializeSOCEngine(bool fullReconstruction)
{
  QString path;
  path = QDir::toNativeSeparators(inputMRCFilePath->text());
  m_MultiResSOC->setInputFile(path.toStdString());

  path = QDir::toNativeSeparators(tempDirPath->text());
  m_MultiResSOC->setTempDir(path.toStdString());

  path = QDir::toNativeSeparators(reconstructedVolumeFileName->text());
  m_MultiResSOC->setOutputFile(path.toStdString());

  path = QDir::toNativeSeparators(inputBrightFieldFilePath->text());
  m_MultiResSOC->setBrightFieldFile(path.toStdString());

  path = QDir::toNativeSeparators(initialReconstructionPath->text());
  m_MultiResSOC->setInitialReconstructionFile(path.toStdString());

  bool ok = false;
  if (fullReconstruction == true)
  {
    m_MultiResSOC->setNumberResolutions(numResolutions->value());
  }
  else
  {
    m_MultiResSOC->setNumberResolutions(1);
  }
  m_MultiResSOC->setFinalResolution(finalResolution->value());
  m_MultiResSOC->setSampleThickness(sampleThickness->text().toFloat(&ok));
  m_MultiResSOC->setTargetGain(targetGain->text().toFloat(&ok));
  m_MultiResSOC->setStopThreshold(stopThreshold->text().toFloat(&ok));
  m_MultiResSOC->setOuterIterations(outerIterations->value());
  m_MultiResSOC->setInnerIterations(innerIterations->value());
  m_MultiResSOC->setSigmaX(sigmaX->text().toFloat(&ok));
  m_MultiResSOC->setMRFShapeParameter(mrf->value());
  m_MultiResSOC->setDefaultOffsetValue(defaultOffset->text().toFloat(&ok));
  m_MultiResSOC->setUseDefaultOffset(useDefaultOffset->isChecked());
  m_MultiResSOC->setTiltSelection(static_cast<SOC::TiltSelection>(tiltSelection->currentIndex()));
  m_MultiResSOC->setExtendObject(extendObject->isChecked());
  m_MultiResSOC->setDefaultVariance(defaultVariance->text().toFloat(&ok));
  m_MultiResSOC->setInitialReconstructionValue(defaultInitialRecon->text().toFloat(&ok));

  m_MultiResSOC->setInterpolateInitialReconstruction(interpolateInitialRecontruction->isChecked());

  AdvancedParametersPtr advParams = AdvancedParametersPtr(new AdvancedParameters);
  SOCEngine::InitializeAdvancedParams(advParams);
  m_MultiResSOC->setAdvParams(advParams);


  std::vector<uint16_t> subvolume(6);
  subvolume[2] = 0;
  subvolume[5] = m_nTilts->text().toUShort(&ok) - 1;
  if (fullReconstruction == true)
  {
    subvolume[0] = xMin->text().toUShort(&ok);
    subvolume[3] = xMax->text().toUShort(&ok);
    subvolume[1] = yMin->text().toUShort(&ok);
    subvolume[4] = yMax->text().toUShort(&ok);
    m_MultiResSOC->setSubvolume(subvolume);
  }
  else
  {
    QLineF line = m_GraphicsView->getXZPlane();
  //  std::cout << "p1: " << line.p1().x() << ", " << line.p1().y()
  //   << "   p2: " << line.p2().x() << ", " << line.p2().y() << std::endl;

    unsigned short x = m_XDim->text().toUShort(&ok);
    // Only reconstruct the middle section of data along the x axis
    subvolume[0] = 0 + x/4;
    subvolume[3] = x - x/4;

    int ySlices = 3 * finalResolution->value();
    // Try and center on the line best we can
    subvolume[1] = line.p1().y() - (ySlices/2);
    subvolume[4] = subvolume[1] + ySlices;

    path = QDir::toNativeSeparators(tempDirPath->text());
    QString tempFolder = QDir::tempPath() + QDir::separator() + QString("EIMTomo");
    m_MultiResSOC->setTempDir(tempFolder.toStdString());
  }
  m_MultiResSOC->setSubvolume(subvolume);

  std::vector<uint8_t> viewMasks;
  if (NULL != m_GainsOffsetsTableModel) {
  QVector<bool> excludedViews = m_GainsOffsetsTableModel->getExcludedTilts();
  for (int i = 0; i < excludedViews.size(); ++i)
  {
      if(excludedViews[i] == true)
      {
        viewMasks.push_back(i);
      }
    }
  }
  m_MultiResSOC->setViewMasks(viewMasks);

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::singleSlicePlaneSet()
{
    m_SingleSliceReconstructionBtn->setEnabled(true);
    QLineF line = m_GraphicsView->getXZPlane();
    int y = line.y1();
    ySingleSliceValue->setText(QString::number(y));
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::singleSliceComplete()
{
  //std::cout << "TomoGui::singleSliceComplete" << std::endl;
  m_SingleSliceReconstructionBtn->setText("Single Slice Reconstruction");
  setWidgetListEnabled(true);
  this->progressBar->setValue(0);
  QString reconVolumeFile = QString::fromStdString(m_MultiResSOC->getTempDir()) + QDir::separator() +
          finalResolution->text() + QString("x") + QDir::separator() + QString::fromStdString(ScaleOffsetCorrection::ReconstructedMrcFile);

  loadSingleSliceReconstruction(reconVolumeFile);


  // Remove all the files that just got created:
  // Remove the Reconstruction.bin file
  QString path = QString::fromStdString(m_MultiResSOC->getTempDir()) + QDir::separator() + finalResolution->text() + QString("x")+ QDir::separator();
  {
    QString filePath = path + QString::fromStdString(m_MultiResSOC->getOutputFile());
    QFile f(filePath);
    f.remove();
  }
  {
    QString filePath = path + ScaleOffsetCorrection::CostFunctionFile.c_str();
    QFile f(filePath);
    f.remove();
  }
  {
    QString filePath = path + ScaleOffsetCorrection::DetectorResponseFile.c_str();
    QFile f(filePath);
    f.remove();
  }
  {
    QString filePath = path + ScaleOffsetCorrection::FinalGainParametersFile.c_str();
    QFile f(filePath);
    f.remove();
  }
  {
    QString filePath = path + ScaleOffsetCorrection::FinalOffsetParametersFile.c_str();
    QFile f(filePath);
    f.remove();
  }
  {
    QString filePath = path + ScaleOffsetCorrection::FinalVariancesFile.c_str();
    QFile f(filePath);
    f.remove();
  }
  {
    QString filePath = path + ScaleOffsetCorrection::ReconstructedBinFile.c_str();
    QFile f(filePath);
    f.remove();
  }
  {
    QString filePath = path + ScaleOffsetCorrection::ReconstructedMrcFile.c_str();
    QFile f(filePath);
    f.remove();
  }
  {
    QString filePath = path + ScaleOffsetCorrection::ReconstructedVtkFile.c_str();
    QFile f(filePath);
    f.remove();
  }
  {
    QString filePath = path + ScaleOffsetCorrection::VoxelProfileFile.c_str();
    QFile f(filePath);
    f.remove();
  }
  // Delete the top level directory
  QDir dir(path);
  dir.rmdir(path);

  m_FullReconstrucionActive = false;
  m_SingleSliceReconstructionActive = false;
  emit pipelineEnded();
  m_MultiResSOC->deleteLater();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::loadSingleSliceReconstruction(QString reconMRCFilePath)
{
  MRCHeader header;
  MRCReader::Pointer reader = MRCReader::New(true);
  int err = reader->readHeader(reconMRCFilePath.toStdString(), &header);
  if(err < 0)
  {
    FREE_FEI_HEADERS( header.feiHeaders )
    return;
  }
  // Read the first image from the file
  int voxelMin[3] =
  { 0, 0, 0 };
  int voxelMax[3] =
  { header.nx, header.ny, header.nz };
  size_t nVoxels = (voxelMax[0] - voxelMin[0]) * (voxelMax[1] - voxelMin[1]) * (voxelMax[2] - voxelMin[2]);

  err = reader->read(reconMRCFilePath.toStdString(), NULL, NULL);
  m_CurrentCorner = 0; // Reset the corner

  char dataType = 0;

  if(err >= 0)
  {
    switch(header.mode)
    {
      case 1:
        dataType = 1;
        break;
      case 2:
        dataType = 2;
        break;
      default:
        std::cout << "Only float and 16 bit signed integers from the MRC file are supported" << std::endl;
        return;
    }
  }

  FREE_FEI_HEADERS( header.feiHeaders )

  QImage image;
  if (dataType == 1) {
    qint16* data = reinterpret_cast<qint16*>(reader->getDataPointer());
    image = xzSigned16CrossSection(data,  nVoxels, voxelMin, voxelMax);
  }
  else if (dataType == 2)
  {
    float* data = reinterpret_cast<float*>(reader->getDataPointer());
    image = xzFloatCrossSection(data, nVoxels, voxelMin, voxelMax);
  }

  image = image.mirrored(false, true);

  m_GraphicsView_2->loadBaseImageFile(image);

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QImage TomoGui::xzFloatCrossSection(float* data, size_t nVoxels, int* voxelMin, int* voxelMax)
{
  float dmax = std::numeric_limits<float>::min();
  float dmin = std::numeric_limits<float>::max();
  for (size_t i = 0; i < nVoxels; ++i)
  {
    if(data[i] > dmax) dmax = data[i];
    if(data[i] < dmin) dmin = data[i];
  }

//  std::cout << "Min float MRC Value:" << dmin << std::endl;
//  std::cout << "Max float MRC Value:" << dmax << std::endl;

  //Scale all the values to 0 and 255 in place over writing the float values with 32 bit ints
  int* iData = reinterpret_cast<int*>(data);
  int imax = std::numeric_limits<int>::min();
  int imin = std::numeric_limits<int>::max();

  for (size_t i = 0; i < nVoxels; ++i)
  {
    iData[i] = (data[i]-dmin) / (dmax - dmin) * 255.0;
    if(iData[i] > imax) imax = iData[i];
    if(iData[i] < imin) imin = iData[i];
  }

//  std::cout << "Min int MRC Value:" << imin << std::endl;
//  std::cout << "Max int MRC Value:" << imax << std::endl;

  QVector<QRgb> colorTable;
  // Generate a Color Table

  int numColors = 256; // We are going to fix this at 256 colors

  // generate the color table

  colorTable.resize(numColors);
  dmin = 0.0;
  dmax = 255;
  float range = 256;
  float r, g, b;
  for (int i = 0; i < numColors; i++)
  {
    int16_t val = static_cast<int16_t>(dmin + ((float)i / numColors) * range);
    getColorCorrespondingTovalue(val, r, g, b, dmax, dmin);
    colorTable[i] = qRgba(r * 255, g * 255, b * 255, 255);
  }

  // Create an RGB Image
  QImage xzImage(voxelMax[0], voxelMax[2], QImage::Format_ARGB32);
//  int xzIdx = 0;
  imax = std::numeric_limits<int>::min();
  imin = std::numeric_limits<int>::max();
  int idx = 0;
  for (int z = voxelMin[2]; z < voxelMax[2]; ++z)
  {
    QImage image(voxelMax[0], voxelMax[1], QImage::Format_ARGB32);

    for (int y = voxelMin[1]; y < voxelMax[1]; ++y)
    {
      for (int x = voxelMin[0]; x < voxelMax[0]; ++x)
      {
       // idx = (voxelMax[0] * voxelMax[1] * z) + (x * voxelMax[1]) + y;
        int iData_idx = iData[idx];
        image.setPixel(x, y, qRgb(iData_idx, iData_idx, iData_idx)); // m_ColorTable[colorIndex]);
        if(iData[idx] > imax) imax = iData[idx];
        if(iData[idx] < imin) imin = iData[idx];
        ++idx;
        if(y == 0)
        {
          xzImage.setPixel(x, z, qRgb(iData[idx], iData[idx], iData[idx]));
        }
      }
    }

    QString fname("/tmp/single_slice_z_");
    fname.append(QString::number(z)).append(".tif");
    image.save(fname);
  }
//  std::cout << "Min int QImage Value:" << imin << std::endl;
//  std::cout << "Max int QImage Value:" << imax << std::endl;
//  xzImage.save("/tmp/xz_image.tif");
  return xzImage;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::loadProgressMRCFile(QString filePath)
{
//  std::cout << "Loading Progress MRC File: " << filePath.toStdString() << std::endl;
  loadSingleSliceReconstruction(filePath);
  QFile fi(filePath);
  fi.remove();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QImage TomoGui::xzSigned16CrossSection(qint16* data, size_t nVoxels, int* voxelMin, int* voxelMax)
{
  qint16 dmax = std::numeric_limits<qint16>::min();
  qint16 dmin = std::numeric_limits<qint16>::max();
  for (size_t i = 0; i < nVoxels; ++i)
  {
    if(data[i] > dmax) dmax = data[i];
    if(data[i] < dmin) dmin = data[i];
  }

  std::cout << "Min MRC Value:" << dmin << std::endl;
  std::cout << "Max MRC Value:" << dmax << std::endl;

  QVector<QRgb> colorTable;
  // Generate a Color Table
  float max = static_cast<float>(dmax);
  float min = static_cast<float>(dmin);
  int numColors = static_cast<int>((max - min) + 1);

  // generate the color table
  if(colorTable.size() != numColors)
  {
    colorTable.resize(numColors);
    float range = max - min;

    float r, g, b;
    for (int i = 0; i < numColors; i++)
    {
      int16_t val = static_cast<int16_t>(min + ((float)i / numColors) * range);
      getColorCorrespondingTovalue(val, r, g, b, max, min);
      colorTable[i] = qRgba(r * 255, g * 255, b * 255, 255);
    }
  }

  // Create an RGB Image
  QImage image(voxelMax[0], voxelMax[2], QImage::Format_ARGB32);

  int idx = 0;
  for (int z = voxelMin[2]; z < voxelMax[2]; ++z)
  {
    for (int x = voxelMin[0]; x < voxelMax[0]; ++x)
    {
      idx = (voxelMin[0] * z) + x;
      int colorIndex = data[idx] - static_cast<int>(dmin);
      image.setPixel(x, z, m_ColorTable[colorIndex]);
    }
  }
  return image;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::pipelineComplete()
{
  // std::cout << "ReconstructionWidget::threadFinished()" << std::endl;
  m_GoBtn->setText("Reconstruct");
  setWidgetListEnabled(true);
  this->progressBar->setValue(0);
  emit pipelineEnded();
  m_MultiResSOC->deleteLater();
  m_FullReconstrucionActive = false;
  m_SingleSliceReconstructionActive = false;

#if 0
  if (m_LayersPalette != NULL)
  {
    m_LayersPalette->getSegmentedImageCheckBox()->setEnabled(true);
    m_LayersPalette->getOriginalImageCheckBox()->setChecked(true);
    m_LayersPalette->getSegmentedImageCheckBox()->setChecked(false);
    m_LayersPalette->getCompositeTypeComboBox()->setCurrentIndex(EmMpm_Constants::Alpha_Blend);
    m_LayersPalette->getOpacitySlider()->setEnabled(true);
    m_LayersPalette->getOpacitySpinBox()->setEnabled(true);
    m_LayersPalette->getCompositeTypeComboBox()->setEnabled(true);
  }
#endif

  setCurrentImageFile(inputMRCFilePath->text());

  QString s = tempDirPath->text();
  s = s.append(QDir::separator()).append(reconstructedVolumeFileName->text());
  setCurrentProcessedFile(s);

//  m_GraphicsView->loadOverlayImageFile(s);


  setWindowTitle(m_CurrentImageFile);
  setWidgetListEnabled(true);

  // Make sure the image manipulating widgets are enabled
  setImageWidgetsEnabled(true);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::pipelineProgress(int val)
{
  this->progressBar->setValue( val );
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::addErrorMessage(QString message)
{
  QString title = "TomoGui Error";
  displayDialogBox(title, message, QMessageBox::Critical);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::addWarningMessage(QString message)
{
  QString title = "TomoGui Warning";
  displayDialogBox(title, message, QMessageBox::Warning);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::addProgressMessage(QString message)
{
  if (NULL != this->statusBar()) {
    this->statusBar()->showMessage(message);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::on_outputDirectoryPathBtn_clicked()
{
  bool canWrite = false;
  QString aDir = QFileDialog::getExistingDirectory(this, tr("Select Output Directory"), m_OpenDialogLastDirectory,
                                            QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

  m_OpenDialogLastDirectory = aDir;
  if (!m_OpenDialogLastDirectory.isNull())
  {
    QFileInfo fi(aDir);
    canWrite = fi.isWritable();
    if (canWrite) {
      this->tempDirPath->setText(m_OpenDialogLastDirectory);
    }
    else
    {
      QMessageBox::critical(this, tr("Output Directory Selection Error"),
                            tr("The Output directory is not writable by your user. Please make it writeable by changing the permissions or select another directory"),
                            QMessageBox::Ok);
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::on_reconstructedVolumeFileNameBtn_clicked()
{
  QString outputFile = m_OpenDialogLastDirectory + QDir::separator() + "Untitled";
  outputFile = QFileDialog::getSaveFileName(this, tr("Save Output File As ..."), outputFile, tr("All files (*.*)"));
  if (outputFile.isEmpty())
  {
    return;
  }

  QFileInfo fi(outputFile);
  QFileInfo fi2(fi.absolutePath());
  if (fi2.isWritable() == true) {
    reconstructedVolumeFileName->setText(fi.fileName());
  }
  else
  {
    QMessageBox::critical(this, tr("Output File Error"),
                          tr("The parent directory of the output file is not writable by your user. Please make it writeable by changing the permissions or select another directory"),
                          QMessageBox::Ok);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::on_inputMRCFilePathBtn_clicked()
{
  //std::cout << "on_actionOpen_triggered" << std::endl;
  QString imageFile =
      QFileDialog::getOpenFileName(this, tr("Select MRC Data file"), m_OpenDialogLastDirectory, tr("MRC Files (*.mrc *.ali)"));

  if (true == imageFile.isEmpty())
  {
    return;
  }
  inputMRCFilePath->setText(imageFile);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::on_inputMRCFilePath_textChanged(const QString & filepath)
{
  if (verifyPathExists(inputMRCFilePath->text(), inputMRCFilePath))
  {
    // Read the header info from the file and populate the GUI with those values
    readMRCHeader(filepath);
    // Now load up the first tilt from the file
    m_CurrentImage = loadMRCTiltImage(filepath, 0);
    // Be sure to properly orient the image which will in turn load the image into
    // the graphics scene
    on_originCB_currentIndexChanged(originCB->currentIndex());

    m_GainsFile = ""; // We are reading a new .mrc file so we probably need a new Gains Offsets File

    setWindowTitle(filepath);
    this->setWindowFilePath(filepath);

#if 0
    m_LayersPalette->getOriginalImageCheckBox()->setChecked(true);
    m_LayersPalette->getSegmentedImageCheckBox()->setChecked(false);
#endif

    // Tell the RecentFileList to update itself then broadcast those changes.
    QRecentFileList::instance()->addFile(filepath);
    setWidgetListEnabled(true);
    setImageWidgetsEnabled(true);
    updateBaseRecentFileList(filepath);
  }


}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::on_inputBrightFieldFilePathBtn_clicked()
{
  //std::cout << "on_actionOpen_triggered" << std::endl;
  QString imageFile =
      QFileDialog::getOpenFileName(this, tr("Select BrightField Data file"), m_OpenDialogLastDirectory, tr("BrightField Files (*.mrc *.ali)"));

  if (true == imageFile.isEmpty())
  {
    return;
  }
  inputBrightFieldFilePath->setText(imageFile);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::on_inputBrightFieldFilePath_textChanged(const QString & filepath)
{
  if (verifyPathExists(inputBrightFieldFilePath->text(), inputBrightFieldFilePath))
  {
    // Read the header info from the file and populate the GUI with those values
//    readBrightFieldHeader(filepath);
    // Now load up the first tilt from the file
//    loadBrightFieldTiltImage(filepath, 0);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::on_reconstructedVolumeFileName_textChanged(const QString & text)
{
  //  verifyPathExists(outputMRCFilePath->text(), movingImageFile);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::on_tempDirPath_textChanged(const QString & text)
{
  verifyPathExists(tempDirPath->text(), tempDirPath);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::on_initialReconstructionPathBtn_clicked()
{
  QString reconFile =
      QFileDialog::getOpenFileName(this, tr("Select Initial Reconstruction file"), m_OpenDialogLastDirectory, tr("All Files (*.*)"));

  if (true == reconFile.isEmpty())
  {
    return;
  }
  initialReconstructionPath->setText(reconFile);
}



// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::on_initialReconstructionPath_textChanged(const QString & text)
{
  verifyPathExists(initialReconstructionPath->text(), initialReconstructionPath);
}




// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::setWidgetListEnabled(bool b)
{
  foreach (QWidget* w, m_WidgetList)
  {
    w->setEnabled(b);
  }

  if (b == true) {

  }

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::setImageWidgetsEnabled(bool b)
{
  foreach (QWidget* w, m_ImageWidgets)
  {
    w->setEnabled(b);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
bool TomoGui::verifyOutputPathParentExists(QString outFilePath, QLineEdit* lineEdit)
{
  QFileInfo fileinfo(outFilePath);
  QDir parent(fileinfo.dir());
  return parent.exists();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
bool TomoGui::verifyPathExists(QString outFilePath, QLineEdit* lineEdit)
{
  QFileInfo fileinfo(outFilePath);
  if (false == fileinfo.exists() && lineEdit != NULL)
  {
    lineEdit->setStyleSheet("border: 1px solid red;");
  }
  else if (lineEdit != NULL)
  {
    lineEdit->setStyleSheet("");
  }
  return fileinfo.exists();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
qint32 TomoGui::checkDirtyDocument()
{
  qint32 err = -1;
  {
    err = 1;
  }
  return err;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::updateBaseRecentFileList(const QString &file)
{
  // Clear the Recent Items Menu
  this->menu_FixedRecentFiles->clear();

  // Get the list from the static object
  QStringList files = QRecentFileList::instance()->fileList();
  foreach (QString file, files)
    {
      QAction* action = new QAction(this->menu_FixedRecentFiles);
      action->setText(QRecentFileList::instance()->parentAndFileName(file));
      action->setData(file);
      action->setVisible(true);
      this->menu_FixedRecentFiles->addAction(action);
      connect(action, SIGNAL(triggered()), this, SLOT(openRecentBaseImageFile()));
    }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::openRecentBaseImageFile()
{
  QAction *action = qobject_cast<QAction *>(sender());
  if (action)
  {
    //std::cout << "Opening Recent file: " << action->data().toString().toStdString() << std::endl;
    QString file = action->data().toString();
    inputMRCFilePath->setText( file );
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::on_actionOpenMRCFile_triggered()
{
  //std::cout << "on_actionOpen_triggered" << std::endl;
  QString file = QFileDialog::getOpenFileName(this, tr("Open Image File"),
    m_OpenDialogLastDirectory,
    tr("Images (*.tif *.tiff *.bmp *.jpg *.jpeg *.png)") );

  if ( true == file.isEmpty() )
  {
    return;
  }
  QFileInfo fi(file);
  m_OpenDialogLastDirectory = fi.absolutePath();
  inputMRCFilePath->setText( file );
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::on_actionOpenOverlayImage_triggered()
{
  //std::cout << "on_actionOpen_triggered" << std::endl;
  QString imageFile = QFileDialog::getOpenFileName(this, tr("Open Segmented Image File"),
    m_OpenDialogLastDirectory,
    tr("Images (*.tif *.tiff *.bmp *.jpg *.jpeg *.png)") );

  if ( true == imageFile.isEmpty() )
  {
    return;
  }
  QFileInfo fi(imageFile);
  m_OpenDialogLastDirectory = fi.absolutePath();
  openOverlayImage(imageFile);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::on_actionSaveCanvas_triggered()
{
  QImage image = m_GraphicsView->getOverlayImage();
#if 0
  if (m_LayersPalette->getOriginalImageCheckBox()->isChecked()
       && m_LayersPalette->getSegmentedImageCheckBox()->isChecked())
  {
    image = m_GraphicsView->getCompositedImage();
  }
  int err = 0;
#endif
  int err = 0;
  QString outputFile = this->m_OpenDialogLastDirectory + QDir::separator() + "Segmented.tif";
  outputFile = QFileDialog::getSaveFileName(this, tr("Save Processed Image As ..."), outputFile, tr("Images (*.tif *.bmp *.jpg *.png)"));
  if (!outputFile.isEmpty())
  {
    bool ok = image.save(outputFile);
    if (ok == true)
    {
      //TODO: Set a window title or something
    }
    else
    {
      //TODO: Add in a GUI dialog to help explain the error or give suggestions.
      err = -1;
    }
  }

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::on_actionAbout_triggered()
{
  ApplicationAboutBoxDialog about(EIMTOMO::LicenseList, this);
  QString an = QCoreApplication::applicationName();
  QString version("");
  version.append(TomoEngine::Version::PackageComplete.c_str());
  about.setApplicationInfo(an, version);
  about.exec();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::on_actionExit_triggered()
{
  this->close();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::openOverlayImage(QString processedImage)
{
  if ( true == processedImage.isEmpty() ) // User cancelled the operation
  {
    return;
  }
  m_GraphicsView->loadOverlayImageFile(processedImage);

#if 0
  if (m_LayersPalette != NULL)
  {
    m_LayersPalette->getOriginalImageCheckBox()->setChecked(true);
    m_LayersPalette->getSegmentedImageCheckBox()->setChecked(true);
    m_LayersPalette->getCompositeTypeComboBox()->setCurrentIndex(EmMpm_Constants::Alpha_Blend);
    m_LayersPalette->getOpacitySlider()->setEnabled(true);
    m_LayersPalette->getOpacitySpinBox()->setEnabled(true);
    m_LayersPalette->getCompositeTypeComboBox()->setEnabled(true);
  }
#endif
  setWidgetListEnabled(true);
  setImageWidgetsEnabled(true);

  updateBaseRecentFileList(processedImage);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::overlayImageFileLoaded(const QString &filename)
{
  // std::cout << "TomoGui::overlayImageFileLoaded" << std::endl;
  reconstructedVolumeFileName->blockSignals(true);
  reconstructedVolumeFileName->setText(filename);
  reconstructedVolumeFileName->blockSignals(false);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::readMRCHeader(QString filepath)
{

  MRCHeader header;
  ::memset(&header, 0, sizeof(header));

  MRCReader::Pointer reader = MRCReader::New(true);
  // Read the header from the file
  int err = reader->readHeader(filepath.toStdString(), &header);
  if(err < 0)
  {
    return;
  }
  int tiltIndex = 0;
  // Transfer the meta data from the MRC Header to the GUI
  m_XDim->setText(QString::number(header.nx));
  m_YDim->setText(QString::number(header.ny));
  m_nTilts->setText(QString::number(header.nz));
  m_XGrid->setText(QString::number(header.mx));
  m_YGrid->setText(QString::number(header.my));
  m_ZGrid->setText(QString::number(header.mz));
  m_XCell->setText(QString::number(header.xlen));
  m_YCell->setText(QString::number(header.ylen));
  m_ZCell->setText(QString::number(header.zlen));
  m_XOrigin->setText(QString::number(header.xorg));
  m_YOrigin->setText(QString::number(header.yorg));
  m_ZOrigin->setText(QString::number(header.zorg));
  // Set the current tilt index range and value
  currentTiltIndex->setRange(0, header.nz - 1);
  currentTiltIndex->setValue(tiltIndex);

  xMin->setText("0");
  yMin->setText("0");
  xMax->setText(QString::number(header.nx-1));
  yMax->setText(QString::number(header.ny-1));

  // If we have the FEI headers get that information
  if(header.feiHeaders != NULL)
  {
    FEIHeader fei = header.feiHeaders[tiltIndex];
    a_tilt->setText(QString::number(fei.a_tilt));
    b_tilt->setText(QString::number(fei.b_tilt));
    x_stage->setText(QString::number(fei.x_stage));
    y_stage->setText(QString::number(fei.y_stage));
    z_stage->setText(QString::number(fei.z_stage));
    x_shift->setText(QString::number(fei.x_shift));
    y_shift->setText(QString::number(fei.y_shift));
    defocus->setText(QString::number(fei.defocus));
    exp_time->setText(QString::number(fei.exp_time));
    mean_int->setText(QString::number(fei.mean_int));
    tiltaxis->setText(QString::number(fei.tiltaxis));
    pixelsize->setText(QString::number(fei.pixelsize));
    magnification->setText(QString::number(fei.magnification));
    voltage->setText(QString::number(fei.voltage));
    QVector<int> indices(header.nz);
    QVector<float> a_tilts(header.nz);
    QVector<float> b_tilts(header.nz);

    QVector<bool>  excludes(header.nz);
    for(int l = 0; l < header.nz; ++l)
    {
      indices[l] = l;
      a_tilts[l] = header.feiHeaders[l].a_tilt;
      b_tilts[l] = header.feiHeaders[l].b_tilt;

      excludes[l] = false;
    }
    if (NULL != m_GainsOffsetsTableModel)
    {
      m_GainsOffsetsTableModel->setTableData(indices, a_tilts, b_tilts, excludes);
    }
  }
  else
  {
    statusBar()->showMessage("FEI Header information was not found in the file and is needed.");
    return;
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QImage TomoGui::loadMRCTiltImage(QString filepath, int tiltIndex)
{

  currentTiltIndex->setValue(tiltIndex);
  QImage image;
  {
    MRCHeader header;
    MRCReader::Pointer reader = MRCReader::New(true);
    // Read the header from the file
    int err = reader->readHeader(filepath.toStdString(), &header);
    if(err < 0)
    {
      FREE_FEI_HEADERS( header.feiHeaders )
      return image;
    }

    // If we have the FEI headers get that information
    if(header.feiHeaders != NULL)
    {
      FEIHeader fei = header.feiHeaders[tiltIndex];
      a_tilt->setText(QString::number(fei.a_tilt));
      b_tilt->setText(QString::number(fei.b_tilt));
      x_stage->setText(QString::number(fei.x_stage));
      y_stage->setText(QString::number(fei.y_stage));
      z_stage->setText(QString::number(fei.z_stage));
      x_shift->setText(QString::number(fei.x_shift));
      y_shift->setText(QString::number(fei.y_shift));
      defocus->setText(QString::number(fei.defocus));
      exp_time->setText(QString::number(fei.exp_time));
      mean_int->setText(QString::number(fei.mean_int));
      tiltaxis->setText(QString::number(fei.tiltaxis));
      pixelsize->setText(QString::number(fei.pixelsize));
      magnification->setText(QString::number(fei.magnification));
      voltage->setText(QString::number(fei.voltage));
    }
    // Read the first image from the file
    int voxelMin[3] =
    { 0, 0, tiltIndex };
    int voxelMax[3] =
    { header.nx - 1, header.ny - 1, tiltIndex };

    err = reader->read(filepath.toStdString(), voxelMin, voxelMax);
    m_CurrentCorner = 0; // Reset the corner

    if(err >= 0)
    {
      switch(header.mode)
      {
        case 0:
          break;
        case 1:
          image = signed16Image(reinterpret_cast<qint16*>(reader->getDataPointer()), header);
          break;
        case 2:
          break;
        default:
          break;
      }
    }
    FREE_FEI_HEADERS( header.feiHeaders )
  }
  return image;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::on_currentTiltIndex_valueChanged(int i)
{
  m_CurrentImage = loadMRCTiltImage(this->windowFilePath(), i);
  // Be sure to properly orient the image which will in turn load the image into
  // the graphics scene
  on_originCB_currentIndexChanged(originCB->currentIndex());
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::drawOrigin(QImage image)
{
  int imageWidth = image.width();
  int imageHeight = image.height();

  int pxHigh = 0;
  int pxWide = 0;

  QFont font("Ariel", 16, QFont::Bold);
  {
    QPainter painter;
    QImage pImage(100, 100, QImage::Format_ARGB32_Premultiplied);
    pImage.fill(0xFFFFFFFF); // All white background
    painter.begin(&pImage);

    painter.setFont(font);
    QFontMetrics metrics = painter.fontMetrics();
    pxHigh = metrics.height();
    pxWide = metrics.width(QString("TD"));
    painter.end();
  }

  int pxOffset = 2 * pxWide;
  int pyOffset = 2 * pxHigh;
  // Get a QPainter object to add some more details to the image

  int pImageWidth = imageWidth + pxOffset * 2;
  int pImageHeight = imageHeight + pyOffset * 2;

  QImage pImage(pImageWidth, pImageHeight, QImage::Format_ARGB32_Premultiplied);
  pImage.fill(0xFFFFFFFF); // All white background

  // Create a Painter backed by a QImage to draw into
  QPainter painter;
  painter.begin(&pImage);
  painter.setRenderHint(QPainter::Antialiasing, true);

  painter.setFont(font);
  QFontMetrics metrics = painter.fontMetrics();
  pxHigh = metrics.height();
  pxWide = metrics.width(QString("TD"));

  QPoint point(pxOffset, pyOffset);
  painter.drawImage(point, image); // Draw the image we just generated into the QPainter's canvas

  qint32 penWidth = 2;
  painter.setPen(QPen(QColor(0, 0, 0, 255), penWidth, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));


  pxWide = metrics.width(QString("(0,0)"));
  painter.drawText(0, pImageHeight - pyOffset + pxHigh + 2, "(0,0)");


  // Draw slightly transparent lines
  painter.setPen(QPen(QColor(0, 0, 0, 180), penWidth, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
//  painter.drawLine(pImageWidth / 2, pImageHeight / 2, pImageWidth - pxOffset, pImageHeight / 2);
//  painter.drawLine(pImageWidth / 2, pImageHeight / 2, pImageWidth / 2, pImageHeight - pyOffset);

  painter.end();

  m_CurrentImage = pImage;

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QImage TomoGui::signed16Image(qint16* data, MRCHeader &header)
{
  qint16 dmax = std::numeric_limits<qint16>::min();
  qint16 dmin = std::numeric_limits<qint16>::max();
  size_t size = header.nx * header.ny;
  for (size_t i = 0; i < size; ++i)
  {
    if(data[i] > dmax) dmax = data[i];
    if(data[i] < dmin) dmin = data[i];
  }

  // Generate a Color Table
  float max = static_cast<float>(dmax);
  float min = static_cast<float>(dmin);
  int numColors = static_cast<int>((max - min) + 1);

  // Only generate the color table if the number of colors does not match
  if(m_ColorTable.size() != numColors)
  {
    m_ColorTable.resize(numColors);
    float range = max - min;

    float r, g, b;
    for (int i = 0; i < numColors; i++)
    {
      int16_t val = static_cast<int16_t>(min + ((float)i / numColors) * range);
      getColorCorrespondingTovalue(val, r, g, b, max, min);
      m_ColorTable[i] = qRgba(r * 255, g * 255, b * 255, 255);
    }
  }


  // Create an RGB Image
  QImage image(header.nx, header.ny, QImage::Format_ARGB32);


  int idx = 0;
  for (int y = 0; y < header.ny; ++y)
  {
    for (int x = 0; x < header.nx; ++x)
    {
      idx = (header.nx * y) + x;
      int colorIndex = data[idx] - static_cast<int>(dmin);
      image.setPixel(x, y, m_ColorTable[colorIndex]);
    }
  }
  return image;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::on_originCB_currentIndexChanged(int corner)
{

//  int corner = 0;
//  if (m_UpperLeftBtn->isChecked()) { corner = 0; }
//  if (m_UpperRightBtn->isChecked()) { corner = 1; }
//  if (m_LowerRightBtn->isChecked()) { corner = 2; }
//  if (m_LowerLeftBtn->isChecked()) { corner = 3; }

  switch(m_CurrentCorner)
  {

    case 0:
      if (corner == 1) {m_CurrentImage = m_CurrentImage.mirrored(true, false);}
      if (corner == 2) {m_CurrentImage = m_CurrentImage.mirrored(true, true);}
      if (corner == 3) {m_CurrentImage = m_CurrentImage.mirrored(false, true);}
      break;
    case 1:
      if (corner == 0) {m_CurrentImage = m_CurrentImage.mirrored(true, false);}
      if (corner == 2) {m_CurrentImage = m_CurrentImage.mirrored(false, true);}
      if (corner == 3) {m_CurrentImage = m_CurrentImage.mirrored(true, true);}
      break;
    case 2:
      if (corner == 0) {m_CurrentImage = m_CurrentImage.mirrored(true, true);}
      if (corner == 1) {m_CurrentImage = m_CurrentImage.mirrored(false, true);}
      if (corner == 3) {m_CurrentImage = m_CurrentImage.mirrored(true, false);}
      break;
    case 3:
      if (corner == 0) {m_CurrentImage = m_CurrentImage.mirrored(false, true);}
      if (corner == 1) {m_CurrentImage = m_CurrentImage.mirrored(true, true);}
      if (corner == 2) {m_CurrentImage = m_CurrentImage.mirrored(true, false);}
      break;
    default:
      break;
  }
  m_CurrentCorner = corner;


  // This will display the image in the graphics scene
  m_GraphicsView->loadBaseImageFile(m_CurrentImage);
}

///getColorCorrespondingToValue ////////////////////////////////////////////////
//
// Assumes you've already generated min and max -- the extrema for the data
// to which you're applying the color map. Then define the number of colorNodes
// and make sure there's a row of three float values (representing r, g, and b
// in a 0.0-1.0 range) for each node. Then call this method for with parameter
// val some float value between min and max inclusive. The corresponding rgb
// values will be returned in the reference-to-float parameters r, g, and b.
//
////////////////////////////////////////////////////////////////////////////////
void TomoGui::getColorCorrespondingTovalue(int16_t val,
                                   float &r, float &g, float &b,
                                   float max, float min)
{
#if GRAY_SCALE
  static const int numColorNodes = 2;
  float color[numColorNodes][3] =
  {
        {0.0f, 0.0f, 0.0f},    // black
        {1.0f, 1.0f, 1.0f}     // white
  };
#else
  static const int numColorNodes = 4;
  float color[numColorNodes][3] =
  {
        {0.25f, 0.2549f, 0.7961f},    // blue
        {0.8274f, 0.8039f, 0.0941f},    // yellow
        {0.1803f, 0.6655f, 0.1490f},    // Green
        {1.0f, 0.0f, 0.0f}     // red
  };
#endif
  float range = max - min;
  for (int i = 0; i < (numColorNodes - 1); i++)
  {
    float currFloor = min + ((float)i / (numColorNodes - 1)) * range;
    float currCeil = min + ((float)(i + 1) / (numColorNodes - 1)) * range;

    if((val >= currFloor) && (val <= currCeil))
    {
      float currFraction = (val - currFloor) / (currCeil - currFloor);
      r = color[i][0] * (1.0f - currFraction) + color[i + 1][0] * currFraction;
      g = color[i][1] * (1.0f - currFraction) + color[i + 1][1] * currFraction;
      b = color[i][2] * (1.0f - currFraction) + color[i + 1][2] * currFraction;
    }
  }
}

#if 0
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::on_importGainsOffsetsBtn_clicked()
{
  QString file =
      QFileDialog::getOpenFileName(this, tr("Select Gains Offsets Data file"), m_OpenDialogLastDirectory, tr("Binary Files (*.bin)"));

  if (true == file.isEmpty())
  {
    return;
  }
  readGainsOffsetsFile (file);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::readGainsOffsetsFile(QString file)
{
  m_GainsFile = file;
  if (verifyPathExists(m_GainsFile, NULL) == false)
  {
    return;
  }
  GainsOffsetsReader::Pointer reader = GainsOffsetsReader::New();
  TomoInputsPtr inputs = TomoInputsPtr(new TomoInputs);
  SOCEngine::InitializeTomoInputs(inputs);
  inputs->fileZSize = m_GainsOffsetsTableModel->rowCount();
  std::vector<int> goodViews(inputs->fileZSize);
  for(int i = 0; i < m_GainsOffsetsTableModel->rowCount(); ++i)
  {
    goodViews[i] = i;
  }
  inputs->goodViews = goodViews;
  inputs->GainsFile = m_GainsFile.toStdString();
  SinogramPtr sinogram = SinogramPtr(new Sinogram);
  SOCEngine::InitializeSinogram(sinogram);
  sinogram->N_theta = m_GainsOffsetsTableModel->rowCount();
  reader->setTomoInputs(inputs);
  reader->setSinogram(sinogram);
  reader->execute();
  if (reader->getErrorCondition() < 0)
  {
    m_GainsFile = ""; // Clear out the variable so we don't try to load it again.
    return;
  }

  RealArrayType::Pointer gains = sinogram->InitialGain;
  RealArrayType::Pointer offsets = sinogram->InitialOffset;
  RealArrayType::Pointer variances = sinogram->InitialVariance;

 // int nDims = gains->getNDims();
  size_t* dims = gains->getDims();
  size_t total = dims[0];
  double* gainsPtr = gains->getPointer();
  double* offsetsPtr = offsets->getPointer();
  double* variancesPtr = variances->getPointer();

  QVector<double> fGains(total);
  QVector<double> fOffsets(total);
  QVector<double> fVariances(total);
  for(size_t i = 0; i < total; ++i)
  {
    fGains[i] = gainsPtr[i];
    fOffsets[i] = offsetsPtr[i];
    fVariances[i] = variancesPtr[i];
  }

  if (m_GainsOffsetsTableModel != NULL)
  {
    m_GainsOffsetsTableModel->setGainsAndOffsets(fGains, fOffsets, fVariances);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::on_exportGainsOffsets_clicked()
{
  QString proposedFile = m_OpenDialogLastDirectory + QDir::separator() + "GainsOffsets.bin";
  QString file = QFileDialog::getSaveFileName(this, tr("Save Gains/Offsets Data"),
                                              proposedFile,
                                              tr("*.bin") );
  if ( true == file.isEmpty() ){ return;  }
  QFileInfo fi(file);
  m_OpenDialogLastDirectory = fi.absolutePath();

  QVector<double> gains;
  QVector<double> offsets;
  QVector<double> variances;
  // Convert to 8 byte double
  m_GainsOffsetsTableModel->getGainsAndOffsets(gains, offsets, variances);

  FILE* f = fopen(file.toStdString().c_str(), "wb");
  if (NULL == f)
  {
    QMessageBox::critical(this, tr("Gains Offsets File Creation"),
                          tr("The Gains & Offsets file could not be created."), QMessageBox::Ok);
    return;
  }
  qint32 written = fwrite(gains.data(), sizeof(double), gains.size(), f);
  if (written != gains.size())
  {
    QMessageBox::critical(this, tr("Gains Offsets File Writing"),
                          tr("The number of gains values did not match the number of rows in the table."), QMessageBox::Ok);
    fclose(f);return;
  }
  written = fwrite(offsets.data(), sizeof(double), offsets.size(), f);
  if (written != gains.size())
  {
    QMessageBox::critical(this, tr("Gains Offsets File Writing"),
                          tr("The number of offset values did not match the number of rows in the table."), QMessageBox::Ok);
    fclose(f);return;
  }
  fclose(f);
}
#endif

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::on_playBtn_clicked()
{
  if (playBtn->text().compare(QString(">") ) == 0)
  {
    playBtn->setText(QString("||") );
 //   qint32 currentIndex = framesPerSecComboBox->currentIndex();
    double rate = 500;
    double update = 1.0/rate * 1000.0;
    this->m_StopAnimation = false;
    m_AnimationTimer->setSingleShot(true);
    m_AnimationTimer->start(static_cast<int>(update) );
  }
  else
  {
    playBtn->setText(">");
    this->m_StopAnimation = true;
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::on_skipEnd_clicked()
{
  currentTiltIndex->setValue(currentTiltIndex->maximum());
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::on_skipStart_clicked()
{
  currentTiltIndex->setValue(currentTiltIndex->minimum());
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::stepForwardFromTimer()
{
  //Stop Playing if the user clicked the play button again
  if (m_StopAnimation)
  {
    this->m_AnimationTimer->stop();
    playBtn->setText(QString(">") );
    return;
  }
  QCoreApplication::processEvents();

  int idx = currentTiltIndex->value();
  if (idx == currentTiltIndex->maximum())
  {
    currentTiltIndex->setValue(0);
  }
  idx = currentTiltIndex->value();

  if (idx < currentTiltIndex->maximum())
  {
    currentTiltIndex->setValue(idx += 1); // This should cause a loading of the image
  }
  else
  {
    m_StopAnimation = true;
  }

  //   qint32 currentIndex = framesPerSecComboBox->currentIndex();
  double rate = 500;
  double update = 1.0/rate * 1000.0;
  m_AnimationTimer->setSingleShot(true);
  m_AnimationTimer->start(static_cast<int>(update) );
  QCoreApplication::processEvents();
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::displayDialogBox(QString title, QString text, QMessageBox::Icon icon)
{

  QMessageBox msgBox;
  msgBox.setText(title);
  msgBox.setInformativeText(text);
  msgBox.setStandardButtons(QMessageBox::Ok);
  msgBox.setDefaultButton(QMessageBox::Ok);
  msgBox.setIcon(icon);
  msgBox.exec();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::on_estimateGainSigma_clicked()
{
 // std::cout << "on_estimateGainSigma_clicked" << std::endl;
  bool ok = false;
  TargetGainSigmaXEstimation::Pointer estimate = TargetGainSigmaXEstimation::New();
  estimate->setInputFile(inputMRCFilePath->text().toStdString());
  estimate->setSampleThickness(sampleThickness->text().toDouble(&ok));
  estimate->setDefaultOffset(defaultOffset->text().toDouble(&ok));
  estimate->setTiltAngles(tiltSelection->currentIndex());
  estimate->addObserver(this);
  estimate->execute();
  this->progressBar->setValue(0);
  targetGain->setText(QString::number(estimate->getTargetGainEstimate()));
  sigmaX->setText(QString::number(estimate->getSigmaXEstimate()));
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::updateProgressAndMessage(const char* message, int progress)
{
  this->progressBar->setValue(progress);
  if (NULL != this->statusBar()) {
    this->statusBar()->showMessage(message);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::updateProgressAndMessage(const std::string &msg, int progress)
{
  this->progressBar->setValue(progress);
  if (NULL != this->statusBar()) {
    this->statusBar()->showMessage(msg.c_str());
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::deleteUserInitArea(ReconstructionArea* recon)
{

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::userInitAreaAdded(ReconstructionArea* recon)
{
  unsigned int xmin, ymin;
  recon->getUpperLeft(xmin, ymin);

  unsigned int xmax, ymax;
  recon->getLowerRight(xmax, ymax);

  xMin->setText(QString::number(xmin));
  xMax->setText(QString::number(xmax));

  yMin->setText(QString::number(ymin));
  yMax->setText(QString::number(ymax));
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::userInitAreaUpdated(ReconstructionArea* recon)
{

  unsigned int xmin, ymin;
  recon->getUpperLeft(xmin, ymin);

  unsigned int xmax, ymax;
  recon->getLowerRight(xmax, ymax);

  xMin->setText(QString::number(xmin));
  xMax->setText(QString::number(xmax));

  yMin->setText(QString::number(ymin));
  yMax->setText(QString::number(ymax));

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::userInitAreaSelected(ReconstructionArea* recon)
{

}



#if 0
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::on_addResolution_clicked()
{
  int count = tomoInputHorzLayout->count() - 1;

  // The user is adding a resolution
  TomoInputWidget* tiw = new TomoInputWidget(this);
  m_TomoInputs.push_back(tiw);
  tomoInputHorzLayout->insertWidget(count, tiw);

  for(int i = 0; i < m_TomoInputs.size(); ++i)
  {
    TomoInputWidget* tiw = qobject_cast<TomoInputWidget*>(m_TomoInputs.at(i));
    if(NULL != tiw)
    {
      tiw->setResolutionMultiple(pow(2, m_TomoInputs.size() - (1+i) ) );
      tiw->setIndexLabel(i);
    }
  }

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::on_removeResolution_clicked()
{
  // int count = tomoInputHorzLayout->count() -1;
  if(m_TomoInputs.count() == 1)
  {
    return;
  }
  // The user is removing a resolution
  QWidget* tiw = m_TomoInputs.at(m_TomoInputs.count() - 1);
  m_TomoInputs.pop_back(); // Remove it from our internal storage
  tomoInputHorzLayout->removeWidget(tiw);
  tiw->deleteLater(); // Delete it later

  for(int i = 0; i < m_TomoInputs.size(); ++i)
  {
    TomoInputWidget* tiw = qobject_cast<TomoInputWidget*>(m_TomoInputs.at(i));
    if(NULL != tiw)
    {
      tiw->setResolutionMultiple(pow(2, m_TomoInputs.size() - (1+i) ) );
      tiw->setIndexLabel(i);
    }
  }

}
#endif

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

void TomoGui::memCalculate()
{
	/*

    float GeomNx,GeomNy,GeomNz;
    float SinoNr,SinoNt,SinoNtheta;
    SinoNr = inputs->xEnd - inputs->xStart+1;
    SinoNt = inputs->yEnd - inputs->yStart+1;
    SinoNtheta = inputs->zEnd - inputs->zStart+1;
    
	AdvancedParametersPtr advancedParams = AdvancedParametersPtr(new AdvancedParameters);
    SOCEngine::InitializeAdvancedParams(advancedParams);
	
	//std::cout<<"Advaced params"<<advancedParams->Z_STRETCH<<std::endl;
	
	
    if(inputs->extendObject == 1)
    {
       		
		      float LengthZ = m_SampleThickness*advancedParams->Z_STRETCH;
				float temp = advancedParams->X_SHRINK_FACTOR * ((SinoN_r * sinogram->delta_r) / cos(maxTilt * M_PI / 180)) + input->LengthZ * tan(max * M_PI / 180);
				temp/= (input->interpolateFactor * sinogram->delta_r);
				float GeomLengthX = floor(temp + 0.5) * inputs->interpolateFactor * sinogram->delta_r;
				GeomN_x = floor(GeomLengthX / inputs->delta_xz);
    }
    else 
	{
        GeomNx = SinoNr/m_FinalResolution;
    }
    
    GeomNy = SinoNt/m_FinalResolution;
    GeomNz = advancedParams->Z_STRETCH*(m_SampleThickness/(m_FinalResolution*sinogram->delta_r));// TODO: need to access Sinogram_deltar and z_stretch. 
	//This is wrong currently. Need to multiply m_FinalResolution by size of voxel in nm 
    
    float dataTypeMem = sizeof(Real_t);
    float ObjectMem = GeomNx*GeomNy*GeomNz*dataTypeMem;
    float SinogramMem = SinoNr*SinoNt*SinoNtheta*dataTypeMem;
    float ErroSinoMem = SinogramMem;
    float WeightMem = SinogramMem; //Weight matrix
    float A_MatrixMem;
    if(0 == inputs->extendObject)
    {
		A_MatrixMem = GeomNx*GeomNz*(m_FinalResolution*3*(dataTypeMem+4)*SinoNtheta);// 4 is the bytes to store the counts 
		//*+4 correspodns to bytes to store a single double and a unsigned into to
		//store the offset. 3*m_FinalRes is the approximate number of detector elements hit per voxel
    }
    else {
        A_MatrixMem = GeomNx*GeomNz*(m_FinalResolution*(dataTypeMem+4)*SinoNtheta); //Since we are reconstructing a larger region there are several voxels with no projection data. so instead of each voxel hitting 3*m_FinalRes det entries we aproximate it by m_FinalRes
    }
    float NuisanceParamMem = SinoNtheta*dataTypeMem*3;//3 is for gains offsets and noise var
	
	if(bf is present) //how to check if file is loaded
		SinogramMem*=2;
	
    float TotalMem = ObjectMem+SinogramMem+ErroSinoMem+WeightMem+A_MatrixMem+NuisanceParamMem;//in bytes
    
    TotalMem/=(1e9);//To get answer in Gb
    
    std::cout<<"Total Max Mem needed ="<<TotalMem<<std::endl;
    
    */
}