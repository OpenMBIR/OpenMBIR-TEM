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

#include "BrightFieldGui.h"

#include <errno.h>

#include <iostream>
#include <sstream>
#include <limits>
#include <fstream>

//-- Qt Includes
#include <QtCore/QPluginLoader>
#include <QtCore/QFileInfo>
#include <QtCore/QDir>
#include <QtCore/QString>
#include <QtCore/QUrl>
#include <QtCore/QThread>
#include <QtCore/QThreadPool>
#include <QtCore/QFileInfoList>
#include <QtCore/QTimer>
#include <QtCore/QRect>

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
#include "MBIRLib/MBIRLib.h"
#include "MBIRLib/MBIRLibVersion.h"
#include "MBIRLib/Common/EIMMath.h"
#include "MBIRLib/Reconstruction/ReconstructionStructures.h"
#include "MBIRLib/Reconstruction/ReconstructionEngine.h"
#include "MBIRLib/GenericFilters/SigmaXEstimation.h"
#include "MBIRLib/IOFilters/MRCHeader.h"
#include "MBIRLib/IOFilters/MRCReader.h"

#include "MBIRLib/BrightField/Filters/GainsOffsetsReader.h"

#include "License/LicenseFiles.h"

#include "CheckBoxDelegate.h"
#include "LayersDockWidget.h"
#include "GainsOffsetsTableModel.h"
#include "ReconstructionArea.h"
#include "MRCInfoWidget.h"
#include "ImageOpenDialog.h"
#include "TomogramTiltLoader.h"

/*
 * Concatenate preprocessor tokens A and B without expanding macro definitions
 * (however, if invoked from a macro, macro arguments are expanded).
 */
#define PPCAT_NX(A, B) A ## B

/*
 * Concatenate preprocessor tokens A and B after macro-expanding them.
 */
#define PPCAT(A, B) PPCAT_NX(A, B)

/*
 * Turn A into a string literal without expanding macro definitions
 * (however, if invoked from a macro, macro arguments are expanded).
 */
#define STRINGIZE_NX(A) #A

/*
 * Turn A into a string literal after macro-expanding it.
 */
#define STRINGIZE(A) STRINGIZE_NX(A)

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

#define WRITE_QRECT_SETTING(prefs, var)\
  prefs.setValue(STRINGIZE(PPCAT(var, _X)), var.x());\
  prefs.setValue(STRINGIZE(PPCAT(var, _Y)), var.y());\
  prefs.setValue(STRINGIZE(PPCAT(var, _W)), var.width());\
  prefs.setValue(STRINGIZE(PPCAT(var, _H)), var.height());

#define READ_QRECT_SETTING(prefs, var)\
  var.setX(prefs.value(STRINGIZE(PPCAT(var, _X))).toInt());\
  var.setY(prefs.value(STRINGIZE(PPCAT(var, _Y))).toInt());\
  var.setWidth(prefs.value(STRINGIZE(PPCAT(var, _W))).toInt());\
  var.setHeight(prefs.value(STRINGIZE(PPCAT(var, _H))).toInt());

#define READ_BOOL_SETTING(prefs, var, emptyValue)\
{ QString s = prefs.value(#var).toString();\
  if (s.isEmpty() == false) {\
  bool bb = prefs.value(#var).toBool();\
  var->setChecked(bb); } else { var->setChecked(emptyValue); } }

#define READ_CHECKBOX_SETTING(prefs, var)\
  bool temp = prefs.value(#var).toBool();\
  var->setChecked(temp);\

#define WRITE_BOOL_SETTING(prefs, var, b)\
  prefs.setValue(#var, (b) );

#define WRITE_CHECKBOX_SETTING(prefs, var)\
  prefs.setValue(#var, var->isChecked() );

#define WRITE_VALUE(prefs, var)\
  prefs.setValue(#var, var);

static int tiffCount;


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
BrightFieldGui::BrightFieldGui(QWidget *parent) :
  QMainWindow(parent),
  m_OutputExistsCheck(false),
  m_LayersPalette(NULL),
  m_WorkerThread(NULL),
  m_MultiResSOC(NULL),
  m_SingleSliceReconstructionActive(false),
  m_FullReconstructionActive(false),
  m_UpdateCachedSigmaX(true)
{
  m_OpenDialogLastDirectory = QDir::homePath();
  setupUi(this);
  setupGui();

  QRecentFileList* recentFileList = QRecentFileList::instance();
  connect(recentFileList, SIGNAL (fileListChanged(const QString &)), this, SLOT(updateBaseRecentFileList(const QString &)));
  // Get out initial Recent File List
  this->updateBaseRecentFileList(QString::null);
  qRegisterMetaType<QVector<double> >("QVector<double>");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
BrightFieldGui::~BrightFieldGui()
{

}

// -----------------------------------------------------------------------------
//  Called when the main window is closed.
// -----------------------------------------------------------------------------
void BrightFieldGui::closeEvent(QCloseEvent *event)
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
void BrightFieldGui::readSettings(QSettings &prefs)
{
  QString val;
  bool ok;
  qint32 i;
  double d;
  prefs.beginGroup("Parameters");

  READ_SETTING(prefs, xWidthFullRecon, ok, i, 100, Int);
  READ_STRING_SETTING(prefs, yMin, "0");
  READ_STRING_SETTING(prefs, yMax, "0");

  // This will auto load the MRC File
  READ_STRING_SETTING(prefs, inputMRCFilePath, "");

  // This will auto load the MRC File
  READ_STRING_SETTING(prefs, inputBrightFieldFilePath, "");

  READ_STRING_SETTING(prefs, initialReconstructionPath, "");
  READ_STRING_SETTING(prefs, reconstructedVolumeFileName, "");


  READ_STRING_SETTING(prefs, sampleThickness, "150");
  READ_STRING_SETTING(prefs, targetGain, "1")
      READ_STRING_SETTING(prefs, smoothness, "1.0");
  READ_STRING_SETTING(prefs, sigma_x, "1.0");
  READ_STRING_SETTING(prefs, bf_offset, "0");
  READ_STRING_SETTING(prefs, bragg_threshold, "3");


  READ_SETTING(prefs, numResolutions, ok, i, 1, Int);
  READ_SETTING(prefs, finalResolution, ok, i, 1, Int);
  READ_SETTING(prefs, outerIterations, ok, i, 1, Int);
  READ_SETTING(prefs, innerIterations, ok, i, 1, Int);

  READ_STRING_SETTING(prefs, defaultOffset, "0");
  READ_STRING_SETTING(prefs, stopThreshold, "0.001");
  READ_SETTING(prefs, mrf, ok, d, 1.2, Double);
  READ_BOOL_SETTING(prefs, extendObject, false);
  READ_BOOL_SETTING(prefs, m_DeleteTempFiles, false);
  READ_BOOL_SETTING(prefs, useDefaultOffset, true);

  QRect rect = QRect(0,0,0,0);
  READ_QRECT_SETTING(prefs, rect);
  m_MRCDisplayWidget->graphicsView()->removeBackgroundSelector();
  m_MRCDisplayWidget->graphicsView()->createBackgroundSelector(rect);

  ok = false;
  i = prefs.value("tiltSelection").toInt(&ok);
  if (false == ok) {i = 0;}
  tiltSelection->setCurrentIndex(i);

  prefs.endGroup();
}

// -----------------------------------------------------------------------------
//  Write our prefs to file
// -----------------------------------------------------------------------------
void BrightFieldGui::writeSettings(QSettings &prefs)
{
  prefs.beginGroup("Parameters");
  WRITE_STRING_SETTING(prefs, inputMRCFilePath);
  WRITE_STRING_SETTING(prefs, inputBrightFieldFilePath);

  WRITE_STRING_SETTING(prefs, initialReconstructionPath);
  //WRITE_STRING_SETTING(prefs, outputDirectoryPath);
  WRITE_STRING_SETTING(prefs, reconstructedVolumeFileName);

  WRITE_STRING_SETTING(prefs, sampleThickness);
  WRITE_STRING_SETTING(prefs, targetGain)
      WRITE_STRING_SETTING(prefs, sigma_x);
  WRITE_STRING_SETTING(prefs, smoothness);
  WRITE_SETTING(prefs, numResolutions);
  WRITE_SETTING(prefs, finalResolution);
  WRITE_SETTING(prefs, outerIterations);
  WRITE_SETTING(prefs, innerIterations);

  WRITE_STRING_SETTING(prefs, bragg_threshold);
  WRITE_STRING_SETTING(prefs, defaultOffset);
  WRITE_STRING_SETTING(prefs, stopThreshold);
  WRITE_STRING_SETTING(prefs, bf_offset);

  WRITE_SETTING(prefs, mrf);
  WRITE_CHECKBOX_SETTING(prefs, extendObject);
  WRITE_CHECKBOX_SETTING(prefs, m_DeleteTempFiles)
      WRITE_CHECKBOX_SETTING(prefs, useDefaultOffset);

  //  WRITE_BOOL_SETTING(prefs, useSubVolume, useSubVolume->isChecked());
  WRITE_SETTING(prefs, xWidthFullRecon);
  WRITE_STRING_SETTING(prefs, yMin);
  WRITE_STRING_SETTING(prefs, yMax);
  //  WRITE_STRING_SETTING(prefs, zMin);
  //  WRITE_STRING_SETTING(prefs, zMax);

  RectangleCreator* rectangle = m_MRCDisplayWidget->graphicsView()->getBackgroundRectangle();
  if (NULL != rectangle)
  {
    QRect rect = m_MRCDisplayWidget->graphicsView()->getBackgroundRectangle()->getMappedRectangleCoordinates();
    WRITE_QRECT_SETTING(prefs, rect);
  }

  prefs.setValue("tiltSelection", tiltSelection->currentIndex());

  prefs.endGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::readWindowSettings(QSettings &prefs)
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
void BrightFieldGui::writeWindowSettings(QSettings &prefs)
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
void BrightFieldGui::on_actionSave_Config_File_triggered()
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

  // Tell the RecentFileList to update itself then broadFt those changes.
  QRecentFileList::instance()->addFile(file);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::on_actionLoad_Config_File_triggered()
{
  QString file = QFileDialog::getOpenFileName(this, tr("Select Configuration File"),
                                              m_OpenDialogLastDirectory,
                                              tr("Configuration Files (*.config *.txt)") );
  if ( true == file.isEmpty() )
  {
    return;
  }
  loadConfigurationFile(file);
  QRecentFileList::instance()->addFile(file);
  updateBaseRecentFileList(file);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::loadConfigurationFile(QString file)
{
  QFileInfo fi(file);
  m_OpenDialogLastDirectory = fi.absolutePath();
  QSettings prefs(file, QSettings::IniFormat, this);
  readSettings(prefs);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::on_actionParameters_triggered()
{
  parametersDockWidget->show();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::setupGui()
{
#ifdef Q_WS_MAC
  // Adjust for the size of the menu bar which is at the top of the screen not in the window
  QSize mySize = size();
  mySize.setHeight(mySize.height() - 30);
  resize(mySize);
#endif

  //  xMin = new QLineEdit(QString("0"),this);
  //  xMin->hide();
  //  xMax = new QLineEdit(QString("0"),this);
  //  xMax->hide();

  reconstructedVolumeFileName->setText("");
  m_ReconstructedDisplayWidget->disableVOISelection();


  connect(m_MRCDisplayWidget->graphicsView(), SIGNAL(fireImageFileLoaded(const QString &)),
          this, SLOT(mrcInputFileLoaded(const QString &)), Qt::QueuedConnection);

  connect(m_MRCDisplayWidget->graphicsView(), SIGNAL(fireReconstructionVOIAdded(ReconstructionArea*)),
          this, SLOT(reconstructionVOIAdded(ReconstructionArea*)), Qt::QueuedConnection);

  connect(m_MRCDisplayWidget->graphicsView(), SIGNAL(fireSingleSliceSelected(int)),
          this, SLOT(singleSlicePlaneSet(int)));

  QFileCompleter* com = new QFileCompleter(this, false);
  inputMRCFilePath->setCompleter(com);
  QObject::connect(com, SIGNAL(activated(const QString &)), this, SLOT(on_inputMRCFilePath_textChanged(const QString &)));

  QFileCompleter* com4 = new QFileCompleter(this, false);
  reconstructedVolumeFileName->setCompleter(com4);
  QObject::connect(com4, SIGNAL(activated(const QString &)), this, SLOT(on_reconstructedVolumeFileName_textChanged(const QString &)));

  // setup the Widget List
  m_WidgetList << inputBrightFieldFilePath << inputBrightFieldFilePathBtn;
  m_WidgetList << reconstructedVolumeFileName << reconstructedVolumeFileNameBtn << initialReconstructionPath << initialReconstructionPathBtn;

  setWidgetListEnabled(false);

  connect(m_MRCDisplayWidget, SIGNAL(memoryCalculationNeedsUpdated()),
          this, SLOT(memCalculate()));

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

  QDoubleValidator* dVal = new QDoubleValidator(this);
  dVal->setDecimals(6);
  smoothness->setValidator(dVal);

  QDoubleValidator* dVal2 = new QDoubleValidator(this);
  dVal2->setDecimals(6);
  sigma_x->setValidator(dVal2);

  // ySingleSliceValue_Label->hide();
  // ySingleSliceValue->hide();
  outputDirectoryPath_OLD->hide();
  outputDirectoryPathBtn->hide();
  outputDirectoryLabel->hide();


  outputTabWidget->removeTab(1);
  initialReconstructionPath->hide();
  initialReconstructionLabel->hide();

  inputBrightFieldFilePath->hide();
  inputBrightFieldFilePathBtn->hide();
  label_37->hide();

  m_MRCInputInfoWidget = new MRCInfoWidget(this);
  m_MRCInputInfoWidget->hide();

  m_MRCOutputInfoWidget = new MRCInfoWidget(this);
  m_MRCOutputInfoWidget->hide();

  // Disable all group boxes except for the Background Selection group box
  singleSliceGroupBox->setEnabled(false);
  fullReconstructionGroupBox->setEnabled(false);
  parametersGroupBox->setEnabled(false);
  advancedParametersGroupBox->setEnabled(false);

  // Hide tabs
  m_ReconstructedDisplayWidget->getControlsTab()->hide();
  m_MRCDisplayWidget->getControlsTab()->hide();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
bool BrightFieldGui::sanityCheckOutputDirectory(QString le, QString msgTitle)
{

  if (le.isEmpty() == true)
  {
    QMessageBox::critical(this, msgTitle,
                          "The output directory has NOT been set. Please set a directory path and try again.",
                          QMessageBox::Ok | QMessageBox::Default);
    return false;
  }

  if (verifyPathExists(le, NULL) == false)
  {
    QString msg("The Output Directory '");
    msg.append(le).append("'\ndoes not exist. Would you like to create it?");
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
      QDir outputDir(le);
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
  verifyPathExists(le, NULL);
  return true;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
bool BrightFieldGui::checkTiltAngles(QVector<float> &tilts)
{
  float sum = 0.0;
  for(int i = 0; i < tilts.count(); ++i)
  {
    sum = sum + abs(tilts[i]);
  }
  if (sum > 1.0)
  {
    return true;
  }
  return false;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::on_m_SingleSliceReconstructionBtn_clicked()
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

  // Get everything up and running in a thred, possibly
  startReconstruction(false);
  m_SingleSliceReconstructionBtn->setEnabled(true);
  m_GoBtn->setEnabled(false);
  m_SingleSliceReconstructionBtn->setText("Cancel");
  m_SingleSliceReconstructionActive = true;
  m_FullReconstructionActive = false;
    // When the QThread finishes, tell this object that it has finished.
  connect(m_WorkerThread, SIGNAL(finished()), this, SLOT( singleSliceComplete() ));
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::on_m_GoBtn_clicked()
{
  tiffCount = 0;
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

  if(reconstructedVolumeFileName->text().isEmpty())
  {
    QString str = QString("The field \"Output Reconstruction File\" is empty.\n\nAborting Reconstruction.");
    QMessageBox::critical(this, tr("Reconstruction Error"), str , QMessageBox::Ok);
    return;
  }

  QFileInfo fileInfo(reconstructedVolumeFileName->text());
  if(fileInfo.isRelative())
  {
    QString str = QString("The field \"Output Reconstruction File\" does not contain a full path name.\n\nAborting Reconstruction.");
    QMessageBox::critical(this, tr("Reconstruction Error"), str , QMessageBox::Ok);
    return;
  }

  m_SingleSliceReconstructionBtn->setEnabled(false);



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

  fi = QFileInfo(reconstructedVolumeFileName->text());
  if (fi.suffix().compare("rec") != 0)
  {
    QString t = reconstructedVolumeFileName->text();
    t.append(".rec");
    reconstructedVolumeFileName->blockSignals(true);
    reconstructedVolumeFileName->setText(t);
    reconstructedVolumeFileName->blockSignals(false);
  }

  // Get the absolute path to the users choice for the output file
  fi = QFileInfo(reconstructedVolumeFileName->text());
  // Get the absolute path to the parent directory of that file
  QString parentPath = fi.absolutePath();

  // Make sure we have an output directory setup and created
  if(false == sanityCheckOutputDirectory(parentPath, QString("OpenMBIR Reconstruction")))
  {
    return;
  }

  // We have a name, make sure the user wants to over write the file
  QFile file(fi.absoluteFilePath());
  if(file.exists() == true)
  {
    int ret = QMessageBox::warning(this, tr("OpenMBIR"), tr("The Output File Already Exists\nDo you want to over write the existing file?"), QMessageBox::No
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
      QString outputFile = m_OpenDialogLastDirectory + QDir::separator() + "Untitled.rec";
      outputFile = QFileDialog::getSaveFileName(this, tr("Save Output File As ..."), outputFile, tr("MRC Reconstruction Files (*.rec)"));
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

  //Write out the current settings so if we crash they can easily start back again
#if defined (Q_OS_MAC)
  QSettings prefs(QSettings::NativeFormat, QSettings::UserScope, QCoreApplication::organizationDomain(), QCoreApplication::applicationName());
#else
  QSettings prefs(QSettings::IniFormat, QSettings::UserScope, QCoreApplication::organizationDomain(), QCoreApplication::applicationName());
#endif
  writeSettings(prefs);
  writeWindowSettings(prefs);

  // Sanity Check the Geometry
  {
    // Sanity Check the Input dimensions
    QImage image = m_MRCDisplayWidget->graphicsView()->getBaseImage();
    QSize size = image.size();

    bool ok = false;
    int x_min = 0;
    int x_max = 0;
    m_MRCDisplayWidget->graphicsView()->reconstructionArea()->getXMinMax(x_min, x_max);


    quint16 y_min = yMin->text().toUShort(&ok);
    quint16 y_max = yMax->text().toUShort(&ok);

    if(x_max <= x_min || x_max > size.width())
    {
      return;
    }
    if(y_max <= y_min || y_max > size.height())
    {
      return;
    }

  }
  GainsOffsetsTableModel* tableModel = qobject_cast<GainsOffsetsTableModel*>(gainsOffsetsTableView->model());
  bool goodTilts = false;
  // Sanity Check the Tilt Selection
  if (tiltSelection->currentIndex() == 0)
  {
    QVector<float> tilts = tableModel->getATilts();
    goodTilts = checkTiltAngles(tilts);
  }
  else
  {
    QVector<float> tilts = tableModel->getBTilts();
    goodTilts = checkTiltAngles(tilts);
  }

  if (false == goodTilts)
  {
    QMessageBox::critical(this, tr("Tilt Selection Error"), tr("The Tilt Selection does not seem to be correct. Please check to make sure you have selected the correct tilts."), QMessageBox::Ok);
    return;
  }

  startReconstruction(true);
  m_SingleSliceReconstructionBtn->setEnabled(false);
  m_GoBtn->setEnabled(true);
  m_GoBtn->setText("Cancel");
  m_SingleSliceReconstructionActive = false;
  m_FullReconstructionActive = true;
    // When the QThread finishes, tell this object that it has finished.
  connect(m_WorkerThread, SIGNAL(finished()), this, SLOT( pipelineComplete() ));
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::startReconstruction(bool fullReconstruction)
{
  // Create a new Reconstruction Object
  m_MultiResSOC = new QMultiResolutionReconstruction(NULL);

  initializeSOCEngine(fullReconstruction);
  if(m_MultiResSOC == NULL)
  {
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

  // Move the Reconstruction object into the thread that we just created.
  m_MultiResSOC->moveToThread(m_WorkerThread);

  reconstructedDisplayGroupBox->setTitle(QString("XZ Reconstruction Plane"));

  /* Connect the signal 'started()' from the QThread to the 'run' slot of the
     * SOCEngine object. Since the SOCEngine object has been moved to another
     * thread of execution and the actual QThread lives in *this* thread then the
     * type of connection will be a Queued connection.
     */
  // When the thread starts its event loop, start the Reconstruction going
  connect(m_WorkerThread, SIGNAL(started()), m_MultiResSOC, SLOT(run()));

  // When the Reconstruction ends then tell the QThread to stop its event loop
  connect(m_MultiResSOC, SIGNAL(finished() ), m_WorkerThread, SLOT(quit()));

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

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::initializeSOCEngine(bool fullReconstruction)
{
  QString path;
  path = QDir::toNativeSeparators(inputMRCFilePath->text());
  m_MultiResSOC->setInputFile(path.toStdString());

  QFileInfo fi(reconstructedVolumeFileName->text());
  path = fi.absolutePath();
  if (fullReconstruction == true)
  {
    QString tempFolder = QDir::tempPath() + QDir::separator() + QString("OpenMBIR");
    m_MultiResSOC->setTempDir(tempFolder.toStdString() /*path.toStdString()*/);
    path = QDir::toNativeSeparators(fi.absoluteFilePath());
    m_MultiResSOC->setOutputFile(path.toStdString());
  }
  else
  {
    QString tempFolder = QDir::tempPath() + QDir::separator() + QString("OpenMBIR");
    m_MultiResSOC->setTempDir(tempFolder.toStdString());
    QString reconVolumeFile = tempFolder + QDir::separator() +
        finalResolution->text() + QString("x") + QDir::separator() + QString::fromStdString(ScaleOffsetCorrection::ReconstructedMrcFile);

    m_MultiResSOC->setOutputFile(reconVolumeFile.toStdString());

    QFileInfo fi(reconVolumeFile);
    QDir dir(fi.absolutePath());
    dir.mkpath(".");
  }

  path = QDir::toNativeSeparators(inputBrightFieldFilePath->text());
  m_MultiResSOC->setBrightFieldFile(path.toStdString());

  path = QDir::toNativeSeparators(initialReconstructionPath->text());
  m_MultiResSOC->setInitialReconstructionFile(path.toStdString());

  bool ok = false;
  if (fullReconstruction == true || singleSliceMultiRes->isChecked() == true)
  {
    m_MultiResSOC->setNumberResolutions(numResolutions->value());
  }
  else
  {
    m_MultiResSOC->setNumberResolutions(1);
  }
  m_MultiResSOC->setFinalResolution(finalResolution->value());
  m_MultiResSOC->setSampleThickness(sampleThickness->text().toFloat(&ok));
  m_MultiResSOC->setTargetGain(1);

  m_MultiResSOC->setBraggThreshold(bragg_threshold->text().toFloat(&ok));
  m_MultiResSOC->setBraggDelta(targetGain->text().toFloat(&ok));
  m_MultiResSOC->setBfOffset(bf_offset->text().toFloat(&ok));

  Real_t true_offset = -log(defaultOffset->text().toDouble()+bf_offset->text().toDouble());
  m_MultiResSOC->setDefaultOffsetValue(true_offset);

  m_MultiResSOC->setUseDefaultOffset(useDefaultOffset->isChecked());
  m_MultiResSOC->setDefaultPixelSize(m_CachedPixelSize);

  m_MultiResSOC->setStopThreshold(stopThreshold->text().toFloat(&ok));
  m_MultiResSOC->setOuterIterations(outerIterations->value());
  m_MultiResSOC->setInnerIterations(innerIterations->value());
  m_MultiResSOC->setSigmaX(sigma_x->text().toFloat(&ok));

  m_MultiResSOC->setMRFShapeParameter(mrf->value());
  if(tiltSelection->currentIndex() == 0)
  {
    m_MultiResSOC->setTilts(m_GainsOffsetsTableModel->getATilts().toStdVector());
  }
  else
  {
    m_MultiResSOC->setTilts(m_GainsOffsetsTableModel->getBTilts().toStdVector());
  }

  if(m_MultiResSOC->getTilts().size() == 0)
  {
    addErrorMessage(QString::fromLatin1("No tilt angles were specified. Either no angles were in the MRC file or no angles were imported or generated"));
    delete m_MultiResSOC;
    m_MultiResSOC = NULL;
    return;
  }

  m_MultiResSOC->setExtendObject(extendObject->isChecked());
  m_MultiResSOC->setDefaultVariance(defaultVariance->text().toFloat(&ok));
  m_MultiResSOC->setInitialReconstructionValue(defaultInitialRecon->text().toFloat(&ok));

  m_MultiResSOC->setInterpolateInitialReconstruction(interpolateInitialRecontruction->isChecked());
  m_MultiResSOC->setDeleteTempFiles(m_DeleteTempFiles->isChecked());
  AdvancedParametersPtr advParams = AdvancedParametersPtr(new AdvancedParameters);
  ReconstructionEngine::InitializeAdvancedParams(advParams);
  m_MultiResSOC->setAdvParams(advParams);


  std::vector<uint16_t> subvolume(6);
  subvolume[2] = 0;
  subvolume[5] = m_nTilts - 1;
  if (fullReconstruction == true)
  {
    // Sanity Check the Input dimensions
    QImage image =  m_MRCDisplayWidget->graphicsView()->getBaseImage();
    //QSize size = image.size();

    int x_min = 0;
    int x_max = 0;
    m_MRCDisplayWidget->graphicsView()->reconstructionArea()->getXMinMax(x_min, x_max);

    quint16 y_min = yMin->text().toUShort(&ok);
    quint16 y_max = yMax->text().toUShort(&ok);

    subvolume[0] = x_min;
    subvolume[1] = y_min;// size.height() - y_max - 1;
    subvolume[3] = x_max;
    subvolume[4] = y_max - 1; // size.height() - y_min - 1;

    m_MultiResSOC->setSubvolume(subvolume);
  }
  else
  {
    QLineF line = m_MRCDisplayWidget->graphicsView()->getXZPlane();
    // std::cout << "p1: " << line.p1().x() << ", " << line.p1().y()
    //  << "   p2: " << line.p2().x() << ", " << line.p2().y() << std::endl;

    QImage image =  m_MRCDisplayWidget->graphicsView()->getBaseImage();
    QSize size = image.size();

    // Only reconstruct the middle section of data along the x axis
    float remWidth = m_XDim * (singleSliceXWidth->value()/100.0f)/2.0f;
    float midWidth = m_XDim/2.0f;

    subvolume[0] = midWidth - remWidth;
    subvolume[3] = midWidth + remWidth;

    //    std::cout << subvolume[0] << ", " << subvolume[3] << std::endl;
    // This is how many slices we are going to reconstruct
    int ySlices = 3 * finalResolution->value();

    // The top of our Reconstruction Volume
    quint16 y_min = size.height() - line.p1().y();
    // The bottom of our Reconstruction Volume
    quint16 y_max = y_min + ySlices - 1;

    if (y_max >= size.height())
    {
      y_min = y_min - (y_max - size.height());
      y_max = size.height() - 1;
    }

    subvolume[1] = y_min;
    subvolume[4] = y_max;

    //    path = QDir::toNativeSeparators(outputDirectoryPath->text());
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
void BrightFieldGui::on_singleSliceXWidth_valueChanged(int value)
{
  m_MRCDisplayWidget->graphicsView()->updateXZLine( static_cast<float>(value/100.0));
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::on_xWidthFullRecon_valueChanged(int value)
{
  m_MRCDisplayWidget->graphicsView()->reconstructionArea()->updateWidth(static_cast<float>(value/100.0));
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::singleSlicePlaneSet(int y)
{
  m_SingleSliceReconstructionBtn->setEnabled(true);
  ySingleSliceValue->setText(QString::number(y));

  singleSliceXWidth->setValue(100);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::singleSliceComplete()
{
  std::cout << "BrightFieldGui::singleSliceComplete" << std::endl;
  m_SingleSliceReconstructionBtn->setText("Single Slice Reconstruction");
  m_GoBtn->setEnabled(true);
  setWidgetListEnabled(true);
  this->progressBar->setValue(0);
  QString reconVolumeFile = QString::fromStdString(m_MultiResSOC->getTempDir()) + QDir::separator() +
      finalResolution->text() + QString("x") + QDir::separator() + QString::fromStdString(ScaleOffsetCorrection::ReconstructedMrcFile);

  m_ReconstructedDisplayWidget->loadXZSliceReconstruction(reconVolumeFile);
  m_ReconstructedDisplayWidget->setMovieWidgetsEnabled(false);

  removeDir(QString::fromStdString(m_MultiResSOC->getTempDir()));

  m_ReconstructedDisplayWidget->getControlsTab()->show();

  // There is no way to hide sub-tabs, so we have to remove it instead
  m_ReconstructedDisplayWidget->getControlsTab()->removeTab(1);   // Removes Advanced Controls tab

  m_FullReconstructionActive = false;
  m_SingleSliceReconstructionActive = false;
  emit pipelineEnded();
  m_MultiResSOC->deleteLater();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::loadProgressMRCFile(QString mrcfilePath)
{
  //  std::cout << "Loading Progress MRC File: " << filePath.toStdString() << std::endl;
  m_ReconstructedDisplayWidget->loadXZSliceReconstruction(mrcfilePath);
  m_ReconstructedDisplayWidget->setImageWidgetsEnabled(true);
  m_ReconstructedDisplayWidget->setMovieWidgetsEnabled(false);
  m_ReconstructedDisplayWidget->setDrawOrigin(false);
#if 0
  QFileInfo fi(mrcfilePath);
  QString parentPath = fi.path();
  QDir dir(parentPath);

  QString filepath("/tmp/");
  filepath = filepath.append(QString::number(tiffCount)).append(".tiff");
  m_ReconstructedDisplayWidget->graphicsView()->getBaseImage().save(filepath);
  tiffCount++;
#endif
}
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
bool BrightFieldGui::removeDir(const QString &dirName)
{
  bool result = true;
  QDir dir(dirName);

  if (dir.exists(dirName)) {
    Q_FOREACH(QFileInfo info, dir.entryInfoList(QDir::NoDotAndDotDot | QDir::System | QDir::Hidden  | QDir::AllDirs | QDir::Files, QDir::DirsFirst)) {
      if (info.isDir()) {
        result = removeDir(info.absoluteFilePath());
      }
      else {
        result = QFile::remove(info.absoluteFilePath());
      }

      if (!result) {
        return result;
      }
    }
    result = dir.rmdir(dirName);
  }

  return result;
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::pipelineComplete()
{
  std::cout << "BrightFieldGui::pipelineComplete()" << std::endl;
  m_GoBtn->setText("Reconstruct");
  m_SingleSliceReconstructionBtn->setEnabled(true);
  setWidgetListEnabled(true);
  this->progressBar->setValue(0);
  emit pipelineEnded();

  m_MultiResSOC->deleteLater();
  m_FullReconstructionActive = false;
  m_SingleSliceReconstructionActive = false;

  setCurrentImageFile(inputMRCFilePath->text());
  setCurrentProcessedFile(reconstructedVolumeFileName->text());

  m_ReconstructedDisplayWidget->setMovieWidgetsEnabled(true);
  m_ReconstructedDisplayWidget->setDrawOrigin(false);
  m_ReconstructedDisplayWidget->loadMRCFile(reconstructedVolumeFileName->text());
  m_ReconstructedDisplayWidget->getControlsTab()->show();

  setWindowTitle(m_CurrentImageFile);
  setWidgetListEnabled(true);

  reconstructedDisplayGroupBox->setTitle(QString("Reconstructed Volume"));
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::pipelineProgress(int val)
{
  this->progressBar->setValue( val );
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::addErrorMessage(QString message)
{
  QString title = "BrightFieldGui Error";
  displayDialogBox(title, message, QMessageBox::Critical);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::addWarningMessage(QString message)
{
  QString title = "BrightFieldGui Warning";
  displayDialogBox(title, message, QMessageBox::Warning);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::addProgressMessage(QString message)
{
  if (NULL != this->statusBar()) {
    this->statusBar()->showMessage(message);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::on_outputDirectoryPathBtn_clicked()
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
      this->outputDirectoryPath_OLD->setText(m_OpenDialogLastDirectory);
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
void BrightFieldGui::on_reconstructedVolumeFileNameBtn_clicked()
{
  QString outputFile = m_OpenDialogLastDirectory + QDir::separator() + "Untitled.rec";
  outputFile = QFileDialog::getSaveFileName(this, tr("Save Output File As ..."), outputFile, tr("MRC Files (*.mrc);;REC Files (*.rec)"));
  if (outputFile.isEmpty())
  {
    return;
  }

  QFileInfo fi(outputFile);
  QFileInfo fi2(fi.absolutePath());
  if (fi2.isWritable() == true) {
    reconstructedVolumeFileName->setText(outputFile);
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
void BrightFieldGui::on_inputMRCFilePathBtn_clicked()
{
  //std::cout << "on_actionOpen_triggered" << std::endl;
  QString imageFile =
      QFileDialog::getOpenFileName(this, tr("Select MRC Data file"), m_OpenDialogLastDirectory, tr("MRC Files (*.mrc);;Aligned MRC Files (*.ali);;Reconstructed MRC Files (*.rec);;All Files (*.*)"));

  if (true == imageFile.isEmpty())
  {
    return;
  }
  inputMRCFilePath->setText(imageFile);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::on_inputMRCFilePath_textChanged(const QString & filepath)
{
  QFileInfo fi(inputMRCFilePath->text());
  if(fi.exists() == true && fi.isDir() == false)
  {
    QSize pictureSize;
    // Read the header info from the file and populate the GUI with those values
    {
      ImageOpenDialog d(this);
      d.show();
      d.activateWindow();
      d.setModal(true);
      readMRCHeader(filepath);
      // Now load up the first tilt from the file

      quint16 halfTilts = m_nTilts/2;

      // Load up the middle tilt which we are assuming is the actual "zero" tilt.
      m_MRCDisplayWidget->resetImageScaling();
      m_MRCDisplayWidget->loadMRCTiltImage(filepath, halfTilts);
      m_MRCDisplayWidget->setImageWidgetsEnabled(true);
      m_MRCDisplayWidget->setMovieWidgetsEnabled(true);
      m_MRCDisplayWidget->on_fitToWindow_clicked();

      QSize imageSize = m_MRCDisplayWidget->currentImage().size();
      pictureSize = imageSize;
      QRectF image(0.0, 0.0, imageSize.width(), imageSize.height() );
      m_MRCDisplayWidget->graphicsView()->createNewReconstructionArea(image);
      m_GainsFile = ""; // We are reading a new .mrc file so we probably need a new Gains Offsets File
      smoothness->setText(QString("1.0"));
      on_estimateSigmaX_clicked();
    }

    setWindowTitle(filepath);
    this->setWindowFilePath(filepath);
#if 0
    m_LayersPalette->getOriginalImageCheckBox()->setChecked(true);
    m_LayersPalette->getSegmentedImageCheckBox()->setChecked(false);
#endif

    // Tell the RecentFileList to update itself then broadcast those changes.
    QRecentFileList::instance()->addFile(filepath);
    setWidgetListEnabled(true);

    updateBaseRecentFileList(filepath);

    m_MRCDisplayWidget->graphicsView()->removeBackgroundSelector();
    m_MRCDisplayWidget->graphicsView()->createBackgroundSelector();

    m_MRCDisplayWidget->getControlsTab()->show();
    m_MRCDisplayWidget->getControlsTab()->setCurrentIndex(0);

  }
  else
  {
    sigma_x->setText(0);
    smoothness->setText(0);
    m_CachedSigmaX = 0.0;
    m_MRCDisplayWidget->loadMRCFile(QString(""));
    setWidgetListEnabled(false);
    m_MRCDisplayWidget->getControlsTab()->hide();
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::on_inputBrightFieldFilePathBtn_clicked()
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
void BrightFieldGui::on_inputBrightFieldFilePath_textChanged(const QString & filepath)
{
  if (verifyPathExists(inputBrightFieldFilePath->text(), inputBrightFieldFilePath))
  {
    targetGain->setText(QString("1"));
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::on_reconstructedVolumeFileName_textChanged(const QString & text)
{

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
//void BrightFieldGui::on_outputDirectoryPath_textChanged(const QString & text)
//{
//  verifyPathExists(outputDirectoryPath_OLD->text(), outputDirectoryPath_OLD);
//}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::on_initialReconstructionPathBtn_clicked()
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
void BrightFieldGui::on_initialReconstructionPath_textChanged(const QString & text)
{
  verifyPathExists(initialReconstructionPath->text(), initialReconstructionPath);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::setWidgetListEnabled(bool b)
{
  foreach (QWidget* w, m_WidgetList)
  {
    w->setEnabled(b);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
bool BrightFieldGui::verifyOutputPathParentExists(QString outFilePath, QLineEdit* lineEdit)
{
  QFileInfo fileinfo(outFilePath);
  QDir parent(fileinfo.dir());
  return parent.exists();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
bool BrightFieldGui::verifyPathExists(QString outFilePath, QLineEdit* lineEdit)
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
qint32 BrightFieldGui::checkDirtyDocument()
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
void BrightFieldGui::updateBaseRecentFileList(const QString &file)
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
void BrightFieldGui::openRecentBaseImageFile()
{
  QAction *action = qobject_cast<QAction *>(sender());
  if (action)
  {
    //std::cout << "Opening Recent file: " << action->data().toString().toStdString() << std::endl;
    QString file = action->data().toString();
    QFileInfo fi(file);
    if (fi.suffix().compare("mrc") == 0 || fi.suffix().compare("ali") == 0)
    {
      inputMRCFilePath->setText( file );
    }
    else if (fi.suffix().compare("txt") == 0 || fi.suffix().compare("config") == 0)
    {
      loadConfigurationFile(file);
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::on_actionLayers_Palette_triggered()
{
  m_LayersPalette->show();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::on_actionOpenMRCFile_triggered()
{
  //std::cout << "on_actionOpen_triggered" << std::endl;
  QString file = QFileDialog::getOpenFileName(this, tr("Open MRC File"), m_OpenDialogLastDirectory, tr("MRC Files (*.mrc *.rec *.ali)"));

  if(true == file.isEmpty())
  {
    return;
  }
  QFileInfo fi(file);
  m_OpenDialogLastDirectory = fi.absolutePath();
  inputMRCFilePath->setText(file);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::on_action_OpenReconstructedMRC_triggered()
{
  QString file = QFileDialog::getOpenFileName(this, tr("Open MRC File"), m_OpenDialogLastDirectory, tr("MRC Files (*.mrc *.rec *.ali)"));

  if(true == file.isEmpty())
  {
    return;
  }
  QFileInfo fi(file);
  m_OpenDialogLastDirectory = fi.absolutePath();
  m_ReconstructedDisplayWidget->setMovieWidgetsEnabled(true);
  m_ReconstructedDisplayWidget->setImageWidgetsEnabled(true);
  m_ReconstructedDisplayWidget->setDrawOrigin(false);
  m_ReconstructedDisplayWidget->loadMRCFile(file);
}

#if 0
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::on_actionOpenOverlayImage_triggered()
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
  openReconstructedMRCFile(imageFile);
}
#endif

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::on_actionAbout_triggered()
{
  ApplicationAboutBoxDialog about(OpenMBIR::LicenseList, this);
  QString an = QCoreApplication::applicationName();
  QString version("");
  version.append(MBIRLib::Version::PackageComplete().c_str());
  about.setApplicationInfo(an, version);
  about.exec();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::on_actionExit_triggered()
{
  this->close();
}

#if 0
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::openReconstructedMRCFile(QString reconMrcFilePath)
{
  if ( true == reconMrcFilePath.isEmpty() ) // User cancelled the operation
  {
    return;
  }


  m_ReconstructedDisplayWidget->graphicsView()->loadOverlayImageFile(reconMrcFilePath);


  setWidgetListEnabled(true);
  setImageWidgetsEnabled(true);

  updateBaseRecentFileList(reconMrcFilePath);

}
#endif

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::mrcInputFileLoaded(const QString &filename)
{
  // std::cout << "BrightFieldGui::overlayImageFileLoaded" << std::endl;
  reconstructedVolumeFileName->blockSignals(true);
  reconstructedVolumeFileName->setText(filename);
  reconstructedVolumeFileName->blockSignals(false);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::readMRCHeader(QString filepath)
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
  m_XDim = header.nx;
  m_nTilts = header.nz;
  /*
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
     */

  {
    QIntValidator* val = new QIntValidator(0, header.ny-2, this);
    yMin->setValidator(val);
    yMin->setText("0");
  }
  {
    QIntValidator* val = new QIntValidator(1, header.ny-1, this);
    yMax->setValidator(val);
    yMax->setText(QString::number(header.ny-1));
  }


  // If we have the FEI headers get that information
  if(header.feiHeaders != NULL)
  {
    FEIHeader fei = header.feiHeaders[tiltIndex];
    /*
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
         */
    QVector<int> indices(header.nz);
    QVector<float> a_tilts(header.nz);
    QVector<float> b_tilts(header.nz);

    QVector<bool>  excludes(header.nz);
    m_CachedLargestAngle = std::numeric_limits<float>::min();
    m_CachedPixelSize = fei.pixelsize;
    for(int l = 0; l < header.nz; ++l)
    {
      indices[l] = l;
      a_tilts[l] = header.feiHeaders[l].a_tilt;
      b_tilts[l] = header.feiHeaders[l].b_tilt;

      if (abs(a_tilts[l]) > m_CachedLargestAngle)
      {
        m_CachedLargestAngle =  abs(a_tilts[l]);
      }
      if (abs(b_tilts[l]) > m_CachedLargestAngle)
      {
        m_CachedLargestAngle =  abs(b_tilts[l]);
      }

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
void BrightFieldGui::displayDialogBox(QString title, QString text, QMessageBox::Icon icon)
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
void BrightFieldGui::on_estimateSigmaX_clicked()
{
  //  std::cout << "on_estimateGainSigma_clicked" << std::endl;
  bool ok = false;
  if (sampleThickness->text().isEmpty() == true)
  {
    return;
  }
  if (defaultOffset->text().isEmpty() == true)
  {
    return;
  }
  if (targetGain->text().isEmpty() == true)
  {
    return;
  }

  int xmin = 0;
  int xmax = 0;
  smoothness->setText(QString::number(1.0));
  ReconstructionArea* reconArea = m_MRCDisplayWidget->graphicsView()->reconstructionArea();
  if (NULL == reconArea)
  {
    return;
  }
  reconArea->getXMinMax(xmin, xmax);

  quint16 ymin = yMin->text().toUShort(&ok);
  quint16 ymax = yMax->text().toUShort(&ok);

  if (verifyPathExists(inputMRCFilePath->text(), inputMRCFilePath))
  {
    SigmaXEstimation::Pointer estimate = SigmaXEstimation::New();
    estimate->setInputFile(inputMRCFilePath->text().toStdString());
    estimate->setSampleThickness(sampleThickness->text().toDouble(&ok));
    estimate->setDefaultOffset(defaultOffset->text().toDouble(&ok));
    estimate->setTargetGain(targetGain->text().toDouble(&ok));
    estimate->setBfOffset(bf_offset->text().toDouble());
    //TODO : Set the Offset from the UI
    if(tiltSelection->currentIndex() == 0)
    {
      estimate->setTiltAngles(m_GainsOffsetsTableModel->getATilts().toStdVector());
    }
    else
    {
      estimate->setTiltAngles(m_GainsOffsetsTableModel->getBTilts().toStdVector());
    }


    estimate->setXDims(xmin, xmax);
    estimate->setYDims(ymin, ymax);
    estimate->addObserver(this);
    estimate->execute();
    this->progressBar->setValue(0);

    //targetGain->setText(QString::number(estimate->getTargetGainEstimate()));

    m_CachedSigmaX = estimate->getSigmaXEstimate();
    //std::cout << "m_CachedSigmaX: " << m_CachedSigmaX << std::endl;
    qreal smth = 1.0/smoothness->text().toDouble(&ok);
    //sigma_x->blockSignals(true);
    sigma_x->setText(QString::number(m_CachedSigmaX * smth));
    //sigma_x->blockSignals(false);

    //sigmaX_ShouldUpdate(false);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::on_sigma_x_textChanged(const QString & text)
{
  //  std::cout << "on_sigma_x_textChanged" << std::endl;
  bool ok = false;
  qreal sigx = sigma_x->text().toDouble(&ok);
  qreal smth = m_CachedSigmaX / sigx;
  smoothness->blockSignals(true);
  smoothness->setText(QString::number(smth));
  smoothness->blockSignals(false);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::on_smoothness_textChanged(const QString & text)
{
  //  std::cout << "on_smoothness_textChanged" << std::endl;
  bool ok = false;
  qreal smth = 1.0/smoothness->text().toDouble(&ok);
  sigma_x->blockSignals(true);
  sigma_x->setText(QString::number(m_CachedSigmaX * smth));
  sigma_x->blockSignals(false);
}

#if 0
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::on_targetGain_editingFinished()
{
  sigmaX_ShouldUpdate(true);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::on_sampleThickness_editingFinished()
{
  sigmaX_ShouldUpdate(true);
  memCalculate();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::on_defaultOffset_editingFinished()
{
  sigmaX_ShouldUpdate(true);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::on_tiltSelection_currentIndexChanged(int index)
{
  sigmaX_ShouldUpdate(true);
}
#endif

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::on_resetSigmaXBtn_clicked()
{
  sigma_x->blockSignals(true);
  sigma_x->setText(QString::number(m_CachedSigmaX));
  sigma_x->blockSignals(false);

  smoothness->blockSignals(true);
  smoothness->setText("1.0");
  smoothness->blockSignals(false);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::updateProgressAndMessage(const char* message, int progress)
{
  this->progressBar->setValue(progress);
  if (NULL != this->statusBar()) {
    this->statusBar()->showMessage(message);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::updateProgressAndMessage(const std::string &msg, int progress)
{
  this->progressBar->setValue(progress);
  if (NULL != this->statusBar()) {
    this->statusBar()->showMessage(msg.c_str());
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::reconstructionVOIAdded(ReconstructionArea* reconVOI)
{
  reconstructionVOIUpdated(reconVOI);

  connect (reconVOI, SIGNAL(fireReconstructionVOIAboutToDelete(ReconstructionArea*)),
           this, SLOT(reconstructionVOIDeleted(ReconstructionArea*)), Qt::DirectConnection);

  connect (reconVOI, SIGNAL (fireReconstructionVOIUpdated(ReconstructionArea*)),
           this, SLOT(reconstructionVOIUpdated(ReconstructionArea*)), Qt::QueuedConnection);

  connect (reconVOI, SIGNAL(fireReconstructionVOISelected(ReconstructionArea*)),
           this, SLOT(reconstructionVOISelected(ReconstructionArea*)), Qt::QueuedConnection);


  //  connect(xMin, SIGNAL(textEdited ( const QString &)),
  //          reconVOI, SLOT(setXMin(const QString &)));
  //  connect(xMax, SIGNAL(textEdited ( const QString &)),
  //          reconVOI, SLOT(setXMax(const QString &)));

  connect(yMin, SIGNAL(textEdited ( const QString &)),
          reconVOI, SLOT(setYMax(const QString &)));
  connect(yMax, SIGNAL(textEdited ( const QString &)),
          reconVOI, SLOT(setYMin(const QString &)));


}

#if 1
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::on_yMin_textChanged(const QString &string)
{
  geometryChanged();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::on_yMax_textChanged(const QString &string)
{
  geometryChanged();
}
#endif


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::geometryChanged()
{
  //   std::cout << "BrightFieldGui::geometryChanged()" << std::endl;
  bool ok = false;

  int x_min = 0;
  int x_max = 0;
  ReconstructionArea* reconArea = m_MRCDisplayWidget->graphicsView()->reconstructionArea();
  if (NULL == reconArea)
  {
    return;
  }
  reconArea->getXMinMax(x_min, x_max);

  //    qint32 x_min = xMin->text().toInt(&ok);
  qint32 y_min = yMin->text().toInt(&ok);
  //    qint32 x_max = xMax->text().toInt(&ok);
  qint32 y_max = yMax->text().toInt(&ok);
  if (y_max < y_min)
  {
    yMin->setStyleSheet("border: 1px solid red;");
    yMax->setStyleSheet("border: 1px solid red;");
    statusBar()->showMessage("The Y End Value is Less than the Y Start Value. Please Correct.");
    m_GoBtn->setEnabled(false);
  }
  else
  {
    yMin->setStyleSheet("");
    yMax->setStyleSheet("");
    emit reconstructionVOIGeometryChanged(x_min, y_min, x_max, y_max);
    m_GoBtn->setEnabled(true);
    statusBar()->showMessage("");
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::reconstructionVOIUpdated(ReconstructionArea* recon)
{
  QImage image =  m_MRCDisplayWidget->graphicsView()->getBaseImage();
  QSize size = image.size();


  int xmin, ymin;
  recon->getUpperLeft(xmin, ymin);

  int xmax, ymax;
  recon->getLowerRight(xmax, ymax);

  //    xMin->setText(QString::number(xmin));
  //    xMax->setText(QString::number(xmax - 1));

  yMin->setText(QString::number(size.height() - ymax));
  yMax->setText(QString::number(size.height() - ymin - 1));

  //sigmaX_ShouldUpdate(true);
  memCalculate();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::reconstructionVOISelected(ReconstructionArea* recon)
{

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::reconstructionVOIDeleted(ReconstructionArea* recon)
{

}

#if 0
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::on_addResolution_clicked()
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
void BrightFieldGui::on_removeResolution_clicked()
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
void BrightFieldGui::on_defaultOffsetUpdateBtn_clicked()
{
  QRect rect = m_MRCDisplayWidget->graphicsView()->getBackgroundRectangle()->getMappedRectangleCoordinates();

  int tilt = m_MRCDisplayWidget->getCurrentTiltIndexBox()->value();

  double mean = BackgroundCalculation::getMeanValue(inputMRCFilePath->text().toStdString(), rect.x(), rect.y(), rect.width(), rect.height(), tilt);

  std::stringstream ss;
  ss << mean;

  defaultOffset->setText(QString::fromStdString(ss.str()));

  // Enable all group boxes
  singleSliceGroupBox->setEnabled(true);
  fullReconstructionGroupBox->setEnabled(true);
  parametersGroupBox->setEnabled(true);
  advancedParametersGroupBox->setEnabled(true);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::on_defaultOffset_textChanged(const QString & text)
{
  if (text == "")
  {
    // Disable all group boxes
    singleSliceGroupBox->setEnabled(false);
    fullReconstructionGroupBox->setEnabled(false);
    parametersGroupBox->setEnabled(false);
    advancedParametersGroupBox->setEnabled(false);
  }
  else
  {
    // Enable all group boxes
    singleSliceGroupBox->setEnabled(true);
    fullReconstructionGroupBox->setEnabled(true);
    parametersGroupBox->setEnabled(true);
    advancedParametersGroupBox->setEnabled(true);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::memCalculate()
{
  bool ok = false;
  float GeomN_x, GeomN_y, GeomN_z;

  int x_min = 0;
  int x_max = 0;
  ReconstructionArea* reconArea = m_MRCDisplayWidget->graphicsView()->reconstructionArea();
  if (NULL == reconArea)
  {
    return;
  }
  reconArea->getXMinMax(x_min, x_max);

  float SinoN_r = x_max - x_min + 1;
  float SinoN_t = yMax->text().toInt(&ok) - yMin->text().toInt(&ok) + 1;
  float SinoNtheta = m_nTilts - 0 + 1;

  float sample_thickness = sampleThickness->text().toFloat(&ok);
  int final_resolution = finalResolution->value();
  int num_resolutions = numResolutions->text().toInt(&ok);
  float interpolate_factor = powf((float)2, (float)num_resolutions-1) * final_resolution;
  float delta_r = m_CachedPixelSize * 1.0e9;
  float delta_xz = delta_r*final_resolution;
  AdvancedParametersPtr advancedParams = AdvancedParametersPtr(new AdvancedParameters);
  ReconstructionEngine::InitializeAdvancedParams(advancedParams);

  //std::cout<<"Advaced params"<<advancedParams->Z_STRETCH<<std::endl;

  if(extendObject->isChecked() == true)
  {
    float maxTilt = m_CachedLargestAngle;
    //    float LengthZ = sample_thickness * advancedParams->Z_STRETCH;
    float temp = advancedParams->X_SHRINK_FACTOR * ((SinoN_r * delta_r) / cos(maxTilt * M_PI / 180)) + sample_thickness * tan(maxTilt * M_PI / 180);
    temp /= (interpolate_factor * delta_r);
    float GeomLengthX = floor(temp + 0.5) * interpolate_factor * delta_r;
    GeomN_x = floor(GeomLengthX / delta_xz);
  }
  else
  {
    GeomN_x = SinoN_r / final_resolution;
  }

  GeomN_y = SinoN_t / final_resolution;
  GeomN_z = advancedParams->Z_STRETCH * (sample_thickness / (final_resolution * delta_r)); // TODO: need to access Sinogram_deltar and z_stretch.
  //This is wrong currently. Need to multiply m_FinalResolution by size of voxel in nm

  float dataTypeMem = sizeof(Real_t);
  float ObjectMem = GeomN_x * GeomN_y * GeomN_z * dataTypeMem;
  float SinogramMem = SinoN_r * SinoN_t * SinoNtheta * dataTypeMem;
  float ErroSinoMem = SinogramMem;
  float WeightMem = SinogramMem; //Weight matrix
  float A_MatrixMem;
  if(extendObject->isChecked() == true)
  {
    A_MatrixMem = GeomN_x * GeomN_z * (final_resolution * 3 * (dataTypeMem + 4) * SinoNtheta); // 4 is the bytes to store the counts
    //*+4 correspodns to bytes to store a single double and a unsigned into to
    //store the offset. 3*m_FinalRes is the approximate number of detector elements hit per voxel
  }
  else
  {
    A_MatrixMem = GeomN_x * GeomN_z * (final_resolution * (dataTypeMem + 4) * SinoNtheta); //Since we are reconstructing a larger region there are several voxels with no projection data. so instead of each voxel hitting 3*m_FinalRes det entries we aproximate it by m_FinalRes
  }
  float NuisanceParamMem = SinoNtheta * dataTypeMem * 3; //3 is for gains offsets and noise var

  if(inputBrightFieldFilePath->text().isEmpty() == false) {
    SinogramMem *= 2;
  }

  float TotalMem = ObjectMem + SinogramMem + ErroSinoMem + WeightMem + A_MatrixMem + NuisanceParamMem; //in bytes

  TotalMem /= (1e9); //To get answer in Gb

  memoryUse->setText(QString::number(TotalMem));

}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::on_actionSaveCanvas_triggered()
{
  m_MRCDisplayWidget->saveCanvas();
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::on_action_InputMRCInfo_triggered()
{
  m_MRCInputInfoWidget->setInfo(inputMRCFilePath->text());
  m_MRCInputInfoWidget->show();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::on_action_OutputMRCInfo_triggered()
{
  m_MRCOutputInfoWidget->setInfo(m_ReconstructedDisplayWidget->getMRCFilePath());
  m_MRCOutputInfoWidget->show();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
//void BrightFieldGui::deleteTempFiles()
//{
//  // Remove Any Temp files that have accumulated
//  for(int i = 0; i < m_TempFilesToDelete.size(); ++i)
//  {
//    QFile f(m_TempFilesToDelete.at(i));
//    f.remove();
//  }
//  m_TempFilesToDelete.clear();
//}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void BrightFieldGui::on_actionLoad_Tilt_Information_triggered()
{
  TomogramTiltLoader loader(this);
  loader.setNumTilts(m_nTilts);
  int ret = loader.exec();
  if(ret == QDialog::Accepted) // the user clicked the OK button, now check what they typed
  {
    QVector<float> a_tilts = loader.getATilts();
    QVector<float> b_tilts = loader.getBTilts();
    if (a_tilts.size() != b_tilts.size() )
    {
      if (a_tilts.size() == 0 && b_tilts.size() != 0)
      {
        a_tilts = QVector<float>(b_tilts.size(), 0.0f);
      }
      else if(b_tilts.size() == 0 && a_tilts.size() != 0)
      {
        b_tilts = QVector<float>(a_tilts.size(), 0.0f);
      }

      //      QMessageBox::critical(this, "Tilt Loading Error", "The A tilts and B Tiles do not have the same number of values.",  QMessageBox::Ok, QMessageBox::Ok);
      //      return;
    }
    QVector<int> indices(a_tilts.size());

    QVector<bool> excludes(a_tilts.size());
    m_CachedLargestAngle = std::numeric_limits<float>::min();
    m_CachedPixelSize = loader.getPixelSize();
    for(int l = 0; l < a_tilts.size(); ++l)
    {
      indices[l] = l;


      if (abs(a_tilts[l]) > m_CachedLargestAngle)
      {
        m_CachedLargestAngle =  abs(a_tilts[l]);
      }
      if (abs(b_tilts[l]) > m_CachedLargestAngle)
      {
        m_CachedLargestAngle =  abs(b_tilts[l]);
      }

      excludes[l] = false;
    }
    if (NULL != m_GainsOffsetsTableModel)
    {
      m_GainsOffsetsTableModel->setTableData(indices, a_tilts, b_tilts, excludes);
    }
  }

  // Automatically estimate the Sigma X value by simulating the cuser clicking the button
  on_estimateSigmaX_clicked();
}




