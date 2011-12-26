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

#include <QtGui/QApplication>
#include <QtGui/QFileDialog>
#include <QtGui/QCloseEvent>
#include <QtGui/QMessageBox>
#include <QtGui/QListWidget>
#include <QtGui/QStringListModel>
#include <QtGui/QLineEdit>
#include <QtGui/QDoubleValidator>

// Our Project wide includes
#include "QtSupport/ApplicationAboutBoxDialog.h"
#include "QtSupport/QRecentFileList.h"
#include "QtSupport/QFileCompleter.h"
#include "QtSupport/ImageGraphicsDelegate.h"
#include "QtSupport/ProcessQueueController.h"
#include "QtSupport/ProcessQueueDialog.h"


//
#include "License/LicenseFiles.h"
#include "TomoEngine/TomoEngineVersion.h"
#include "TomoEngineTask.h"
#include "LayersDockWidget.h"

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
m_QueueController(NULL),
m_LayersPalette(NULL),
#if defined(Q_WS_WIN)
m_OpenDialogLastDirectory("C:\\")
#else
m_OpenDialogLastDirectory("~/")
#endif
{
  setupUi(this);
  setupGui();

#if defined (Q_OS_MAC)
  QSettings prefs(QSettings::NativeFormat, QSettings::UserScope, QCoreApplication::organizationDomain(), QCoreApplication::applicationName());
#else
  QSettings prefs(QSettings::IniFormat, QSettings::UserScope, QCoreApplication::organizationDomain(), QCoreApplication::applicationName());
#endif
  readIOSettings(prefs);
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
    writeIOSettings(prefs);
    writeSettings(prefs);
    writeWindowSettings(prefs);
    event->accept();
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::writeIOSettings(QSettings &prefs)
{
  prefs.beginGroup("Input_Output");
  WRITE_STRING_SETTING(prefs, inputMRCFilePath);
  WRITE_STRING_SETTING(prefs, outputMRCFilePath);

  prefs.endGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::readIOSettings(QSettings &prefs)
{
  prefs.beginGroup("Input_Output");
  READ_STRING_SETTING(prefs, outputMRCFilePath, "");


  prefs.endGroup();

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


  prefs.endGroup();
}

// -----------------------------------------------------------------------------
//  Write our prefs to file
// -----------------------------------------------------------------------------
void TomoGui::writeSettings(QSettings &prefs)
{
  prefs.beginGroup("Parameters");

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
  QString proposedFile = m_OpenDialogLastDirectory + QDir::separator() + "EMMPM-Config.txt";
  QString file = QFileDialog::getSaveFileName(this, tr("Save EM/MPM Configuration"),
                                              proposedFile,
                                              tr("*.txt") );
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
                                                 tr("Configuration File (*.txt)") );
  if ( true == file.isEmpty() ){return;  }
  QFileInfo fi(file);
  m_OpenDialogLastDirectory = fi.absolutePath();
  QSettings prefs(file, QSettings::IniFormat, this);
  readSettings(prefs);
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
  mySize.setHeight( mySize.height() -30);
  resize(mySize);
#endif

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

  connect (m_GraphicsView, SIGNAL(fireOverlayImageFileLoaded(const QString &)),
           this, SLOT(overlayImageFileLoaded(const QString &)), Qt::QueuedConnection);

  connect (zoomIn, SIGNAL(clicked()),
           m_GraphicsView, SLOT(zoomIn()), Qt::QueuedConnection);
  connect(zoomOut, SIGNAL(clicked()),
          m_GraphicsView, SLOT(zoomOut()), Qt::QueuedConnection);

  QFileCompleter* com = new QFileCompleter(this, false);
  inputMRCFilePath->setCompleter(com);
  QObject::connect(com, SIGNAL(activated(const QString &)), this, SLOT(on_inputMRCFilePath_textChanged(const QString &)));

  QFileCompleter* com4 = new QFileCompleter(this, false);
  outputMRCFilePath->setCompleter(com4);
  QObject::connect(com4, SIGNAL(activated(const QString &)), this, SLOT(on_outputMRCFile_textChanged(const QString &)));


//  m_QueueDialog->setVisible(false);
  cancelBtn->setVisible(false);


  // setup the Widget List
 // m_WidgetList << m_NumClasses << m_EmIterations << m_MpmIterations << m_Beta << m_MinVariance;

  setWidgetListEnabled(false);

  m_ImageWidgets << zoomIn << zoomOut << fitToWindow << layersPalette;
  setImageWidgetsEnabled(false);

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
void TomoGui::on_processBtn_clicked()
{

  {
    QFileInfo fi(inputMRCFilePath->text());
    if (fi.exists() == false)
    {
      QMessageBox::critical(this, tr("Fixed Image File Error"), tr("Input Image does not exist. Please check the path."), QMessageBox::Ok);
      return;
    }

    if (outputMRCFilePath->text().isEmpty() == true)
    {
      QMessageBox::critical(this, tr("Output Image File Error"), tr("Please select a file name for the segmented image to be saved as."), QMessageBox::Ok);
      return;
    }
    QFile file(outputMRCFilePath->text());
    if (file.exists() == true)
    {
      int ret = QMessageBox::warning(this, tr("EIM Tomo GUI"),
                                     tr("The Output File Already Exists\nDo you want to over write the existing file?"),
                                     QMessageBox::No | QMessageBox::Default,
                                     QMessageBox::Yes,
                                     QMessageBox::Cancel);
      if (ret == QMessageBox::Cancel)
      {
        return;
      }
      else if (ret == QMessageBox::Yes)
      {
        setOutputExistsCheck(true);
      }
      else
      {
        QString outputFile = getOpenDialogLastDirectory() + QDir::separator() + "Untitled.tif";
        outputFile = QFileDialog::getSaveFileName(this, tr("Save Output File As ..."), outputFile, tr("TIF (*.tif)"));
        if (!outputFile.isNull())
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

    fi = QFileInfo(outputMRCFilePath->text());
    QDir outputDir(fi.absolutePath());
    if (outputDir.exists() == false)
    {
      bool ok = outputDir.mkpath(".");
      if (ok == false)
      {
        QMessageBox::critical(this, tr("Output Directory Creation"), tr("The output directory could not be created."), QMessageBox::Ok);
        return;
      }
    }
  }


  m_QueueDialog->clearTable();
  if (getQueueController() != NULL)
  {
    getQueueController()->deleteLater();
  }
  ProcessQueueController* queueController = new ProcessQueueController(this);
  setQueueController(queueController);
 // bool ok;

  InputOutputFilePairList filepairs;

#if 0
  m_LayersPalette->getSegmentedImageCheckBox()->setEnabled(true);
  m_LayersPalette->getSegmentedImageCheckBox()->setChecked(true);
  m_LayersPalette->getCompositeTypeComboBox()->setEnabled(true);
  m_LayersPalette->getUseColorTable()->setEnabled(enableUserDefinedAreas->isChecked());
#endif

  QString inputFile = inputMRCFilePath->text();
  QString outputFile = outputMRCFilePath->text();
  TomoEngineTask* task = newTomoEngineTask(inputFile, outputFile, queueController);

  queueController->addTask(static_cast<QThread*> (task));
  connect(cancelBtn, SIGNAL(clicked()), task, SLOT(cancel()));

  connect(task, SIGNAL(progressTextChanged(QString )), this, SLOT(processingMessage(QString )), Qt::QueuedConnection);
  connect(task, SIGNAL(updateImageAvailable(QImage)), m_GraphicsView, SLOT(setOverlayImage(QImage)));
  connect(task, SIGNAL(histogramsAboutToBeUpdated()), this, SLOT(clearProcessHistograms()));
  connect(task, SIGNAL(updateHistogramAvailable(QVector<double>)), this, SLOT(addProcessHistogram(QVector<double>)));
  this->addProcess(task);

  setInputOutputFilePairList(filepairs);

  // When the event loop of the controller starts it will signal the ProcessQueue to run
  connect(queueController, SIGNAL(started()), queueController, SLOT(processTask()));
  // When the QueueController finishes it will signal the QueueController to 'quit', thus stopping the thread
  connect(queueController, SIGNAL(finished()), this, SLOT(queueControllerFinished()));

  connect(queueController, SIGNAL(started()), this, SLOT(processingStarted()));
  connect(queueController, SIGNAL(finished()), this, SLOT(processingFinished()));

//  getQueueDialog()->setParent(this);
//  m_QueueDialog->setVisible(true);
  processBtn->setVisible(false);
  cancelBtn->setVisible(true);

  setWidgetListEnabled(false);
  setImageWidgetsEnabled(true);

  queueController->start();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
TomoEngineTask* TomoGui::newTomoEngineTask( QString inputFile, QString outputFile, ProcessQueueController* queueController)
{
  bool ok = false;
  TomoEngineTask* task = new TomoEngineTask(NULL);
  //FIXME: Setup all the input/output parameters
  return task;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::addProcess(TomoEngineTask* task)
{
  m_QueueDialog->addProcess(task);
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::on_cancelBtn_clicked()
{
  std::cout << "on_cancelBtn_clicked" << std::endl;
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::processingStarted()
{
//  std::cout << "TomoGui::processingStarted()" << std::endl;
  processBtn->setText("Cancel");
  processBtn->setVisible(false);
  this->statusBar()->showMessage("Processing Images...");
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::processingFinished()
{
//  std::cout << "IPHelper::processingFinished()" << std::endl;
  /* Code that cleans up anything from the processing */
  processBtn->setText("Segment");
  processBtn->setVisible(true);
  cancelBtn->setVisible(false);
//  this->statusBar()->showMessage("Processing Complete");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::processingMessage(QString str)
{
  this->statusBar()->showMessage(str);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::queueControllerFinished()
{

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

    setCurrentImageFile (inputMRCFilePath->text());
    setCurrentProcessedFile(outputMRCFilePath->text());
    m_GraphicsView->loadOverlayImageFile(outputMRCFilePath->text());
//    m_LayersPalette->getSegmentedImageCheckBox()->setChecked(true);

  setWindowTitle(m_CurrentImageFile);
  setWidgetListEnabled(true);

  getQueueController()->deleteLater();
  setQueueController(NULL);

  // Make sure the image manipulating widgets are enabled
  setImageWidgetsEnabled(true);
}



// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::on_outputDirectoryBtn_clicked()
{
  bool canWrite = false;
  QString aDir = QFileDialog::getExistingDirectory(this, tr("Select Output Directory"), getOpenDialogLastDirectory(),
                                            QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

  setOpenDialogLastDirectory(aDir);
  if (!getOpenDialogLastDirectory().isNull())
  {
    QFileInfo fi(aDir);
    canWrite = fi.isWritable();
    if (canWrite) {
      this->outputDirectoryPath->setText(getOpenDialogLastDirectory());
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
void TomoGui::on_inputMRCFilePathBtn_clicked()
{
  //std::cout << "on_actionOpen_triggered" << std::endl;
  QString imageFile =
      QFileDialog::getOpenFileName(this, tr("Select Fixed Image"), getOpenDialogLastDirectory(), tr("Images (*.tif *.tiff *.bmp *.jpg *.jpeg *.png)"));

  if (true == imageFile.isEmpty())
  {
    return;
  }
  inputMRCFilePath->setText(imageFile);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::on_outputMRCButton_clicked()
{
  QString outputFile = getOpenDialogLastDirectory() + QDir::separator() + "Untitled.tif";
  outputFile = QFileDialog::getSaveFileName(this, tr("Save Output File As ..."), outputFile, tr("Images (*.tif *.tiff *.bmp *.jpg *.jpeg *.png)"));
  if (outputFile.isEmpty())
  {
    return;
  }

  QFileInfo fi(outputFile);
  QFileInfo fi2(fi.absolutePath());
  if (fi2.isWritable() == true) {
    outputMRCFilePath->setText(outputFile);
  }
  else
  {
    QMessageBox::critical(this, tr("Output Image File Error"),
                          tr("The parent directory of the output image is not writable by your user. Please make it writeable by changing the permissions or select another directory"),
                          QMessageBox::Ok);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::on_inputMRCFilePath_textChanged(const QString & text)
{
  if (verifyPathExists(inputMRCFilePath->text(), inputMRCFilePath))
  {
    openBaseImageFile(text);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::on_outputMRCFilePath_textChanged(const QString & text)
{
  //  verifyPathExists(outputMRCFilePath->text(), movingImageFile);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::on_outputDirectory_textChanged(const QString & text)
{
  verifyPathExists(outputDirectoryPath->text(), outputDirectoryPath);
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
  if (false == fileinfo.exists())
  {
    lineEdit->setStyleSheet("border: 1px solid red;");
  }
  else
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
    openBaseImageFile( file );
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::on_actionOpenMRCFile_triggered()
{
  //std::cout << "on_actionOpen_triggered" << std::endl;
  QString imageFile = QFileDialog::getOpenFileName(this, tr("Open Image File"),
    m_OpenDialogLastDirectory,
    tr("Images (*.tif *.tiff *.bmp *.jpg *.jpeg *.png)") );

  if ( true == imageFile.isEmpty() )
  {
    return;
  }
  QFileInfo fi(imageFile);
  m_OpenDialogLastDirectory = fi.absolutePath();
  openBaseImageFile(imageFile);
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
void TomoGui::openBaseImageFile(QString imageFile)
{
  if ( true == imageFile.isEmpty() ) // User cancelled the operation
  {
    return;
  }

  inputMRCFilePath->blockSignals(true);
  inputMRCFilePath->setText(imageFile);
  inputMRCFilePath->blockSignals(false);

  setWindowTitle(imageFile);
  this->setWindowFilePath(imageFile);

#if 0
  m_LayersPalette->getOriginalImageCheckBox()->setChecked(true);
  m_LayersPalette->getSegmentedImageCheckBox()->setChecked(false);
#endif

  // Tell the RecentFileList to update itself then broadcast those changes.
  QRecentFileList::instance()->addFile(imageFile);
  setWidgetListEnabled(true);
  setImageWidgetsEnabled(true);
  updateBaseRecentFileList(imageFile);

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
  outputMRCFilePath->blockSignals(true);
  outputMRCFilePath->setText(filename);
  outputMRCFilePath->blockSignals(false);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::on_actionParameters_triggered()
{
  ParametersDockWidget->show();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGui::on_actionLayers_Palette_triggered()
{
  m_LayersPalette->show();
}
