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

#ifndef EMMPMGUI_H_
#define EMMPMGUI_H_

//-- Qt Includes
#include <QtCore/QObject>
#include <QtCore/QString>
#include <QtCore/QSettings>
#include <QtGui/QCloseEvent>
#include <QtGui/QMainWindow>
#include <QtGui/QWidget>
#include <QtGui/QGraphicsScene>
#include <QtGui/QMessageBox>



class LayersDockWidget;
class GainsOffsetsTableModel;
class ReconstructionArea;
class MRCInfoWidget;

//-- UIC generated Header
#include "ui_TomoGui.h"
#include "QMultiResolutionSOC.h"

#include "TomoEngine/IO/MRCHeader.h"
#include "TomoEngine/Common/Observer.h"


/**
 * @class TomoGui TomoGui.h Code/Applications/TomoGui/TomoGui.h
 * @brief This is the implementation of the Main Window for the Tomo Gui application.
 *
 * @author Michael A. Jackson for BlueQuartz Software
 * @date Dec 26, 2011
 * @version 1.0
 */
class TomoGui :  public QMainWindow, private Ui::TomoGui, public Observer
{
    Q_OBJECT;

  public:
    TomoGui(QWidget *parent = 0);
    virtual ~TomoGui();

    MXA_INSTANCE_PROPERTY(QString, CurrentImageFile)
    MXA_INSTANCE_PROPERTY(QString, CurrentProcessedFile)

    MXA_INSTANCE_PROPERTY(bool, OutputExistsCheck)

    /**
     * @brief Reads the Preferences from the users pref file
     */
    void readSettings(QSettings &prefs);

    /**
     * @brief Writes the preferences to the users pref file
     */
    void writeSettings(QSettings &prefs);

    void readWindowSettings(QSettings &prefs);
    void writeWindowSettings(QSettings &prefs);

    void readMRCHeader(QString filepath);

    QImage loadMRCTiltImage(QString filepath, int tiltIndex);

    /**
     * @brief Loads the Single Slice Reconstruction in the 2nd Graphics View
     * @param filepath The path to the single slice reconstruction
     */
    void loadSingleSliceReconstruction(QString filepath);

   // void openReconstructedMRCFile(QString reconMrcFilePath);

    void readGainsOffsetsFile(QString file);

    bool sanityCheckOutputDirectory(QLineEdit* le, QString msgTitle);

    /**
      * @brief Either prints a message or sends the message to the User Interface
      * @param message The message to print
      * @param progress The progress of the GrainGenerator normalized to a value between 0 and 100
      */
     virtual void updateProgressAndMessage(const char* message, int progress);
     virtual void updateProgressAndMessage(const std::string &msg, int progress);

  public slots:

  // Manual hookup slots to get signals from the graphics view
    void mrcInputFileLoaded(const QString &filename);


    void reconstructionVOIAdded(ReconstructionArea* recon);
    void reconstructionVOIUpdated(ReconstructionArea* recon);
    void reconstructionVOISelected(ReconstructionArea* recon);
    void reconstructionVOIDeleted(ReconstructionArea* recon);


  protected slots:
  //Manual Hookup Menu Actions
  // File Menu
    void on_actionOpenMRCFile_triggered(); // Open a Data File
 //   void on_actionOpenOverlayImage_triggered();
    void on_actionAbout_triggered();
    void on_actionExit_triggered();
    void on_actionSave_Config_File_triggered();
    void on_actionLoad_Config_File_triggered();
    void on_actionSaveCanvas_triggered();
    void on_actionMRC_Info_triggered();



//Window Menu
    void on_actionParameters_triggered();
    void on_actionLayers_Palette_triggered();

    void memCalculate();

    /* slots for the buttons in the GUI */
    void on_m_GoBtn_clicked();

    void on_estimateGainSigma_clicked();

    void on_m_SingleSliceReconstructionBtn_clicked();


    /**
     * @brief Qt Slot that fires in response to a click on a "Recent File' Menu entry.
     */
    void openRecentBaseImageFile();

    /**
     * @brief Updates the QMenu 'Recent Files' with the latest list of files. This
     * should be connected to the Signal QRecentFileList->fileListChanged
     * @param file The newly added file.
     */
    void updateBaseRecentFileList(const QString &file);

    // -----------------------------------------------------------------------------
    //  Input Tab Widgets
    void on_inputMRCFilePathBtn_clicked();
    void on_inputMRCFilePath_textChanged(const QString &string);

    void on_inputBrightFieldFilePathBtn_clicked();
    void on_inputBrightFieldFilePath_textChanged(const QString &string);

    void on_reconstructedVolumeFileNameBtn_clicked();
    void on_reconstructedVolumeFileName_textChanged(const QString & text);

    void on_outputDirectoryPathBtn_clicked();
    void on_tempDirPath_textChanged(const QString & text);

    void on_initialReconstructionPathBtn_clicked();
    void on_initialReconstructionPath_textChanged(const QString & text);

    void singleSlicePlaneSet(int y);

  protected:

    /**
    * @brief Implements the CloseEvent to Quit the application and write settings
    * to the preference file
    */
    void closeEvent(QCloseEvent *event);


    /**
    * @brief Initializes some of the GUI elements with selections or other GUI related items
    */
    void setupGui();

    /**
    * @brief Checks the currently open file for changes that need to be saved
    * @return
    */
    qint32 checkDirtyDocument();

    /**
    * @brief Enables or Disables all the widgets in a list
    * @param b
    */
    void setWidgetListEnabled(bool b);

    /**
    * @brief Verifies that a path exists on the file system.
    * @param outFilePath The file path to check
    * @param lineEdit The QLineEdit object to modify visuals of (Usually by placing a red line around the QLineEdit widget)
    */
    bool verifyPathExists(QString outFilePath, QLineEdit* lineEdit);

    /**
    * @brief Verifies that a parent path exists on the file system.
    * @param outFilePath The parent file path to check
    * @param lineEdit The QLineEdit object to modify visuals of (Usually by placing a red line around the QLineEdit widget)
    */
    bool verifyOutputPathParentExists(QString outFilePath, QLineEdit* lineEdit);

    void initializeSOCEngine(bool fullReconstruction);

    qint32 initImageViews();

    void setImageWidgetsEnabled(bool b);

    void drawOrigin(QImage image);

    void displayDialogBox(QString title, QString text, QMessageBox::Icon icon);

    bool checkTiltAngles(QVector<float> &tilts);

  signals:

    /**
    * @brief Signal emitted when a process is started
    */
    void pipelineStarted();

    /**
    * @brief Signal Emitted when a process has ended.
    */
    void pipelineEnded();

    /**
    * @brief Signal Emitted when the process has been canceled.
    */
    void cancelPipeline();

    /**
    * @brief Signal emitted when a message is available for display to the user
    * @param
    */
    void pipelineWarningMessage(const QString &);
    void pipelineErrorMessage(const QString &);

  private slots:
    // slots for our worker thread to communicate
    virtual void addErrorMessage(QString message);
    virtual void addWarningMessage(QString message);
    virtual void addProgressMessage(QString message);
    virtual void loadProgressMRCFile(QString filePath);

    virtual void pipelineComplete();
    virtual void pipelineProgress(int value);
    virtual void singleSliceComplete();

  private:
    LayersDockWidget*  m_LayersPalette;
    QMap<QObject*, QWidget*> m_TasksMap;

    QList<QWidget*>       m_WidgetList;

    float                 m_CachedLargestAngle;
    float                 m_CachedPixelSize;

    QString               m_GainsFile;
    QString               m_OffsetsFile;
    QString               m_VarianceFile;
    QThread*              m_WorkerThread;
    QMultiResolutionSOC*   m_MultiResSOC;
    GainsOffsetsTableModel*  m_GainsOffsetsTableModel;

    bool                  m_SingleSliceReconstructionActive;
    bool                  m_FullReconstrucionActive;
    MRCInfoWidget*        m_MRCInfoWidget;
    quint16               m_XDim;
    quint16               m_nTilts;
    QString               m_OpenDialogLastDirectory;




    TomoGui(const TomoGui&); // Copy Constructor Not Implemented
    void operator=(const TomoGui&); // Operator '=' Not Implemented
};

#endif /* EMMPMGUI_H_ */
