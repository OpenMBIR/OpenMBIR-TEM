/* ============================================================================
 * Copyright (c) 2012 Michael A. Jackson (BlueQuartz Software)
 * Copyright (c) 2012 Singanallur Venkatakrishnan (Purdue University)
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
 * Neither the name of Singanallur Venkatakrishnan, Michael A. Jackson, the Pudue
 * Univeristy, BlueQuartz Software nor the names of its contributors may be used
 * to endorse or promote products derived from this software without specific
 * prior written permission.
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
 *
 *  This code was written under United States Air Force Contract number
 *                           FA8650-07-D-5800
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */



#ifndef QMRCDISPLAYWIDGET_H_
#define QMRCDISPLAYWIDGET_H_


//-- Qt Includes
#include <QtCore/QObject>
#include <QtGui/QWidget>


#include "ui_QMRCDisplayWidget.h"
#include "TomoEngine/TomoEngine.h"
#include "TomoEngine/IO/MRCHeader.h"


class QTimer;



/*
 *
 */
class QMRCDisplayWidget : public QWidget, private Ui::QMRCDisplayWidget
{
  Q_OBJECT;

  public:
    QMRCDisplayWidget(QWidget *parent = 0);
    virtual ~QMRCDisplayWidget();

    static QImage xzSigned16CrossSection(qint16* data, size_t nVoxels, int* voxelMin, int* voxelMax);
    static QImage xzFloatCrossSection(float* data, size_t nVoxels, int* voxelMin, int* voxelMax);
    static void getColorCorrespondingTovalue(int16_t val,
                                       float &r, float &g, float &b,
                                       float max, float min);


    MRCGraphicsView* graphicsView();
    QImage currentImage();
  //  int currentCorner();

    void loadMRCFile(QString mrcFilePath);

    void loadMRCTiltImage(QString filepath, int tiltIndex);
    void fireOriginCB_Changed();
    void saveCanvas();
    void loadXZSliceReconstruction(QString reconMRCFilePath);

  protected:
    void setupGui();

    void showWidgets(bool b, QList<QWidget*> &list);
    void setImageWidgetsEnabled(bool b);

    QImage signed16Image(qint16* data, MRCHeader &header);
    void drawOrigin(QImage image);



  protected slots:
    void on_playBtn_clicked();
    void on_skipStart_clicked();
    void on_skipEnd_clicked();
    void stepForwardFromTimer();
    void z10_triggered();
    void z25_triggered();
    void z50_triggered();
    void z100_triggered();
    void z125_triggered();
    void z150_triggered();
    void z200_triggered();
    void z400_triggered();
    void z600_triggered();
    void on_fitToWindow_clicked();

//    void on_originCB_currentIndexChanged(int i);

    void on_currentTiltIndex_valueChanged(int i);



  signals:
    void memoryCalculationNeedsUpdated();


  private:
    bool                  m_StopAnimation;     // Trigger to stop a running animation
    QTimer*               m_AnimationTimer;
    QList<QWidget*>       m_ImageWidgets;
    QList<QWidget*>       m_MovieWidgets;

 //   int                   m_CurrentCorner;
    QImage                m_CurrentImage;
    QString               m_CurrentMRCFilePath;
    QString               m_OpenDialogLastDirectory;

    QMRCDisplayWidget(const QMRCDisplayWidget&); // Copy Constructor Not Implemented
    void operator=(const QMRCDisplayWidget&); // Operator '=' Not Implemented
};

#endif /* QMRCDISPLAYWIDGET_H_ */
