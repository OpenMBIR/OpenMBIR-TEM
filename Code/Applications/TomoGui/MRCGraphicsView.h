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

#ifndef _TomoGui_GRAPHICS_VIEW_H_
#define _TomoGui_GRAPHICS_VIEW_H_

#include <iostream>

#include <QtCore/QVector>
#include <QtGui/QGraphicsView>
#include <QtGui/QGraphicsItem>
#include <QtGui/QRubberBand>

class TomoGui;
class ReconstructionArea;

namespace TomoGui_Constants {

enum ImageDisplayType {
  OriginalImage = 0,
  SegmentedImage,
  CompositedImage,
  UnknownDisplayType
};

enum CompositeType
{
  Exclusion,
  Difference,
  Alpha_Blend,
  UnknownCompositeType
};

}

/**
 * @class MRCGraphicsView MRCGraphicsView.h APplications/TomoGui/MRCGraphicsView.h
 * @brief
 * @author Michael A. Jackson for BlueQuartz Software
 * @date March 13, 2012
 * @version 1.0
 */
class MRCGraphicsView : public QGraphicsView
{
    Q_OBJECT;

  public:

    MRCGraphicsView( QWidget *parent = NULL);

    void setWidget(QWidget* gui);

    void setAddUserArea(bool b);

    void disableVOISelection(bool disable);


    /**
    * @brief Over-riding implementation from base class
    * @param event QDragEnterEvent Event fired when dragging enters the QGraphicsView
    */
    void dragEnterEvent(QDragEnterEvent *event);

    /**
    * @brief Over-riding implementation from base class
    * @param event QDropEvent Event fired when object is dropped on QGraphicsView
    */
    void dropEvent(QDropEvent *event);

    /**
    * @brief Over-riding implementation from base class
    * @param event QDragLeaveEvent Event fired when dragging leaves QGraphicsView
    */
    void dragLeaveEvent(QDragLeaveEvent *event);

    void mousePressEvent( QMouseEvent* event );
    void mouseMoveEvent( QMouseEvent* event );
    void mouseReleaseEvent( QMouseEvent* event );


    QImage getBaseImage();

    void loadBaseImageFile(QImage image);

    void addNewInitArea(ReconstructionArea* userInitArea);
    void createNewUserInitArea(const QRectF brect);

    ReconstructionArea* reconstructionArea();

    QLineF getXZPlane();

  public slots:
    void zoomIn();

    void zoomOut();

    void setOverlayTransparency(float f);

    void setZoomIndex(int index);

    void setImageDisplayType(TomoGui_Constants::ImageDisplayType displayType);

    void setCompositeMode(TomoGui_Constants::CompositeType mode);

    void addUserInitArea(bool b);

    void updateDisplay();

  //  void userInitAreaUpdated(ReconstructionArea* uia);

  signals:
    void fireBaseMRCFileLoaded();

    void fireImageFileLoaded(const QString &filename);

    void fireSingleSliceSelected(int y);

    void fireReconstructionVOIAdded(ReconstructionArea* reconVOI);

   /*
   void fireReconstructionVOIUpdated(ReconstructionArea* reconVOI);
   void fireReconstructionVOISelected(ReconstructionArea* reconVOI);

   void fireUserInitAreaAboutToDelete(ReconstructionArea* reconVOI);
    */
     void fireReconstructionVOILostFocus();

  protected:


  private:
   QGraphicsItem* m_ImageGraphicsItem;
   QImage         m_BaseImage;
   bool           m_DisableVOISelection;
   bool           m_AddUserInitArea;
   QRubberBand*   m_RubberBand;
   QPoint         m_MouseClickOrigin;
   float          m_ZoomFactors[10];
   QGraphicsPolygonItem m_XZLine;

   QWidget*       m_MainGui;
   TomoGui_Constants::ImageDisplayType  m_ImageDisplayType;
   bool                                 m_ShowOverlayImage;
   QPainter::CompositionMode            m_composition_mode;
   float                                m_OverlayTransparency;

   ReconstructionArea*                  m_ReconstructionArea;

   MRCGraphicsView(const MRCGraphicsView&); // Copy Constructor Not Implemented
   void operator=(const MRCGraphicsView&); // Operator '=' Not Implemented
};

#endif /* _TomoGui_GRAPHICS_VIEW_H_ */
