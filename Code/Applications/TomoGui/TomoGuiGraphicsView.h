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


class TomoGuiGraphicsView : public QGraphicsView
{
    Q_OBJECT;

  public:

    TomoGuiGraphicsView( QWidget *parent = NULL);



    void setTomoGui(TomoGui* gui);

    void setAddUserArea(bool b);

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

    //void setBaseImage(QImage image);
    QImage getBaseImage();

    //void setOverlayImage(QImage image);
    QImage getOverlayImage();

    QImage getCompositedImage();

    void loadBaseImageFile(QImage image);
    void loadOverlayImageFile(const QString &filename);

    QImage& blend(QImage& src, QImage& dst, float opacity);

    void addNewInitArea(ReconstructionArea* userInitArea);
    void createNewUserInitArea(const QRectF brect);


    QLineF getXZPlane();

  public slots:
    void zoomIn();

    void zoomOut();

    void setOverlayTransparency(float f);

    void setZoomIndex(int index);

    void setImageDisplayType(TomoGui_Constants::ImageDisplayType displayType);

    void setCompositeMode(TomoGui_Constants::CompositeType mode);

    void addUserInitArea(bool b);

    void setOverlayImage(QImage image);

    void updateDisplay();

    void userInitAreaUpdated(ReconstructionArea* uia);

  signals:
   void fireBaseMRCFileLoaded();
   void fireOverlayImageFileLoaded(const QString &filename);

   void fireUserInitAreaAdded(ReconstructionArea* uia);
   void fireUserInitAreaLostFocus();
  protected:


  private:
   QGraphicsItem* m_ImageGraphicsItem;
   QImage         m_BaseImage;
   QImage         m_OverlayImage;
   QImage         m_CompositedImage;

   bool           m_AddUserInitArea;
   QRubberBand*   m_RubberBand;
   QPoint         m_MouseClickOrigin;
   float          m_ZoomFactors[10];
   QGraphicsPolygonItem m_XZLine;

   TomoGui*      m_MainGui;
   TomoGui_Constants::ImageDisplayType  m_ImageDisplayType;
   bool                                 m_ShowOverlayImage;
   QPainter::CompositionMode            m_composition_mode;
   float                                m_OverlayTransparency;

   ReconstructionArea*                  m_ReconstructionArea;

   TomoGuiGraphicsView(const TomoGuiGraphicsView&); // Copy Constructor Not Implemented
   void operator=(const TomoGuiGraphicsView&); // Operator '=' Not Implemented
};

#endif /* _TomoGui_GRAPHICS_VIEW_H_ */
