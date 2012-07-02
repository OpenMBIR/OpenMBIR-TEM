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


#include "TomoGuiGraphicsView.h"

#include <QtCore/QFileInfo>
#include <QtCore/QUrl>
#include <QtGui/QWidget>
#include <QDragEnterEvent>
#include <QDragLeaveEvent>
#include <QDropEvent>
#include <QtDebug>
#include <QtGui/QPixmap>
#include <QtGui/QGraphicsPolygonItem>

#include "TomoGui.h"
#include "ReconstructionArea.h"

namespace UIA
{
  const static int Alpha = 155;
}



// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
TomoGuiGraphicsView::TomoGuiGraphicsView(QWidget *parent)
: QGraphicsView(parent),
  m_ImageGraphicsItem(NULL),
  m_ReconstructionArea(NULL)
{
  setAcceptDrops(true);
  setDragMode(RubberBandDrag);
  m_AddUserInitArea = true;

  m_ZoomFactors[0] = 0.1f;
  m_ZoomFactors[1] = 0.25f;
  m_ZoomFactors[2] = 0.5f;
  m_ZoomFactors[3] = 1.0f;
  m_ZoomFactors[4] = 1.250f;
  m_ZoomFactors[5] = 1.500f;
  m_ZoomFactors[6] = 2.000f;
  m_ZoomFactors[7] = 4.000f;
  m_ZoomFactors[8] = 6.000f;
  m_ZoomFactors[9] = -1.0f;
  m_MainGui = NULL;
  m_RubberBand = NULL;
  m_ImageDisplayType = TomoGui_Constants::OriginalImage;
  m_composition_mode = QPainter::CompositionMode_SourceOver;
  m_OverlayTransparency = 1.0f; // Fully opaque

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGuiGraphicsView::setOverlayTransparency(float f)
{
  m_OverlayTransparency = f;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGuiGraphicsView::addUserInitArea(bool b)
{
  m_AddUserInitArea = b;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGuiGraphicsView::zoomIn()
{
  scale(1.1, 1.1);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGuiGraphicsView::zoomOut()
{
  scale(1.0 / 1.1, 1.0 / 1.1);
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGuiGraphicsView::setZoomIndex(int index)
{
  if (index == 9)
  {
    QRectF r = scene()->sceneRect();
    fitInView(r, Qt::KeepAspectRatio);
  }
  else
  {
    QTransform transform;
    transform.scale(m_ZoomFactors[index], m_ZoomFactors[index]);
    setTransform(transform);
  }

}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGuiGraphicsView::userInitAreaUpdated(ReconstructionArea* uia)
{
  updateDisplay();
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGuiGraphicsView::dragEnterEvent(QDragEnterEvent *event)
{
 // qWarning("QFSDroppableGraphicsView::dragEnterEvent(QDragEnterEvent *event)");
  // accept just text/uri-list mime format
  if (event->mimeData()->hasFormat("text/uri-list"))
  {
    event->acceptProposedAction();
  }
  this->setStyleSheet("border: 1px solid green;");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGuiGraphicsView::dragLeaveEvent(QDragLeaveEvent *event)
{
//  qWarning("QFSDroppableGraphicsView::dragLeaveEvent(QDragLeaveEvent *event)");
  this->setStyleSheet("");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGuiGraphicsView::dropEvent(QDropEvent *event)
{
  this->setStyleSheet("");
//  qWarning("QFSDroppableGraphicsView::dropEvent(QDropEvent *event)");
  QList<QUrl> urlList;
  QString fName;
  QFileInfo info;
#if 0
  if (event->mimeData()->hasUrls())
  {
    urlList = event->mimeData()->urls(); // returns list of QUrls
    // if just text was dropped, urlList is empty (size == 0)

    if ( urlList.size() > 0) // if at least one QUrl is present in list
    {
      fName = urlList[0].toLocalFile(); // convert first QUrl to local path
      info.setFile( fName ); // information about file
      QString ext = info.suffix().toLower();
      if (ext.compare("tif") == 0
          || ext.compare("tiff") == 0
          || ext.compare("jpg") == 0
          || ext.compare("jpeg") == 0
          || ext.compare("png") == 0
          || ext.compare("bmp") == 0)
      {
        m_MainGui->on_inputMRCFilePath_textChanged(fName);
      }
    }
  }
#endif
  event->acceptProposedAction();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QImage TomoGuiGraphicsView::getCompositedImage()
{
  return m_CompositedImage;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QImage& TomoGuiGraphicsView::blend(QImage& src, QImage& dst, float opacity)
{
    if (src.width() <= 0 || src.height() <= 0)
        return dst;
    if (dst.width() <= 0 || dst.height() <= 0)
        return dst;

    if (src.width() != dst.width() || src.height() != dst.height()) {
#ifndef NDEBUG
        std::cerr << "WARNING: ImageEffect::blend : src and destination images are not the same size\n";
#endif
        return dst;
    }

    if (opacity < 0.0 || opacity > 1.0) {
#ifndef NDEBUG
        std::cerr << "WARNING: ImageEffect::blend : invalid opacity. Range [0, 1]\n";
#endif
        return dst;
    }

    if (src.depth() != 32) src = src.convertToFormat(QImage::Format_ARGB32);
    if (dst.depth() != 32) dst = dst.convertToFormat(QImage::Format_ARGB32);

    int pixels = src.width() * src.height();
    {
#ifdef WORDS_BIGENDIAN   // ARGB (skip alpha)
        register unsigned char *data1 = (unsigned char *)dst.bits() + 1;
        register unsigned char *data2 = (unsigned char *)src.bits() + 1;
#else                    // BGRA
        register unsigned char *data1 = (unsigned char *)dst.bits();
        register unsigned char *data2 = (unsigned char *)src.bits();
#endif

        for (register int i=0; i<pixels; i++)
        {
#ifdef WORDS_BIGENDIAN
            *data1 += (unsigned char)((*(data2++) - *data1) * opacity);
            data1++;
            *data1 += (unsigned char)((*(data2++) - *data1) * opacity);
            data1++;
            *data1 += (unsigned char)((*(data2++) - *data1) * opacity);
            data1++;
#else
            *data1 += (unsigned char)((*(data2++) - *data1) * opacity);
            data1++;
            *data1 += (unsigned char)((*(data2++) - *data1) * opacity);
            data1++;
            *data1 += (unsigned char)((*(data2++) - *data1) * opacity);
            data1++;
#endif
            data1++; // skip alpha
            data2++;
        }
    }

    return dst;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGuiGraphicsView::updateDisplay()
{
 // QVector<QRgb > colorTable(256);
  if (m_OverlayImage.isNull() == false)
  {

  }
  else
  {
    return;
  }

//  std::cout << "TomoGuiGraphicsView::updateDisplay()" << std::endl;
  QPainter painter;
  QSize pSize(0, 0);
  if (m_BaseImage.isNull() == false)
  {
   pSize = m_BaseImage.size();
  }
  else
  {
    return;
  }

  QImage paintImage(pSize, QImage::Format_ARGB32_Premultiplied);
  QPoint point(0, 0);
  painter.begin(&paintImage);
  painter.setPen(Qt::NoPen);
  if (m_ImageDisplayType == TomoGui_Constants::OriginalImage)
  {
    painter.drawImage(point, m_BaseImage);
  }
  else if (m_ImageDisplayType == TomoGui_Constants::SegmentedImage)
  {
    painter.drawImage(point, m_OverlayImage);
  }
  else if (m_ImageDisplayType == TomoGui_Constants::CompositedImage)
  {

    if (m_composition_mode == QPainter::CompositionMode_SourceOver)
    {
      QImage img = m_OverlayImage.copy(0, 0, m_OverlayImage.width(), m_OverlayImage.height());
      img = blend(m_BaseImage, img, m_OverlayTransparency);
      painter.drawImage(point, img);
    }
    else
    {
      painter.drawImage(point, m_BaseImage);
      if (m_OverlayImage.isNull() == false) {
      painter.setCompositionMode(m_composition_mode);
      painter.drawImage(point, m_OverlayImage);
      }
    }
  }
  painter.end();
  m_CompositedImage = paintImage;

  if (paintImage.isNull() == true)
  {
    return;
  }
  QGraphicsPixmapItem *pixItem = qgraphicsitem_cast<QGraphicsPixmapItem*> (m_ImageGraphicsItem);
  pixItem->setPixmap(QPixmap::fromImage(paintImage));

  this->update();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGuiGraphicsView::setImageDisplayType(TomoGui_Constants::ImageDisplayType displayType)
{
  m_ImageDisplayType = displayType;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGuiGraphicsView::loadBaseImageFile(QImage image)
{
  m_BaseImage = image;
  if (m_BaseImage.isNull() == true)
  {
    return;
  }
  QSize pSize(0, 0);
  m_OverlayImage = m_BaseImage;
  m_CompositedImage = m_BaseImage;

  m_BaseImage = m_BaseImage.convertToFormat(QImage::Format_ARGB32_Premultiplied);

  QGraphicsScene* gScene = scene();
  if (gScene == NULL)
  {
    gScene = new QGraphicsScene(this);
    setScene(gScene);
  }
  else
  {
    gScene->removeItem(m_ImageGraphicsItem);
  }

  if (NULL != m_ReconstructionArea)
  {
    m_ReconstructionArea->setParentItem(NULL);
  }


//  if (NULL == m_ImageGraphicsItem)
  {
    m_ImageGraphicsItem = gScene->addPixmap(QPixmap::fromImage(m_BaseImage));
    if (NULL != m_ReconstructionArea)
        {m_ReconstructionArea->setParentItem(m_ImageGraphicsItem);}
  }


  m_ImageGraphicsItem->setAcceptDrops(true);
  m_ImageGraphicsItem->setZValue(-1);
  QRectF rect = m_ImageGraphicsItem->boundingRect();
  gScene->setSceneRect(rect);

  // Line Color
  m_XZLine.setPen(QPen(QColor(255, 255, 0, UIA::Alpha)));
  // Fill Color
  m_XZLine.setBrush(QBrush(QColor(255, 255, 0, UIA::Alpha)));
  m_XZLine.setZValue(1);
  m_XZLine.setCacheMode(QGraphicsItem::DeviceCoordinateCache);
  m_XZLine.setVisible(m_XZLine.isVisible());
  if (scene()->items().contains(&m_XZLine) == false) {
    scene()->addItem(&m_XZLine);
  }

  this->updateDisplay();
  emit fireBaseMRCFileLoaded();
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGuiGraphicsView::loadOverlayImageFile(const QString &filename)
{

  m_OverlayImage = QImage(filename);
  if (m_OverlayImage.isNull() == true)
  {
    std::cout << "Error Loading Overlay Image: " << filename.toStdString() << std::endl;
    return;
  }

  m_ImageDisplayType = TomoGui_Constants::CompositedImage;
  setOverlayImage(m_OverlayImage);
  emit fireOverlayImageFileLoaded(filename);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGuiGraphicsView::setOverlayImage(QImage image)
{
  m_OverlayImage = image;

  QSize size = m_OverlayImage.size();
 // std::cout << "Overlay Image Size: " << size.width() << " x " << size.height() << std::endl;
  if (size.width() == 0 || size.height() == 0)
  {
    return;
  }

// Convert to an Pre multiplied Image for faster rendering
  m_OverlayImage = m_OverlayImage.convertToFormat(QImage::Format_ARGB32_Premultiplied);

  QGraphicsScene* gScene = scene();
  if (gScene == NULL)
  {
    gScene = new QGraphicsScene(this);
    setScene(gScene);
  }

  // If the GraphicsScene Item does not exist yet lets make one. This would happen
  // if the user loads a segmented image first.
  if (NULL == m_ImageGraphicsItem) {
    m_ImageGraphicsItem = gScene->addPixmap(QPixmap::fromImage(m_OverlayImage));
  }
  m_ImageGraphicsItem->setAcceptDrops(true);
  m_ImageGraphicsItem->setZValue(-1);
  QRectF rect = m_ImageGraphicsItem->boundingRect();
  gScene->setSceneRect(rect);
  updateDisplay();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QImage TomoGuiGraphicsView::getBaseImage()
{
  return m_BaseImage;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
//void TomoGuiGraphicsView::setBaseImage(QImage image)
//{
//  m_BaseImage = image;
//}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QImage TomoGuiGraphicsView::getOverlayImage()
{
  return m_OverlayImage;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGuiGraphicsView::setAddUserArea(bool b)
{
  m_AddUserInitArea = b;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGuiGraphicsView::mousePressEvent(QMouseEvent *event)
{
  // std::cout << "TomoGuiGraphicsView::mousePressEvent accepted:" << (int)(event->isAccepted()) << std::endl;
  m_AddUserInitArea = true;
  QGraphicsView::mousePressEvent(event);

  // std::cout << "    event->accepted() == false" << std::endl;
  if(m_AddUserInitArea == true)
  {
    m_MouseClickOrigin = event->pos();
    if(!m_RubberBand) m_RubberBand = new QRubberBand(QRubberBand::Rectangle, this);
    m_RubberBand->setGeometry(QRect(m_MouseClickOrigin, QSize()));
    m_RubberBand->show();
  }

}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGuiGraphicsView::mouseMoveEvent(QMouseEvent *event)
{
  if(m_AddUserInitArea == true && m_RubberBand != NULL)
  {
    m_RubberBand->setGeometry(QRect(m_MouseClickOrigin, event->pos()).normalized());
  }
  else
  {
    QGraphicsView::mouseMoveEvent(event);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGuiGraphicsView::mouseReleaseEvent(QMouseEvent *event)
{

  if (event->modifiers() == Qt::ShiftModifier)
  {
    m_RubberBand->hide();
    QRectF sr = sceneRect();

    QPointF mappedPoint = mapToScene(QPoint(0, m_MouseClickOrigin.y()));
    QPointF p0(0, mappedPoint.y());
    QPointF p1(sr.width(), mappedPoint.y());

  //  std::cout << "SceneRect: " << sr.x() << ", " << sr.y() << "  " << sr.width() << ", " << sr.height() << std::endl;
    QVector<QPointF> line;
    line.push_back(p0);
    line.push_back(p1);
    QPolygonF polygon(line);

    m_XZLine.setPolygon(polygon);
    m_XZLine.setVisible(true);

    emit fireSingleSliceSelected();
  }
  else if (m_AddUserInitArea == true)
  {
    m_RubberBand->hide();
    QRect box = QRect(m_MouseClickOrigin, event->pos()).normalized();
    QPolygonF sceneBox = mapToScene(box);
    createNewUserInitArea(sceneBox.boundingRect());
    m_AddUserInitArea = false;
  }
  else
  {
    QGraphicsView::mouseReleaseEvent(event);
    if (scene())
    {
      QList<QGraphicsItem *> selected;
      selected = scene()->selectedItems();
      if (selected.count() == 0)
      {
        emit fireUserInitAreaLostFocus();
      }
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QLineF TomoGuiGraphicsView::getXZPlane()
{
 QPointF p0 = m_XZLine.polygon().at(0);
 QPointF p1 = m_XZLine.polygon().at(1);

 QLineF line(p0, p1);
 return line;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGuiGraphicsView::createNewUserInitArea(const QRectF brect)
{
  ReconstructionArea* userInitArea = new ReconstructionArea( brect);
  userInitArea->setGraphicsView(this);
  // Line Color
  userInitArea->setPen(QPen(QColor(225, 225, 225, UIA::Alpha)));
  // Fill Color
  userInitArea->setBrush(QBrush(QColor(28, 28, 200, UIA::Alpha)));
  userInitArea->setZValue(1);
  userInitArea->setCacheMode(QGraphicsItem::DeviceCoordinateCache);

  addNewInitArea(userInitArea);
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGuiGraphicsView::addNewInitArea(ReconstructionArea* userInitArea)
{
 // std::cout << "EMMPMGraphicsView::addNewInitArea()" << std::endl;


  // Set the Parent Item
  userInitArea->setParentItem(m_ImageGraphicsItem);
  // Add it to the vector of UserInitAreas
  if (m_ReconstructionArea != NULL)
  {
    m_ReconstructionArea->deleteLater();
  }
  m_ReconstructionArea = userInitArea;

  // Hook up the signals and slots
  connect (userInitArea, SIGNAL(fireUserInitAreaAboutToDelete(ReconstructionArea*)),
           m_MainGui, SLOT(deleteUserInitArea(ReconstructionArea*)), Qt::DirectConnection);

  connect (userInitArea, SIGNAL (fireUserInitAreaUpdated(ReconstructionArea*)),
           m_MainGui, SLOT(userInitAreaUpdated(ReconstructionArea*)), Qt::QueuedConnection);

  connect (userInitArea, SIGNAL(fireUserInitAreaSelected(ReconstructionArea*)),
           m_MainGui, SLOT(userInitAreaSelected(ReconstructionArea*)), Qt::QueuedConnection);

  connect (userInitArea, SIGNAL (fireUserInitAreaUpdated(ReconstructionArea*)),
           this, SLOT(userInitAreaUpdated(ReconstructionArea*)), Qt::QueuedConnection);


  emit fireUserInitAreaAdded(userInitArea);
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGuiGraphicsView::setTomoGui(TomoGui* gui)
{
  m_MainGui = gui;
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TomoGuiGraphicsView::setCompositeMode(TomoGui_Constants::CompositeType mode)
{
  m_ImageDisplayType = TomoGui_Constants::CompositedImage;
  switch(mode)
  {
    case TomoGui_Constants::Exclusion: m_composition_mode = QPainter::CompositionMode_Exclusion; break;
    case TomoGui_Constants::Difference: m_composition_mode = QPainter::CompositionMode_Difference; break;
    case TomoGui_Constants::Alpha_Blend:
      m_composition_mode = QPainter::CompositionMode_SourceOver;
      break;
#if 0
    case 2: m_composition_mode = QPainter::CompositionMode_Plus; break;
    case 3: m_composition_mode = QPainter::CompositionMode_Multiply; break;
    case 4: m_composition_mode = QPainter::CompositionMode_Screen; break;
    case 5: m_composition_mode = QPainter::CompositionMode_Darken; break;
    case 6: m_composition_mode = QPainter::CompositionMode_Lighten; break;
    case 7: m_composition_mode = QPainter::CompositionMode_ColorDodge; break;
    case 8: m_composition_mode = QPainter::CompositionMode_ColorBurn; break;
    case 9: m_composition_mode = QPainter::CompositionMode_HardLight; break;
    case 10: m_composition_mode = QPainter::CompositionMode_SoftLight; break;

    case 12: m_composition_mode = QPainter::CompositionMode_Destination; break;
    case 13: m_composition_mode = QPainter::CompositionMode_Source; break;
    case 14: m_composition_mode = QPainter::CompositionMode_DestinationOver; break;
    case 15: m_composition_mode = QPainter::CompositionMode_SourceIn; break;
    case 16: m_composition_mode = QPainter::CompositionMode_DestinationIn; break;
    case 17: m_composition_mode = QPainter::CompositionMode_DestinationOut; break;
    case 18: m_composition_mode = QPainter::CompositionMode_SourceAtop; break;
    case 19: m_composition_mode = QPainter::CompositionMode_DestinationAtop; break;
    case 20: m_composition_mode = QPainter::CompositionMode_Overlay; break;
    case 21: m_composition_mode = QPainter::CompositionMode_Clear; break;
#endif
  default:
    m_composition_mode = QPainter::CompositionMode_Exclusion; break;
  }

  this->setImageDisplayType(m_ImageDisplayType);
}
