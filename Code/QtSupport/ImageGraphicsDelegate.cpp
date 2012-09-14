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
#include "ImageGraphicsDelegate.h"

//-- STL includes
#include <iostream>

//-- Qt Includes
#include <QtCore/QPoint>
#include <QtGui/QMainWindow>
#include <QtGui/QTableWidget>
#include <QtGui/QPixmap>
#include <QtGui/QGraphicsItem>
#include <QtGui/QGraphicsScene>
#include <QtGui/QGraphicsView>
#include <QtGui/QHeaderView>
#include <QtGui/QGraphicsPixmapItem>

#define ZOOM_INDEX_MAX 9
#define ZOOM_INDEX_MIN 0

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ImageGraphicsDelegate::ImageGraphicsDelegate(QObject* parent) :
  QObject(parent),
  m_MainWindow(NULL),
  m_GraphicsView(NULL),
  m_GraphicsScene(NULL),
  m_CompositeImages(false),
  m_CurrentGraphicsItem(NULL),
  _shouldFitToWindow(false)
{

  _zoomFactors[0] = 0.1;
  _zoomFactors[1] = 0.25;
  _zoomFactors[2] = 0.5;
  _zoomFactors[3] = 1.0;
  _zoomFactors[4] = 1.25;
  _zoomFactors[5] = 1.5;
  _zoomFactors[6] = 2.0;
  _zoomFactors[7] = 4.0;
  _zoomFactors[8] = 6.0;
  _zoomFactors[9] = -1.0;
  _zoomIndex = 3;

  m_composition_mode = QPainter::CompositionMode_Exclusion;
  this->m_CachedImage = QImage();
  this->m_OverlayImage = QImage();
  this->m_CompositedImage = QImage();
  m_DelegateName = "Default ImageGraphicsDelegate";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ImageGraphicsDelegate::~ImageGraphicsDelegate()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImageGraphicsDelegate::resetCaches()
{
  this->m_CachedImage = QImage();
  this->m_OverlayImage = QImage();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImageGraphicsDelegate::displayTextMessage(QString message)
{
  _displayTextMessage(message);
}

#if 0
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImageGraphicsDelegate::increaseZoom()
{
  _shouldFitToWindow = false;
  if (m_CachedImage.isNull() == true)
  {
    return;
  }
  // Find the next scaling factor up from where we are currently
  QSize imageSize = this->m_CachedImage.size();
  int gvWidth = m_GraphicsView->size().width();
  int gvHeight = m_GraphicsView->size().height();
  gvWidth -= 4;
  gvHeight -= 4;
  if (imageSize.width() > imageSize.height())
  {
    for (int i = 0; i < ZOOM_INDEX_MAX - 1; ++i)
    {
      if (_zoomFactor < this->_zoomFactors[i] && i > 0)
      {
        this->_zoomIndex = i;
        this->_zoomFactor = this->_zoomFactors[this->_zoomIndex];
        break;
      }
    }
  }

  for (int i = 0; i < ZOOM_INDEX_MAX - 1; ++i)
  {
    if (_zoomFactor < this->_zoomFactors[i] && i > 0)
    {
      this->_zoomIndex = i;
      this->_zoomFactor = this->_zoomFactors[this->_zoomIndex];
      break;
    }
  }

  updateGraphicsView();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImageGraphicsDelegate::decreaseZoom()
{
  _shouldFitToWindow = false;
  // Find the next scaling factor down
  if (m_CachedImage.isNull() == true)
  {
    return;
  }
  QSize imageSize = this->m_CachedImage.size();
  int gvWidth = m_GraphicsView->size().width();
  int gvHeight = m_GraphicsView->size().height();
  gvWidth -= 4;
  gvHeight -= 4;
  if (imageSize.width() > imageSize.height())
  {
    for (int i = 0; i < ZOOM_INDEX_MAX - 1; ++i)
    {
      if (_zoomFactor < this->_zoomFactors[i] && i > 0)
      {
        this->_zoomIndex = i - 1;
        this->_zoomFactor = this->_zoomFactors[this->_zoomIndex];
        break;
      }
    }
  }

  for (int i = ZOOM_INDEX_MAX - 1; i >= 0; --i)
  {
    if (_zoomFactor > this->_zoomFactors[i] && i > 0)
    {
      this->_zoomIndex = i;
      this->_zoomFactor = this->_zoomFactors[this->_zoomIndex];
      break;
    }
  }
  updateGraphicsView();
}
#endif


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
#if 0
void ImageGraphicsDelegate::fitToWindow(int checkbox_state)
{
  if (m_CachedImage.isNull() == true)
  {
    std::cout << "ImageGraphicsDelegate::fitToWindow  m_CachedImage.isNull() == true" << std::endl;
    return;
  }
  if (checkbox_state == Qt::Checked)
  {
    _shouldFitToWindow = true;
  }
  else
  {
    _shouldFitToWindow = false;
    return;
  }
  //  std::cout << m_DelegateName.toStdString() << " fitToWindow." << std::endl;
  _zoomIndex = ZOOM_INDEX_MAX;
  this->setZoomIndex(_zoomFactors[_zoomIndex]);

  QSize imageSize = this->m_CachedImage.size();
  int gvWidth = m_GraphicsView->size().width();
  int gvHeight = m_GraphicsView->size().height();
  double zfW = (double)(gvWidth) / (double)(imageSize.width());
  double zfH = (double)(gvHeight) / (double)(imageSize.height());
  if (zfW < zfH)
  {
    this->setZoomIndex(zfW);
  }
  else
  {
    this->setZoomIndex(zfH);
  }
}
#endif


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImageGraphicsDelegate::on_parentResized()
{
  if (_shouldFitToWindow == true)
  {
    setZoomIndex(9); // Force a rescaling of the image
  }
  updateGraphicsView();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int ImageGraphicsDelegate::getZoomIndex()
{
  return this->_zoomIndex;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImageGraphicsDelegate::setZoomIndex(int zoomIndex)
{
  this->_zoomIndex = zoomIndex;
  if (this->_zoomIndex == 9)
  {
    _shouldFitToWindow = true;
  }
  else
  {
    _shouldFitToWindow = false;
  }

  m_ScaledCachedImage = _scaleImage();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QImage ImageGraphicsDelegate::_scaleImage()
{
  QSize imageSize = this->m_CachedImage.size();
  double _zoomFactor = _zoomFactors[_zoomIndex];

  if (_shouldFitToWindow)
  {
    QSize imageSize = this->m_CachedImage.size();
    int gvWidth = m_GraphicsView->size().width();
    int gvHeight = m_GraphicsView->size().height();
    double zfW = (double)(gvWidth) / (double)(imageSize.width());
    double zfH = (double)(gvHeight) / (double)(imageSize.height());
    if (zfW < zfH)
    {
      imageSize *= zfW;
    }
    else
    {
      imageSize *= zfH;
    }
    return this->m_CachedImage.scaled(imageSize, Qt::KeepAspectRatio, Qt::SmoothTransformation);
  }
  else if (_zoomIndex != 3)
  {
    imageSize *= _zoomFactor;
    return this->m_CachedImage.scaled(imageSize, Qt::KeepAspectRatio, Qt::SmoothTransformation);
  }
  return this->m_CachedImage;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QImage ImageGraphicsDelegate::_scaleImage(QImage image)
{
  QSize imageSize = image.size();
  double _zoomFactor = _zoomFactors[_zoomIndex];

  if (_shouldFitToWindow)
  {
    QSize imageSize = this->m_CachedImage.size();
    int gvWidth = m_GraphicsView->size().width();
    int gvHeight = m_GraphicsView->size().height();
    double zfW = (double)(gvWidth) / (double)(imageSize.width());
    double zfH = (double)(gvHeight) / (double)(imageSize.height());
    if (zfW < zfH)
    {
      imageSize *= zfW;
    }
    else
    {
      imageSize *= zfH;
    }
    return image.scaled(imageSize, Qt::KeepAspectRatio, Qt::SmoothTransformation);
  }
  else if (_zoomIndex != 3)
  {
    imageSize *= _zoomFactor;
    return image.scaled(imageSize, Qt::KeepAspectRatio, Qt::SmoothTransformation);
  }
  return image;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImageGraphicsDelegate::updateGraphicsView(bool updateGraphicsScene)
{

  if (this->m_CachedImage.isNull() == true)
  {
    return;
  }
  if (NULL != m_CurrentGraphicsItem)
  {
    m_GraphicsScene->removeItem(m_CurrentGraphicsItem); //Remove the image that is displaying
    m_CurrentGraphicsItem->setParentItem(NULL); // Set the parent to NULL
    delete m_CurrentGraphicsItem; // Delete the object
  }

  QImage dataImage;
  if (_zoomIndex != 3) // There was scaling and a scaled image already exists
  {
    dataImage = m_ScaledCachedImage;
  }
  else // 100% view, ie, no scaling just use the existing cached image
  {
    dataImage = m_CachedImage;
  }

  QPixmap imagePixmap;
  if (m_CompositeImages == true && m_OverlayImage.isNull() == false)
  {
    QImage topImage = _scaleImage(m_OverlayImage);
    QPainter painter;
    QImage paintImage(dataImage.size(), QImage::Format_ARGB32_Premultiplied);
    QPoint point(0, 0);
    painter.begin(&paintImage);
    // Draw the fixed Image first
    painter.setPen(Qt::NoPen);
    painter.drawImage(point, topImage);
    // Draw the moving image next
    painter.setCompositionMode(m_composition_mode);
    painter.drawImage(point, dataImage);
    painter.end();
    imagePixmap = QPixmap::fromImage(paintImage);
    m_CompositedImage = paintImage;
  }
  else
  {
    m_CompositedImage = QImage();
    imagePixmap = QPixmap::fromImage(dataImage);
  }

  m_CurrentGraphicsItem = m_GraphicsScene->addPixmap(imagePixmap); // Add the new image into the display
  m_CurrentGraphicsItem->setAcceptDrops(true);
  QRectF rect = m_CurrentGraphicsItem->boundingRect();
  m_GraphicsScene->setSceneRect(rect);
  m_GraphicsView->setScene(m_GraphicsScene);
  m_GraphicsView->centerOn(m_CurrentGraphicsItem);
  if (updateGraphicsScene)
  {
    m_GraphicsScene->update(rect);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImageGraphicsDelegate::_displayTextMessage(QString message)
{
  if (NULL != m_CurrentGraphicsItem)
  {
    m_GraphicsScene->removeItem(m_CurrentGraphicsItem); //Remove the image that is displaying
    m_CurrentGraphicsItem->setParentItem(NULL); // Set the parent to NULL
    delete m_CurrentGraphicsItem; // Delete the object
  }
  m_CurrentGraphicsItem = NULL;
  QGraphicsTextItem* tItem = m_GraphicsScene->addText(message); // Add the new image into the display

  QRectF rect = tItem->boundingRect();
  m_GraphicsScene->setSceneRect(rect);
  m_GraphicsView->setScene(m_GraphicsScene);
  m_GraphicsView->centerOn(tItem);
  m_GraphicsScene->update(rect);
}

