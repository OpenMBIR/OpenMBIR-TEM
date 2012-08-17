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


#include "QMRCDisplayWidget.h"

// STL includes
#include <limits>

//-- Qt Core Includes
#include <QtCore/QTimer>
#include <QtCore/QDir>

//-- Qt GUI Includes
#include <QtGui/QFileDialog>
#include <QtGui/QImage>
#include <QtGui/QMenu>
#include <QtGui/QAction>

//-- TomoEngine Includes
#include "TomoEngine/TomoEngine.h"
#include "TomoEngine/TomoEngineVersion.h"
#include "TomoEngine/IO/MRCHeader.h"
#include "TomoEngine/IO/MRCReader.h"

#define GRAY_SCALE 1

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QMRCDisplayWidget::QMRCDisplayWidget(QWidget *parent) :
QWidget(parent),
m_StopAnimation(true),
m_ImageWidgetsEnabled(false),
m_MovieWidgetsEnabled(false)
{
  m_OpenDialogLastDirectory = QDir::homePath();
  setupUi(this);
  setupGui();

  m_AnimationTimer = new QTimer(this);
  connect(m_AnimationTimer, SIGNAL(timeout() ), this, SLOT(stepForwardFromTimer() ));

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QMRCDisplayWidget::~QMRCDisplayWidget()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
MRCGraphicsView* QMRCDisplayWidget::graphicsView()
{
  return m_GraphicsView;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QImage QMRCDisplayWidget::currentImage()
{
  return m_CurrentImage;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
//int QMRCDisplayWidget::currentCorner()
//{
//  return m_CurrentCorner;
//}


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
void QMRCDisplayWidget::z##var##_triggered() {\
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
void QMRCDisplayWidget::on_fitToWindow_clicked()
{
  m_GraphicsView->setZoomIndex(9);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void QMRCDisplayWidget::disableVOISelection()
{
  if (m_GraphicsView != NULL) { m_GraphicsView->disableVOISelection(true); }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void QMRCDisplayWidget::enableVOISelection()
{
  if (m_GraphicsView != NULL) { m_GraphicsView->disableVOISelection(false); }
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void QMRCDisplayWidget::setupGui()
{
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

  m_GraphicsView->setWidget(this);

  // Just place a really big white image to get our GUI to layout properly
  QImage image(1000, 1000, QImage::Format_ARGB32_Premultiplied);
  image.fill(0);
  m_GraphicsView->loadBaseImageFile(image);

  connect(zoomIn, SIGNAL(clicked()), m_GraphicsView, SLOT(zoomIn()), Qt::QueuedConnection);
  connect(zoomOut, SIGNAL(clicked()), m_GraphicsView, SLOT(zoomOut()), Qt::QueuedConnection);

  m_ImageWidgets << m_SaveCanvasBtn << zoomButton << zoomIn << zoomOut << fitToWindow;

  m_MovieWidgets << indexLabel << currentTiltIndex << skipEnd << skipStart << playBtn << stopBtn;

  showWidgets(m_ImageWidgetsEnabled, m_ImageWidgets);
  showWidgets(m_MovieWidgetsEnabled, m_MovieWidgets);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void QMRCDisplayWidget::on_playBtn_clicked()
{
  playBtn->setEnabled(false);
  stopBtn->setEnabled(true);
  double rate = 500;
  double update = 1.0 / rate * 1000.0;
  this->m_StopAnimation = false;
  m_AnimationTimer->setSingleShot(true);
  m_AnimationTimer->start(static_cast<int>(update));

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void QMRCDisplayWidget::on_stopBtn_clicked()
{
  playBtn->setEnabled(true);
  stopBtn->setEnabled(false);
  this->m_StopAnimation = true;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void QMRCDisplayWidget::on_skipEnd_clicked()
{
  currentTiltIndex->setValue(currentTiltIndex->maximum());
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void QMRCDisplayWidget::on_skipStart_clicked()
{
  currentTiltIndex->setValue(currentTiltIndex->minimum());
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void QMRCDisplayWidget::stepForwardFromTimer()
{
  //Stop Playing if the user clicked the stop button
  if (m_StopAnimation)
  {
    this->m_AnimationTimer->stop();
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
void QMRCDisplayWidget::on_currentTiltIndex_valueChanged(int i)
{
  loadMRCTiltImage(m_CurrentMRCFilePath, i);
  // Be sure to properly orient the image which will in turn load the image into
  // the graphics scene
//  on_originCB_currentIndexChanged(originCB->currentIndex());
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void QMRCDisplayWidget::on_m_SaveCanvasBtn_clicked()
{
  saveCanvas();
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void QMRCDisplayWidget::saveCanvas()
{
  QImage image = m_GraphicsView->getBaseImage();
#if 0
  if (m_LayersPalette->getOriginalImageCheckBox()->isChecked()
       && m_LayersPalette->getSegmentedImageCheckBox()->isChecked())
  {
    image = m_GraphicsView->getCompositedImage();
  }
  int err = 0;
#endif
  int err = 0;
  QString outputFile = this->m_OpenDialogLastDirectory + QDir::separator() + "Untitled.tif";
  outputFile = QFileDialog::getSaveFileName(this, tr("Save Image As ..."), outputFile, tr("Images (*.tif *.bmp *.jpg *.png)"));
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
void QMRCDisplayWidget::loadMRCFile(QString mrcFilePath)
{
    m_CurrentMRCFilePath = mrcFilePath;
    loadMRCTiltImage(m_CurrentMRCFilePath, 0);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void QMRCDisplayWidget::loadMRCTiltImage(QString mrcFilePath, int tiltIndex)
{

  m_CurrentMRCFilePath = mrcFilePath;
  if (m_CurrentMRCFilePath.isEmpty() == true)
  {
    m_GraphicsView->clearContent();
  }

  QImage image;

  MRCHeader header;
  ::memset(&header, 0, sizeof(MRCHeader));
  MRCReader::Pointer reader = MRCReader::New(true);
  // Read the header from the file
  int err = reader->readHeader(mrcFilePath.toStdString(), &header);
  if(err < 0)
  {
    FREE_FEI_HEADERS( header.feiHeaders)
    return;
  }
  currentTiltIndex->setRange(0, header.nz - 1);
  currentTiltIndex->setValue(tiltIndex);

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
  }
  // Read the first image from the file
  int voxelMin[3] =
  { 0, 0, tiltIndex };
  int voxelMax[3] =
  { header.nx - 1, header.ny - 1, tiltIndex };

  err = reader->read(mrcFilePath.toStdString(), voxelMin, voxelMax);

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
        image = floatImage(reinterpret_cast<float*>(reader->getDataPointer()), header);
        break;
      default:
        break;
    }
  }

  FREE_FEI_HEADERS( header.feiHeaders)


  // put the origin in the lower left corner
  image = image.mirrored(false, true);

  drawOrigin(image);
  //m_CurrentImage = image;

  // This will display the image in the graphics scene
  m_GraphicsView->loadBaseImageFile(m_CurrentImage);

  // Calculate the approx memory usage
  emit memoryCalculationNeedsUpdated();

  if (m_ImageWidgetsEnabled == true) { showWidgets(true, m_ImageWidgets); }
  if (m_MovieWidgetsEnabled == true) { showWidgets(true, m_MovieWidgets); }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void QMRCDisplayWidget::setImageWidgetsEnabled(bool b)
{
  m_ImageWidgetsEnabled = b;
  showWidgets(b, m_ImageWidgets);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void QMRCDisplayWidget::setMovieWidgetsEnabled(bool b)
{
  m_MovieWidgetsEnabled = b;
  showWidgets(b, m_MovieWidgets);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void QMRCDisplayWidget::showWidgets(bool b, QList<QWidget*> &list)
{
  foreach (QWidget* w, list)
  {
    w->setVisible(b);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QImage QMRCDisplayWidget::signed16Image(qint16* data, MRCHeader &header)
{
  qint16 dmax = std::numeric_limits<qint16>::min();
  qint16 dmin = std::numeric_limits<qint16>::max();
  size_t nVoxels = header.nx * header.ny;
  for (size_t i = 0; i < nVoxels; ++i)
  {
    if(data[i] > dmax) dmax = data[i];
    if(data[i] < dmin) dmin = data[i];
  }

  // Generate a Color Table
  float max = static_cast<float>(dmax);
  float min = static_cast<float>(dmin);
  int numColors = static_cast<int>((max - min) + 1);

  QVector<QRgb>         colorTable;


  // Only generate the color table if the number of colors does not match
  if(colorTable.size() != numColors)
  {
    colorTable.resize(numColors);
    float range = max - min;

    float r, g, b;
    for (int i = 0; i < numColors; i++)
    {
      int16_t val = static_cast<int16_t>(min + ((float)i / numColors) * range);
      QMRCDisplayWidget::getColorCorrespondingTovalue(val, r, g, b, max, min);
      colorTable[i] = qRgba(r * 255, g * 255, b * 255, 255);
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
      image.setPixel(x, y, colorTable[colorIndex]);
    }
  }
  return image;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QImage QMRCDisplayWidget::floatImage(float* data, MRCHeader &header)
{

  float dmax = std::numeric_limits<float>::min();
  float dmin = std::numeric_limits<float>::max();
  size_t nVoxels = header.nx * header.ny;
  for (size_t i = 0; i < nVoxels; ++i)
  {
    if(data[i] > dmax) dmax = data[i];
    if(data[i] < dmin) dmin = data[i];
  }

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

  // Generate a Color Table
  int numColors = 256; // We are going to fix this at 256 colors
  QVector<QRgb> colorTable(numColors);
  // generate the color table
  dmin = 0.0;
  dmax = 255;
  float range = 256;
  float r, g, b;
  for (int i = 0; i < numColors; i++)
  {
    int16_t val = static_cast<int16_t>(dmin + ((float)i / numColors) * range);
    QMRCDisplayWidget::getColorCorrespondingTovalue(val, r, g, b, dmax, dmin);
    colorTable[i] = qRgba(r * 255, g * 255, b * 255, 255);
  }

  // Create an RGB Image
  QImage image(header.nx, header.ny, QImage::Format_ARGB32);

  int idx = 0;
  for (int y = 0; y < header.ny; ++y)
  {
    for (int x = 0; x < header.nx; ++x)
    {
      idx = (header.nx * y) + x;
      int colorIndex = iData[idx] - static_cast<int>(imin);
      image.setPixel(x, y, colorTable[colorIndex]);
    }
  }

//  std::cout << "Min int QImage Value:" << imin << std::endl;
//  std::cout << "Max int QImage Value:" << imax << std::endl;
//  xzImage.save("/tmp/xz_image.tif");
  return image;

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QImage QMRCDisplayWidget::xzSigned16CrossSection(qint16* data, size_t nVoxels, int* voxelMin, int* voxelMax)
{
  qint16 dmax = std::numeric_limits<qint16>::min();
  qint16 dmin = std::numeric_limits<qint16>::max();
  for (size_t i = 0; i < nVoxels; ++i)
  {
    if(data[i] > dmax) dmax = data[i];
    if(data[i] < dmin) dmin = data[i];
  }

//  std::cout << "Min MRC Value:" << dmin << std::endl;
//  std::cout << "Max MRC Value:" << dmax << std::endl;


  // Generate a Color Table
  float max = static_cast<float>(dmax);
  float min = static_cast<float>(dmin);
  int numColors = static_cast<int>((max - min) + 1);
  QVector<QRgb> colorTable(numColors);

  float range = max - min;

  float r, g, b;
  for (int i = 0; i < numColors; i++)
  {
    int16_t val = static_cast<int16_t>(min + ((float)i / numColors) * range);
    QMRCDisplayWidget::getColorCorrespondingTovalue(val, r, g, b, max, min);
    colorTable[i] = qRgba(r * 255, g * 255, b * 255, 255);
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
      image.setPixel(x, z, colorTable[colorIndex]);
    }
  }
  return image;
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QImage QMRCDisplayWidget::xzFloatCrossSection(float* data, size_t nVoxels, int* voxelMin, int* voxelMax)
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

  // Generate a Color Table
  int numColors = 256; // We are going to fix this at 256 colors
  QVector<QRgb> colorTable(numColors);
  // generate the color table
  dmin = 0.0;
  dmax = 255;
  float range = 256;
  float r, g, b;
  for (int i = 0; i < numColors; i++)
  {
    int16_t val = static_cast<int16_t>(dmin + ((float)i / numColors) * range);
    QMRCDisplayWidget::getColorCorrespondingTovalue(val, r, g, b, dmax, dmin);
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
#if 0
    // This code section writes all the individual images to files in the /tmp directory
    // this is useful for debugging
    QString fname("/tmp/single_slice_z_");
    fname.append(QString::number(z)).append(".tif");
    image.save(fname);
#endif
  }
//  std::cout << "Min int QImage Value:" << imin << std::endl;
//  std::cout << "Max int QImage Value:" << imax << std::endl;
//  xzImage.save("/tmp/xz_image.tif");
  return xzImage;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void QMRCDisplayWidget::drawOrigin(QImage image)
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

  int pImageWidth = imageWidth;// + pxOffset * 2;
  int pImageHeight = imageHeight;// + pyOffset * 2;

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

  QPoint point(0, 0);
  painter.drawImage(point, image); // Draw the image we just generated into the QPainter's canvas

  qint32 penWidth = 2;
  painter.setPen(QPen(QColor(255, 0, 0, 255), penWidth, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));


  pxWide = metrics.width(QString("(0,0)"));
  painter.drawText(5, pImageHeight - 5, "(0,0)");


  // Draw slightly transparent lines
  painter.setPen(QPen(QColor(255, 255, 0, 180), penWidth, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
//  painter.drawLine(pImageWidth / 2, pImageHeight / 2, pImageWidth - pxOffset, pImageHeight / 2);
//  painter.drawLine(pImageWidth / 2, pImageHeight / 2, pImageWidth / 2, pImageHeight - pyOffset);

  painter.end();

  m_CurrentImage = pImage;

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
//void QMRCDisplayWidget::fireOriginCB_Changed()
//{
//  on_originCB_currentIndexChanged(originCB->currentIndex());
//}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
#if 0
void QMRCDisplayWidget::on_originCB_currentIndexChanged(int corner)
{

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

  // Calculate the approx memory usage
  emit memoryCalculationNeedsUpdated();

  // This will display the image in the graphics scene
  m_GraphicsView->loadBaseImageFile(m_CurrentImage);
}
#endif

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
void QMRCDisplayWidget::getColorCorrespondingTovalue(int16_t val,
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




// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void QMRCDisplayWidget::loadXZSliceReconstruction(QString reconMRCFilePath)
{
 // std::cout << "QMRCDisplayWidget::loadXZSliceReconstruction(" << reconMRCFilePath.toStdString() << ")" << std::endl;
  MRCHeader header;
  header.feiHeaders = NULL;
  ::memset(&header, 0, 1024); // Splat zeros across the entire structure
  MRCReader::Pointer reader = MRCReader::New(true);
  int err = reader->readHeader(reconMRCFilePath.toStdString(), &header);
  if(err < 0)
  {
    FREE_FEI_HEADERS( header.feiHeaders )
    return;
  }
  //Read the Entire File
  int voxelMin[3] =
  { 0, 0, 0 };
  int voxelMax[3] =
  { header.nx, header.ny, header.nz };
  size_t nVoxels = (voxelMax[0] - voxelMin[0]) * (voxelMax[1] - voxelMin[1]) * (voxelMax[2] - voxelMin[2]);

  err = reader->read(reconMRCFilePath.toStdString(), NULL, NULL);

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
  if(dataType == 1)
  {
    qint16* data = reinterpret_cast<qint16*>(reader->getDataPointer());
    image = QMRCDisplayWidget::xzSigned16CrossSection(data, nVoxels, voxelMin, voxelMax);
  }
  else if(dataType == 2)
  {
    float* data = reinterpret_cast<float*>(reader->getDataPointer());
    image = QMRCDisplayWidget::xzFloatCrossSection(data, nVoxels, voxelMin, voxelMax);
  }

  m_CurrentImage = image.mirrored(false, true);

  if (m_ImageWidgetsEnabled == true) { showWidgets(true, m_ImageWidgets); }
  if (m_MovieWidgetsEnabled == true) { showWidgets(true, m_MovieWidgets); }

  m_GraphicsView->loadBaseImageFile(m_CurrentImage);
 // std::cout << "QMRCDisplayWidget::loadXZSliceReconstruction( COMPLETE )" << std::endl;

}

