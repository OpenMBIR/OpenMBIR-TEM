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
#ifndef IMAGEGRAPHICSDELEGATE_H_
#define IMAGEGRAPHICSDELEGATE_H_

#include "MXA/Common/MXASetGetMacros.h"

//-- Qt Includes
#include <QtCore/QObject>
#include <QtGui/QGraphicsPixmapItem>
#include <QtGui/QPainter>

//-- STL Includes
#include <string>

// --- Forward declarations of Classes ----
class QMainWindow;
class QGraphicsItem;
class QGraphicsView;
class QGraphicsScene;

/**
 * @class ImageGraphicsDelegate ImageGraphicsDelegate.h ImageGraphicsDelegate.h
 * @brief A Delegate class that renders Images using a QGraphicsView as its output display.
 * @author Mike Jackson
 * @date May 29, 2007
 * @version $Revision: 1.5 $
 */
class ImageGraphicsDelegate : public QObject
{
    Q_OBJECT

  public:
    ImageGraphicsDelegate(QObject* parent = 0);
    virtual ~ImageGraphicsDelegate();

    MXA_INSTANCE_PROPERTY(QString, DelegateName)
    MXA_INSTANCE_PROPERTY(QMainWindow*, MainWindow)
    MXA_INSTANCE_PROPERTY(QGraphicsView*, GraphicsView)
    MXA_INSTANCE_PROPERTY(QGraphicsScene*, GraphicsScene)
    MXA_INSTANCE_PROPERTY(QImage, CachedImage)
    MXA_INSTANCE_PROPERTY(QImage, ScaledCachedImage)
    MXA_INSTANCE_PROPERTY(QImage, OverlayImage)
    MXA_INSTANCE_PROPERTY(QImage, CompositedImage)
    MXA_INSTANCE_PROPERTY(bool, CompositeImages)

    /**
     * @brief Displays a Text message in the graphics view. This is typically used
     * when there is an error. You should NOT use this to display general String
     * data. Use an QHDFStringDataWindow instead
     * @param message The message to display
     */
    void displayTextMessage(QString message);

    /**
     * @brief sets all cached values to NULL or empty
     */
    void resetCaches();

  public slots:

    /**
     * @brief This slot can be called when the parent Widget gets resized so the
     * Graphics view can also be resized.
     */
    void on_parentResized();

    /**
     * @brief Sets the zoom value as a floating point number
     * where 1.0 is 100% or NO ZOOM and -1.0 represents fit to current window.
     * @param zoomFactor The value of the zoom Factor
     */
    void setZoomIndex(int zoomIndex);

    int getZoomIndex();

#if 0
    /**
     * @brief Increase the zoom level of the view which is the same as zooming in.
     */
    void increaseZoom();

    /**
     * @brief Decrease the zoom level of the view which is the same as zooming out.
     */
    void decreaseZoom();
#endif


    /**
     * @brief If checkbox_state is "Qt::Checked" then ensures the Image fits into the current size of the QGraphicsView object.
     * @param checkbox_state Use eithet Qt::Checked or Qt::UnChecked.
     */
    //  void fitToWindow(int checkbox_state);

    /**
     * @brief Forcibly updates the embedded QGraphicsView with an option to update the
     * QGraphicsScene also
     * @param updateGraphicsScene Boolean flag to optionally update the graphics Scene. Defaults to true.
     */
    void updateGraphicsView(bool updateGraphicsScene = true);

    /**
     * @brief Sets the QPainter composition mode to Clear
     */
    void setClearMode()
    {
      m_composition_mode = QPainter::CompositionMode_Clear;
    }

    /**
     * @brief Sets the QPainter composition mode to Source
     */
    void setSourceMode()
    {
      m_composition_mode = QPainter::CompositionMode_Source;
    }

    /**
     * @brief Sets the QPainter composition mode to Destination
     */
    void setDestinationMode()
    {
      m_composition_mode = QPainter::CompositionMode_Destination;
    }

    /**
     * @brief Sets the QPainter composition mode to SourceOver
     */
    void setSourceOverMode()
    {
      m_composition_mode = QPainter::CompositionMode_SourceOver;
    }

    /**
     * @brief Sets the QPainter composition mode to Destination Over
     */
    void setDestinationOverMode()
    {
      m_composition_mode = QPainter::CompositionMode_DestinationOver;
    }

    /**
     * @brief Sets the QPainter composition mode to SourceIn
     */
    void setSourceInMode()
    {
      m_composition_mode = QPainter::CompositionMode_SourceIn;
    }

    /**
     * @brief Sets the QPainter composition mode to DestIn
     */
    void setDestInMode()
    {
      m_composition_mode = QPainter::CompositionMode_DestinationIn;
    }

    /**
     * @brief Sets the QPainter composition mode to SourceOut
     */
    void setSourceOutMode()
    {
      m_composition_mode = QPainter::CompositionMode_SourceOut;
    }

    /**
     * @brief Sets the QPainter composition mode to DestOut
     */
    void setDestOutMode()
    {
      m_composition_mode = QPainter::CompositionMode_DestinationOut;
    }

    /**
     * @brief Sets the QPainter composition mode to SourceAtop
     */
    void setSourceAtopMode()
    {
      m_composition_mode = QPainter::CompositionMode_SourceAtop;
    }

    /**
     * @brief Sets the QPainter composition mode to DestAtop
     */
    void setDestAtopMode()
    {
      m_composition_mode = QPainter::CompositionMode_DestinationAtop;
    }

    /**
     * @brief Sets the QPainter composition mode to XOr
     */
    void setXorMode()
    {
      m_composition_mode = QPainter::CompositionMode_Xor;
    }

    /**
     * @brief Sets the QPainter composition mode to Plus
     */
    void setPlusMode()
    {
      m_composition_mode = QPainter::CompositionMode_Plus;
    }

    /**
     * @brief Sets the QPainter composition mode to Mutliply
     */
    void setMultiplyMode()
    {
      m_composition_mode = QPainter::CompositionMode_Multiply;
    }

    /**
     * @brief Sets the QPainter composition mode to Screen
     */
    void setScreenMode()
    {
      m_composition_mode = QPainter::CompositionMode_Screen;
    }

    /**
     * @brief Sets the QPainter composition mode to Overlay
     */
    void setOverlayMode()
    {
      m_composition_mode = QPainter::CompositionMode_Overlay;
    }

    /**
     * @brief Sets the QPainter composition mode to Darken
     */
    void setDarkenMode()
    {
      m_composition_mode = QPainter::CompositionMode_Darken;
    }

    /**
     * @brief Sets the QPainter composition mode to Lighten
     */
    void setLightenMode()
    {
      m_composition_mode = QPainter::CompositionMode_Lighten;
    }

    /**
     * @brief Sets the QPainter composition mode to ColorDodge
     */
    void setColorDodgeMode()
    {
      m_composition_mode = QPainter::CompositionMode_ColorDodge;
    }

    /**
     * @brief Sets the QPainter composition mode to ColorBurn
     */
    void setColorBurnMode()
    {
      m_composition_mode = QPainter::CompositionMode_ColorBurn;
    }

    /**
     * @brief Sets the QPainter composition mode to Hard Light
     */
    void setHardLightMode()
    {
      m_composition_mode = QPainter::CompositionMode_HardLight;
    }

    /**
     * @brief Sets the QPainter composition mode to SoftLight
     */
    void setSoftLightMode()
    {
      m_composition_mode = QPainter::CompositionMode_SoftLight;
    }

    /**
     * @brief Sets the QPainter composition mode to Difference
     */
    void setDifferenceMode()
    {
      m_composition_mode = QPainter::CompositionMode_Difference;
    }

    /**
     * @brief Sets the QPainter composition mode to Exclusion
     */
    void setExclusionMode()
    {
      m_composition_mode = QPainter::CompositionMode_Exclusion;
    }

  private:
    QGraphicsPixmapItem* m_CurrentGraphicsItem;

    double _zoomFactors[10];
    int _zoomIndex;
    bool _shouldFitToWindow;
    QPainter::CompositionMode m_composition_mode;

    /**
     * @brief Displays a Text message in the QGraphicsScene
     * @param message The message to display
     */
    void _displayTextMessage(QString message);

    /**
     * @brief Scales the cached image by a _zoomFactor
     * @return Returns a QImage that is scaled copy of the cached image
     */
    QImage _scaleImage();

    /**
     * @brief Scales a QImage by the current scaling factor
     * @param image The input QImage to scale
     * @return The Scaled Image
     */
    QImage _scaleImage(QImage image);

    ImageGraphicsDelegate(const ImageGraphicsDelegate&); //Copy Constructor Not Implemented
    void operator=(const ImageGraphicsDelegate&); //Copy Assignment Not Implemented
};

#endif /* IMAGEGRAPHICSDELEGATE_H_ */
