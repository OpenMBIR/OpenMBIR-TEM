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

#ifndef _UIA_H_
#define _UIA_H_

#include <QtCore/QObject>
#include <QtCore/QSettings>
#include <QtGui/QGraphicsPolygonItem>
#include <QtGui/QGraphicsRectItem>
#include <QtGui/QBrush>

namespace UIA
{
    const static int Alpha = 155;
}


class RectangleCreator : public QObject, public QGraphicsPolygonItem
{
    Q_OBJECT;
    
public:
    enum { Type = UserType + 1 };
    enum CTRL_POINTS
    {
        NO_CTRL_POINT, UPPER_LEFT_CTRL_POINT, UPPER_RIGHT_CTRL_POINT, LOWER_RIGHT_CTRL_POINT, LOWER_LEFT_CTRL_POINT
    };
    
    
    RectangleCreator(const QPolygonF &polygon,
                 QGraphicsItem *parent = 0);
    virtual ~RectangleCreator();
    
    
    
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget = 0);
    
    static void deleteAllRectangleCreators(QGraphicsScene* scene);
    
    static RectangleCreator* NewRectangleCreator(QSettings &prefs, int index, QGraphicsItem *parent = 0);
    
    public slots:
    
    void setColor(QColor color);
    QColor getColor() { return brush().color(); }
    
    //  void setUpperLeft(unsigned int x, unsigned int y);
    void getUpperLeft(unsigned int &x, unsigned int &y);
    
    //  void setLowerRight(unsigned int x, unsigned int y);
    void getLowerRight(unsigned int &x, unsigned int &y);
    
    void setLineWidth(qreal w);
    qreal getLineWidth();
    
    // void setVisible(bool visible);
    
    
    
signals:
    
    void fireRectangleCreatorUpdated(RectangleCreator*);
    void fireRectangleCreatorAboutToDelete(RectangleCreator*);
    void fireRectangleCreatorDeleted(RectangleCreator*);
    void fireRectangleCreatorSelected(RectangleCreator*);
    void fireRectangleChanged(const QRect &);
    
protected:
    virtual void mouseMoveEvent(QGraphicsSceneMouseEvent *event);
    virtual void hoverMoveEvent(QGraphicsSceneHoverEvent *event);
    virtual void mousePressEvent(QGraphicsSceneMouseEvent *event);
    virtual void mouseReleaseEvent(QGraphicsSceneMouseEvent *event);
    virtual void keyPressEvent(QKeyEvent *event);
    
    virtual int type() const;
    CTRL_POINTS isInResizeArea(const QPointF &pos);
    
    
    static void deleteSelectedItems(QGraphicsScene *scene);
    
    static void duplicateSelectedItems(QGraphicsScene *scene);
    static void growSelectedItems(QGraphicsScene *scene);
    static void shrinkSelectedItems(QGraphicsScene *scene);
    
    
private:
    QPoint m_rectangleLeftCorner;
    bool m_isResizing;
    CTRL_POINTS m_CurrentResizeHandle;
    float ctrlPointSize;
    int m_UpperLeft[2];
    int m_LowerRight[2];
    qreal m_LineWidth;
};


#endif // _UIA_H_