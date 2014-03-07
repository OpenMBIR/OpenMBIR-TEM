#ifndef _TomogramTiltLoader_H_
#define _TomogramTiltLoader_H_

#include <QtCore/QVector>

#include <QtGui/QDialog>


#include "ui_TomogramTiltLoader.h"

class TomogramTiltLoader : public QDialog, private Ui::TomogramTiltLoader
{

    Q_OBJECT

  public:
    TomogramTiltLoader(QWidget* parent = NULL);

    virtual ~TomogramTiltLoader();

    QVector<float> getATilts();
    QVector<float> getBTilts();
    float getPixelSize();


  protected:
    void setupGui();

  private:

    TomogramTiltLoader(const TomogramTiltLoader&); // Copy Constructor Not Implemented
    void operator=(const TomogramTiltLoader&); // Operator '=' Not Implemented



};


#endif /* _TomogramTiltLoader_H_ */
