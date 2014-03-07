#ifndef _TiltAngleGroupBox_H_
#define _TiltAngleGroupBox_H_


#include <QtCore/QVector>
#include <QtGui/QGroupBox>

#include "ui_TiltAngleGroupBox.h"

class QAbstractButton;

class TiltAngleGroupBox : public QGroupBox, private Ui::TiltAngleGroupBox
{
    Q_OBJECT

  public:
    TiltAngleGroupBox(QWidget* parent = NULL);

    virtual ~TiltAngleGroupBox();
    QVector<float> getAngles();
    void setNumTilts(int nTilts);

  public slots:
    void on_selectRawTlt_clicked();
    void on_rawTltFile_textChanged(const QString& text);
    void on_increment_textChanged(const QString& text);
    void on_startingAngle_textChanged(const QString& text);
    void on_startIncrementBtn_toggled(bool b);
    void on_rawTltBtn_toggled(bool b);

  protected:
    void setupGui();
    void toggleWidgetsForIncrement();
    void toggleWidgetsForRawTltFile();

    bool verifyPathExists(QString filePath, QLineEdit* lineEdit);
    void loadAngles();

  private:
    QVector<float> m_Angles;
    int m_NumTilts;

    static QString    m_OpenDialogLastDirectory;

    TiltAngleGroupBox(const TiltAngleGroupBox&); // Copy Constructor Not Implemented
    void operator=(const TiltAngleGroupBox&); // Operator '=' Not Implemented
};


#endif /* _TiltAngleGroupBox_H_ */
