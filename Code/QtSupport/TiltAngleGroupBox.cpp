#include "TiltAngleGroupBox.h"

#include <QtCore/QDebug>
#include <QtCore/QString>
#include <QtCore/QDir>
#include <QtCore/QFile>
#include <QtCore/QFileInfo>
#include <QtCore/QStringList>
#include <QtCore/QStringListIterator>
#include <QtGui/QMessageBox>
#include <QtGui/QFileDialog>
#include <QtGui/QAbstractButton>

#include "QtSupport/QFileCompleter.h"


// Initialize private static member variable
QString TiltAngleGroupBox::m_OpenDialogLastDirectory = "";

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
TiltAngleGroupBox::TiltAngleGroupBox(QWidget* parent) :
  QGroupBox(parent),
  m_NumTilts(0)
{
  setupUi(this);
  setupGui();
  if ( m_OpenDialogLastDirectory.isEmpty() )
  {
    m_OpenDialogLastDirectory = QDir::homePath();
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
TiltAngleGroupBox::~TiltAngleGroupBox()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TiltAngleGroupBox::setupGui()
{
  // Put a file path completer to help out the user to select a valid file
  QFileCompleter* com = new QFileCompleter(this, false);
  rawTltFile->setCompleter(com);
  QObject::connect( com, SIGNAL(activated(const QString&)),
                    this, SLOT(on_rawTltFile_textChanged(const QString&)));

  {
    QDoubleValidator* dVal = new QDoubleValidator(this);
    dVal->setDecimals(6);
    startingAngle->setValidator(dVal);
  }

  {
    QDoubleValidator* dVal = new QDoubleValidator(this);
    dVal->setDecimals(6);
    increment->setValidator(dVal);
  }

  rawTltBtn->toggle();
  rawTltBtn->toggle();
  noTilts->setChecked(true);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TiltAngleGroupBox::on_startingAngle_textChanged(const QString& text)
{
  // Generate tilts if this button is ON
  if (startIncrementBtn->isChecked() == true)
  {
    m_Angles.clear();
    bool ok = false;
    float start = startingAngle->text().toFloat(&ok);
    float incr = increment->text().toFloat(&ok);
    for(int i = 0; i < m_NumTilts; i++)
    {
      m_Angles.push_back(start + (float)(i) * incr);
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TiltAngleGroupBox::on_increment_textChanged(const QString& text)
{
  // Generate tilts if this button is ON
  if (startIncrementBtn->isChecked() == true)
  {
    m_Angles.clear();
    bool ok = false;
    float start = startingAngle->text().toFloat(&ok);
    float incr = increment->text().toFloat(&ok);
    for(int i = 0; i < m_NumTilts; i++)
    {
      m_Angles.push_back(start + (float)(i) * incr);
    }
  }
}
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TiltAngleGroupBox::on_startIncrementBtn_toggled(bool b)
{
  startLabel->setEnabled(startIncrementBtn->isChecked());
  startingAngle->setEnabled(startIncrementBtn->isChecked());
  incrementLabel->setEnabled(startIncrementBtn->isChecked());
  increment->setEnabled(startIncrementBtn->isChecked());

  // Generate tilts if this button is ON
  if (startIncrementBtn->isChecked() == true)
  {
    m_Angles.clear();
    bool ok = false;
    float start = startingAngle->text().toFloat(&ok);
    float incr = increment->text().toFloat(&ok);
    for(int i = 0; i < m_NumTilts; i++)
    {
      m_Angles.push_back(start + (float)(i) * incr);
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TiltAngleGroupBox::on_rawTltBtn_toggled(bool b)
{
  m_Angles.clear();
  rawTltFile->setEnabled(rawTltBtn->isChecked());
  selectRawTlt->setEnabled(rawTltBtn->isChecked());
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TiltAngleGroupBox::on_selectRawTlt_clicked()
{

  QString ext("*.rawtlt");
  QString s = QString("RawTlt Files (") + ext + QString(");;All Files(*.*)");
  QString defaultName = m_OpenDialogLastDirectory + QDir::separator() + "Untitled";
  QString file = QFileDialog::getOpenFileName(this, tr("Select Input File"), defaultName, s);

  if(true == file.isEmpty())
  {
    return;
  }
  file = QDir::toNativeSeparators(file);
  // Store the last used directory into the private instance variable
  QFileInfo fi(file);
  m_OpenDialogLastDirectory = fi.path();
  rawTltFile->setText(file);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TiltAngleGroupBox::on_rawTltFile_textChanged(const QString& text)
{
  if(verifyPathExists(text, rawTltFile) == true)
  {
    loadAngles();
  }
}



// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
bool TiltAngleGroupBox::verifyPathExists(QString filePath, QLineEdit* lineEdit)
{
  QFileInfo fileinfo(filePath);
  if (false == fileinfo.exists())
  {
    lineEdit->setStyleSheet("border: 1px solid red;");
  }
  else
  {
    lineEdit->setStyleSheet("");
  }
  return fileinfo.exists();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TiltAngleGroupBox::loadAngles()
{

  //clear the current set of angles
  m_Angles.clear();

  // Open the file and read everything into a text buffer
  QFileInfo fi(rawTltFile->text());
  QFile source(fi.absoluteFilePath());
  source.open(QFile::ReadOnly);
  bool hasErrors = false;
  bool ok = false;

  while(source.atEnd() == false)
  {
    QString text = source.readLine();
//    qDebug() << text;
//    qDebug() << source.error();
//    qDebug() << source.errorString();

    m_Angles.push_back(text.toFloat(&ok));
    if(ok == false)
    {
      hasErrors = true;
    }
  }

  if(hasErrors == true)
  {
    QString ss = QString("Error occurred reading in the angle values from the file. Some of the angle values could not be parsed correctly");
    QMessageBox::StandardButton reply;
    reply = QMessageBox::critical(NULL, QString("Tilt Angle Error"), ss, QMessageBox::Ok);
  }

  source.close();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QVector<float> TiltAngleGroupBox::getAngles()
{
  return m_Angles;
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TiltAngleGroupBox::setNumTilts(int nTilts)
{
  m_NumTilts = nTilts;
}
