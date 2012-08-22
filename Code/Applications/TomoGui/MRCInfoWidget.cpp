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


#include "MRCInfoWidget.h"

#include <QtCore/QFileInfo>


#include "TomoEngine/IO/MRCHeader.h"
#include "TomoEngine/IO/MRCReader.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
MRCInfoWidget::MRCInfoWidget(QWidget *parent) :
QWidget(parent)
{
  setupUi(this);
  this->setWindowFlags(Qt::Tool);

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
MRCInfoWidget::~MRCInfoWidget()
{}

#define SET_MRC_INFO(var)\
  var->setText(QString::number(header.var));


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void MRCInfoWidget::setInfo(QString mrcFilePath)
{

  QFileInfo fi(mrcFilePath);

  this->setWindowTitle(fi.fileName());
  MRCHeader header;
  ::memset(&header, 0, sizeof(header));

  MRCReader::Pointer reader = MRCReader::New(true);
  // Read the header from the file
  int err = reader->readHeader(mrcFilePath.toStdString(), &header);
  if(err < 0)
  {
      std::cout << "Failed to open MRC File: " << mrcFilePath.toStdString() << std::endl;
    return;
  }
  int tiltIndex = 0;
  // Transfer the meta data from the MRC Header to the GUI
  SET_MRC_INFO(nx)
  SET_MRC_INFO(ny)
  SET_MRC_INFO(nz)
  SET_MRC_INFO(mode)
  SET_MRC_INFO(nxstart)
  SET_MRC_INFO(nystart)
  SET_MRC_INFO(nzstart)
  SET_MRC_INFO(mx)
  SET_MRC_INFO(my)
  SET_MRC_INFO(mz)
  SET_MRC_INFO(xlen)
  SET_MRC_INFO(ylen)
  SET_MRC_INFO(zlen)
  SET_MRC_INFO(alpha)
  SET_MRC_INFO(beta)
  SET_MRC_INFO(gamma)
  SET_MRC_INFO(mapc)
  SET_MRC_INFO(mapr)
  SET_MRC_INFO(maps)
  SET_MRC_INFO(amin)
  SET_MRC_INFO(amax)
  SET_MRC_INFO(amean)
  SET_MRC_INFO(idtype)
  SET_MRC_INFO(lens)
  SET_MRC_INFO(nd1)
  SET_MRC_INFO(nd2)
  SET_MRC_INFO(vd1)
  SET_MRC_INFO(vd2)
  SET_MRC_INFO(xorg)
  SET_MRC_INFO(yorg)
  SET_MRC_INFO(zorg)

  // If we have the FEI headers get that information
  if(header.feiHeaders != NULL)
  {
    FEIHeader fei = header.feiHeaders[tiltIndex];
//    a_tilt->setText(QString::number(fei.a_tilt));
//    b_tilt->setText(QString::number(fei.b_tilt));
//    x_stage->setText(QString::number(fei.x_stage));
//    y_stage->setText(QString::number(fei.y_stage));
//    z_stage->setText(QString::number(fei.z_stage));
//    x_shift->setText(QString::number(fei.x_shift));
//    y_shift->setText(QString::number(fei.y_shift));
//    defocus->setText(QString::number(fei.defocus));
//    exp_time->setText(QString::number(fei.exp_time));
//    mean_int->setText(QString::number(fei.mean_int));
//    tiltaxis->setText(QString::number(fei.tiltaxis));
//    pixelsize->setText(QString::number(fei.pixelsize));
//    magnification->setText(QString::number(fei.magnification));
//    voltage->setText(QString::number(fei.voltage));
  }
  else
  {
    return;
  }
}



