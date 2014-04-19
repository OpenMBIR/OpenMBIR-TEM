/* ============================================================================
 * Copyright (c) 2010, Michael A. Jackson (BlueQuartz Software)
 * Copyright (c) 2010, Dr. Michael A. Groeber (US Air Force Research Laboratories
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
 * Neither the name of Michael A. Groeber, Michael A. Jackson, the US Air Force,
 * BlueQuartz Software nor the names of its contributors may be used to endorse
 * or promote products derived from this software without specific prior written
 * permission.
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
#ifndef SGMDFITEMDELEGATE_H_
#define SGMDFITEMDELEGATE_H_

#include <iostream>

#include <QtCore/QModelIndex>
#include <QtGui/QComboBox>
#include <QtGui/QPainter>
#include <QtGui/QStyleOptionViewItemV4>
#include <QtGui/QLineEdit>
#include <QtGui/QCheckBox>
#include <QtGui/QDoubleValidator>
#include <QtGui/QIntValidator>
#include <QtGui/QStyledItemDelegate>


#include "GainsOffsetsTableModel.h"

/**
 * @class GainsOffsetsItemDelegate GainsOffsetsItemDelegate.h StatsGenerator/GainsOffsetsItemDelegate.h
 * @brief This class creates the appropriate Editor Widget for the Tables
 * @author Michael A. Jackson for BlueQuartz Software
 * @date Dec 28, 2010
 * @version 1.0
 */
class GainsOffsetsItemDelegate : public QStyledItemDelegate
{
    Q_OBJECT

  public:
    explicit GainsOffsetsItemDelegate(QObject* parent = 0) :
      QStyledItemDelegate(parent)
    {
    }

    // -----------------------------------------------------------------------------
    //
    // -----------------------------------------------------------------------------
    void paint(QPainter* painter, const QStyleOptionViewItem& option, const QModelIndex& index) const
    {
      QStyledItemDelegate::paint(painter, option, index);
    }

    // -----------------------------------------------------------------------------
    //
    // -----------------------------------------------------------------------------
    QWidget* createEditor(QWidget* parent, const QStyleOptionViewItem& option, const QModelIndex& index) const
    {
//      QLineEdit* editor;
      //  QCheckBox* checkbox;
      //  QIntValidator* intValidator;
//      QDoubleValidator* dblValidator;

//      qint32 col = index.column();
//      switch(col)
//      {
//        case GainsOffsetsTableModel::TiltIndex:
//          editor = new QLineEdit(parent);
//          editor->setFrame(false);
//          intValidator = new QIntValidator(editor);
//          editor->setValidator(intValidator);
//          return editor;
//        case GainsOffsetsTableModel::A_Tilt:
//          editor = new QLineEdit(parent);
//          editor->setFrame(false);
//          dblValidator = new QDoubleValidator(editor);
//          dblValidator->setDecimals(4);
//          editor->setValidator(dblValidator);
//          return editor;
//        case GainsOffsetsTableModel::Gains:
//          editor = new QLineEdit(parent);
//          editor->setFrame(false);
//          dblValidator = new QDoubleValidator(editor);
//          dblValidator->setDecimals(4);
//          editor->setValidator(dblValidator);
//          return editor;
//        case GainsOffsetsTableModel::Offsets:
//          editor = new QLineEdit(parent);
//          editor->setFrame(false);
//          dblValidator = new QDoubleValidator(editor);
//          dblValidator->setDecimals(4);
//          editor->setValidator(dblValidator);
//          return editor;
//        case GainsOffsetsTableModel::Exclude:
//          checkbox = new QCheckBox(parent);
//          checkbox->setText("");
//          checkbox->setTristate(false);
//          return checkbox;
//        default:
//          break;
//      }
      return QStyledItemDelegate::createEditor(parent, option, index);
    }

    // -----------------------------------------------------------------------------
    //
    // -----------------------------------------------------------------------------
    void setEditorData(QWidget* editor, const QModelIndex& index) const
    {
      qint32 col = index.column();
      if (col == GainsOffsetsTableModel::TiltIndex
          || col == GainsOffsetsTableModel::A_Tilt)
      {
        QLineEdit* lineEdit = qobject_cast<QLineEdit* > (editor);
        Q_ASSERT(lineEdit);
        lineEdit->setText(index.model()->data(index).toString());
      }
//      else if (col == GainsOffsetsTableModel::Exclude)
//      {
//        QCheckBox* checkbox = qobject_cast<QCheckBox*>(editor);
//        Q_ASSERT(checkbox);
//        checkbox->setChecked(index.model()->data(index).toBool());
//      }
      else
      {
        QStyledItemDelegate::setEditorData(editor, index);
      }
    }

    // -----------------------------------------------------------------------------
    //
    // -----------------------------------------------------------------------------
    void setModelData(QWidget* editor, QAbstractItemModel* model, const QModelIndex& index) const
    {
      //  std::cout << "GainsOffsetsItemDelegate::setModelData" << std::endl;
      qint32 col = index.column();
      //  bool ok = false;
      if (col == GainsOffsetsTableModel::TiltIndex)
      {
        QLineEdit* lineEdit = qobject_cast<QLineEdit* > (editor);
        Q_ASSERT(lineEdit);
        bool ok = false;
        int v = lineEdit->text().toInt(&ok);
        model->setData(index, v);
      }
      else if (col == GainsOffsetsTableModel::A_Tilt )
      {
        QLineEdit* lineEdit = qobject_cast<QLineEdit* > (editor);
        Q_ASSERT(lineEdit);
        bool ok = false;
        double v = lineEdit->text().toFloat(&ok);
        model->setData(index, v);
      }
//      else if (col == GainsOffsetsTableModel::Exclude)
//      {
//        QCheckBox* checkbox = qobject_cast<QCheckBox*>(editor);
//        Q_ASSERT(checkbox);
//        model->setData(index, checkbox->isChecked());
//      }
      else
      {
        QStyledItemDelegate::setModelData(editor, model, index);
      }
    }

  private:
    QModelIndex m_Index;

};


#endif /* SGMDFITEMDELEGATE_H_ */
