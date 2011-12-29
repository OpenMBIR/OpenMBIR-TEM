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

#include "GainsOffsetsTableModel.h"

#include <iostream>

#include <QApplication>
#include <QtGui/QStyleOptionComboBox>
#include <QtGui/QAbstractItemDelegate>

#include "GainsOffsetsItemDelegate.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
GainsOffsetsTableModel::GainsOffsetsTableModel(QObject* parent) :
QAbstractTableModel(parent),
m_RowCount(0)
{
  m_ColumnCount = ColumnCount;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
GainsOffsetsTableModel::~GainsOffsetsTableModel()
{
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
Qt::ItemFlags GainsOffsetsTableModel::flags(const QModelIndex &index) const
{
  //  std::cout << "GainsOffsetsTableModel::flags" << std::endl;
  if (!index.isValid())
  {
    return Qt::NoItemFlags;
  }
  Qt::ItemFlags theFlags = QAbstractTableModel::flags(index);
  if (index.isValid())
  {
  //  theFlags |= Qt::ItemIsEditable | Qt::ItemIsSelectable | Qt::ItemIsEnabled;

    int col = index.column();
    if (col == Gains)
    {
      theFlags = Qt::ItemIsEditable | Qt::ItemIsSelectable | Qt::ItemIsEnabled;
    }
    else if (col == Offsets)
    {
      theFlags = Qt::ItemIsEditable | Qt::ItemIsSelectable | Qt::ItemIsEnabled;
    }
    else if (col == Exclude)
    {
      theFlags = Qt::ItemIsSelectable | Qt::ItemIsEnabled;
    }

  }
  return theFlags;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QVariant GainsOffsetsTableModel::data(const QModelIndex &index, qint32 role) const
{

  if (!index.isValid())
  {
    return QVariant();
  }

  if (role == Qt::SizeHintRole)
  {
    QStyleOptionFrame comboBox;
    QLineEdit contents("1");
    switch(index.column())
    {
      case TiltIndex:
      {
        contents.setText( QString("11.") );
        const QString header = headerData(TiltIndex, Qt::Horizontal, Qt::DisplayRole).toString();
        if (header.length() > contents.text().length()) contents.text() = header;
        break;
      }
      case A_Tilt:
      {
        contents.setText( QString("11.") );
        const QString header = headerData(A_Tilt, Qt::Horizontal, Qt::DisplayRole).toString();
        if (header.length() > contents.text().length()) contents.text() = header;
        break;
      }
      case B_Tilt:
      {
        contents.setText( QString("11.") );
        const QString header = headerData(B_Tilt, Qt::Horizontal, Qt::DisplayRole).toString();
        if (header.length() > contents.text().length()) contents.text() = header;
        break;
      }
      case Gains:
      {
        contents.setText( QString("1,1,1") );
        const QString header = headerData(Gains, Qt::Horizontal, Qt::DisplayRole).toString();
        if (header.length() > contents.text().length()) contents.text() = header;
        break;
      }
      case Offsets:
      {
        contents.setText( QString("11.") );
        const QString header = headerData(Offsets, Qt::Horizontal, Qt::DisplayRole).toString();
        if (header.length() > contents.text().length()) contents.text() = header;
        break;
      }
      case Exclude:
      {
        contents.setText( QString("11.") );
        const QString header = headerData(Exclude, Qt::Horizontal, Qt::DisplayRole).toString();
        if (header.length() > contents.text().length()) contents.text() = header;
        break;
      }
      default:
        Q_ASSERT(false);
    }
    QFontMetrics fontMetrics(data(index, Qt::FontRole) .value<QFont > ());
    comboBox.fontMetrics = fontMetrics;
    QSize size(fontMetrics.width(contents.text()), fontMetrics.height());
    return qApp->style()->sizeFromContents(QStyle::CT_ComboBox, &comboBox, size);
  }
  else if (role == Qt::TextAlignmentRole)
  {
    return int(Qt::AlignRight | Qt::AlignVCenter);
  }
  else if (role == Qt::DisplayRole || role == Qt::EditRole)
  {
    int col = index.column();
    if (col == TiltIndex)
    {
      return QVariant(m_AngleIndexes[index.row()]);
    }
    else if (col == A_Tilt)
    {
      return QVariant(m_ATilts[index.row()]);
    }
    else if (col == B_Tilt)
    {
      return QVariant(m_BTilts[index.row()]);
    }
    else if (col == Gains)
    {
      return QVariant(m_Gains[index.row()]);
    }
    else if (col == Offsets)
    {
      return QVariant(m_Offsets[index.row()]);
    }
    else if(col == Exclude)
    {
      return QVariant(m_Excludes[index.row()]);
    }
  }
  else if(role == Qt::CheckStateRole)
  {
    int col = index.column();
    if(col == Exclude)
    {
      return QVariant(m_Excludes[index.row()]);
    }
  }

  return QVariant();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QVariant GainsOffsetsTableModel::headerData(int section, Qt::Orientation orientation, int role) const
{
  if (orientation == Qt::Horizontal && role == Qt::DisplayRole)
  {
    switch(section)
    {
      case TiltIndex:
        return QVariant(QString("Index"));
      case A_Tilt:
        return QVariant(QString("A Tilts"));
      case B_Tilt:
        return QVariant(QString("B Tilts"));
      case Gains:
        return QVariant(QString("Gains"));
        break;
      case Offsets:
        return QVariant(QString("Offsets"));
        break;
      case Exclude:
        return QVariant(QString("Exclude Tilt"));
        break;
      default:
        break;
    }

  }
  return QVariant();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int GainsOffsetsTableModel::rowCount(const QModelIndex &index) const
{
  return index.isValid() ? 0 : m_RowCount;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int GainsOffsetsTableModel::columnCount(const QModelIndex &index) const
{
  return index.isValid() ? 0 : m_ColumnCount;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
bool GainsOffsetsTableModel::setHeaderData(int col, Qt::Orientation o, const QVariant& var, int role)
{
  return false;
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
bool GainsOffsetsTableModel::setData(const QModelIndex & index, const QVariant & value, int role)
{
  // std::cout << "GainsOffsetsTableModel::setData " << value.toString().toStdString() << std::endl;
  if(!index.isValid() || index.row() < 0 || index.row() >= m_ATilts.count() || index.column() < 0 || index.column() >= m_ColumnCount)
  {
    return false;
  }

  if(role == Qt::EditRole)
  {
    bool ok;
    qint32 row = index.row();
    qint32 col = index.column();
    switch(col)
    {
      case TiltIndex:
        m_AngleIndexes[row] = value.toInt(&ok);
        break;
      case A_Tilt:
        m_ATilts[row] = value.toFloat(&ok);
        break;
      case B_Tilt:
        m_BTilts[row] = value.toFloat(&ok);
        break;
      case Gains:
        m_Gains[row] = value.toFloat(&ok);
        break;
      case Offsets:
        m_Offsets[row] = value.toFloat(&ok);
        break;
      case Exclude:
        m_Excludes[row] = value.toBool();
        break;
      default:
        Q_ASSERT(false);

    }
    emit dataChanged(index, index);
  }


  return true;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
bool GainsOffsetsTableModel::insertRows(int row, int count, const QModelIndex& index)
{
  QString axis("<0,0,1>");
  int angleIndx = 0;
  float angle = 0.0;
  float gain = 0.0;
  float offset = 0.0f;
  bool exclude = false;

  beginInsertRows(QModelIndex(), row, row + count - 1);
  for (int i = 0; i < count; ++i)
  {
    m_AngleIndexes.append(angleIndx);
    m_ATilts.append(angle);
    m_BTilts.append(angle);
    m_Offsets.append(offset);
    m_Gains.append(gain);
    m_Excludes.append(exclude);
    m_RowCount = m_ATilts.count();
  }
  endInsertRows();
  emit dataChanged(index, index);
  return true;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
bool GainsOffsetsTableModel::removeRows(int row, int count, const QModelIndex& index)
{
  if (count < 1)
  {
    return true;
  } // No Rows to remove
  beginRemoveRows(QModelIndex(), row, row + count - 1);
  for (int i = 0; i < count; ++i)
  {
    m_AngleIndexes.remove(row);
    m_ATilts.remove(row);
    m_BTilts.remove(row);
    m_Offsets.remove(row);
    m_Gains.remove(row);
    m_Excludes.remove(row);
    m_RowCount = m_ATilts.count();
  }
  endRemoveRows();
  emit dataChanged(index, index);
  return true;
}

#if 0

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void GainsOffsetsTableModel::setColumnData(int col, QVector<float> &data)
{
  switch(col)
  {
    case A_Tilt:
      m_ATilts = data; break;
    case Offsets:
      m_Offsets = data;break;
    case Gains:
//      m_Gains = data; break;
    default:
      Q_ASSERT(false);
  }
}
#endif
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void GainsOffsetsTableModel::setTableData(QVector<int> angleIndexes,
                                          QVector<float> a_tilts,
                                          QVector<float> b_tilts,
                                          QVector<float> gains,
                                          QVector<float> offsets,
                                          QVector<bool> excludes)
{
  qint32 count = angleIndexes.count();
  qint32 row = 0;
  // Remove all the current rows in the table model
  removeRows(0, rowCount());
  QModelIndex topLeft;
  QModelIndex botRight;

  if (count > 0) {
    // Now mass insert the data to the table then emit that the data has changed
    beginInsertRows(QModelIndex(), row, row + count - 1);
    m_AngleIndexes = angleIndexes;
    m_ATilts = a_tilts;
    m_BTilts = b_tilts;
    m_Gains = gains;
    m_Offsets = offsets;
    m_Excludes = excludes;

    m_RowCount = count;
    endInsertRows();
    QModelIndex topLeft = createIndex(0, 0);
    QModelIndex botRight = createIndex(count-1, ColumnCount);
  }
  emit dataChanged(topLeft, botRight);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void GainsOffsetsTableModel::setGainsAndOffsets(QVector<float> gains,
                  QVector<float> offsets)
{
  qint32 count = gains.count();
  qint32 row = 0;

  QVector<int> angleIndexes(m_AngleIndexes);
  QVector<float> aTilts(m_ATilts);
  QVector<float> bTilts(m_BTilts);
  QVector<bool> excludes(m_Excludes);

  // Remove all the current rows in the table model
  removeRows(0, rowCount());
  QModelIndex topLeft;
  QModelIndex botRight;


  if (count > 0) {
    // Now mass insert the data to the table then emit that the data has changed
    beginInsertRows(QModelIndex(), row, row + count - 1);
    m_Gains = gains;
    m_Offsets = offsets;

    m_RowCount = count;
    endInsertRows();
    QModelIndex topLeft = createIndex(0, 0);
    QModelIndex botRight = createIndex(count-1, ColumnCount);
  }

  m_AngleIndexes = angleIndexes;
  m_ATilts = aTilts;
  m_BTilts = bTilts;
  m_Excludes = excludes;

  emit dataChanged(topLeft, botRight);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QAbstractItemDelegate* GainsOffsetsTableModel::getItemDelegate()
{
  return new GainsOffsetsItemDelegate;
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void GainsOffsetsTableModel::setInitialValues()
{

}

