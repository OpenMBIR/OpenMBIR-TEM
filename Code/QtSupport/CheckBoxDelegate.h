/*
 * This code was taken directly from
 * http://stackoverflow.com/questions/3363190/qt-qtableview-how-to-have-a-checkbox-only-column
 *
 * There was no license attached and the code was on a public web server so I will
 * assume the author put the code in the public domain.
 */


#ifndef CHECKBOXDELEGATE_H_
#define CHECKBOXDELEGATE_H_

#include <QtGui/QStyledItemDelegate>

/**
 * @brief
 *
 */
class CheckBoxDelegate : public QStyledItemDelegate
{
    Q_OBJECT;

  public:
    explicit CheckBoxDelegate(QObject* parent = NULL);
    virtual ~CheckBoxDelegate();


    void paint ( QPainter * painter, const QStyleOptionViewItem & option, const QModelIndex & index ) const;

    bool editorEvent(QEvent *event,
                     QAbstractItemModel *model,
                     const QStyleOptionViewItem &option,
                     const QModelIndex &index);
};

#endif /* CHECKBOXDELEGATE_H_ */
