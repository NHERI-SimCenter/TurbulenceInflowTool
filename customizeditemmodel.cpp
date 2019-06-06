#include "customizeditemmodel.h"

// this is a wrapper class to customize the appearance of the treeView

CustomizedItemModel::CustomizedItemModel()
{

}

QVariant
CustomizedItemModel::data(const QModelIndex &index, int role) const
{
    if (role == Qt::TextAlignmentRole) {
        return Qt::AlignCenter;
    }
    else if (role == Qt::SizeHintRole) {
        return QSize(20, 50);
    }
    else {
        return QStandardItemModel::data(index, role);
    }
}
