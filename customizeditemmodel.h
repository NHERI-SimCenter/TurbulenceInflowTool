#ifndef CUSTOMIZEDITEMMODEL_H
#define CUSTOMIZEDITEMMODEL_H

#include <QStandardItemModel>

class CustomizedItemModel : public QStandardItemModel
{
    Q_OBJECT
public:
    CustomizedItemModel();
    QVariant data(const QModelIndex &index, int role) const;

signals:

public slots:
};

#endif // CUSTOMIZEDITEMMODEL_H
