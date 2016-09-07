/****************************************************************************
** Meta object code from reading C++ file 'selectedpointwidget.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.6.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../source/selectedpointwidget.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'selectedpointwidget.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.6.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_SelectedPointWidget_t {
    QByteArrayData data[11];
    char stringdata0[161];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_SelectedPointWidget_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_SelectedPointWidget_t qt_meta_stringdata_SelectedPointWidget = {
    {
QT_MOC_LITERAL(0, 0, 19), // "SelectedPointWidget"
QT_MOC_LITERAL(1, 20, 19), // "show_new_slider_pos"
QT_MOC_LITERAL(2, 40, 0), // ""
QT_MOC_LITERAL(3, 41, 3), // "pos"
QT_MOC_LITERAL(4, 45, 11), // "update_view"
QT_MOC_LITERAL(5, 57, 18), // "set_new_slider_pos"
QT_MOC_LITERAL(6, 76, 6), // "accept"
QT_MOC_LITERAL(7, 83, 6), // "reject"
QT_MOC_LITERAL(8, 90, 23), // "set_transform_direction"
QT_MOC_LITERAL(9, 114, 22), // "set_interpolation_type"
QT_MOC_LITERAL(10, 137, 23) // "set_transform_selection"

    },
    "SelectedPointWidget\0show_new_slider_pos\0"
    "\0pos\0update_view\0set_new_slider_pos\0"
    "accept\0reject\0set_transform_direction\0"
    "set_interpolation_type\0set_transform_selection"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_SelectedPointWidget[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       8,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    1,   54,    2, 0x08 /* Private */,
       4,    1,   57,    2, 0x08 /* Private */,
       5,    0,   60,    2, 0x08 /* Private */,
       6,    0,   61,    2, 0x08 /* Private */,
       7,    0,   62,    2, 0x08 /* Private */,
       8,    1,   63,    2, 0x08 /* Private */,
       9,    1,   66,    2, 0x08 /* Private */,
      10,    1,   69,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void, QMetaType::Int,    3,
    QMetaType::Void, QMetaType::Int,    3,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Int,    2,

       0        // eod
};

void SelectedPointWidget::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        SelectedPointWidget *_t = static_cast<SelectedPointWidget *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->show_new_slider_pos((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 1: _t->update_view((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 2: _t->set_new_slider_pos(); break;
        case 3: _t->accept(); break;
        case 4: _t->reject(); break;
        case 5: _t->set_transform_direction((*reinterpret_cast< const int(*)>(_a[1]))); break;
        case 6: _t->set_interpolation_type((*reinterpret_cast< const int(*)>(_a[1]))); break;
        case 7: _t->set_transform_selection((*reinterpret_cast< const int(*)>(_a[1]))); break;
        default: ;
        }
    }
}

const QMetaObject SelectedPointWidget::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_SelectedPointWidget.data,
      qt_meta_data_SelectedPointWidget,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *SelectedPointWidget::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *SelectedPointWidget::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_SelectedPointWidget.stringdata0))
        return static_cast<void*>(const_cast< SelectedPointWidget*>(this));
    return QDialog::qt_metacast(_clname);
}

int SelectedPointWidget::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 8)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 8;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 8)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 8;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
