#include <QtWidgets/QTableWidget>
#include <QKeyEvent>
#include <QColor>

#ifndef _QBUFFERBAR_H_
#define _QBUFFERBAR_H_

/*
class QBufferBarItem : public QTableWidgetItem
{
 public:
    QBufferBarItem(const QString &text);
    ~QBufferBarItem();
 protected:
};
*/

class QBufferBar : public QTableWidget
{
 public: 
  QBufferBar(QWidget* parent=0);
  ~QBufferBar();
  void addBuffer(QColor color,QString comment,int colorIndex=0);
  enum DataRole
  {
      COLOR_INDEX = 32
  };

 protected:
  void keyPressEvent ( QKeyEvent * event );

 private:
   int id;

};

#endif
