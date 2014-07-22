#include <QTableWidgetItem>
#include <QHeaderView>
#include "qbufferbar.h"

QBufferBar::QBufferBar(QWidget* parent):QTableWidget(parent)
{
    setColumnCount(1);
    id=1;
    verticalHeader()->hide();
    horizontalHeader()->hide();
}

QBufferBar::~QBufferBar()
{

}

void QBufferBar::addBuffer(QColor color, QString comment, int colorIndex){
    QTableWidgetItem* item = new QTableWidgetItem(QString::number(id++));
    item->setFlags(Qt::ItemIsSelectable|Qt::ItemIsEnabled);
    item->setBackground(color);
    item->setForeground(QColor("#d0d0d0"));
    item->setToolTip(comment);
    item->setData(COLOR_INDEX,colorIndex);

    insertRow(rowCount());
    setItem(rowCount()-1,0,item);
    setRowHeight(rowCount()-1,18);
}

void QBufferBar::keyPressEvent(QKeyEvent* event)
{
  if(event->key()==Qt::Key_Escape)
  {
      clearSelection();
  }
  if(event->key()==Qt::Key_Down)
  {

      if(currentRow()<rowCount()-1)
      {
        //clearSelection();
        setCurrentCell(currentRow()+1,currentColumn());        
      }
  }
  if(event->key()==Qt::Key_Up)
  {
      if(currentRow()>0)
      {
        //clearSelection();
        setCurrentCell(currentRow()-1,currentColumn());        
      }
        
  }

}
