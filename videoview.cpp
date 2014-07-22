#include "videoview.h"
#include "ui_videoview.h"
#include <QPixmap>
#include <QMouseEvent>
#include <QPainter>


VideoView::VideoView(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::VideoView)
{
    ui->setupUi(this);
    connect(ui->videoLabel,SIGNAL(rectSelected(QRect)),this,SIGNAL(rectSelected(QRect)));
}

void VideoView::setImage(QImage* image){
    ui->videoLabel->resize(image->size());
    ui->videoLabel->setPixmap(QPixmap::fromImage(*image));
}

VideoView::~VideoView()
{
    delete ui;
}


void VideoView::on_toolButton_clicked()
{
    ui->videoLabel->startSelection();
}

VideoLabel::VideoLabel(QWidget *parent) :
    QLabel(parent)
{

}

VideoLabel::~VideoLabel(){

}

void VideoLabel::mousePressEvent(QMouseEvent *event){
    if(selecting)
        if(event->button()==Qt::LeftButton){
            first = event->pos();
        }
}

void VideoLabel::mouseMoveEvent(QMouseEvent *event){
    if(selecting)
            second = event->pos();
}

void VideoLabel::mouseReleaseEvent(QMouseEvent *event){
    if(selecting)
        if(event->button()==Qt::LeftButton){
            second = event->pos();
            selecting = false;
            emit(rectSelected(QRect(first,second).normalized()));
        }
}

void VideoLabel::paintEvent(QPaintEvent *event){
    QLabel::paintEvent(event);
    QPainter p;
    p.begin(this);
    p.setPen(Qt::blue);
    p.drawRect(QRect(first,second));
    p.end();
}

void VideoLabel::startSelection(){
    selecting = true;
}
