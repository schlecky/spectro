#ifndef VIDEOVIEW_H
#define VIDEOVIEW_H

#include <QDialog>
#include <QImage>
#include <QLabel>

namespace Ui {
class VideoView;
}

class VideoLabel : public QLabel
{
    Q_OBJECT;
public:
    explicit VideoLabel(QWidget *parent = 0);
    ~VideoLabel();

private:
    bool selecting;
    QPoint first,second;

    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);

    void paintEvent(QPaintEvent *event);

public slots:
    void startSelection();

signals:
    void rectSelected(QRect rect);
};



class VideoView : public QDialog
{
    Q_OBJECT
    
public:
    explicit VideoView(QWidget *parent = 0);
    ~VideoView();

    
private:
    Ui::VideoView *ui;
signals:
    void rectSelected(QRect rect);
public slots:
    void setImage(QImage* image);
private slots:
    void on_toolButton_clicked();
};

#endif // VIDEOVIEW_H
