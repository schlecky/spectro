#ifndef SPECTRO_H
#define SPECTRO_H

#include <QMainWindow>
#include <QTimer>
#include <QPoint>
#include <QVector>
#include <QMap>
#include <opencv2/opencv.hpp>

#include "videoview.h"
#include "spectrometer.h"
#include "qcustomplot.h"
#include "jcampdx.h"


namespace Ui {
class Spectro;
}

class Calibration
{
public:
    QString comment;
    QString name;
    double focal;
    double sin_i;
    double pas;
    QRect spectrumRect;
};



class Spectro : public QMainWindow
{
    Q_OBJECT
    
public:
    explicit Spectro(QWidget *parent = 0);
    ~Spectro();
    
private:

    enum State {
        NONE = 0,
        SELECT_POINT1 = 1,
        SELECT_POINT2 = 2,
    };

    State currentState;

    Ui::Spectro *ui;

    VideoView video;

    Spectrometer* spectro;

    QVector<QColor> colors;
    QVector<Tjcampdx> spectra;
    Tjcampdx lastSpectrum;

    int currentColor;

    QString lastDir;

    bool divData,addData,addConst,multConst,subData,calibration;

    QCPItemTracer *crosshairTracer;
    QCPItemStraightLine *wlLine;

    QCPItemText *wlLabel;

    QList<QTableWidgetItem *> selectedForOp;
    QVector<QCPGraph*> graphs;
    QCPGraph* liveSpec,*bbSpec;

    QList<Calibration> calibrationList;


    void resizeEvent(QResizeEvent *);
    void changeScale(double xmin, double xmax, double ymin, double ymax);
    void readCalibrationList();

private slots:
    void on_actionHideSelectedItems_triggered();
    void on_actionShowSelectedItems_triggered();
    void on_actionShowScale_triggered();
    void on_tabWidget_tabCloseRequested(int index);


    void uncheckMathButtons();
    void on_btnDivBuffer_clicked();
    void on_btnAddBuffer_clicked();
    void on_btnSubBuffer_clicked();

    void on_bufferBar_cellClicked(int row, int column);
    void on_actionShowMath_triggered();

    void on_btnOkScale_clicked();

    void on_btnAutoScale_clicked();

    void xScaleChanged(QCPRange range);
    void yScaleChanged(QCPRange range);

    void plotClicked(QMouseEvent *event);
    void plotMouseMoved(QMouseEvent *event);
    void focusChanged(QWidget* old,QWidget* now);
    void on_btnSetRef_clicked();

    void on_actionSpectro_triggered();
    void on_btnUnsetRef_clicked();
    //void on_btnCalibrate_clicked();
    void on_btnCalibrateX_clicked();
    void on_actionCalibrate_Wavelength_triggered();
    void on_btnSetDark_clicked();
    void on_btnMeasureGo_clicked();
    void on_chkLiveSpec_toggled(bool checked);
    void on_sldSaturation_valueChanged(int value);
    void on_sldGain_valueChanged(int value);

    void setCalibEnabled(bool enable);




    void on_actionSaveFile_triggered();

    void on_actionOpenFile_triggered();

    void on_actionDeleteSelected_triggered();

    void on_btnUnsetDark_clicked();

    void on_spinNMeas_valueChanged(int arg1);

    void on_spinRepeat_valueChanged(int arg1);

    void on_chkShowBB_toggled(bool checked);

    void on_edtTempBB_editingFinished();

    void on_pushButton_clicked();




    void on_btnMultConst_clicked();

    void on_btnAddConst_clicked();

    void on_btnMathOk_clicked();

    void on_btnSaveCalX_clicked();

    void on_btnDeleteCal_clicked();

    void on_btnLoadCal_clicked();

public slots:
    void newSpectrum(Tjcampdx spectrum);
    //void newSpectrum_rgb(Tjcampdx specR,Tjcampdx specG,Tjcampdx specB);
    void addSpectrum(Tjcampdx spec);
    void bufferChange(bool vis);
    void setSpectrumRect(QRect r);
    void spectroSaturation();
    void newMeasure(Tjcampdx spec);
};

#endif // SPECTRO_H
