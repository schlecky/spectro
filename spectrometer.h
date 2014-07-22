#ifndef SPECTROMETER_H
#define SPECTROMETER_H

#include <QObject>
#include <QRect>
#include <QTimer>
#include <QImage>
#include <QVector>
#include <QPainter>

#include "jcampdx.h"

#include <opencv2/opencv.hpp>

#include <gsl/gsl_roots.h>

typedef struct{
    int pixelNum;
    double wavelength;
} CalibrationPoint;

class Spectrometer : public QObject
{
    Q_OBJECT
public:
    explicit Spectrometer(int videoSrc=1);
    Tjcampdx getLastSpectrum();
    Tjcampdx getBlackBodySpec(double max,double T);
    void resetCalibration();
    void setCalibration(QVector<CalibrationPoint> cal) {
        calibration = cal;
        calibrate();
    }

    bool isWlCalib(){
        return isXCalibrated;
    }

    double getFocal() {return focal;}
    double getSin() {return sin_angle_i;}
    double getPas() {return pas;}

    void addCalibrationPoint(CalibrationPoint);
    bool isOk() {return ok;}


    
private slots:
    void captureFrame();

public slots:
    void setSpectrumRect(QRect rect);
    QRect getSpectrumRect(){return spectrumRect;}

    /*
    void setReferenceSpec(Tjcampdx spec);
    void unsetReference();

    void setDarkSpec(Tjcampdx spec);
    void resetDark(){useDark = false;}
    */

    void setSaturation(double sat);
    void setGain(double gain);

    /*
    void calibrateBB(Tjcampdx spectra, double T, double max);
    void resetYCalibration(){
        isYCalibrated = false;
    }
    */

    void setFrameRate(int frameRate);
    double getFrameRate(){
        if(timerCapture.interval()>0)
            return 1000/timerCapture.interval();
        else
            return 15;
    }

    void setIntegNum(int num){
        iFrame = 0;
        nFramesInteg = num;
    }

    void startMeasure(int integNum=1);

    void setParameters(double f,double s,double d){
      focal = f;
      sin_angle_i = s;
      pas = d;
      isXCalibrated = true;
    }


    
signals:
    void newSpectrum(Tjcampdx spectrum);
    void measureFinished(Tjcampdx spectrum);
    void newSpectrum_rgb(Tjcampdx spectrum_r,Tjcampdx spectrum_g,Tjcampdx spectrum_b);
    void newFrame(QImage* frame);
    void saturation();
    void measureProgressed(int progress);

private:
    QRect spectrumRect;

    QTimer timerCapture;
    cv::VideoCapture cap;

    QImage lastFrame;
    //Frame integration
    int nFramesInteg;               // frames to integrate for one spectrum
    int iFrame;                     // index of current frame
    int nMeasure;
    int iMeasure;

    // Calibration
    // parametres pour le calcul de f
    struct calcFParams {double l1; double l2; double x1; double x2;double a;};
    static double calcFFunc(double f,void* params);
    void calibrate();
    QVector<CalibrationPoint> calibration;
    double focal;   // Distance focale image de la lentille
    double sin_angle_i; // Angle d'incidence sur le réseau
    double pas;     // Pas du réseau

    double imageWidth=640;

    bool measuring;
    //bool useReference;
    //bool useDark;
    bool isXCalibrated;
    //bool isYCalibrated;
    bool ok;

    //Tjcampdx referenceSpec;
    //Tjcampdx darkSpec;
    //Tjcampdx calibSpectrum;

    QVector<double> pixDataMeasure;
    QVector<double> pixData,pixDataf;
    QVector<double> pixDataRed,pixDataRedf;
    QVector<double> pixDataGreen,pixDataGreenf;
    QVector<double> pixDataBlue,pixDataBluef;
    QVector<double> pixNum,pixNumf;

    QVector<double> pixToWavelength(QVector<double> pix);
    double pixNumToWaveLength(double pixnum);
    void getSpectrumFromImage(QImage* image);
    double blackBodyIntens(double lambda, double T);
    double blackBodyIntensNorm(double lambda, double T);
};

#endif // SPECTROMETER_H
