#include <QDebug>
#include "spectrometer.h"

Spectrometer::Spectrometer(int videoSrc)
{
    iFrame = 0;
    nFramesInteg = 1;

    measuring = false;
    //useReference = false;
    isXCalibrated = false;
    //isYCalibrated = false;

    spectrumRect = QRect(QPoint(10,100),QPoint(300,130)).normalized();
    pixData.resize(spectrumRect.width()-1);
    pixDataRed.resize(spectrumRect.width()-1);
    pixDataGreen.resize(spectrumRect.width()-1);
    pixDataBlue.resize(spectrumRect.width()-1);
    pixNum.resize(spectrumRect.width()-1);

    cap=cv::VideoCapture(videoSrc);
    if(cap.isOpened()){
        ok=true;
        connect(&timerCapture,SIGNAL(timeout()),this,SLOT(captureFrame()));
        cap.set(CV_CAP_PROP_GAIN,0.5);
        //cap.set(CV_CAP_PROP_AUTO_EXPOSURE,0);
        //cap.set(CV_CAP_PROP_FPS,10);
        cap.set(CV_CAP_PROP_SATURATION,0);
    }
    else {
        ok=false;
    }
}

void Spectrometer::setSaturation(double sat){
    cap.set(CV_CAP_PROP_SATURATION,sat);
}

void Spectrometer::setGain(double gain){
    cap.set(CV_CAP_PROP_GAIN,gain);
}

void Spectrometer::setFrameRate(int frameRate){
    timerCapture.start(1000.0/frameRate);
}

void Spectrometer::captureFrame(){
    cv::Mat frame;
    cap >> frame; // get a new frame from camera
    lastFrame = QImage((uchar*)frame.data,frame.cols,frame.rows,frame.step,QImage::Format_RGB888);
    lastFrame = lastFrame.rgbSwapped();
    getSpectrumFromImage(&lastFrame);

    emit(newFrame(&lastFrame));
    //video.setImage(qframe);
}
/*
void Spectrometer::calibrateBB(Tjcampdx spectra, double T, double max){
    Tjcampdx bbSpec = getBlackBodySpec(max,T);
    calibSpectrum = bbSpec/spectra;
    isYCalibrated=true;
}
*/
double Spectrometer::calcFFunc(double f, void *p){
    struct calcFParams * params = (struct calcFParams*)p;
    double l1 = params->l1;
    double l2 = params->l2;
    double x1 = params->x1;
    double x2 = params->x2;
    double a= params->a;
    //qDebug()<<f;
    return l2-l1-a*(x1/sqrt(f*f+x1*x1)-x2/sqrt(f*f+x2*x2));
    //return f*f-5;
}

void Spectrometer::calibrate(){
    if(calibration.size()!=2)
        return;

    double a=1e3; // grating step in nm
    double l1 = calibration[0].wavelength;
    double x1 = calibration[0].pixelNum-imageWidth/2.;
    double l2 = calibration[1].wavelength;
    double x2 = calibration[1].pixelNum-imageWidth/2.;

    struct calcFParams params = {l1,l2,x1,x2,a};

    double f_lo= 200;
    double f_hi= 3000;

    gsl_function F;
    F.function = &Spectrometer::calcFFunc;
    F.params = &params;

    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    //T = gsl_root_fsolver_brent;
    T = gsl_root_fsolver_bisection;
    s = gsl_root_fsolver_alloc(T);
    gsl_root_fsolver_set (s, &F, f_lo, f_hi); //f in pixels guessed to be between 300 and 5000

    double r=500;
    int iter=0;
    int max_iter=100;
    int status;
    do{
          iter++;
          status = gsl_root_fsolver_iterate(s);
          r = gsl_root_fsolver_root(s);
          f_lo = gsl_root_fsolver_x_lower(s);
          f_hi = gsl_root_fsolver_x_upper(s);
          status = gsl_root_test_interval (f_lo, f_hi,
                                           0, 0.001);

    } while (status == GSL_CONTINUE && iter < max_iter);

    if(status==GSL_SUCCESS) {
        focal = r;
        pas = a;
        sin_angle_i = l1/pas + x1/sqrt(r*r+x1*x1);
        isXCalibrated = true;
        qDebug()<<"Calibration finished !";
        qDebug()<<"focal : "<< focal;
        qDebug()<<"angle i : "<< 180/3.1415*asin(sin_angle_i);
        qDebug()<<"angle i2 : "<< 180/3.1415*asin( l2/pas + x2/sqrt(r*r+x2*x2));
    }

    gsl_root_fsolver_free (s);

}


double Spectrometer::blackBodyIntens(double lambda, double T){
    double hc=6.63e-34 * 3.0e8;
    double k=1.38e-23;
    double pi=3.1416;
    return 8*pi*hc/pow(lambda,5)/(exp(hc/(lambda*k*T)-1));
}

double Spectrometer::blackBodyIntensNorm(double lambda, double T){
    double max=blackBodyIntens(2.9e-3/T,T);
    return blackBodyIntens(lambda,T)/max;
}

Tjcampdx Spectrometer::getBlackBodySpec(double max, double T){
    QVector<double> lambda = pixToWavelength(pixNum);
    QVector<double> intens;
    intens.resize(lambda.size());
    for(int i=0;i<lambda.size();i++){
        intens[i]=blackBodyIntensNorm(lambda[i]*1e-9,T)*max;
    }
    Tjcampdx spectrum;
    spectrum.LoadData(lambda.toStdVector(),intens.toStdVector());
    return spectrum;
}

/*
void Spectrometer::setReferenceSpec(Tjcampdx spec){
    referenceSpec = spec;
    useReference = true;
}

void Spectrometer::setDarkSpec(Tjcampdx spec){
    darkSpec = spec;
    useDark = true;
}

void Spectrometer::unsetReference(){
    useReference = false;
}
*/

void Spectrometer::setSpectrumRect(QRect rect){
    spectrumRect = rect;
    pixData.resize(spectrumRect.width()-1);
    pixDataRed.resize(spectrumRect.width()-1);
    pixDataGreen.resize(spectrumRect.width()-1);
    pixDataBlue.resize(spectrumRect.width()-1);
    pixNum.resize(spectrumRect.width()-1);
}

Tjcampdx Spectrometer::getLastSpectrum(){
    Tjcampdx spectrum;
    spectrum.LoadData(pixNumf.toStdVector(),pixDataf.toStdVector());
    /*
    if(useDark)
        spectrum = spectrum - darkSpec;
    if(useReference)
        spectrum = spectrum/referenceSpec;
    */
    return spectrum;
}

// Start a measurement
void Spectrometer::startMeasure(int integNum){
    qDebug()<<"Starting measure "<<integNum;
    measuring = true;
    iMeasure = 0;
    nMeasure = integNum;
    pixDataMeasure.resize(pixData.size());
    pixDataMeasure.fill(0);
}

void Spectrometer::getSpectrumFromImage(QImage* image){
    if(iFrame==nFramesInteg){
        iFrame = 0;
        pixDataf = pixData;
        pixNumf = pixNum;
        pixData.fill(0);
        pixNum.fill(0);

        Tjcampdx spec,spec_r,spec_g,spec_b;
        spec.LoadData(pixToWavelength(pixNumf).toStdVector(),pixDataf.toStdVector());
        spec_r.LoadData(pixToWavelength(pixNumf).toStdVector(),pixDataRed.toStdVector());
        spec_g.LoadData(pixToWavelength(pixNumf).toStdVector(),pixDataGreen.toStdVector());
        spec_b.LoadData(pixToWavelength(pixNumf).toStdVector(),pixDataBlue.toStdVector());
        /*
        if(useDark)
            spec = spec-darkSpec;
        if(isYCalibrated)
            spec = spec*calibSpectrum;
        if(useReference)
            spec = spec/referenceSpec;
        */
        //emit(newSpectrum_rgb(spec_r,spec_g,spec_b));
        emit(newSpectrum(spec));

        pixDataRed.fill(0);
        pixDataGreen.fill(0);
        pixDataBlue.fill(0);
    }

    if(measuring && (iMeasure == nMeasure)){
        Tjcampdx spec;
        spec.LoadData(pixToWavelength(pixNumf).toStdVector(),pixDataMeasure.toStdVector());
        emit(measureFinished(spec));
        measuring = false;
    }

    bool sat=false;
    for(int iLine=spectrumRect.left();iLine<spectrumRect.right();iLine++){
        double total = 0;
        double total_r = 0;
        double total_g = 0;
        double total_b = 0;

        for(int iHeight=spectrumRect.top();iHeight<spectrumRect.bottom();iHeight++){
            QRgb pixel = image->pixel(iLine,iHeight);
            total_r+=qRed(pixel);
            total_g+=qGreen(pixel);
            total_b+=qBlue(pixel);
            total+=qBlue(pixel) + qRed(pixel) + qGreen(pixel);
            sat |= ((qRed(pixel)==255) || (qGreen(pixel)==255) || (qBlue(pixel)==255));
        }
        if(measuring){
            pixDataMeasure[iLine-spectrumRect.left()]+=total/spectrumRect.height()/nMeasure;
            emit(measureProgressed((iMeasure+1)*100/nMeasure));
        }
        pixData[iLine-spectrumRect.left()]+=total/spectrumRect.height()/nFramesInteg;
        pixDataRed[iLine-spectrumRect.left()]+=total_r/spectrumRect.height()/nFramesInteg;
        pixDataGreen[iLine-spectrumRect.left()]+=total_g/spectrumRect.height()/nFramesInteg;
        pixDataBlue[iLine-spectrumRect.left()]+=total_b/spectrumRect.height()/nFramesInteg;
        pixNum[iLine-spectrumRect.left()]=iLine;
    }
    if(sat){
        emit(saturation());
    }
    iFrame++;
    if(measuring)
        iMeasure++;
}

// Convert one pixel to wavelength
double Spectrometer::pixNumToWaveLength(double pixnum){
    if(!isXCalibrated)
        return pixnum;
    else {
        double x = pixnum-imageWidth/2.;
        return pas*(sin_angle_i - x/sqrt(focal*focal+x*x));
    }
}

// Convert pixels to wavelength
QVector<double> Spectrometer::pixToWavelength(QVector<double> pix){
    QVector<double> wl;
    wl.resize(pix.size());
    for(int i=0;i<pix.count();i++){
        wl[i] = pixNumToWaveLength(pix[i]);
    }
    return wl;
}

void Spectrometer::resetCalibration(){
    isXCalibrated = false;
    calibration.clear();
}

void Spectrometer::addCalibrationPoint(CalibrationPoint point){
    calibration.append(point);
}
