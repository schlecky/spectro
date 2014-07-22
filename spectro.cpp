#include "spectro.h"
#include "ui_spectro.h"
#include "videoview.h"
#include <opencv2/opencv.hpp>
#include <QDebug>
#include <QRgb>
#include <QImage>
#include <QFileDialog>
//#include <opencv2/highgui/highgui.hpp>

//using namespace cv;
Spectro::Spectro(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::Spectro)
{
    ui->setupUi(this);

    divData = false;
    calibration = false;

    colors.append(Qt::black);
    colors.append(Qt::red);
    colors.append(Qt::darkGreen);
    colors.append(Qt::darkRed);
    colors.append(Qt::blue);
    colors.append(Qt::darkCyan);
    colors.append(Qt::darkMagenta);
    colors.append(Qt::darkYellow);
    colors.append(Qt::darkBlue);
    currentColor = 0;

    ui->tabWidget->removeTab(ui->tabWidget->indexOf(ui->tabWidgetScale));
    ui->tabWidget->removeTab(ui->tabWidget->indexOf(ui->tabWidgetMaths));
    ui->tabWidget->removeTab(ui->tabWidget->indexOf(ui->tabWidgetCalib));
    ui->tabWidget->hide();

    spectro = new Spectrometer(1);


    QSettings settings;
    setCalibEnabled(true);

   /*
    if(!settings.value("calibration/focal").isNull()){
        double focal = settings.value("calibration/focal").toDouble();
        double sin_i = settings.value("calibration/sin_i").toDouble();
        double pas = settings.value("calibration/pas").toDouble();
        spectro->setParameters(focal,sin_i,pas);
        ui->tabWidgetCalib->setEnabled(false);
        if(spectro->isWlCalib())
            setCalibEnabled(false);
    }
    */

    readCalibrationList();
    if(!settings.value("calibration_num",0).isNull()){
        int num = settings.value("calibration_num").toInt();
        qDebug()<<"num : "<<num;
        if((num>=0) && (num <calibrationList.size())){
            qDebug()<<"calibrating";
            Calibration cal = calibrationList[num];
            spectro->setParameters(cal.focal,cal.sin_i,1e3);
            ui->tabWidgetCalib->setEnabled(false);
            if(spectro->isWlCalib())
                setCalibEnabled(false);
        }
    }


    lastDir = settings.value("application/fileDir",QApplication::applicationDirPath()).toString();

    if(!settings.value("spectro/x1").isNull()){
        double x1 = settings.value("spectro/x1").toInt();
        double y1 = settings.value("spectro/y1").toInt();
        double x2 = settings.value("spectro/x2").toInt();
        double y2 = settings.value("spectro/y2").toInt();
        setSpectrumRect(QRect(QPoint(x1,y1),QPoint(x2,y2)).normalized());
    }

    if(!spectro->isOk())
        close();
    //connect(&video,SIGNAL(rectSelected(QRect)),spectro,SLOT(setSpectrumRect(QRect)));
    connect(&video,SIGNAL(rectSelected(QRect)),this,SLOT(setSpectrumRect(QRect)));

    connect(spectro,SIGNAL(newSpectrum(Tjcampdx)),this,SLOT(newSpectrum(Tjcampdx)));
    //connect(spectro,SIGNAL(newSpectrum_rgb(Tjcampdx,Tjcampdx,Tjcampdx)),this,SLOT(newSpectrum_rgb(Tjcampdx,Tjcampdx,Tjcampdx)));
    connect(spectro,SIGNAL(saturation()),this,SLOT(spectroSaturation()));
    connect(spectro,SIGNAL(measureFinished(Tjcampdx)),this,SLOT(newMeasure(Tjcampdx)));

    connect(spectro,SIGNAL(measureProgressed(int)),ui->progressMeasure,SLOT(setValue(int)));

    connect(spectro,SIGNAL(newFrame(QImage*)),&video,SLOT(setImage(QImage*)));
    connect(ui->plot->xAxis,SIGNAL(rangeChanged(QCPRange)),this,SLOT(xScaleChanged(QCPRange)));
    connect(ui->plot->yAxis,SIGNAL(rangeChanged(QCPRange)),this,SLOT(yScaleChanged(QCPRange)));
    connect(ui->actionShowVideo,SIGNAL(toggled(bool)),&video,SLOT(setVisible(bool)));
    connect(ui->plot,SIGNAL(mousePress(QMouseEvent*)),this,SLOT(plotClicked(QMouseEvent*)));
    connect(qApp,SIGNAL(focusChanged(QWidget*,QWidget*)),this,SLOT(focusChanged(QWidget*,QWidget*)));

    // graph for live spectrum
    liveSpec = ui->plot->addGraph();

    /*
    liveSpecR = ui->plot->addGraph();
    liveSpecR->setPen(QPen(Qt::red));
    liveSpecR->setVisible(false);

    liveSpecG = ui->plot->addGraph();
    liveSpecG->setPen(QPen(Qt::green));
    liveSpecG->setVisible(false);


    liveSpecB = ui->plot->addGraph();
    liveSpecB->setPen(QPen(Qt::blue));
    liveSpecB->setVisible(false);
    */
    bbSpec = ui->plot->addGraph();
    bbSpec->setPen(QPen(Qt::black));
    bbSpec->setVisible(false);


    wlLine = new QCPItemStraightLine(ui->plot);

    ui->plot->addItem(wlLine);
    wlLine->setPen(QPen(Qt::red));


    // add the text label at the top:
    wlLabel = new QCPItemText(ui->plot);
    ui->plot->addItem(wlLabel);
    wlLabel->setPositionAlignment(Qt::AlignTop|Qt::AlignRight);
    wlLabel->setClipToAxisRect(false);
    wlLabel->position->setType(QCPItemPosition::ptAxisRectRatio);
    wlLabel->position->setCoords(1, 0); // place position at right/top of axis rect
    wlLabel->setText("\u03BB = ??? nm");

    connect(ui->plot,SIGNAL(mouseMove(QMouseEvent*)),this,SLOT(plotMouseMoved(QMouseEvent*)));


    if(!settings.value("scale/xmin").isNull()){
        double xmin = settings.value("scale/xmin",400).toDouble();
        double xmax = settings.value("scale/xmax",1000).toDouble();
        double ymin = settings.value("scale/ymin",0).toDouble();
        double ymax = settings.value("scale/ymax",700).toDouble();
        changeScale(xmin,xmax,ymin,ymax);
    }

    ui->spinNMeas->setValue(settings.value("application/nAverage",1).toInt());
    ui->spinRepeat->setValue(settings.value("application/nRepeat",1).toInt());

    ui->sldGain->setValue(settings.value("spectro/gain",0).toInt());
    spectro->setGain(ui->sldGain->value()/100.);

    /*
    QFileInfo gainFile = QFileInfo(QFileInfo(settings.fileName()).absolutePath()+
                                   QString("dark_g%1.dx").arg(ui->sldGain->value()));
    if(gainFile.exists())
    {
        Tjcampdx darkSpec;
        darkSpec.Open(gainFile.absoluteFilePath().toStdString());
        spectro->setDarkSpec(darkSpec);
    }
    */

    spectro->setFrameRate(15);
    video.hide();
}

Spectro::~Spectro(){
    QSettings settings;
    settings.setValue("scale/xmin",ui->plot->xAxis->range().lower);
    settings.setValue("scale/xmax",ui->plot->xAxis->range().upper);
    settings.setValue("scale/ymin",ui->plot->yAxis->range().lower);
    settings.setValue("scale/ymax",ui->plot->yAxis->range().upper);

    settings.setValue("application/fileDir",lastDir);

    int nMeas=ui->spinNMeas->value();
    if(nMeas!=0)
        settings.setValue("application/nAverage",nMeas);

    int r=ui->spinRepeat->value();
    if(r!=0)
        settings.setValue("application/nRepeat",r);

    settings.setValue("spectro/gain",ui->sldGain->value());
    delete ui;
}


void Spectro::setSpectrumRect(QRect r){
    spectro->setSpectrumRect(r);
    QSettings settings;
    settings.setValue("spectro/x1",r.left());
    settings.setValue("spectro/y1",r.top());
    settings.setValue("spectro/x2",r.right());
    settings.setValue("spectro/y2",r.bottom());
}

void Spectro::plotMouseMoved(QMouseEvent *event){
    double wl = ui->plot->graph()->keyAxis()->pixelToCoord(event->pos().x());
    wlLine->point1->setCoords(wl,0);
    wlLine->point2->setCoords(wl,1);
    wlLabel->setText(QString("\u03BB = %1 nm").arg(QString::number(wl,'f',2)));
    ui->plot->replot();
}

void Spectro::focusChanged(QWidget *old, QWidget *now){
    if(now == ui->edtPt1Pix)
        currentState = SELECT_POINT1;
    else if(now == ui->edtPt2Pix)
            currentState = SELECT_POINT2;
    else
            currentState = NONE;
    qDebug()<<currentState;
}

void Spectro::resizeEvent(QResizeEvent *event){

}

void Spectro::changeScale(double xmin, double xmax, double ymin, double ymax){
    ui->plot->xAxis->setRange(xmin,xmax);
    ui->plot->yAxis->setRange(ymin,ymax);
}

void Spectro::addSpectrum(Tjcampdx spec){
    QColor color = colors[currentColor];
    ui->bufferBar->addBuffer(color,QString("Spectrum"),currentColor);

    spectra.append(spec);
    graphs.append(ui->plot->addGraph());
    graphs.back()->setData(QVector<double>::fromStdVector(spec.DataX()),
                               QVector<double>::fromStdVector(spec.DataY()));
    graphs.back()->setPen(colors[currentColor]);

    ui->plot->replot();
    currentColor = (currentColor + 1) % colors.size();
}

void Spectro::bufferChange(bool vis){

}

void Spectro::spectroSaturation(){
    ui->statusBar->showMessage("Saturation !",100);
}

void Spectro::newMeasure(Tjcampdx spec){
    addSpectrum(spec);
    ui->progressMeasure->setEnabled(false);
}

void Spectro::newSpectrum(Tjcampdx spectrum){
    lastSpectrum = spectrum;
    liveSpec->setData(QVector<double>::fromStdVector(spectrum.DataX()),
                                QVector<double>::fromStdVector(spectrum.DataY()));
    ui->plot->replot();
}

/*
void Spectro::newSpectrum_rgb(Tjcampdx specR,Tjcampdx specG,Tjcampdx specB){
    liveSpecR->setData(QVector<double>::fromStdVector(specR.DataX()),
                       QVector<double>::fromStdVector(specR.DataY()));
    liveSpecG->setData(QVector<double>::fromStdVector(specG.DataX()),
                       QVector<double>::fromStdVector(specG.DataY()));
    liveSpecB->setData(QVector<double>::fromStdVector(specB.DataX()),
                       QVector<double>::fromStdVector(specB.DataY()));
    ui->plot->replot();
}
*/

void Spectro::xScaleChanged(QCPRange range){
    ui->edtXmin->setText(QString("%1").arg(range.lower));
    ui->edtXmax->setText(QString("%1").arg(range.upper));
}

void Spectro::yScaleChanged(QCPRange range){
    ui->edtYmin->setText(QString("%1").arg(range.lower));
    ui->edtYmax->setText(QString("%1").arg(range.upper));
}

void Spectro::on_actionHideSelectedItems_triggered()
{
    QList<QTableWidgetItem *> selected = ui->bufferBar->selectedItems();
    for(int i=0; i<selected.size();i++){
      QTableWidgetItem* item=selected[i];
      int row=item->row();
      int column = item->column();
      graphs[row]->setVisible(false);
      ui->bufferBar->item(row,column)->setBackground(Qt::white);
    }
    ui->plot->replot();
}

void Spectro::on_actionShowSelectedItems_triggered()
{
    QList<QTableWidgetItem *> selected = ui->bufferBar->selectedItems();
    for(int i=0; i<selected.size();i++){
      QTableWidgetItem* item=selected[i];
      int row=item->row();
      int column = item->column();
      graphs[row]->setVisible(true);
      ui->bufferBar->item(row,column)->setBackground(colors[item->data(QBufferBar::COLOR_INDEX).toInt()]);
    }
    ui->plot->replot();
}



void Spectro::on_tabWidget_tabCloseRequested(int index)
{
    ui->tabWidget->removeTab(index);
    if(ui->tabWidget->count()==0)
        ui->tabWidget->hide();
}


void Spectro::on_bufferBar_cellClicked(int row, int column)
{
    if(!selectedForOp.isEmpty())
    {
      Tjcampdx other = Tjcampdx(spectra[row]);
      for(int i=0; i<selectedForOp.size();i++)
      {
          QTableWidgetItem* item=selectedForOp[i];
          Tjcampdx current = Tjcampdx(spectra[item->row()]);
          if(divData)
            current = current / other;
          else if(addData)
            current = current + other;
          else if(subData)
            current = current - other;
          addSpectrum(current);
      }

      ui->btnDivBuffer->setChecked(false);
      ui->btnAddBuffer->setChecked(false);
      ui->btnSubBuffer->setChecked(false);
      ui->plot->replot();
      selectedForOp.clear();
    }
}


void Spectro::on_btnMathOk_clicked()
{
    if(!selectedForOp.isEmpty())
    {
      for(int i=0; i<selectedForOp.size();i++)
      {
          QTableWidgetItem* item=selectedForOp[i];
          Tjcampdx current = Tjcampdx(spectra[item->row()]);
          if(addConst)
              current = current + ui->edtMathConstant->text().toDouble();
          else if(multConst)
              current = current * ui->edtMathConstant->text().toDouble();
          addSpectrum(current);
      }
      ui->btnAddConst->setChecked(false);
      ui->btnMultConst->setChecked(false);

      ui->lblMathConst->setEnabled(false);
      ui->btnMathOk->setEnabled(false);
      ui->edtMathConstant->setEnabled(false);

      ui->plot->replot();
      selectedForOp.clear();
    }
}

void Spectro::on_btnDivBuffer_clicked()
{
    divData = !divData;
    if(divData)
    {
      selectedForOp = ui->bufferBar->selectedItems();
      if(!selectedForOp.isEmpty()){
        divData = true;
        uncheckMathButtons();
      }
      else
        divData = false;
    }
    ui->btnDivBuffer->setChecked(divData);
    subData = addConst = multConst = addData = false;
}

void Spectro::on_btnAddBuffer_clicked()
{
    addData = ui->btnAddBuffer->isChecked();
    if(addData)
    {
      selectedForOp = ui->bufferBar->selectedItems();
      if(!selectedForOp.isEmpty())
      {
        addData = true;
        uncheckMathButtons();
      }
      else{
            addData = false;
      }
    }
    subData = addConst = multConst = divData = false;
    ui->btnAddBuffer->setChecked(addData);
}

void Spectro::on_btnSubBuffer_clicked()
{
    subData = ui->btnSubBuffer->isChecked();
    if(subData)
    {
      selectedForOp = ui->bufferBar->selectedItems();
      if(!selectedForOp.isEmpty())
      {
        subData = true;
        uncheckMathButtons();
      }
      else{
            subData = false;
      }
    }
    addData = addConst = multConst = divData = false;
    ui->btnSubBuffer->setChecked(subData);
}

void Spectro::on_btnMultConst_clicked()
{
    multConst = ui->btnMultConst->isChecked();
    if(multConst)
    {
      selectedForOp = ui->bufferBar->selectedItems();
      if(!selectedForOp.isEmpty())
      {
        multConst = true;
        uncheckMathButtons();
      }
      else{
            multConst = false;
      }
    }
    ui->edtMathConstant->setEnabled(multConst);
    ui->lblMathConst->setEnabled(multConst);
    ui->btnMathOk->setEnabled(multConst);
    addData = addConst = subData = divData = false;
    ui->btnMultConst->setChecked(multConst);
}

void Spectro::on_btnAddConst_clicked()
{
    addConst = ui->btnAddConst->isChecked();
    if(addConst)
    {
      selectedForOp = ui->bufferBar->selectedItems();
      if(!selectedForOp.isEmpty())
      {
        addConst = true;
        uncheckMathButtons();
      }
      else{
            addConst = false;
      }
    }
    ui->edtMathConstant->setEnabled(addConst);
    ui->lblMathConst->setEnabled(addConst);
    ui->btnMathOk->setEnabled(addConst);
    addData = multConst = subData = divData = false;
    ui->btnAddConst->setChecked(addConst);
}


void Spectro::uncheckMathButtons(){
    ui->btnAddBuffer->setChecked(false);
    ui->btnDivBuffer->setChecked(false);
    ui->btnAddConst->setChecked(false);
    ui->btnSubBuffer->setChecked(false);
    ui->btnMultConst->setChecked(false);
}





/*
void Spectro::on_actionSubBuffer_triggered()
{
    if(subData)
    {
      subData=false;
    }
    else
    {
      selectedForOp = ui->bufferBar->selectedItems();
      if(!selectedForOp.isEmpty())
        subData = true;
      else
        subData = false;
    }
    ui->actionSubBuffer->setChecked(subData);
    ui->btnSubBuffer->setChecked(subData);
}
*/
void Spectro::on_actionShowScale_triggered()
{
    ui->tabWidget->show();
    if(ui->tabWidget->indexOf(ui->tabWidgetScale)==-1){
        ui->tabWidget->insertTab(0,ui->tabWidgetScale,QString("Scale"));
    }
    ui->tabWidget->setCurrentWidget(ui->tabWidgetScale);
}

void Spectro::on_actionShowMath_triggered()
{
    ui->tabWidget->show();
    if(ui->tabWidget->indexOf(ui->tabWidgetMaths)==-1){
        ui->tabWidget->insertTab(0,ui->tabWidgetMaths,QString("Maths"));
    }
    ui->tabWidget->setCurrentWidget(ui->tabWidgetMaths);
}

void Spectro::on_btnOkScale_clicked()
{
    changeScale(ui->edtXmin->text().toFloat(),
                ui->edtXmax->text().toFloat(),
                ui->edtYmin->text().toFloat(),
                ui->edtYmax->text().toFloat());
}

void Spectro::on_btnAutoScale_clicked()
{
    ui->plot->rescaleAxes(true);
}


void Spectro::on_btnSetRef_clicked()
{
    QList<QTableWidgetItem*> selection = ui->bufferBar->selectedItems();
    if(selection.count()!=1){
        ui->statusBar->showMessage(QString("Select exactly one spectrum as reference"));
    }
    else {
        ui->lblRef->setText(QString("Reference : yes"));
        //spectro->setReferenceSpec(spectra[selection[0]->row()]);
    }

}


void Spectro::on_btnUnsetRef_clicked()
{
    //spectro->unsetReference();
    ui->lblRef->setText(QString("Reference : no"));
}


void Spectro::on_actionSpectro_triggered()
{
    ui->tabWidget->show();
    if(ui->tabWidget->indexOf(ui->tabWidgetSpectro)==-1){
        ui->tabWidget->insertTab(0,ui->tabWidgetSpectro,QString("Spectro"));
    }
    ui->tabWidget->setCurrentWidget(ui->tabWidgetSpectro);

}

/*
void Spectro::on_btnCalibrate_clicked()
{
    if(calibration){
        calibration = false;
        ui->btnCalibrateX->setText("Calibrate");
    }

    else{
        calibration = true;
        ui->btnCalibrateX->setText("End Calibr.");
    }

}
*/

void Spectro::plotClicked(QMouseEvent *event){
    if(currentState != NONE){
        int pix;
        pix = ui->plot->graph()->keyAxis()->pixelToCoord(event->pos().x());
        if(currentState == SELECT_POINT1){
            ui->edtPt1Pix->setText(QString("%1").arg(pix));
        }
        else if(currentState == SELECT_POINT2){
                ui->edtPt2Pix->setText(QString("%1").arg(pix));
        }
    }

}


void Spectro::on_btnCalibrateX_clicked()
{
    // are we calibrating
    if(ui->edtPt1Pix->isEnabled()){
        bool ok;
        bool ready=true;
        QVector<CalibrationPoint> cal;
        CalibrationPoint point;
        point.pixelNum = ui->edtPt1Pix->text().toInt(&ok);
        ready = ready && ok;
        point.wavelength = ui->edtPt1Wl->text().toDouble(&ok);
        ready = ready && ok;
        cal.append(point);

        point.pixelNum = ui->edtPt2Pix->text().toInt(&ok);
        ready = ready && ok;
        point.wavelength = ui->edtPt2Wl->text().toDouble(&ok);
        ready = ready && ok;
        cal.append(point);

        if(ready){
            spectro->setCalibration(cal);
            if(spectro->isWlCalib()){
                QSettings settings;
                settings.setValue("calibration/focal",spectro->getFocal());
                settings.setValue("calibration/sin_i",spectro->getSin());
                settings.setValue("calibration/pas",spectro->getPas());
                setCalibEnabled(true);
                //ui->tabWidgetCalib->setEnabled(false);
                setCalibEnabled(false);
                ui->statusBar->showMessage("Calibration ok.");
            }
            else{
                ui->statusBar->showMessage("Calibration failed !");
            }
        }
        else {
            ui->statusBar->showMessage("More calibration data needed.");
        }
    }
    else{
        setCalibEnabled(true);
        spectro->resetCalibration();
    }


}

void Spectro::on_actionCalibrate_Wavelength_triggered()
{
    ui->tabWidget->show();
    if(ui->tabWidget->indexOf(ui->tabWidgetCalib)==-1){
        ui->tabWidget->insertTab(0,ui->tabWidgetCalib,QString("Calibration"));
    }
    ui->tabWidget->setCurrentWidget(ui->tabWidgetCalib);
    ui->tabWidgetCalib->setEnabled(true);
}


void Spectro::on_btnSetDark_clicked()
{
    QList<QTableWidgetItem*> selection = ui->bufferBar->selectedItems();
    if(selection.count()!=1){
        ui->statusBar->showMessage(QString("Select exactly one spectrum as dark"));
    }
    else {
        QSettings settings;
        ui->lblDark->setText(QString("Dark : yes"));
        //spectro->setDarkSpec(spectra[selection[0]->row()]);
        //spectra[selection[0]->row()].Save((QFileInfo(settings.fileName()).absolutePath()+
        //                                   QString("dark_g%1.dx").arg(ui->sldGain->value())).toStdString());
    }
}


void Spectro::on_btnMeasureGo_clicked()
{
    int nMeasure = ui->spinNMeas->value();
    if(nMeasure>0){
        spectro->startMeasure(nMeasure);
        ui->progressMeasure->setEnabled(true);
        ui->progressMeasure->setValue(0);
    }
}


void Spectro::on_chkLiveSpec_toggled(bool checked)
{
    liveSpec->setVisible(checked);
}


void Spectro::on_sldSaturation_valueChanged(int value)
{
    double sat = (double)value / ui->sldSaturation->maximum();
    spectro->setSaturation(sat);
}


void Spectro::setCalibEnabled(bool enable){
    ui->edtPt1Pix->setEnabled(enable);
    ui->edtPt2Pix->setEnabled(enable);
    ui->edtPt1Wl->setEnabled(enable);
    ui->edtPt2Wl->setEnabled(enable);
    ui->btnSaveCalX->setEnabled(!enable);
    if(enable)
        ui->btnCalibrateX->setText("Accept Calib.");
    else
        ui->btnCalibrateX->setText("Start Calib.");
}


void Spectro::on_sldGain_valueChanged(int value)
{
    double gain = (double)value / ui->sldGain->maximum();
    spectro->setGain(gain);
}


void Spectro::on_actionSaveFile_triggered()
{
    QList<QTableWidgetItem*> selection = ui->bufferBar->selectedItems();
    if(selection.count()!=1){
        ui->statusBar->showMessage(tr("Select exactly one spectrum to save"));
    }
    else{
        QString filename = QFileDialog::getSaveFileName(this,tr("Save spectrum as"),lastDir,tr("JCampDX (*.dx)"));
        if(!filename.isEmpty()){
            spectra[selection[0]->row()].Save(filename.toStdString());
            lastDir=QFileInfo(filename).absolutePath();
        }
    }
}

void Spectro::on_actionOpenFile_triggered()
{
    QStringList filenames = QFileDialog::getOpenFileNames(this,"Open spectrum file",lastDir,tr("JCampDX (*.dx)"));
    if(!filenames.isEmpty())
    {
        for(int i=0; i<filenames.length();i++){
            Tjcampdx spectrum;
            spectrum.Open(filenames[i].toStdString());
            addSpectrum(spectrum);
        }
        lastDir=QFileInfo(filenames.back()).absolutePath();
    }
}

void Spectro::on_actionDeleteSelected_triggered()
{
    QList<QTableWidgetItem*> selection = ui->bufferBar->selectedItems();
    for(int i=0; i<selection.size();i++)
    {
        QTableWidgetItem* item=selection[i];
        int iSpectrum = item->row();
        ui->bufferBar->removeRow(iSpectrum);
        spectra.takeAt(iSpectrum);
        ui->plot->removeGraph(graphs[iSpectrum]);
        graphs.remove(iSpectrum);
        //ui->plot->removeGraph(ui->plot->griSpectrum);
    }
    ui->plot->replot();
    qDebug()<<"item deleted";
}

void Spectro::on_btnUnsetDark_clicked()
{
    //spectro->resetDark();
    ui->lblDark->setText(tr("Dark : no"));
}

void Spectro::on_spinNMeas_valueChanged(int arg1)
{
    ui->lblEstTime->setText(QString::number(ui->spinRepeat->value()*arg1/spectro->getFrameRate(),'f',2)+QString(" s"));
}

void Spectro::on_spinRepeat_valueChanged(int arg1)
{
    ui->lblEstTime->setText(QString::number(ui->spinNMeas->value()*arg1/spectro->getFrameRate(),'f',2)+QString(" s"));
}

/*
void Spectro::on_btnCalibrateY_clicked()
{
    QList<QTableWidgetItem*> selection = ui->bufferBar->selectedItems();
    if(selection.count()!=1){
        ui->statusBar->showMessage(tr("Select exactly one spectrum for calibration"));
    }
    else{
        spectro->calibrateBB(spectra[selection[0]->row()],
                ui->edtTempBB->text().toDouble(),
                ui->edtMaxBB->text().toDouble());
    }
}
*/

void Spectro::on_chkShowBB_toggled(bool checked)
{
    bbSpec->setVisible(checked);
}

void Spectro::on_edtTempBB_editingFinished()
{
    bool ok1,ok2;
    //bbSpec->setVisible(true);
    double T = ui->edtTempBB->text().toDouble(&ok1);
    double max= ui->edtMaxBB->text().toDouble(&ok2);
    if(ok1 && ok2){
        Tjcampdx spec = spectro->getBlackBodySpec(max,T);
        bbSpec->setData(QVector<double>::fromStdVector(spec.DataX()),
                        QVector<double>::fromStdVector(spec.DataY()));
        ui->plot->replot();
    }
}

void Spectro::on_pushButton_clicked()
{
    bool ok1,ok2;
    //bbSpec->setVisible(true);
    double T = ui->edtTempBB->text().toDouble(&ok1);
    double max= ui->edtMaxBB->text().toDouble(&ok2);
    if(ok1 && ok2){
        Tjcampdx spec = spectro->getBlackBodySpec(max,T);
        addSpectrum(spec);
    }
    ui->plot->replot();
}


void Spectro::readCalibrationList(){
    QSettings settings;
    int size = settings.beginReadArray("calibrations");
    for(int i=0;i<size; i++){
        settings.setArrayIndex(i);
        Calibration calib;
        calib.sin_i = settings.value("sin_angle_i",40).toDouble();
        calib.comment = settings.value("comment","RIEN").toString();
        calib.focal = settings.value("focal",10).toDouble();
        calib.spectrumRect = settings.value("specRect").toRect();
        calib.name = settings.value("name").toString();
        calibrationList.append(calib);
        ui->cmbCalibrations->addItem(calib.name,i);
    }
    settings.endArray();
}


void Spectro::on_btnSaveCalX_clicked()
{
    if(spectro->isWlCalib()){
        QSettings settings;
        settings.beginWriteArray("calibrations");
        settings.setArrayIndex(calibrationList.size());
        settings.setValue("sin_angle_i",spectro->getSin());
        settings.setValue("focal",spectro->getFocal());
        QString name = QInputDialog::getText(this,"Calibration name","name");
        settings.setValue("name",name);
        QString comment = QInputDialog::getText(this,"Set comment","comment");
        settings.setValue("comment",comment);
        settings.setValue("specRect",spectro->getSpectrumRect());
        settings.endArray();
    }
}



void Spectro::on_btnDeleteCal_clicked()
{

}

void Spectro::on_btnLoadCal_clicked()
{
    if(!calibrationList.isEmpty()){
        if(ui->cmbCalibrations->currentIndex()<calibrationList.size()){
            Calibration cal = calibrationList[ui->cmbCalibrations->currentIndex()];
            spectro->setParameters(cal.focal,cal.sin_i,cal.pas);
        }
        if(!spectro->isWlCalib()){
            ui->statusBar->showMessage("Calibration failed");
        }
    }
}
