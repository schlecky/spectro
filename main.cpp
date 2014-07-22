#include "spectro.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    QCoreApplication::setOrganizationName("Renaud Schleck");
    QCoreApplication::setOrganizationDomain("renaud.schleck.free.fr");
    QCoreApplication::setApplicationName("Spectrometer");
    Spectro w;
    w.show();
    
    return a.exec();
}
