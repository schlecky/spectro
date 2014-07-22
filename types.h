//----------------------------------------------------------------------------------------
// types -- Define COMPLEX, TPointDouble and TRectDouble types
// (c) 1998-2009 R.P.S.M. Lobo - lobo@espci.fr
// File created: 1998/01/01
// Last update: 2010/09/27
// namespace "lobo"
//
// Redefines type COMPLEX;
//
// Class lobo::TPointDouble
//   public:
//     TPointDouble();
//     TPointDouble(double x, double y);
//     TPointDouble(const TPointDouble &P);
//
//     inline double &x();
//     inline const double &x() const;
//     inline double &y();
//     inline const double &y() const;
//     void Assign(double x, double y);
//
//     const TPointDouble &operator=(const TPointDouble &P);
//
//     friend TPointDouble operator+(const TPointDouble &P1, const TPointDouble &P2);
//     friend TPointDouble operator-(const TPointDouble &P1, const TPointDouble &P2);
//     friend std::ostream &operator<<(std::ostream &os, const TPointDouble &P);
//     friend std::istream &operator>>(std::istream &is, TPointDouble &P);
//
// Class lobo::TRectDouble
//   public:
//     TRectDouble();
//     TRectDouble(double x1, double y1, double x2, double y2);
//     TRectDouble(const TRectDouble &R);
// 
//     inline double &x1();
//     inline const double &x1() const;
//
//     inline double &x2();
//     inline const double &x2() const;
// 
//     inline double &y1();
//     inline const double &y1() const;
//
//     inline double &y2();
//     inline const double &y2() const;
// 
//     void Assign(double x1, double y1, double x2, double y2);
// 
//     const TRectDouble &operator=(const TRectDouble &R);
//
//     friend std::ostream &operator<<(std::ostream &os, const TRectDouble &R);
//     friend std::istream &operator>>(std::istream &is, TRectDouble &R);
//----------------------------------------------------------------------------------------
#ifndef typesH
#define typesH
//----------------------------------------------------------------------------------------
#include <complex>
//----------------------------------------------------------------------------------------
// Type re-definitions
//----------------------------------------------------------------------------------------
typedef std::complex<double> COMPLEX;    // Short definition for complex numbers
//----------------------------------------------------------------------------------------
// New classes
//----------------------------------------------------------------------------------------
namespace lobo
{
//----------------------------------------------------------------------------------------
// TPointDouble class
//----------------------------------------------------------------------------------------
class TPointDouble
{
  public:
    TPointDouble();
    TPointDouble(double x, double y);
    TPointDouble(const TPointDouble &P);

    ~TPointDouble() {}

    inline double &x()
      { return _x; }
    inline const double &x() const
      { return _x; }

    inline double &y()
      { return _y; }
    inline const double &y() const
      { return _y; }

    void Assign(double x, double y);

    const TPointDouble &operator=(const TPointDouble &P);

    friend TPointDouble operator+(const TPointDouble &P1, const TPointDouble &P2)
    {
      TPointDouble Q(P1._x + P2._x, P1._y + P2._y);
      return Q;
    }
    friend TPointDouble operator-(const TPointDouble &P1, const TPointDouble &P2)
    {
      TPointDouble Q(P1._x - P2._x, P1._y - P2._y);
      return Q;
    }
    friend bool operator<(const TPointDouble &P1, const TPointDouble &P2)
    {
      return (P1._x < P2._x);
    }

    friend std::ostream &operator<<(std::ostream &os, const TPointDouble &P);
    friend std::istream &operator>>(std::istream &is, TPointDouble &P);

  private:
    double _x, _y;
};        // End of TPointDouble class

//----------------------------------------------------------------------------------------
// TRectDouble class
//----------------------------------------------------------------------------------------
class TRectDouble
{
  public:
    TRectDouble();
    TRectDouble(double x1, double y1, double x2, double y2);
    TRectDouble(const TRectDouble &R);

    ~TRectDouble() {}

    inline double &x1()
      { return _x1; }
    inline const double &x1() const
      { return _x1; }

    inline double &x2()
      { return _x2; }
    inline const double &x2() const
      { return _x2; }

    inline double &y1()
      { return _y1; }
    inline const double &y1() const
      { return _y1; }

    inline double &y2()
      { return _y2; }
    inline const double &y2() const
      { return _y2; }

    void Assign(double x1, double y1, double x2, double y2);

    const TRectDouble &operator=(const TRectDouble &R);
    friend std::ostream &operator<<(std::ostream &os, const TRectDouble &R);
    friend std::istream &operator>>(std::istream &is, TRectDouble &R);

  private:
    double _x1, _y1, _x2, _y2;
};        // End of TRectDouble class

}         // End of namespace "lobo"
//----------------------------------------------------------------------------------------
#endif

