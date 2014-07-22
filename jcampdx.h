//----------------------------------------------------------------------------------------
// TJCampDX - Class for saving, reading and storing files in J-Camp DX format
// R.P.S.M. Lobo - lobo@espci.fr
// File created: 2009/02/16
// Last modified: 2012/02/20
// This class is utterly incomplete and written with my IR data in mind. Data can be
// uniformly spaced or simple (x,y) pairs. The tags recognized are (R means used in
// reading and W in writing the file):
// ##TITLE=<title>  (RW -- used as main comment for the file)
// ##JCAMP-DX=4.24  (RW -- used to check file format)
// ##SAMPLING PROCEDURE=<Comment> (R -- optional -- added to TITLE variable for writing)
// ##XUNITS=<units> (RW -- known values 1/CM, THz, eV, meV, nm, um, Generic)
// ##YUNITS=<units> (RW -- anything is valid)
// ##NPOINTS=<value>  (RW for the 4 variables -- At least 3 are required for equally 
// ##FIRSTX=<value>    spaced data. If all four are defined, then LASTX is ignored or 
// ##LASTX=<value>     used for checking only.)
// ##DELTAX=<value>
// ##YFACTOR=<value>  (RW -- value by which Y must be multiplied -- if absent assumes 1)
// ##RESOLUTION=<value> (R -- optional -- added to TITLE variable for writing)
// ##XYDATA=(X++(Y..Y)) or (XY..XY) (RW -- used to indicate equal or unequal spacing)
// "##END=" (RW)
// Some comments: the data is sorted by ascending x values.
//----------------------------------------------------------------------------------------
#ifndef TJCAMPDX_H
#define TJCAMPDX_H
//----------------------------------------------------------------------------------------
#include <vector>
#include <string>
#include <gsl/gsl_sort_vector_double.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_vector.h>
//----------------------------------------------------------------------------------------
//#include "types.h"

//----------------------------------------------------------------------------------------
// Excpetion class to handle errors in Tjcampdx
//----------------------------------------------------------------------------------------
class Ejcampdx
{
  public:
    // Known errors
    enum ErrorType
		{
      BAD_QNG_FORMAT,       // qng file has unknown format
      CANNOT_CREATE_FILE,   // File cannot be created
      CANNOT_MERGE,         // Data ranges do not allow merging
      DATA_NOT_COMPATIBLE,  // The data are not compatible
      EMPTY_DATA,           // Final data size is empty
      END_NOT_FOUND,        // ##END= tag not found
      FILE_NOT_FOUND,       // File does not exist
      GENERIC_MESSAGE,      // Generic message, for example from a caught exception
      HEADER_FAILURE,       // Failed reading JCamp-dx header
      INCOMPATIBLE_XF,      // XF incompatible with xi, dx, and npt
      NO_NEGATIVE_LOG,      // x value requested for log scale is negative
      NO_VERSION_OR_TYPE,   // No information on version or data type
      OUT_OF_BOUNDS,        // Trying to access data outside the range
      RANGE_TOO_SMALL,      // Generated range is too small
      UNKNOWN_COMPAT,       // Compatibility tests are incompetent
      XY_NO_SIZE_MATCH,     // Size of x and y vectors passed to dx are different
      ZERO_CONVERSION       // Conversion factor cannot be zero
    };

    // Contructor
    Ejcampdx(ErrorType error) :
      Error(error) {}

    // Sets the generic message
    void SetMessage(const std::string &msg)
    {
      GenericMessage = msg;
    }

    // Returns the current error message
    std::string Message()
    {
      switch (Error)
      {
        case BAD_QNG_FORMAT:
          return "qng file has unknown format.";
        case CANNOT_CREATE_FILE:
          return "File cannot be created.";
        case CANNOT_MERGE:
          return "Range of data sets do not allow merging them together";
        case DATA_NOT_COMPATIBLE:
          return "The data are not compatible for this operation.";
        case EMPTY_DATA:
          return "Final data size is empty";
        case END_NOT_FOUND:
          return "##END tag not found in file.";
        case FILE_NOT_FOUND:
          return "File does not exist";
        case GENERIC_MESSAGE:
          return GenericMessage;
        case HEADER_FAILURE:
          return "Failed reading JCamp-dx header.";
        case INCOMPATIBLE_XF:
          return "Last x is incompatible with information on xi, dx, and npt.";
        case NO_NEGATIVE_LOG:
          return "Range requested for log scale contains negative values.";
        case NO_VERSION_OR_TYPE:
          return "File does not contain information on version or data type.";
        case OUT_OF_BOUNDS:
          return "Trying to access data outside the defined range.";
        case RANGE_TOO_SMALL:
          return "The data range is too small.";
        case UNKNOWN_COMPAT:
          return "Compatibility tests are incompetent.";
        case XY_NO_SIZE_MATCH:
          return "The sizes of x and y vectors are different";
        case ZERO_CONVERSION:
          return "Conversion factor cannot be zero";
        default :
          return "Unknown error. How did you get here?";
      }
    }

  private:
    ErrorType Error;  // Contains the error type
    std::string GenericMessage;
};
//----------------------------------------------------------------------------------------
// Structure to hold extreme limits of data. Used as a pointer, allows to have a const
// class and, yet, change the values of the extremes.
//----------------------------------------------------------------------------------------
struct TExtremes
{
  double x_min,x_max,y_min,y_max;
};

//----------------------------------------------------------------------------------------
class Tjcampdx;   // Forward declaration of Tjcampdx class need for namespace jcamp
//----------------------------------------------------------------------------------------
// namespace jcamp to define functions acting on Tjcampdx objects
//----------------------------------------------------------------------------------------
namespace jcamp
{
  const double ZERO = 5E-5;      // Constant defining a value to be interpreted as zero

  // Averages a list of Tjcampdx data which addresses are in jc. jc[0] is the master data.
  Tjcampdx Average(const std::vector<const Tjcampdx *> &jc);

  // (FRIEND) Function to make two spectra compatible. If it succeeds, AT LEAST ONE OF  
  // THE DATA BLOCKS WILL BE OVERWRITTEN. Operation will change files so that the smallest 
  // range and the data spacing of jc1 is preserved. 
  void MakeRangeEquivalent(Tjcampdx &jc1, Tjcampdx &jc2);

  // (FRIEND) Merges two data sets. Keeps the resolution of the lowest frequency data.
  Tjcampdx Merge(const Tjcampdx &jc1, const Tjcampdx &jc2, int Mode);

  // Math functions -- do not overwrite data (overwriting versions are in the class definition)
  Tjcampdx ABS(const Tjcampdx &jc);           // Modulus
  Tjcampdx ACOS(const Tjcampdx &jc);          // Arc-cosinus
  Tjcampdx ASIN(const Tjcampdx &jc);          // Arc-sinus
  Tjcampdx ATAN(const Tjcampdx &jc);          // Arc-tangent
  Tjcampdx COS(const Tjcampdx &jc);           // Cosinus
  Tjcampdx Differentiate(const Tjcampdx &jc); // Derivative
  Tjcampdx EXP(const Tjcampdx &jc);           // Exponential
  Tjcampdx EXP10(const Tjcampdx &jc);         // 10^y
  Tjcampdx Integrate(const Tjcampdx &jc);     // Integral
  Tjcampdx Inv(const Tjcampdx &jc);           // 1 / y
  Tjcampdx LOG(const Tjcampdx &jc);           // Logarithm
  Tjcampdx LOG10(const Tjcampdx &jc);         // Base 10 logarithm
  Tjcampdx SIN(const Tjcampdx &jc);           // Sinus
  Tjcampdx Sqrt(const Tjcampdx &jc);          // Square root
  Tjcampdx Square(const Tjcampdx &jc);        // y^2
  Tjcampdx TAN(const Tjcampdx &jc);           // Tangent

  // Data modifying functions -- do not overwrite data (overwriting versions are in the class definition)
  Tjcampdx ArbSpacing(const Tjcampdx &jc);       
  Tjcampdx Cut(const Tjcampdx &jc, double x0, double x1);
  Tjcampdx Filter(const Tjcampdx &jc, double Low, double High);
  Tjcampdx Morph(const Tjcampdx &jc, const Tjcampdx &jc_Orig);                 
  Tjcampdx Noisify(const Tjcampdx &jc, double Amplitude, bool Prop = true);
  Tjcampdx Resize(const Tjcampdx &jc, double x_0, double x_1, unsigned int n, bool LogScale = false);
  Tjcampdx Resize(const Tjcampdx &jc, unsigned int n);
  Tjcampdx Smooth(const Tjcampdx &jc, unsigned WinSize, unsigned Order);
  Tjcampdx UnitConv(const Tjcampdx &jc, double ConvFactor);
}
//----------------------------------------------------------------------------------------
// JCamp-DX class definition
//----------------------------------------------------------------------------------------
class Tjcampdx
{
  public:
    enum      // Enumeration for compatibility of data
    {
      EQUIVALENT,     // Evenly spaced and equivalent
      SUB_SET,        // Range is a sub set of this block
      SUPER_SET,      // Range is a super set of this block
      ABOVE_OVERLAP,  // Range is partially overlapping but goes above range
      BELOW_OVERLAP,  // Range is partially overlapping but goes below range
			SAME_RANGE,     // Range is the same but do not have the same resolution or spacing type
      NO_OVERLAP      // Ranges are not compatible - no overlap
    };

    // Enumeration to define merged spectra resolution between two overlapping buffers
    enum
    {
      mgHighRes,      // Produces a final spectrum with the resolution equals the highest
      mgLowRes,       // Produces a final spectrum with the resolution equals the lowest
      mgIndependent   // Don't change resolution of both parts
    };

    enum                  // Known units
    {
      un_Wavenumber,
      un_eV,
      un_meV,
      un_THz,
      un_Microns,
      un_nm,
      un_NotDefined
    };

    // Constructors / Destructor
    Tjcampdx();                                     // Default
    Tjcampdx(const Tjcampdx &jc);                   // Copy
    Tjcampdx(const std::string &FileName);          // Loading data from a file
    Tjcampdx(unsigned n, double x0, double dx);     // n Data [= 0] points
    Tjcampdx(unsigned n, double x0, double dx, double offset);      // n Data [y = x + offset] points
    ~Tjcampdx();                                    // Destructor

    // Member functions
    void Clear() { *this = Tjcampdx(); }            // Erases all data
    int EquiRange(const Tjcampdx &jc) const;        // Is range equivalent to this? Return enum.
    void Export_asc(const std::string &FileName,    // Exports file in ascii format
                    unsigned Skip = 0,              // ...# of points to skip between two exported values 
                    char Spacer = ' ') const;       // ...spacer character for the file.
    void GenerateXY(double x1, double x2,                           // Generates a pair of x,y vectors with range [x1,x2]. If at least
                    std::vector<double> &x, std::vector<double> &y, // 5 points cannot be generated throws a RANGE_TOO_SMALL exception.
                    int n = -1) const;                              // n < 0 uses all points. Otherwise interpolates to have n points.
    void Import_asc(const std::string &FileName,    // Imports file containing ascii data
                    unsigned Skip = 0);             // ...# of line to skip in the beginning of the file
    void Import_csv(const std::string &FileName);   // Import simple two column csv format
    void Import_qng(const std::string &FileName);   // Import from Qenga format
    int Index(double x);                            // Returns the index corresponding value <= x. 
    double Interpolate(double x) const;             // Finds a value at arbritrary x
    bool IsEmpty() { return jcDataY.empty(); }       // True if the jcampdx object has no data
    void LoadData(double x0, double dx, const std::vector<double> &y);      // Loads data starting at x0, spaced by dx and placed in y
    void LoadData(const std::vector<double> &x, const std::vector<double> &y);  // Loads data from two vectors
    double MatchFactor(const Tjcampdx &jc, double x0, double x1); // Factor = <jc>/<this>; where <...> is the average y from x0 to x1
    void Open(const std::string &FileName);         // Opens a JCamp-DX file
    void Save(const std::string &FileName) const;   // Saves a JCamp-DX file
    unsigned Search(double x, double y, double dx, double dy);  // Returns index of closest point to P within the box defined by +/dP. No point found, returns the size of the data.

    // Member functions -- overwrite data -- see functions in namespace jcamp to data preserving versions
    void ArbSpacing();                              // Forces the data to become arbitrarily spaced
    void Cut(double x0, double x1);                 // Changes the range (keeps x0 and x1) without changing anything else
    void Filter(double Low, double High);           // Filters the data set
    void Morph(const Tjcampdx &jc);                 // Makes data of the same range and resolution as jc 
    void Noisify(double Amplitude,                  // Adds noise. Amplitude is in %. 
                     bool Prop = true);             // Prop = (true,false) gives noise ~ [y, max(y)]
    void Resize(double x_0, double x_1,             // Resizes, but does not extrapolate, the data. In linear scale, interpolates
                unsigned int n,                     // to place n points between x0 and x1. The output data is evenly spaced. 
                bool LogScale = false);             // In log scale the output is arbitrarily spaced.
    void Resize(unsigned int n);                    // Keeps limits and changes (evenly spaced) number of points 
    void Smooth(unsigned WinSize, unsigned Order);  // Savitzky-Golay smoothing
    void UnitConv(double ConvFactor);               // Converts the x units using a factor (> 0 -> k * x; < 0 -> k / x) -- See PhysConst.h The function does not change the type or title of the unit

    //  Data access functions
    inline const std::string &Comment() const { return jcComment; }
    inline void Comment(const std::string &Com) { jcComment = Com; }

    inline bool IsEvenSpaced() const { return jcConstSpacing; }

    inline const TExtremes &Extremes() const { return *Extr; }

    const char *Units() const;                // Known units are: "1/cm"; "eV";
    const char *UnitsHTML() const;
    void Units(const std::string &Line);            // ..."meV"; "THz"; "um"; and "nm".
    void Units(int un) { jcUnits = un; }            // If user wants to use one of the enumerations

    inline const std::string &XTitle() const { return jcXTitle; }
    inline void XTitle(const std::string Title) { jcXTitle = Title; }

    inline const std::string &YTitle() const { return jcYTitle; }
    inline void YTitle(const std::string Title) { jcYTitle = Title; }

    inline double x0() const { return jcx0; }
    inline double dx() const { return jcdx; }
    inline double xf() const { return jcxf; }
    inline unsigned int Size() const { return jcDataY.size(); }

    inline const std::vector<double> &DataX() const { return jcDataX; }
    inline const std::vector<double> &DataY() const { return jcDataY; }

    // Math functions -- overwrite data -- see functions in namespace jcamp to data preserving versions
    void ABS();           // Modulus
    void ACOS();          // Arc-cosinus
    void ASIN();          // Arc-sinus
    void ATAN();          // Arc-tangent
    void COS();           // Cosinus
    void Differentiate(); // Derivative  
    void EXP();           // Exponential
    void EXP10();         // 10^y
    void Integrate();     // Integral
    void Inv();           // 1 / y
    void LOG();           // Logarithm
    void LOG10();         // Base 10 logarithm
    void SIN();           // Sinus
    void Sqrt();          // Square root
    void Square();        // y^2
    void TAN();           // Tangent

    // Friend functions declared in jcamp namespace (see beginning of file)
    friend void jcamp::MakeRangeEquivalent(Tjcampdx &jc1, Tjcampdx &jc2);             // Make two data sets compatible
    friend Tjcampdx jcamp::Merge(const Tjcampdx &jc1, const Tjcampdx &jc2, int Mode = mgHighRes); // Merges two data sets

    // Friend operators and streams
    friend Tjcampdx operator+(const Tjcampdx &jc1, const Tjcampdx &jc2);
    friend Tjcampdx operator-(const Tjcampdx &jc1, const Tjcampdx &jc2);
    friend Tjcampdx operator/(const Tjcampdx &jc1, const Tjcampdx &jc2);
    friend Tjcampdx operator*(const Tjcampdx &jc1, const Tjcampdx &jc2);

    friend Tjcampdx operator+(const Tjcampdx &jc, double a);
    friend Tjcampdx operator+(double a, const Tjcampdx &jc);

    friend Tjcampdx operator-(const Tjcampdx &jc, double a);
    friend Tjcampdx operator-(double a, const Tjcampdx &jc);

    friend Tjcampdx operator*(const Tjcampdx &jc, double a);
    friend Tjcampdx operator*(double a, const Tjcampdx &jc);

    friend Tjcampdx operator/(const Tjcampdx &jc, double a);
    friend Tjcampdx operator/(double a, const Tjcampdx &jc);

    friend std::ostream &operator<<(std::ostream &dataout, const Tjcampdx &jc);
    friend std::istream &operator>>(std::istream &datain, Tjcampdx &jc);

    // Operators
    double x(unsigned int i) const;       // Returns the i-th x. Not an operator but here is its place
    double &operator[](unsigned int Index) { return jcDataY[Index]; }
    const double &operator[](unsigned int Index) const { return jcDataY[Index]; }
    const Tjcampdx &operator=(const Tjcampdx &jc);

    Tjcampdx operator+=(const Tjcampdx &jc) { *this = (*this) + jc; return *this; }
    Tjcampdx operator-=(const Tjcampdx &jc) { *this = (*this) - jc; return *this; }
    Tjcampdx operator*=(const Tjcampdx &jc) { *this = (*this) * jc; return *this; }
    Tjcampdx operator/=(const Tjcampdx &jc) { *this = (*this) / jc; return *this; }

    Tjcampdx operator+=(double a) { *this = (*this) + a; return *this; }
    Tjcampdx operator-=(double a) { *this = (*this) - a; return *this; }
    Tjcampdx operator*=(double a) { *this = (*this) * a; return *this; }
    Tjcampdx operator/=(double a) { *this = (*this) / a; return *this; }

  private:
    static const unsigned POLY_INT = 4;     // Number of points for polynomial interpolation
    static const unsigned MIN_XY_SIZE = 5;  // Minimum number of points for XY generation

    // Mandatory JCamp-DX fields. Version and XY data type are always mandatory.
    // In equally spaced data, three out of the other four fields are mandatory.
    enum      
    {
      mf_Version = 0x01,
      mf_FirstX = 0x02,
      mf_LastX = 0x04,
      mf_DeltaX = 0x08,
      mf_NPoints = 0x10,
      mf_XYData = 0x20
    };

    // Private functions
    int Locate(const std::vector<double> &xx, double x);
    int Hunt(const std::vector<double> &xx, double x, int guess) const;
    double Random(int i=1);
    template<class T> void VectorKeep(std::vector<T> &v, unsigned i0, unsigned i1);
    void ReorderVector(const gsl_permutation *p, std::vector<double> &v);
    std::vector<double> vector_gsl2std(const gsl_vector *vgsl);
    gsl_vector *vector_std2gsl(const std::vector<double> &vstd);
    void SortTwoVectors(std::vector<double> &v1, std::vector<double> &v2);
    void FindBounds(double x_0, double x_1, unsigned &LB, unsigned &UB) const;  // Lower and upper indexes for sub-set in the range [x_0,x_1]
    void FindExtremes() const;                                        // Find extreme values for x and y
    void FindOverlap(const Tjcampdx &jc, double &x0, double &x1) const; // Loads x0 and x1 with the overlapping regions between jc and *this
    std::string GetKey(std::string &Line);                            // Gets a key and replaces Line by key value
    void ParseJCampDX(const std::vector<std::string> &Contents);      // Contents to Tjcampdx
    void ParseJCampDXHeader(const std::vector<std::string> &Contents, // Handles the file header
                            unsigned int DataStart);

    typedef std::vector<double>::iterator Iter;   // Shortcut for vector iterator

    // Data for JCamp-dx format definition
    bool jcConstSpacing;            // true if the data has a constant spacing

    double jcx0, jcxf, jcdx,        // First x, last x and increment
           jcYFactor;               // Factor by which y must be multiplied

    int jcUnits;                    // Units for the x axis

    Tjcampdx *jcUnevenCompat;       // Pointer to unevenly spaced data made compatible to this

    std::string jcComment,          // Title, comment, procedure, etc...
                jcXTitle,           // Title for the x axis
                jcYTitle;           // Title for the y axis

    std::vector<double> jcDataX,    // Vector to hold X points if jcConstSpacing == false
                        jcDataY;    // Vector to hold the data

    TExtremes *Extr;                // Pointer to struct holding max and min x & y values
};
//----------------------------------------------------------------------------------------
#endif // TJCAMPDX_H
