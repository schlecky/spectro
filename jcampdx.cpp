//----------------------------------------------------------------------------------------
// TJCampDX - Class for saving, reading and storing files in J-Camp DX format
// R.P.S.M. Lobo - lobo@espci.fr
// File created: 2009/02/16
// Last modified: 2012/02/20
//----------------------------------------------------------------------------------------
#include "jcampdx.h"
//----------------------------------------------------------------------------------------
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <ctime>
#include <string>
#include <vector>
#include <gsl/gsl_sort_vector_double.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_vector.h>
//----------------------------------------------------------------------------------------
//#include "numbers.h"
//#include "letters.h"
//#include "numerics.h"
//----------------------------------------------------------------------------------------
using std::endl;


// trim from start
static inline std::string &ltrim(std::string &s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
        return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) {
        s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
        return s;
}

// trim from both ends
static inline std::string &trim(std::string &s) {
        return ltrim(rtrim(s));
}


//****************************************************************************************
// NAMESPACE jcamp
//****************************************************************************************
//----------------------------------------------------------------------------------------
// Averages a list of Tjcampdx data which addresses are in jc. jc[0] is the master data.
//----------------------------------------------------------------------------------------
Tjcampdx jcamp::Average(const std::vector<const Tjcampdx *> &jc)
{
  if (jc.empty())                           // There is no data in jc...
    throw Ejcampdx(Ejcampdx::EMPTY_DATA);   // ...throws an exception.

  Tjcampdx retval = *jc[0];       // Copies first data into return value
  for (unsigned i = 1; i < jc.size(); i++)
  {
    try { retval += *jc[i]; }   // Tries to add the next data to retval
    catch(Ejcampdx e) { throw; }          // On error stops the whole process
  }
  
  return ( retval / double(jc.size()) );  // Divides the sum by the number of spectra
}
//----------------------------------------------------------------------------------------
// namespace jcamp function -- friend to Tjcampdx
// Function to make two spectra compatible. If it succeeds, AT LEAST ONE OF THE DATA
// BLOCKS WILL BE OVERWRITTEN. Operation will change files so that the smallest range
// and the data spacing of jc1 is preserved.
//----------------------------------------------------------------------------------------
void jcamp::MakeRangeEquivalent(Tjcampdx &jc1, Tjcampdx &jc2)
{
  if (jc1.EquiRange(jc2) == Tjcampdx::EQUIVALENT)    // The ranges are already compatible
    return;

  double x0, x1;
  
  try { jc1.FindOverlap(jc2, x0, x1); }   // Looks for overlapping region between jc1 and jc2
  catch(Ejcampdx e) { throw; }            // On error passes the buck

  if (jc1.IsEvenSpaced())               // jc1 data is evenly spaced
  {
    unsigned n = 1 + int((x1 - x0) / jc1.dx()); // Calculates the number of points that preserves jc1's dx
    try { jc1.Resize(x0, x1, n); }              // Tries to resize jc1 to new limits
    catch(Ejcampdx e) { throw; }                // On error passes the buck
  }
  else                                // jc1 has an arbitrary data spacing
  {
    try { jc1.Cut(x0, x1); }                    // Cuts jc1 to fit the limits
    catch(Ejcampdx e) { throw; }                // On error passes the buck
  }

  try { jc2.Morph(jc1); }             // Morphs jc2 to have the same properties as jc1
  catch(Ejcampdx e) { throw; }        // On error passes the buck

  jc1.jcUnevenCompat = &jc2;          // Sets the compatibility pointers between jc1 and jc2 that...
  jc2.jcUnevenCompat = &jc1;          // ...are used in arbitrary spaced data compatibility.
}
//----------------------------------------------------------------------------------------
// namespace jcamp function -- friend to Tjcampdx
// Merges two data sets. When the data overlap, keeps the resolution defined by Mode.
//----------------------------------------------------------------------------------------
Tjcampdx jcamp::Merge(const Tjcampdx &jc1, const Tjcampdx &jc2, int Mode)
{
  Tjcampdx jclow, jchigh,       // Spectra containing lower and higher frequecies data
           Midlow, Midhigh;     // Spectra containing overlapping regions

  bool Overlap = true,          // Initially assumes that the data have an overlapping range...
       Resamp = true;           // ...and that data will be resampled for a single resolution

  int i0, i1, n;                // Indexes defining the overlapping range and number of points

  double x0, x1,                // Defines the overlapping range
         w, width;              // Weight and width for overlapping interval average

  switch (jc1.EquiRange(jc2))      // Checks which kind of overlapping exists
  {
    case Tjcampdx::ABOVE_OVERLAP: // jc2 overlaps jc1 and jc2 frequencies are higher
      jclow = jc1;                // Sets the low frequency data
      jchigh = jc2;               // Sets the high frequency data
      break;

    case Tjcampdx::BELOW_OVERLAP: // jc2 overlaps jc1 and jc2 frequencies are lower
      jclow = jc2;                // Sets the low frequency data
      jchigh = jc1;               // Sets the high frequency data
      break;

    case Tjcampdx::NO_OVERLAP:                    // jc1 and jc2 do not overlap
      jclow = jc1.x0() < jc2.x0() ? jc1 : jc2;    // Sets the low frequency data
      jchigh = jc1.x0() > jc2.x0() ? jc1 : jc2;   // Sets the high frequency data
      Overlap = false;                            // And sets the overlapping flag to false
      break;

    default:                    // Covers results COMPAT, SUB_SET, SUPER_SET and SAME_RANGE. All unmergeable.
      throw Ejcampdx(Ejcampdx::CANNOT_MERGE);   // Throws exception
  }

  // Decides about merging spacing and resolution
  if ( (Mode == Tjcampdx::mgIndependent) ||                     // User wants independent merging
       (!Overlap) ||                                            // Data does not overlap
       ((!jclow.IsEvenSpaced()) && (!jchigh.IsEvenSpaced())) )  // No data set is even spaced
  {
    Resamp = false;                                           // Changes resample flag to force appending data
  }
  else    // Data will be resampled to a single resolution throughout the data set
  {
    // Resamples, if necessary, jclow to contain final resolution
    if (jchigh.IsEvenSpaced())  // If this is not true, the LF already has the final resolution
    {
      if ( ((Mode == Tjcampdx::mgHighRes) && (jchigh.dx() < jclow.dx())) ||   // Wants to keep highest resolution and it is in jchigh or...
           ((Mode == Tjcampdx::mgLowRes) && (jchigh.dx() > jclow.dx())) ||    // ...wants to keep lowest resolution and it is again in jchigh.
           (!jclow.IsEvenSpaced()) )                                          // Only jchigh is even spaced
      {
        n = 1 + int((jclow.xf() - jclow.x0()) / jchigh.dx());   // Calculates the number of points to create in jclow to match jchigh resolution
        jclow.Resize(jclow.x0(), jclow.xf(), n);                // Resizes jclow to match higher resolution of jchigh
      }
    }
  }


  if (Resamp)       // Resample data for a single resolution based on jclow (as redefined above)
  {
    jclow.FindOverlap(jchigh, x0, x1);      // Find overlapping range (only x0 is important, x1 will be interpolated)

    i0 = jclow.Index(x0);           // Calculates the position just below x0 in jclow
    x0 = jclow.x(i0);               // Makes sure that x0 belongs to jclow

    n = 1 + int((jchigh.xf() - x0) / jclow.dx());   // Calculates the number of points to create in jchigh to match jclow resolution
    jchigh.Resize(x0, jchigh.xf(), n);              // Changes jchigh to constant spacing matching jclow properties

    Midlow = jclow;                   // Copies LF and...
    Midhigh = jchigh;                 // ...HF data to their respective mid portions.

    jclow.Cut(jclow.x0(), x0);        // Cuts LF data below overlapping part [xa,xb]

    x0 += 1.1 * jclow.dx();   // Adds a step to define mid portion initial x. 1.1 avoids truncation error.
    x1 = Midlow.xf();         // And goes to the end of the LF data
    Midlow.Cut(x0, x1);       // Gets LF overlapping part (xb,xc]
    Midhigh.Cut(x0, x1);      // Makes MidHigh the same range as MidLow (Morph also works but it is slower)

    x1 += 1.1 * jclow.dx();   // Adds a step to define HF initial x. 1.1 avoids truncation error.

    jchigh.Cut(x1, jchigh.xf());  // Cuts HF data above overlapping part (xc,xd]

    width = x1 - x0;              // Interval width
    x1 = Midlow.x0();             // x1 becomes the running x
    for (unsigned int i = 0; i < Midlow.Size(); i++)  // Weighted average of the two midbits - Puts it in Midlow
    {
      w = (x1 - x0) / width;                                                      // Weight and...
      Midlow.jcDataY[i] = Midlow.jcDataY[i] * (1.0 - w) + Midhigh.jcDataY[i] * w; // ...average
      x1 += Midlow.dx();          // Increments x1
    }

    // Paste all y data bits together
    jclow.jcDataY.insert(jclow.jcDataY.end(), Midlow.jcDataY.begin(), Midlow.jcDataY.end());
    jclow.jcDataY.insert(jclow.jcDataY.end(), jchigh.jcDataY.begin(), jchigh.jcDataY.end());

    jclow.jcxf = jclow.x0() + (jclow.Size() - 1) * jclow.dx();    // Updates last x
  }
  else              // Append data
  {
    jclow.ArbSpacing();       // Converts, if necessary, both vectors...
    jchigh.ArbSpacing();      // ...to arbitrary data spacing.

    if (Overlap)    // Data overlap -- create new jclow with overlapping portion and erases it from jchigh
    {
      x0 = jchigh.x0();       // Beginning of overlapping range
      x1 = jclow.xf();        // End of overlapping range

      //i0 = numerics::Locate(jclow.jcDataX, x0) + 1;     // Finds index just above x0 in jclow
      i0=0;
      while((i0<jclow.DataX().size()) &&   (jclow.DataX()[i0]<x0))
          i0++;

      //i1 = numerics::Locate(jchigh.jcDataX, x1) + 1;    // Finds index just above x1 in jchigh
      i1=0;
      while((i1<jchigh.DataX().size()) &&   (jchigh.DataX()[i1]<x1))
          i1++;

      Midlow.jcDataX.assign(jclow.jcDataX.begin() + i0, jclow.jcDataX.end()); // Loads mid portion x's from LF data
      Midlow.jcDataY.assign(jclow.jcDataY.begin() + i0, jclow.jcDataY.end()); // Loads mid portion y's from LF data

      Midhigh.jcDataX.assign(jchigh.jcDataX.begin(), jchigh.jcDataX.begin() + i1); // Loads mid portion x's from HF data
      Midhigh.jcDataY.assign(jchigh.jcDataY.begin(), jchigh.jcDataY.begin() + i1); // Loads mid portion y's from HF data

      Midlow.jcConstSpacing = Midhigh.jcConstSpacing = false;   // Makes sure that data is set to arbitrary spacing
      Midhigh.Morph(Midlow);                                    // Morphs Midhigh with Midlow properties

      jclow.jcDataX.erase(jclow.jcDataX.begin() + i0, jclow.jcDataX.end());   // Erases overlapping x range from LF data
      jclow.jcDataY.erase(jclow.jcDataY.begin() + i0, jclow.jcDataY.end());   // Erases overlapping y range from LF data

      jchigh.jcDataX.erase(jchigh.jcDataX.begin(), jchigh.jcDataX.begin() + i1);  // Erases overlapping x range from HF data
      jchigh.jcDataY.erase(jchigh.jcDataY.begin(), jchigh.jcDataY.begin() + i1);  // Erases overlapping y range from HF data

      x0 = Midlow.jcDataX.front();                      // Updates beginning of overlapping interval
      width = Midlow.jcDataX.back() - x0;               // Interval width
      for (unsigned int i = 0; i < Midlow.Size(); i++)  // Weighted average of the two midbits - Puts it in Midlow
      {
        w = (Midlow.jcDataX[i] - x0) / width;                                       // Weight and...
        Midlow.jcDataY[i] = Midlow.jcDataY[i] * (1.0 - w) + Midhigh.jcDataY[i] * w; // ...average
      }

      // Paste overlapping block into jclow
      jclow.jcDataX.insert(jclow.jcDataX.end(), Midlow.jcDataX.begin(), Midlow.jcDataX.end());
      jclow.jcDataY.insert(jclow.jcDataY.end(), Midlow.jcDataY.begin(), Midlow.jcDataY.end());
    }

    // Appends high frequency data to low frequency data
    jclow.jcDataX.insert(jclow.jcDataX.end(), jchigh.jcDataX.begin(), jchigh.jcDataX.end());
    jclow.jcDataY.insert(jclow.jcDataY.end(), jchigh.jcDataY.begin(), jchigh.jcDataY.end());

    jclow.jcConstSpacing = false;         // Sets arbitrary data spacing
    jclow.jcx0 = jclow.jcDataX.front();   // Sets first x
    jclow.jcxf = jclow.jcDataX.back();    // Sets last x
    jclow.jcdx = 0;                       // dx has no meaning
  }

  jclow.FindExtremes();       // Looks for extreme values of data

  return jclow;               // Returns merged data
}
//----------------------------------------------------------------------------------------
// namespace jcamp function -- Modulus
//----------------------------------------------------------------------------------------
Tjcampdx jcamp::ABS(const Tjcampdx &jc)         
{
  Tjcampdx retval = jc;
  retval.ABS();
  return retval;
}
//----------------------------------------------------------------------------------------
// namespace jcamp function -- Arc-cosinus
//----------------------------------------------------------------------------------------
Tjcampdx jcamp::ACOS(const Tjcampdx &jc)
{
  Tjcampdx retval = jc;
  retval.ACOS();
  return retval;
}
//----------------------------------------------------------------------------------------
// namespace jcamp function -- Arc-sinus
//----------------------------------------------------------------------------------------
Tjcampdx jcamp::ASIN(const Tjcampdx &jc)
{
  Tjcampdx retval = jc;
  retval.ASIN();
  return retval;
}
//----------------------------------------------------------------------------------------
// namespace jcamp function -- Arc-tangent
//----------------------------------------------------------------------------------------
Tjcampdx jcamp::ATAN(const Tjcampdx &jc)
{
  Tjcampdx retval = jc;
  retval.ATAN();
  return retval;
}
//----------------------------------------------------------------------------------------
// namespace jcamp function -- Cosinus
//----------------------------------------------------------------------------------------
Tjcampdx jcamp::COS(const Tjcampdx &jc)
{
  Tjcampdx retval = jc;
  retval.COS();
  return retval;
}
//----------------------------------------------------------------------------------------
// namespace jcamp function -- Derivative
//----------------------------------------------------------------------------------------
Tjcampdx jcamp::Differentiate(const Tjcampdx &jc)
{
  Tjcampdx retval = jc;
  retval.Differentiate();
  return retval;
}
//----------------------------------------------------------------------------------------
// namespace jcamp function -- Exponential
//----------------------------------------------------------------------------------------
Tjcampdx jcamp::EXP(const Tjcampdx &jc) 
{
  Tjcampdx retval = jc;
  retval.EXP();
  return retval;
}
//----------------------------------------------------------------------------------------
// namespace jcamp function -- 10^y
//----------------------------------------------------------------------------------------
Tjcampdx jcamp::EXP10(const Tjcampdx &jc)      
{
  Tjcampdx retval = jc;
  retval.EXP10();
  return retval;
}
//----------------------------------------------------------------------------------------
// namespace jcamp function -- Integral
//----------------------------------------------------------------------------------------
Tjcampdx jcamp::Integrate(const Tjcampdx &jc)
{
  Tjcampdx retval = jc;
  retval.Integrate();
  return retval;
}
//----------------------------------------------------------------------------------------
// namespace jcamp function -- 1 / y
//----------------------------------------------------------------------------------------
Tjcampdx jcamp::Inv(const Tjcampdx &jc)        
{
  Tjcampdx retval = jc;
  retval.Inv();
  return retval;
}
//----------------------------------------------------------------------------------------
// namespace jcamp function -- Logarithm
//----------------------------------------------------------------------------------------
Tjcampdx jcamp::LOG(const Tjcampdx &jc) 
{
  Tjcampdx retval = jc;
  retval.LOG();
  return retval;
}
//----------------------------------------------------------------------------------------
// namespace jcamp function -- Base 10 logarithm
//----------------------------------------------------------------------------------------
Tjcampdx jcamp::LOG10(const Tjcampdx &jc)       
{
  Tjcampdx retval = jc;
  retval.LOG10();
  return retval;
}
//----------------------------------------------------------------------------------------
// namespace jcamp function -- Sinus
//----------------------------------------------------------------------------------------
Tjcampdx jcamp::SIN(const Tjcampdx &jc)
{
  Tjcampdx retval = jc;
  retval.SIN();
  return retval;
}
//----------------------------------------------------------------------------------------
// namespace jcamp function -- Square root
//----------------------------------------------------------------------------------------
Tjcampdx jcamp::Sqrt(const Tjcampdx &jc)   
{
  Tjcampdx retval = jc;
  retval.Sqrt();
  return retval;
}
//----------------------------------------------------------------------------------------
// namespace jcamp function -- y^2
//----------------------------------------------------------------------------------------
Tjcampdx jcamp::Square(const Tjcampdx &jc)  
{
  Tjcampdx retval = jc;
  retval.Square();
  return retval;
}
//----------------------------------------------------------------------------------------
// namespace jcamp function -- Tangent
//----------------------------------------------------------------------------------------
Tjcampdx jcamp::TAN(const Tjcampdx &jc)
{
  Tjcampdx retval = jc;
  retval.TAN();
  return retval;
}
//----------------------------------------------------------------------------------------
// namespace jcamp function -- Forces the data to become arbitrarily spaced
//----------------------------------------------------------------------------------------
Tjcampdx jcamp::ArbSpacing(const Tjcampdx &jc)
{
  Tjcampdx retval = jc;
  retval.ArbSpacing();
  return retval;
}
//----------------------------------------------------------------------------------------
// namespace jcamp function -- Changes the range (keeps x0 and x1) without changing 
// anything else
//----------------------------------------------------------------------------------------
Tjcampdx jcamp::Cut(const Tjcampdx &jc, double x0, double x1)
{
  Tjcampdx retval = jc;
  retval.Cut(x0,x1);
  return retval;
}
//----------------------------------------------------------------------------------------
// namespace jcamp function -- Filters the data set
//----------------------------------------------------------------------------------------
Tjcampdx jcamp::Filter(const Tjcampdx &jc, double Low, double High)
{
  Tjcampdx retval = jc;
  retval.Filter(Low,High);
  return retval;
}
//----------------------------------------------------------------------------------------
// namespace jcamp function -- Makes data of the same range and resolution as jc 
//----------------------------------------------------------------------------------------
Tjcampdx jcamp::Morph(const Tjcampdx &jc, const Tjcampdx &jc_Orig)
{
  Tjcampdx retval = jc;
  retval.Morph(jc_Orig);
  return retval;
}
//----------------------------------------------------------------------------------------
// namespace jcamp function -- Adds noise
//----------------------------------------------------------------------------------------
Tjcampdx jcamp::Noisify(const Tjcampdx &jc, double Amplitude, bool Prop)
{
  Tjcampdx retval = jc;
  retval.Noisify(Amplitude, Prop);
  return retval;
}
//----------------------------------------------------------------------------------------
// namespace jcamp function -- Resizes
//----------------------------------------------------------------------------------------
Tjcampdx jcamp::Resize(const Tjcampdx &jc, double x_0, double x_1, unsigned int n, bool LogScale)
{
  Tjcampdx retval = jc;
  retval.Resize(x_0, x_1, n, LogScale);
  return retval;
}
//----------------------------------------------------------------------------------------
// namespace jcamp function -- Keeps limits and changes (evenly spaced) number of points 
//----------------------------------------------------------------------------------------
Tjcampdx jcamp::Resize(const Tjcampdx &jc, unsigned int n)
{
  Tjcampdx retval = jc;
  retval.Resize(n);
  return retval;
}
//----------------------------------------------------------------------------------------
// namespace jcamp function -- Savitzky-Golay smoothing
//----------------------------------------------------------------------------------------
Tjcampdx jcamp::Smooth(const Tjcampdx &jc, unsigned WinSize, unsigned Order)
{
  Tjcampdx retval = jc;
  retval.Smooth(WinSize, Order);
  return retval;
}
//----------------------------------------------------------------------------------------
// namespace jcamp function -- x units conversion smoothing
//----------------------------------------------------------------------------------------
Tjcampdx jcamp::UnitConv(const Tjcampdx &jc, double ConvFactor)
{
  Tjcampdx retval = jc;
  retval.UnitConv(ConvFactor);
  return retval;
}
//****************************************************************************************
// CLASS DEFINITION Tjcampdx
//****************************************************************************************
//----------------------------------------------------------------------------------------
// Constructor with no arguments
//----------------------------------------------------------------------------------------
Tjcampdx::Tjcampdx() : jcConstSpacing(true), jcx0(0), jcxf(0), jcdx(0), jcYFactor(1),
  jcUnits(un_NotDefined), jcUnevenCompat(NULL), jcComment(""), jcXTitle(""), jcYTitle("")
{
  jcDataX.clear();
  jcDataY.clear();

  Extr = new TExtremes;
  Extr->x_min = Extr->x_max = Extr->y_max = Extr->y_min = 0;
}
//----------------------------------------------------------------------------------------
// Copy constructor
//----------------------------------------------------------------------------------------
Tjcampdx::Tjcampdx(const Tjcampdx &jc) : jcConstSpacing(jc.jcConstSpacing),
  jcx0(jc.jcx0), jcxf(jc.jcxf), jcdx(jc.jcdx), jcYFactor(jc.jcYFactor),
  jcUnits(jc.jcUnits), jcUnevenCompat(jc.jcUnevenCompat), 
  jcComment(jc.jcComment), jcXTitle(jc.jcXTitle), jcYTitle(jc.jcYTitle)
{
  if (jcConstSpacing)
    jcDataX.resize(0);
  else
    jcDataX = jc.jcDataX;

  jcDataY = jc.jcDataY;

  Extr = new TExtremes;
  *Extr = *(jc.Extr);
}
//----------------------------------------------------------------------------------------
// Constructor loading data from a file
//----------------------------------------------------------------------------------------
Tjcampdx::Tjcampdx(const std::string &FileName) : jcUnevenCompat(NULL)
{
  Extr = new TExtremes;
  Open(FileName);
}
//----------------------------------------------------------------------------------------
// Constructor for n Data = 0 points
//----------------------------------------------------------------------------------------
Tjcampdx::Tjcampdx(unsigned n, double x0, double dx) :  jcConstSpacing(true),
  jcx0(x0), jcdx(dx), jcYFactor(1), jcUnits(un_NotDefined), jcUnevenCompat(NULL),
  jcComment(""), jcXTitle(""), jcYTitle(""),
  jcDataY(n,0)
{
  jcDataX.clear();
  jcxf = jcx0 + (n - 1) * jcdx;
  Extr = new TExtremes;
  Extr->x_min = jcx0;
  Extr->y_min = -0.005;
  Extr->x_max = jcxf;
  Extr->y_max = 0.005;
}
//----------------------------------------------------------------------------------------
// Constructor for n Data [y = x + offset] points
//----------------------------------------------------------------------------------------
Tjcampdx::Tjcampdx(unsigned n, double x0, double dx, double offset) : jcConstSpacing(true),
  jcx0(x0), jcdx(dx), jcYFactor(1), jcUnits(un_NotDefined), jcUnevenCompat(NULL),
  jcComment(""), jcXTitle(""), jcYTitle("")
{
  jcDataX.clear();
  jcxf = jcx0 + (n - 1) * jcdx;

  jcDataY.resize(n);
  Iter Item = jcDataY.begin();
  double y = offset + jcx0;               
  while (Item != jcDataY.end()) 
  {
    *Item++ = y;            // Sets value and increments pointer
    y += jcdx;              // Next value
  }

  Extr = new TExtremes;
  FindExtremes();
}
//----------------------------------------------------------------------------------------
// Destructor
//----------------------------------------------------------------------------------------
Tjcampdx::~Tjcampdx()
{
  delete Extr;
}
//----------------------------------------------------------------------------------------
// Check whether jc is compatible with current data. The return value is one of the enum 
// values defined at the beginning of the Tjcampdx class.
//----------------------------------------------------------------------------------------
int Tjcampdx::EquiRange(const Tjcampdx &jc) const
{
  if ( jcConstSpacing && jc.jcConstSpacing &&                              
       fabs(jcx0 - jc.jcx0) < jcamp::ZERO && fabs(jcdx - jc.jcdx) < jcamp::ZERO && 
       jcDataY.size() == jc.Size() )  
    return EQUIVALENT;            // Data are evenly spaced and compatible

  if ( !jcConstSpacing && jc.jcUnevenCompat == this ) // Unevenly spaced data was...
    return EQUIVALENT;                                // ...made compatible to this.

  if ( (jc.jcx0 > jcx0 && jc.jcxf <= jcxf) ||
       (jc.jcx0 >= jcx0 && jc.jcxf < jcxf) )
    return SUB_SET;               // Data range is a sub set of this block

  if ( (jc.jcx0 < jcx0 && jc.jcxf >= jcxf) ||
       (jc.jcx0 <= jcx0 && jc.jcxf > jcxf) )
    return SUPER_SET;             // Data range is a super set of this block

  if ( jc.jcx0 > jcx0 && jc.jcx0 <= jcxf && jc.jcxf > jcxf )
    return ABOVE_OVERLAP;         // Data partially overlapping but goes above this range

  if ( jc.jcx0 < jcx0 && jc.jcxf >= jcx0 && jc.jcxf < jcxf )
    return BELOW_OVERLAP;         // Data partially overlapping but goes below this range

  if ( fabs(jcx0 - jc.jcx0) < jcamp::ZERO && 
       fabs(jc.jcxf - jcxf) < jcamp::ZERO )
    return SAME_RANGE;            // Data have same range but not same resolution or spacing type

  if ( jcx0 >= jc.jcxf || jc.jcx0 >= jcxf )
    return NO_OVERLAP;            // Data do not overlap

  throw Ejcampdx(Ejcampdx::UNKNOWN_COMPAT);   // Should never get here
}
//----------------------------------------------------------------------------------------
// Exports the dx file in a two column x,y ascii format. If 'skip' is not zero makes the 
// export file shrink the data size by skipping points. 'spacer' is the character used to
// separate x and y values. 
//----------------------------------------------------------------------------------------
void Tjcampdx::Export_asc(const std::string &FileName, unsigned Skip, char Spacer) const
{
  std::ofstream os(FileName.c_str(), std::ios::out);    // Creates the file

  if (!os)
    throw Ejcampdx(Ejcampdx::CANNOT_CREATE_FILE);       // Cannot create file exception

  os.precision(15);                   // Max precision for double
  unsigned s = Skip + 1;

  if (jcConstSpacing)                 // Evenly spaced data
  {
    double x = jcx0,                  // x to save
           incr = s * jcdx;           // Actual increment taking into account skiped lines

    for (unsigned i = 0; i < jcDataY.size(); i += s)    // Gets every s points
    {
      os << x << Spacer << jcDataY[i] << std::endl;     // Save pair of values...
      x += incr;                                        // ...and next x
    }
  }
  else                                // Arbitrarily spaced data
  {
    for (unsigned i = 0; i < jcDataY.size(); i += s)            // Gets every s points
      os << jcDataX[i] << Spacer << jcDataY[i] << std::endl;    // Save pair of x,y values
  }

  os.close();                                         // Closes the file
}
//----------------------------------------------------------------------------------------
// Generates a pair of x,y vectors with range [x1,x2]. If at least 5 points cannot 
// be generated throws a RANGE_TOO_SMALL exception. If n is negative all points in
// the range are in the output vector. Otherwise interpolates to have n points.
//----------------------------------------------------------------------------------------
void Tjcampdx::GenerateXY(double x1, double x2,
                std::vector<double> &x, std::vector<double> &y, int n) const
{
  unsigned LB, UB;                // Indexes for sub set past high limit

  try { FindBounds(x1, x2, LB, UB); } // Looks for the indexes. If there is an error...
  catch(Ejcampdx e) { throw; }        // ...passes the buck.

  // Defines first and last values of x
  double xi = jcConstSpacing ? jcx0 + LB * jcdx : jcDataX[LB],
         xf = jcConstSpacing ? jcx0 + UB * jcdx : jcDataX[UB];

  // Calculates the number of points to put in vectors.
  int npt = jcConstSpacing ? int(double(1.0 + (xf - xi) / jcdx)) : UB - LB; 
  if (n > 0) npt = n;         // User wants to specify the number of points

  if (npt < MIN_XY_SIZE)      // Not enough points for any serious fitting...
    throw Ejcampdx(Ejcampdx::RANGE_TOO_SMALL);  // ...throws an exception

  x.resize(npt);        // Resizes x...
  y.resize(npt);        // ...and y vectors.

  if (n > 0)      // User defined own number of points
  {
    double DX = (xf - xi) / double(npt - 1);  // Calculates data increment
    for (int i = 0; i < npt; i++)    
    {
      x[i] = xi;                          // Loads the current x value.
      try { y[i] = Interpolate(x[i]); }   // Tries to interpolate the data. On error...
      catch (Ejcampdx e) { throw; }       // ...passes the buck.
      xi += DX;                           // Increments x to its next value
    }
  }
  else            // Uses original points
  {
    if (jcConstSpacing)     // Handles x's for envenly spaced data
      for (int i = 0; i < npt; i++)
      {
        x[i] = xi;              // Loads the current x value.
        xi += jcdx;             // Increments x to its next value
      }
    else                    // Handles x's for uneven spacing
      x.assign(jcDataX.begin() + LB, jcDataX.begin() + UB + 1); // Copies relevant x data

    y.assign(jcDataY.begin() + LB, jcDataY.begin() + UB + 1);   // Copies relevant y data
  }
}
//----------------------------------------------------------------------------------------
// Import a 2 column ascii file. Ignores all lines beginning with '#'. Skip the first
// 'skip' lines
//----------------------------------------------------------------------------------------
void Tjcampdx::Import_asc(const std::string &FileName, unsigned Skip)
{
  std::ifstream is(FileName.c_str(), std::ios::in); // Opens the file

  if (!is)                                          // File not found...
    throw Ejcampdx(Ejcampdx::FILE_NOT_FOUND);       // ...throws an exception and quits.

  std::vector<double> x, y;       // Vectors for data
  std::string Line;               // Each line of the file
  std::istringstream iss;         // std::string stream to read the data from

  unsigned int skipped = 0;                 // Counting the lines
  while (!is.eof() && skipped++ < Skip)     // Checks stop cond's and increments skipped
    getline(is, Line);                      // Reads a line from the file

  while (!is.eof())               // Go through the file
  {
    getline(is, Line);            // Read one line

    if (Line[0] != '#' && !Line.empty()) // First character is not '#' and line not empty
    {
      iss.clear();                  // Clear the std::string stream
      iss.str(Line);                // Load the line into the std::string stream
      double xi, yi;                // Temp's individual x's and y's
      iss >> xi >> yi;              // Read xi and yi
      x.push_back(xi);              // Load xi and...
      y.push_back(yi);              // ...yi into their containers.
    }
	}

  is.close();                       // Closes the file

  // Sorting: the initial assumption is that the file is sorted in ascending or descending frequencies.
  // So, we spend some time checking for this before attempting to sort the vectors.
	// If the last frequency is smaller than the first, there is a chance that the file was sorted in descending order.
  // We reverse the vectors.
	if ( x[0] > x[x.size()-1] )
	{
		std::reverse(x.begin(), x.end());
		std::reverse(y.begin(), y.end());
	}

  bool Sorted = true;       // Assumes that the file is now sorted but let's make sure.
  for (unsigned i = 1; i < x.size(); i++)
    if (x[i] < x[i-1])      // The file was not sorted
    {
      Sorted = false;
      break;
    }

  if (!Sorted)              // If needed sorts x and changes y accordingly
    SortTwoVectors(x, y);

  try { LoadData(x,y); }          // Tries to load a (x,y) pair of vectors into dx object.
  catch (Ejcampdx e) { throw; }   // On error passes the buck.

  jcUnits = un_Wavenumber;
  jcXTitle = "Wavenumbers <sup>-1</sup>";
  jcYTitle = "Optical response";
}

//----------------------------------------------------------------------------------------
// Convert gsl vectors to std
//----------------------------------------------------------------------------------------
std::vector<double> Tjcampdx::vector_gsl2std(const gsl_vector *vgsl)
{
  std::vector<double> v(vgsl->size);              // Create a vector of the proper size
  v.assign(vgsl->data, vgsl->data + vgsl->size);  // Copy gsl vector into v
}

//----------------------------------------------------------------------------------------
// Convert std vectors to gsl. The created vector must be destroyed later in the code.
//----------------------------------------------------------------------------------------
gsl_vector *Tjcampdx::vector_std2gsl(const std::vector<double> &vstd)
{
  gsl_vector *v = gsl_vector_alloc(vstd.size());     // Create gsl vector
  memcpy(v->data, &vstd[0], vstd.size() * v->stride); // Copy std vector into v
}

//----------------------------------------------------------------------------------------
// Rearranges v according to the permutations defined in Indexes 'p'
//----------------------------------------------------------------------------------------
void Tjcampdx::ReorderVector(const gsl_permutation *p, std::vector<double> &v)
{
  std::vector<double> vtmp = v;           // Copy v to a temporary variable

  for (int i = 0; i < v.size(); i++)
    v[i] = vtmp[p->data[i]];              // Reorder v from the permutation in p
}

//----------------------------------------------------------------------------------------
// Given an array xx[0..n-1] return the index so that xx[i] < x < xx[i+1]. xx must be
// monotonic (increasing or decreasing). The returned index is -1 or n-1 when x is out
// of range. Guess is a first guess of the value. At worst this routine is a factor 2
// slower than Locate and at best a factor log_2(n) faster.
//----------------------------------------------------------------------------------------
int Tjcampdx::Hunt(const std::vector<double> &xx, double x, int guess) const
{
  int n = xx.size();

  int l = guess,          // Estimation of the lower bound
      u,                  // Upper bound
      m, inc;

  bool ascnd(xx[n - 1] >= xx[0]);   // True if ascending table

  if (l < 0 || l > n - 1) // Guess is not useful
  {
    l = -1;
    u = n;
  }
  else
  {
    inc = 1;
    if (x >= xx[l] == ascnd)    // Hunt up
    {
      if (l == n - 1)
        return l;
      u = l + 1;
      while (x >= xx[u] == ascnd)
      {
        l = u;
        inc += inc;
        u = l + inc;
        if (u > n - 1)
        {
          u = n;
          break;
        }
      }
    }
    else                        // Hunt down
    {
      if (l == 0)
        return -1;
      u = l--;
      while (x < xx[l] == ascnd)
      {
        u = l;
        inc *= 2;          // Weird multiplication by 2
        if (inc >= u)
        {
          l = -1;
          break;
        }
        else
          l = u - inc;
      }
    }
  }

  while (u - l != 1)
  {
    m = (u + l) / 2;
    if (x >= xx[m] == ascnd)
      l = m;
    else
      u = m;
  }

  if (x == xx[n - 1])
    return n - 2;
  else if (x == xx[0])
    return 0;
  else
    return l;
}


//----------------------------------------------------------------------------------------
// For a vector xx[0..n-1], return the index 'i' so that xx[i] < x < xx[i+1]. xx must be
// monotonic (increasing or decreasing). The returned index is -1 or n-1 when x is out
// of range.
//----------------------------------------------------------------------------------------
int Tjcampdx::Locate(const std::vector<double> &xx, double x)
{
  int n = xx.size();

  int l(-1),            // Lower and...
      u(n),             // ...upper limits.
      m;

  bool ascnd(xx[n - 1] >= xx[0]);    // True is ascending table.

  while (u - l > 1)               // Iterates through indexes that are "far" appart
  {
    m = (u + l) / 2;              // Gets half way index
    if (x >= xx[m] == ascnd)      // x is bigger that x_m in ascending table...
      l = m;                      // ...so the new lower limit is m.
    else
      u = m;                      // In every other case, m corresponds to the new upper limit
  }

  if (x == xx[0])                 // Handle specific cases of x being the first...
    return 0;
  else if (x == xx[n-1])          // ...or the last element in the vector.
    return n - 2;
  else
    return l;
}


//----------------------------------------------------------------------------------------
// Sort v1 in ascending order and re-arranges v2 in the same manner.
//----------------------------------------------------------------------------------------
void Tjcampdx::SortTwoVectors(std::vector<double> &v1, std::vector<double> &v2)
{
  gsl_vector *vgsl = vector_std2gsl(v1);    // Create gsl vector for v1

  gsl_permutation *p = gsl_permutation_alloc(vgsl->size); // Allocate and...
  gsl_sort_vector_index(p,vgsl);                          // ...find permutation.

  ReorderVector(p, v1);         // Sort the first vector and...
  ReorderVector(p, v2);         // ...reorder the second accordingly.

  // Memory clean-up
  gsl_permutation_free(p);
  gsl_vector_free(vgsl);
}

//----------------------------------------------------------------------------------------
// Import simple two column csv format
//----------------------------------------------------------------------------------------
void Tjcampdx::Import_csv(const std::string &FileName)
{
  std::ifstream is(FileName.c_str(), std::ios::in); // Opens the file

  if (!is)                                          // File not found...
    throw Ejcampdx(Ejcampdx::FILE_NOT_FOUND);       // ...throws an exception and quits.

  std::vector<double> x, y;       // Vectors for data
  std::string Line;               // Each line of the file
  std::istringstream iss;         // std::string stream to read the data from

  while (!is.eof())               // Go through the file
  {
    getline(is, Line);            // Read one line

    if (!Line.empty())            // First character is not '#' and line not empty
    {
      iss.clear();                // Clear the std::string stream
      iss.str(Line);              // Load the line into the std::string stream

      std::string vs;
      double vn;
      std::getline(iss, vs, ',');
      vn = std::strtod(vs.c_str(),NULL);
      x.push_back(vn);

      std::getline(iss, vs, ',');
      vn = std::strtod(vs.c_str(),NULL);
      y.push_back(vn);
    }
  }

  is.close();                       // Closes the file

  // Sorting: the initial assumption is that the file is sorted in ascending or descending frequencies.
  // So, we spend some time checking for this before attempting to sort the vectors.
  // If the last frequency is smaller than the first, there is a chance that the file was sorted in descending order.
  // We reverse the vectors.
  if ( x[0] > x[x.size()-1] )
  {
    std::reverse(x.begin(), x.end());
    std::reverse(y.begin(), y.end());
  }

  bool Sorted = true;       // Assumes that the file is now sorted but let's make sure.
  for (unsigned i = 1; i < x.size(); i++)
    if (x[i] < x[i-1])      // The file was not sorted
    {
      Sorted = false;
      break;
    }

  if (!Sorted)              // If needed sorts x and changes y accordingly
    SortTwoVectors(x, y);

  try { LoadData(x,y); }          // Tries to load a (x,y) pair of vectors into dx object.
  catch (Ejcampdx e) { throw; }   // On error passes the buck.

  jcUnits = un_Wavenumber;
  jcXTitle = "Wavenumbers <sup>-1</sup>";
  jcYTitle = "Optical response";
}
//----------------------------------------------------------------------------------------
// Import from Qenga format
//----------------------------------------------------------------------------------------
void Tjcampdx::Import_qng(const std::string &FileName)     
{
  std::ifstream is(FileName.c_str(), std::ios::in); // Opens the file

  if (!is)                                      // File not found...
    throw Ejcampdx(Ejcampdx::FILE_NOT_FOUND);   // ...throws an exception and quits.

  std::string buffer;
  unsigned int n, unit;
  double v;

  getline(is, buffer);                        // Read the first line...
  if (GetKey(buffer) != "COMMENT")            // ...if it is not a comment...
    throw Ejcampdx(Ejcampdx::BAD_QNG_FORMAT); // ...throws an exception.
  jcComment = buffer;                         // Copies the comment to the appropriate location

  getline(is, buffer);                        // Read the first line...
  if (GetKey(buffer) != "X AXIS")             // ...if it is not the x axis title...
    throw Ejcampdx(Ejcampdx::BAD_QNG_FORMAT); // ...throws an exception.
  jcXTitle = buffer;                          // Copies the comment to the appropriate location

  getline(is, buffer);                        // Read the first line...
  if (GetKey(buffer) != "Y AXIS")             // ...if it is not the x axis title...
    throw Ejcampdx(Ejcampdx::BAD_QNG_FORMAT); // ...throws an exception.
  jcYTitle = buffer;                          // Copies the comment to the appropriate location

  is >> jcx0 >> jcdx >> unit >> n;            // Read the x0, dx, units and size.

  // Converts the units from qng to dx format
  switch (unit)
  {
    case 0x00:
    case 0x01:
      jcUnits = un_Wavenumber; break;
    case 0x02:
      jcUnits = un_eV; break;
    case 0x04:
      jcUnits = un_meV; break;
    case 0x08:
      jcUnits = un_THz; break;
    case 0x10:
    default:
      jcUnits = un_NotDefined; break;
  }

  jcDataX.resize(0);                      // Makes sure that (unused) x vector is empty
  jcDataY.clear();                        // Clears the y container
  jcDataY.reserve(n);                     // Reserves space to avoid memory re-allocation
  for (unsigned int i = 0; i < n; i++)
  {
    is >> v;                              // Reads each value...
    jcDataY.push_back(v);                 // ...and puts it into the container
  }

  is.close();                             // Closes the file

  jcConstSpacing = true;                  // qng format is always constant spacing
  FindExtremes();                         // Looks for data extreme values
  jcxf = Extr->x_max;                     // Sets the last x value
  jcYFactor = 1;                          // Sets the multiplying factor for y
}
//----------------------------------------------------------------------------------------
// Returns the index of value x (or index just below)
//--------- -------------------------------------------------------------------------------
int Tjcampdx::Index(double x)
{
  if (jcConstSpacing)
  {
    double ix = jcamp::ZERO + (x - jcx0) / jcdx;  // Avoids truncation error -- only works for ix > 0
    return int(round(ix - 0.5));                  // Calculates the position <= x
  }
  else
    return Locate(jcDataX, x);
//  return NumRecipes::Locate(jcDataX, x);
}
//----------------------------------------------------------------------------------------
// Interpolates (linearly) current data to find value at position x
//----------------------------------------------------------------------------------------
double Tjcampdx::Interpolate(double x) const
{
  double dydx, yy, x1;
  int Index;

  if (jcConstSpacing)
  {
    Index = int(double( (x - jcx0) / jcdx ));  // Index of point x
    if (Index >= jcDataY.size() - 1) Index = jcDataY.size() - 2;
    dydx = (jcDataY[Index + 1] - jcDataY[Index]) / jcdx;
    x1 = jcx0 + Index * jcdx;
  }
  else
  {
    Index = int((x - jcx0) * jcDataX.size() / (jcxf- jcx0));  // Guesses the location of x assuming that the array is approximatelly evenly spaced
    Index = Hunt(jcDataX, x, Index);                          // Finds the real location
    if (Index < 0) Index = 0;
    if (Index >= jcDataY.size() - 1) Index = jcDataY.size() - 2;
    dydx = (jcDataY[Index+1] - jcDataY[Index]) / (jcDataX[Index+1] - jcDataX[Index]);
    x1 = jcDataX[Index];
  }

  yy = dydx * (x - x1) + jcDataY[Index];
  return yy;
}
//----------------------------------------------------------------------------------------
// Loads a y vector having evenly dx spaced abscissas starting at x0
//----------------------------------------------------------------------------------------
void Tjcampdx::LoadData(double x0, double dx, const std::vector<double> &y)
{
  jcConstSpacing = true;    // Data has constant x spacing
  
  jcDataX.resize(0);        // Makes sure that unused x vector does not occupy the memory
  jcDataY.assign(y.begin(), y.end());     // Copies the data to y vector

  jcx0 = x0;                // Sets the first x value...
  jcdx = dx;                // ...and abscissas constant increment.

  FindExtremes();           // Looks for data extreme values

  jcxf = Extr->x_max;       // Sets the value of the last x
  jcYFactor = 1;              // Sets the mulitplying factor for y
}
//----------------------------------------------------------------------------------------
// Loads a pair of (x,y) vectors into the dx object
//----------------------------------------------------------------------------------------
void Tjcampdx::LoadData(const std::vector<double> &x, const std::vector<double> &y)
{
  if (x.size() != y.size())                       // If vectors have different sizes...
    throw Ejcampdx(Ejcampdx::XY_NO_SIZE_MATCH);   // ...throws an exception.

  jcConstSpacing = false;     // Data will be treated as unevenly spaced
  
  jcDataX.assign(x.begin(), x.end());   // Copies the x values to x vector
  jcDataY.assign(y.begin(), y.end());   // Copies the y values to y vector

  FindExtremes();     // Looks for data extreme values

  jcx0 = Extr->x_min;   // Sets the first x value
  jcxf = Extr->x_max;   // Sets the last x value
  jcdx = 0;                     // dx has no meaning in arbitrarily spaced data
  jcYFactor = 1;                // Sets the multiplying factor for y
}
//----------------------------------------------------------------------------------------
// Calculates factor: <jc> / <this>; where <...> is the average y from x0 to x1
//----------------------------------------------------------------------------------------
double Tjcampdx::MatchFactor(const Tjcampdx &jc, double x0, double x1)
{
  Tjcampdx jc1, jc2;          // Variables to hold data between x0 and x1
  double ave1, ave2;          // Average of data in the x0 to x1 interval

  try { jc1 = jcamp::Cut(*this, x0, x1);}   // Loads local data into jc1
  catch(Ejcampdx e) { throw; }              // The range is not good, throws an exception

  try { jc2 = jcamp::Cut(jc, x0, x1);}      // Loads jc data into jc2
  catch(Ejcampdx e) { throw; }              // The range is not good, throws an exception

  // Calculates average of local data in the x0 to x1 range
  ave1 = 0;
  for (Iter it = jc1.jcDataY.begin(); it != jc1.jcDataY.end(); it++)
    ave1 += (*it);
  ave1 /= jc1.jcDataY.size();

  // Calculates average of jc in the x0 to x1 range
  ave2 = 0;
  for (Iter it = jc2.jcDataY.begin(); it != jc2.jcDataY.end(); it++)
    ave2 += (*it);
  ave2 /= jc2.jcDataY.size();

  return (ave2 / ave1);     // Returns the matching factor
}
//----------------------------------------------------------------------------------------
// Open one data file
//----------------------------------------------------------------------------------------
void Tjcampdx::Open(const std::string &FileName)
{
  std::ifstream FileIn(FileName.c_str());         // Defines and opens the file.
  if (!FileIn)                                    // File not found...
    throw Ejcampdx(Ejcampdx::FILE_NOT_FOUND);     // ...throws an exception and quits.

  try { FileIn >> *this; }                        // Tries to load the data. On error...
  catch (Ejcampdx e) { throw; }                   // ...passes the buck.

  FileIn.close();                                 // Closes the file

  FindExtremes();                                 // Finds data extreme values
}
//----------------------------------------------------------------------------------------
// Saves one data file
//----------------------------------------------------------------------------------------
void Tjcampdx::Save(const std::string &FileName) const
{
  std::ofstream FileOut(FileName.c_str());        // Defines and opens the file.
  if (!FileOut)                                   // Impossible to create file...
    throw Ejcampdx(Ejcampdx::CANNOT_CREATE_FILE); // ...throws an exception and quits.

  try { FileOut << (*this); }                     // Tries to save the data. On error...
  catch (Ejcampdx e) { throw; }                   // ...passes the buck.    

  FileOut.close();                                // Closes the file
}

//----------------------------------------------------------------------------------------
// Keeps the values ranging from index i0 to i1 (included) in vector
//----------------------------------------------------------------------------------------
template<class T>
void Tjcampdx::VectorKeep(std::vector<T> &v, unsigned i0, unsigned i1)
{
  unsigned i_0 = i1 > i0 ? i0 : i1,
           i_1 = i0 < i1 ? i1 : i0;

  v.erase(v.begin() + i_1 + 1, v.end());       // Erases the end of the vector -- i1 stays
  v.erase(v.begin(), v.begin() + i_0);         // Erases the beginning of the vector -- i0 stays
}


//----------------------------------------------------------------------------------------
// Returns index of closest point to P within the box defined by +/dP. If no point is
// found returns the size of the Y vector (which is above the index of the last point).
//----------------------------------------------------------------------------------------
unsigned Tjcampdx::Search(double x, double y, double dx, double dy)
{
  int nX0 = Index(x - dx),      // Index of point at or below x - dx
      nX1 = Index(x + dx) + 1;  // Index of point at or above x + dx

  if (nX0 > nX1) std::swap(nX0,nX1);

  if (nX1 <= 0 || nX0 >= Size())        // The search range is out of data bounds
    return Size();                      // Returns the no point found value

  if (nX0 < 0) nX0 = 0;                 // We know that nX1 is inside data range. Makes sure that nX0 is also.
  if (nX1 >= Size()) nX1 = Size();      // We know that nX0 is inside data range. Makes sure that nX1 is also.

  double yMax = y + dy,         // Upper bound of y box
         yMin = y - dy;         // Lower bound of y box

  if (yMin > yMax) std::swap(yMin,yMax);   // Puts y bounds in ascending order


  if (nX0 == nX1)                                                       // Upper and lower boundaries are the same then...
    return (jcDataY[nX0] > yMax || jcDataY[nX0] < yMin) ? Size() : nX0; // ...the point is found if y is also within boundaries.

  unsigned iRet;                    // Index to be returned by the function.
  double ddx, ddy, mind2, d2;         // Variables to calculate distances.

  while ((jcDataY[nX0] > yMax || jcDataY[nX0] < yMin) && nX0++ <= nX1); // Looks for first y within boundaries
  if (nX0 > nX1) return Size();      // No point found within boundaries

  iRet = nX0;                 // Starts with the first point in the box
  ddx = jcDataX[iRet] - x;       // Difference in the x direction
  ddy = jcDataY[iRet] - y; // Difference in the y direction
  mind2 = ddx * ddx + ddy * ddy;  // Distance (squared) from data to clicked point

  for (int i = nX0 + 1; i < nX1; i++)      // Browses through the other points
    if (jcDataY[i] < yMax && jcDataY[i] > yMin)   // Only proceeds with points having y within boundaries
    {
      ddx = jcDataX[i] - x;            // Difference in the x direction
      ddy = jcDataY[i] - y;             // Difference in the y direction
      d2 = ddx * ddx + ddy * ddy;       // Distance (squared) from data to clicked point

      if (d2 < mind2)               // If the distance decreased, updates the return index.
      {
        mind2 = d2;
        iRet = i;
      }
    }

  return iRet; // The point is found
}
//----------------------------------------------------------------------------------------
// Forces the data to become arbitrarily spaced
//----------------------------------------------------------------------------------------
void Tjcampdx::ArbSpacing()
{
  if (!jcConstSpacing)    // Nothing do do
    return;

  double x = jcx0;        // Running x initialized to first x

  jcDataX.resize(jcDataY.size());       // Makes room for x values

  // Loads x with appropriate values
  for (Iter Item = jcDataX.begin(); Item != jcDataX.end(); Item++)
  {
    *Item = x;
    x += jcdx;
  }

  jcConstSpacing = false;   // Sets arbitrary space flag
  jcdx = 0;                 // dx has no meaning in evenly spaced data
  FindExtremes();           // Finds data extreme values
}
//----------------------------------------------------------------------------------------
// Changes the range of the data without changing any other data property
//----------------------------------------------------------------------------------------
void Tjcampdx::Cut(double x0, double x1)
{
  unsigned LB, UB;                // Indexes for sub set past high limit

  try { FindBounds(x0, x1, LB, UB); } // Looks for the indexes. If there is an error...
  catch(Ejcampdx e) { throw; }          // ...passes the buck.

  if (jcConstSpacing)       // Evenly spaced data
  {
    jcxf = x(UB);             // Updates last x first as it depends on first x.
    jcx0 = x(LB);             // Updates first x only after updating last x.
    VectorKeep<double>(jcDataY, LB, UB);    // Keeps data between LB and UB (included)
  }
  else                      // Arbitrarily spaced data
  {
    VectorKeep<double>(jcDataX, LB, UB);  // Keeps data between LB and UB (included) in both X...
    VectorKeep<double>(jcDataY, LB, UB);  // ...and Y vectors.
    jcx0 = jcDataX.front();             // Updates first x...
    jcxf = jcDataX.back();              // ...and last x.
  }

  FindExtremes();         // Find extreme values
}
//----------------------------------------------------------------------------------------
// Filters the data set to fit y values between boundaries defined by High and Low.
//----------------------------------------------------------------------------------------
void Tjcampdx::Filter(double Low, double High)
{
  for (unsigned int i = 0; i < jcDataY.size(); i++)
  {
    if (jcDataY[i] > High) jcDataY[i] = High;
    if (jcDataY[i] < Low) jcDataY[i] = Low;
  }
  FindExtremes();           // Finds data extreme values
}
//----------------------------------------------------------------------------------------
// Makes data of the same range and resolution as jc 
//----------------------------------------------------------------------------------------
void Tjcampdx::Morph(const Tjcampdx &jc)
{
  int Continue = EquiRange(jc);

  if (Continue == EQUIVALENT) return;       // Nothing to do
  
  if (Continue == SUPER_SET || Continue == ABOVE_OVERLAP || 
      Continue == BELOW_OVERLAP || Continue == NO_OVERLAP)      // Cannot morph data
      throw Ejcampdx(Ejcampdx::DATA_NOT_COMPATIBLE);            // Throws an exception

  // Getting here, all we have to do is convert the data to have the same properties of jc
  if (jc.IsEvenSpaced())          // jc is evenly spaced
  {
    try { Resize(jc.jcx0, jc.jcxf, jc.Size()); }  // Tries to resize the data to match jc...
    catch(Ejcampdx e) { throw; }                  // ...on error passes the buck.
  }
  else                          // jc has arbitrary spacing
  {
    std::vector<double> YY(jc.Size());
    
    for (unsigned i = 0; i < YY.size(); i++)
      YY[i] = Interpolate(jc.jcDataX[i]);       // Calculates y values for new x's

    jcDataX = jc.jcDataX;         // Loads new vectors into object
    jcDataY = YY;

    jcConstSpacing = false;       // Not constant spacing 
    jcx0 = jcDataX.front();       // Sets the first x
    jcxf = jcDataX.back();        // Sets the last x
    jcdx = 0;                     // dx has no meaning
  }

  FindExtremes();       // Looks for data extreme values
}

//----------------------------------------------------------------------------------------
// Generates a pseudo random number in the -1 to 1 range. Call with i = 0 to initialize
// the generator. Call with i < 0 to provide a user defined seed. Call with no parameters
// (or i positive) to get the next number in the series.
//----------------------------------------------------------------------------------------
double Tjcampdx::Random(int i)
{
  static unsigned long j;         // j is made static to be preserved between calls

  if (i == 0)                     // Initializes random number generator
    j = time(NULL) % 65536;       // Uses time mod 2^16
  else if (i < 0)                 // Allows for "debug" initialization...
    j = -i;                       // ...with a user provided value.

  j = 1664525L * j + 1013904223L;       // Gets a "random" number between 0 and 2^32
  return double(j)/2147483648.0 - 1.0;  // Converts it to a number in the interval (-1,1).
}

//----------------------------------------------------------------------------------------
// Adds noise. Amplitude is in %. Prop = (true,false) gives noise ~ [y, max(y)]
//----------------------------------------------------------------------------------------
void Tjcampdx::Noisify(double Amplitude, bool Prop)
{
  double A = Amplitude / 100;
  Iter Item = jcDataY.begin();
  Random(0);      // Initialize random number generator

  if (Prop)
    while (Item != jcDataY.end()) *Item++ = (*Item) * (1 + A * Random());
  else
  {
    A *= std::max(std::abs(Extr->y_min),std::abs(Extr->y_max));
    while (Item != jcDataY.end()) *Item++ = (*Item) + A * Random();
  }
  FindExtremes();           // Finds data extreme values
}
//---------------------------------------------------------------------------------------- 
// Resizes, but does not extrapolate, the data. In linear scale, interpolates to place n 
// points between x0 and x1. The output data is evenly spaced. In log scale the output is 
// arbitrarily spaced.
//----------------------------------------------------------------------------------------
void Tjcampdx::Resize(double x_0, double x_1, unsigned int n, bool LogScale)
{
  if (n == 0)                               // No points in extrapolation...
    throw Ejcampdx(Ejcampdx::EMPTY_DATA);   // ...throws an exception.

  double x0, x1, dx, x;                 // Variables for resizing range
	std::vector<double> YY, XX;           // Will hold new values

  x0 = x_0 < jcx0 ? jcx0 : x_0;         // Makes sure that lower bound is within data range
  x1 = x_1 > jcxf ? jcxf : x_1;         // Makes sure that upper bound is within data range

  if (LogScale)     // Log scale resizing
  {
    if (x0 < 0 || x1 < 0)                         // If one of the values is negative...
      throw Ejcampdx(Ejcampdx::NO_NEGATIVE_LOG);  // ...aborts the operation.

    if (n < 2)                                    // No points in extrapolation...
      throw Ejcampdx(Ejcampdx::EMPTY_DATA);       // ...throws an exception.

    x0 = (x0 == 0) ? 0.01 * x1 : x0;    // If first x is zero, sets it to 1% of last x

    dx = exp(log(x1 / x0) / double(n - 1) );  // Calculate log increment [x(i) = k * x(i-1)]

    XX.clear();
    x = x0;                             // Running variable
    for (unsigned i = 0; i < n; i++)  
    {
      XX.push_back(x);
      YY.push_back(Interpolate(x));    // Interpolates the data at each point
      x *= dx;                         // Log increment of the running variable
    }

    jcDataX.assign(XX.begin(),XX.end());
    jcConstSpacing = false; // Makes sure that the data is flagged as arbitrary spacing
    jcdx = 0;               // dx has no meaning

  }
  else              // Linear scale resizing
  {
    dx = (x1 - x0) / (n - 1);           // Sets dx value

    if (n == 0)                               // No points in extrapolation...
      throw Ejcampdx(Ejcampdx::EMPTY_DATA);   // ...throws an exception.

    x = x0;                             // Running variable
    for (unsigned i = 0; i < n; i++)  
    {
      YY.push_back(Interpolate(x));     // Interpolates the data at each point
      x += dx;                          // Increments the running variable
    }

    jcConstSpacing = true;  // Makes sure that the data is flagged as constant spacing
    jcDataX.resize(0);      // Erases unused x vector.
    jcdx = dx;              // Updates the increment.
  }

  jcx0 = x0;              // Updates the initial x.
  jcxf = x1;              // Updates the final x.
  jcDataY.assign(YY.begin(), YY.end());      // Loads tmp into y
  FindExtremes();         // Finds extreme values
}
//----------------------------------------------------------------------------------------
// Resizes the data block. Keeps x0 and xf and changes the number of points. Result is
// evenly spaced data.
//----------------------------------------------------------------------------------------
void Tjcampdx::Resize(unsigned int n)
{
  // Checks that the number of points is enough
  if (n < 2)
    throw Ejcampdx(Ejcampdx::EMPTY_DATA);

  try { Resize(jcx0, jcxf, n); }              // Resizes vector. On error...
  catch(Ejcampdx e) { throw; }                // ...passes the buck.
}
//----------------------------------------------------------------------------------------
// Savitzky-Golay smoothing
//----------------------------------------------------------------------------------------
void Tjcampdx::Smooth(unsigned WinSize, unsigned Order)
{
  /*std::vector<double> Smoothed;

  try 
    { NumRecipes::SavGolaySmoothing(jcDataY, Smoothed, WinSize, Order); }   // Tries to smooth the data
  catch (ENumRecipes e)
  { 
    Ejcampdx ej(Ejcampdx::GENERIC_MESSAGE);               // Creates a JCamp generic exception
    ej.SetMessage(e.Message());                           // Loads the numerical recipes message into the JCamp's
    throw ej;                                             // Throws a JCamp exception based on the numerical
  }

  jcDataY = Smoothed;    // Loads smoothed data into return variable                        
  FindExtremes();           // Finds data extreme values
  */
}
//----------------------------------------------------------------------------------------
// Converts the x units using a factor (> 0 -> k * x; < 0 -> k / x) -- See PhysConst.h
// The function does not change the type or the name of the unit. This has to be done
// separately using the Units(int) and XTitle(std::string) functions.
//----------------------------------------------------------------------------------------
void Tjcampdx::UnitConv(double ConvFactor)
{
  double k = fabs(ConvFactor);

  if (ConvFactor > 0)           // Linear conversion in x axis
  {
    if (jcConstSpacing)
    {
      jcx0 *= k;                            // Renormalizes first x
      jcdx *= k;                            // Renormalizes x increment
      jcxf = jcx0 + (Size() - 1) * jcdx;    // Calculates last x
    }
    else
    {
      for (Iter it = jcDataX.begin(); it != jcDataX.end(); it++)
        *it *= k;                                       // Calculates k / x

      jcx0 *= k;
      jcxf *= k;
    }
  }
  else if (ConvFactor < 0)      // Inverse conversion in x axis
  {
    ArbSpacing();                                     // If necessary makes the data arbitrarily spaced

    for (Iter it = jcDataX.begin(); it != jcDataX.end(); it++)
      *it = k / (*it);                                // Calculates k / x

    std::reverse(jcDataX.begin(), jcDataX.end());     // Reverses x...
    std::reverse(jcDataY.begin(), jcDataY.end());     // ...and y vectors.

    jcx0 = jcDataX.front();
    jcxf = jcDataX.back();
  }
  else
    throw Ejcampdx(Ejcampdx::ZERO_CONVERSION);        // Conversion factor error

  FindExtremes();         // Find new extreme values
}
//----------------------------------------------------------------------------------------
// Sets the abcissa units of the data
//----------------------------------------------------------------------------------------
void Tjcampdx::Units(const std::string &Line)
{
  std::string s = Line;
  std::transform(s.begin(), s.end(), s.begin(), ::tolower);
  if (s == "1/cm") jcUnits = un_Wavenumber;
  else if (s == "ev") jcUnits = un_eV;
  else if (s == "mev") jcUnits = un_meV;
  else if (s == "thz") jcUnits = un_THz;
  else if (s == "um") jcUnits = un_Microns;
  else if (s == "nm") jcUnits = un_nm;
  else jcUnits = un_NotDefined;
}
//----------------------------------------------------------------------------------------
// Gets the abcissa units of the data
//----------------------------------------------------------------------------------------
const char *Tjcampdx::Units() const
{
   if (jcUnits == un_Wavenumber) return "1/cm";
   else if (jcUnits == un_eV) return "eV";
   else if (jcUnits == un_meV) return "meV";
   else if (jcUnits == un_THz) return "THz";
   else if (jcUnits == un_Microns) return "um";
   else if (jcUnits == un_nm) return "nm";
   else return "NotDef";
}
//----------------------------------------------------------------------------------------
// Gets the abcissa units of the data
//----------------------------------------------------------------------------------------
const char *Tjcampdx::UnitsHTML() const
{
   if (jcUnits == un_Wavenumber) return "Wavenumber [cm<sup>-1</sub>]";
   else if (jcUnits == un_eV) return "Energy [eV]";
   else if (jcUnits == un_meV) return "Energy [meV]";
   else if (jcUnits == un_THz) return "Frequency [THz]";
   else if (jcUnits == un_Microns) return "Wavelength [&mu;m]";
   else if (jcUnits == un_nm) return "Wavelength [nm]";
   else return "NotDef";
}
//----------------------------------------------------------------------------------------
// Modulus -- Overwrites data
//----------------------------------------------------------------------------------------
void Tjcampdx::ABS()         
{
  Iter Item = jcDataY.begin();
  while (Item != jcDataY.end()) *Item++ = std::abs(*Item);
  FindExtremes();           // Finds data extreme values
}
//----------------------------------------------------------------------------------------
// Arc-cosinus -- Overwrites data
//----------------------------------------------------------------------------------------
void Tjcampdx::ACOS()
{
  Iter Item = jcDataY.begin();
  while (Item != jcDataY.end()) *Item++ = std::acos(*Item);
  FindExtremes();           // Finds data extreme values
}
//----------------------------------------------------------------------------------------
// Arc-sinus -- Overwrites data
//----------------------------------------------------------------------------------------
void Tjcampdx::ASIN()
{
  Iter Item = jcDataY.begin();
  while (Item != jcDataY.end()) *Item++ = std::asin(*Item);
  FindExtremes();           // Finds data extreme values
}
//----------------------------------------------------------------------------------------
// Arc-tangent -- Overwrites data
//----------------------------------------------------------------------------------------
void Tjcampdx::ATAN()
{
  Iter Item = jcDataY.begin();
  while (Item != jcDataY.end()) *Item++ = std::atan(*Item);
  FindExtremes();           // Finds data extreme values
}
//----------------------------------------------------------------------------------------
// Cosinus -- Overwrites data
//----------------------------------------------------------------------------------------
void Tjcampdx::COS()
{
  Iter Item = jcDataY.begin();
  while (Item != jcDataY.end()) *Item++ = std::cos(*Item);
  FindExtremes();           // Finds data extreme values
}
//----------------------------------------------------------------------------------------
// Derivative -- Overwrites data
//----------------------------------------------------------------------------------------
void Tjcampdx::Differentiate()
{
  double Y0 = jcDataY[0],     // Storage for previous Y value
         Y1;                  // Storage for next Y value 

  if (jcConstSpacing)         // Evenly spaced data
  {
    for (unsigned i = 1; i < Size(); i++) 
    {
      Y1 = jcDataY[i];                // Stores current Y
      jcDataY[i] = (Y1 - Y0) / jcdx;  // Derivative
      Y0 = Y1;                        // Preserves current Y for next step
    }
  }
  else                            // Arbitrarily spaced data
  {
    double dx;                              // Spacing between two consecutive x's
    for (unsigned i = 1; i < Size(); i++)   
    {
      Y1 = jcDataY[i];                  // Stores current Y
      dx = jcDataX[i] - jcDataX[i-1];   // Spacing between x(i) & x(i-1)
      jcDataY[i] = (Y1 - Y0) / dx;      // Derivative
      Y0 = Y1;                          // Preserves current Y for next step
    }
  }
  jcDataY[0] = jcDataY[1];      // Derivative at 1st point set to be the same as 2nd
  FindExtremes();           // Finds data extreme values
}
//----------------------------------------------------------------------------------------
// Exponential -- Overwrites data
//----------------------------------------------------------------------------------------
void Tjcampdx::EXP() 
{
  Iter Item = jcDataY.begin();
  while (Item != jcDataY.end()) *Item++ = std::exp(*Item);
  FindExtremes();           // Finds data extreme values
}
//----------------------------------------------------------------------------------------
// 10^y -- Overwrites data
//----------------------------------------------------------------------------------------
void Tjcampdx::EXP10()      
{
  Iter Item = jcDataY.begin();
  while (Item != jcDataY.end()) *Item++ = std::exp((*Item) * 2.302585092994);
  FindExtremes();           // Finds data extreme values
}
//----------------------------------------------------------------------------------------
// Integral -- Overwrites data
//----------------------------------------------------------------------------------------
void Tjcampdx::Integrate()
{
  double Y0 = jcDataY[0],     // Storage for previous Y value
         Y1;                  // Storage for next Y value 
  jcDataY[0] = 0.0;           // Imposes boudary condition

  if (jcConstSpacing)           // Evenly spaced data
    for (unsigned i = 1; i < Size(); i++)
    {
      Y1 = jcDataY[i];                                      // Stores current Y  
      jcDataY[i] = jcDataY[i - 1] + 0.5 * jcdx * (Y0 + Y1); // Trapezoidal integration
      Y0 = Y1;                                              // Preserves current Y for next step
    }
  else                          // Arbitrarily spaced data
  {
    double dx;                              // Spacing between two consecutive x's
    for (unsigned i = 1; i < Size(); i++)   
    {
      Y1 = jcDataY[i];                                      // Stores current Y  
      dx = jcDataX[i] - jcDataX[i-1];                       // Spacing between x(i) & x(i-1)
      jcDataY[i] = jcDataY[i - 1] + 0.5 * dx * (Y0 + Y1);   // Trapezoidal integration
      Y0 = Y1;                                              // Preserves current Y for next step
    }
  }
  FindExtremes();           // Finds data extreme values
}
//----------------------------------------------------------------------------------------
// 1 / y -- Overwrites data
//----------------------------------------------------------------------------------------
void Tjcampdx::Inv()        
{
  Iter Item = jcDataY.begin();
  while (Item != jcDataY.end()) *Item++ = ((*Item) != 0) ? 1.0 / (*Item) : 1E100;
  FindExtremes();           // Finds data extreme values
}
//----------------------------------------------------------------------------------------
// Logarithm -- Overwrites data
//----------------------------------------------------------------------------------------
void Tjcampdx::LOG() 
{
  Iter Item = jcDataY.begin();
  while (Item != jcDataY.end()) *Item++ = log(fabs(*Item));
  FindExtremes();           // Finds data extreme values
}
//----------------------------------------------------------------------------------------
// Base 10 logarithm -- Overwrites data
//----------------------------------------------------------------------------------------
void Tjcampdx::LOG10()       
{
  Iter Item = jcDataY.begin();
  while (Item != jcDataY.end()) *Item++ = std::log10(fabs(*Item));
  FindExtremes();           // Finds data extreme values
}
//----------------------------------------------------------------------------------------
// Sinus -- Overwrites data
//----------------------------------------------------------------------------------------
void Tjcampdx::SIN()
{
  Iter Item = jcDataY.begin();
  while (Item != jcDataY.end()) *Item++ = std::sin(*Item);
  FindExtremes();           // Finds data extreme values
}
//----------------------------------------------------------------------------------------
// Square root -- Overwrites data
//----------------------------------------------------------------------------------------
void Tjcampdx::Sqrt()   
{
  Iter Item = jcDataY.begin();
  while (Item != jcDataY.end()) *Item++ = std::sqrt(*Item);
  FindExtremes();           // Finds data extreme values
}
//----------------------------------------------------------------------------------------
// y^2 -- Overwrites data
//----------------------------------------------------------------------------------------
void Tjcampdx::Square()  
{
  Iter Item = jcDataY.begin();
  while (Item != jcDataY.end()) *Item++ = (*Item) * (*Item);
  FindExtremes();           // Finds data extreme values
}
//----------------------------------------------------------------------------------------
// Tangent -- Overwrites data
//----------------------------------------------------------------------------------------
void Tjcampdx::TAN()
{
  Iter Item = jcDataY.begin();
  while (Item != jcDataY.end()) *Item++ = std::tan(*Item);
  FindExtremes();           // Finds data extreme values
}
//----------------------------------------------------------------------------------------
// Returns the x value of index i
//----------------------------------------------------------------------------------------
double Tjcampdx::x(unsigned int i) const
{
  if (i > Size())
    throw Ejcampdx::OUT_OF_BOUNDS;

  if (jcConstSpacing)
    return (jcx0 + i * jcdx);
  else
    return jcDataX[i];
}
//----------------------------------------------------------------------------------------
// Assignment operator
//----------------------------------------------------------------------------------------
const Tjcampdx &Tjcampdx::operator=(const Tjcampdx &jc)
{
  if (this == &jc) return *this;      // Avoids self assignement

  jcConstSpacing = jc.jcConstSpacing;
  jcx0 = jc.jcx0;
  jcxf = jc.jcxf;
  jcdx = jc.jcdx;
  jcYFactor = jc.jcYFactor;
  jcUnits = jc.jcUnits;
  jcUnevenCompat = jc.jcUnevenCompat;
  jcComment = jc.jcComment;
  jcXTitle = jc.jcXTitle;
  jcYTitle = jc.jcYTitle;

  if (jcConstSpacing)
    jcDataX.resize(0);
  else
    jcDataX = jc.jcDataX;

  jcDataY = jc.jcDataY;

  Extr = new TExtremes;
  *Extr = *(jc.Extr);

  return *this;
}

//----------------------------------------------------------------------------------------
// Operator jc1 + jc2 - friend
//----------------------------------------------------------------------------------------
Tjcampdx operator+(const Tjcampdx &jc1, const Tjcampdx &jc2)
{
  Tjcampdx jca = jc1,         // Copies data to...
           jcb = jc2;         // local variables.

  try { jcamp::MakeRangeEquivalent(jca,jcb); }        // Checks/tries to make data equivalent
  catch(Ejcampdx e) { throw; }                        // On error throws an exception

  for (unsigned i = 0; i < jca.jcDataY.size(); i++)
    jca.jcDataY[i] += jcb.jcDataY[i];                 // Adds the data points

  jca.FindExtremes();           // Finds data extreme values
  return jca;
}
//----------------------------------------------------------------------------------------
// Operator jc1 - jc2 - friend
//----------------------------------------------------------------------------------------
Tjcampdx operator-(const Tjcampdx &jc1, const Tjcampdx &jc2)
{
  Tjcampdx jca = jc1,         // Copies data to...
           jcb = jc2;         // local variables.

  try { jcamp::MakeRangeEquivalent(jca,jcb); }        // Checks/tries to make data equivalent
  catch(Ejcampdx e) { throw; }                        // On error throws an exception

  for (unsigned i = 0; i < jca.jcDataY.size(); i++)
    jca.jcDataY[i] -= jcb.jcDataY[i];                 // Subtracts the data points

  jca.FindExtremes();           // Finds data extreme values
  return jca;
}
//----------------------------------------------------------------------------------------
// Operator jc1 / jc2 - friend
//----------------------------------------------------------------------------------------
Tjcampdx operator/(const Tjcampdx &jc1, const Tjcampdx &jc2)
{
  Tjcampdx jca = jc1,         // Copies data to...
           jcb = jc2;         // local variables.

  try { jcamp::MakeRangeEquivalent(jca,jcb); }        // Checks/tries to make data equivalent
  catch(Ejcampdx e) { throw; }                        // On error throws an exception

  double den;
  for (unsigned i = 0; i < jca.jcDataY.size(); i++)
  {
    den = jcb.jcDataY[i] != 0.0 ? jcb.jcDataY[i] : -1.0E-8;   // Avoids div by 0
    jca.jcDataY[i] /= den;                                    // Makes the division
  }

  jca.FindExtremes();           // Finds data extreme values
  return jca;
}
//----------------------------------------------------------------------------------------
// Operator jc1 * jc2 - friend
//----------------------------------------------------------------------------------------
Tjcampdx operator*(const Tjcampdx &jc1, const Tjcampdx &jc2)
{
  Tjcampdx jca = jc1,         // Copies data to...
           jcb = jc2;         // local variables.

  try { jcamp::MakeRangeEquivalent(jca,jcb); }        // Checks/tries to make data equivalent
  catch(Ejcampdx e) { throw; }                        // On error throws an exception

  for (unsigned i = 0; i < jca.jcDataY.size(); i++)
    jca.jcDataY[i] *= jcb.jcDataY[i];                 // Multiplies the data points

  jca.FindExtremes();           // Finds data extreme values
  return jca;
}
//----------------------------------------------------------------------------------------
// Operator JCamp + constant - friend
//----------------------------------------------------------------------------------------
Tjcampdx operator+(const Tjcampdx &jc, double a)
{
  Tjcampdx jc1 = jc;
  for (Tjcampdx::Iter it = jc1.jcDataY.begin(); it != jc1.jcDataY.end(); it++)
    *it += a;
  jc1.FindExtremes();           // Finds data extreme values
  return jc1;
}
//----------------------------------------------------------------------------------------
// Operator constant + JCamp - friend
//----------------------------------------------------------------------------------------
Tjcampdx operator+(double a, const Tjcampdx &jc)
{
  Tjcampdx jc1 = jc;
  for (Tjcampdx::Iter it = jc1.jcDataY.begin(); it != jc1.jcDataY.end(); it++)
    *it += a;
  jc1.FindExtremes();           // Finds data extreme values
  return jc1;
}
//----------------------------------------------------------------------------------------
// Operator JCamp - constant - friend
//----------------------------------------------------------------------------------------
Tjcampdx operator-(const Tjcampdx &jc, double a)
{
  Tjcampdx jc1 = jc;
  for (Tjcampdx::Iter it = jc1.jcDataY.begin(); it != jc1.jcDataY.end(); it++)
    *it -= a;
  jc1.FindExtremes();           // Finds data extreme values
  return jc1;
}
//----------------------------------------------------------------------------------------
// Operator constant - JCamp - friend
//----------------------------------------------------------------------------------------
Tjcampdx operator-(double a, const Tjcampdx &jc)
{
  Tjcampdx jc1 = jc;
  for (Tjcampdx::Iter it = jc1.jcDataY.begin(); it != jc1.jcDataY.end(); it++)
    *it = a - (*it);
  jc1.FindExtremes();           // Finds data extreme values
  return jc1;
}
//----------------------------------------------------------------------------------------
// Operator JCamp * constant - friend
//----------------------------------------------------------------------------------------
Tjcampdx operator*(const Tjcampdx &jc, double a)
{
  Tjcampdx jc1 = jc;
  for (Tjcampdx::Iter it = jc1.jcDataY.begin(); it != jc1.jcDataY.end(); it++)
    *it *= a;
  jc1.FindExtremes();           // Finds data extreme values
  return jc1;
}
//----------------------------------------------------------------------------------------
// Operator constant * JCamp - friend
//----------------------------------------------------------------------------------------
Tjcampdx operator*(double a, const Tjcampdx &jc)
{
  Tjcampdx jc1 = jc;
  for (Tjcampdx::Iter it = jc1.jcDataY.begin(); it != jc1.jcDataY.end(); it++)
    *it *= a;
  jc1.FindExtremes();           // Finds data extreme values
  return jc1;
}
//----------------------------------------------------------------------------------------
// Operator JCamp * constant - friend
//----------------------------------------------------------------------------------------
Tjcampdx operator/(const Tjcampdx &jc, double a)
{
  Tjcampdx jc1 = jc;
  double den = a != 0 ? a : -1.0E-8;    // Avoids div by 0
  for (Tjcampdx::Iter it = jc1.jcDataY.begin(); it != jc1.jcDataY.end(); it++)
    *it /= den;
  jc1.FindExtremes();           // Finds data extreme values
  return jc1;
}
//----------------------------------------------------------------------------------------
// Operator constant * JCamp - friend
//----------------------------------------------------------------------------------------
Tjcampdx operator/(double a, const Tjcampdx &jc)
{
  Tjcampdx jc1 = jc;
  double den;
  for (Tjcampdx::Iter it = jc1.jcDataY.begin(); it != jc1.jcDataY.end(); it++)
  {
    den = (*it) != 0.0 ? *it : -1.0E-8;   // Avoids div by 0
    *it = a / den;
  }
  jc1.FindExtremes();           // Finds data extreme values
  return jc1;
}
//----------------------------------------------------------------------------------------
// Stream out - friend
//----------------------------------------------------------------------------------------
std::ostream &operator<<(std::ostream &dataout, const Tjcampdx &jc)
{
  dataout.precision(15);               // Max precision for double

  long int long_y;            // Data will be exported as long int
  double y1, y2, yfactor;     // Variables to calculate multiplying y factor

  y1 = fabs(jc.Extremes().y_min);            // Modulus of lowest y extreme
  y2 = fabs(jc.Extremes().y_max);            // Modulus of highest y extreme
  yfactor = (y1 > y2) ? 2.0E9 / y1 : 2.0E9 / y2;  // yfactor is such that largest |y| = 2E9
  if (yfactor == 0) yfactor = 1;            // All numbers are zero, so yfactor is meaningless

  // Fields common for all data formats
  dataout << "##TITLE=" << jc.jcComment << endl
          << "##JCAMP-DX=4.24" << endl
          << "##XUNITS=" << jc.Units() << endl
          << "##YUNITS=" << jc.jcYTitle << endl
          << "##YFACTOR=" << (jc.jcConstSpacing ? 1.0/yfactor : 1) << endl
          << "##FIRSTX=" << jc.jcx0 << endl
          << "##LASTX=" << jc.jcxf << endl
          << "##NPOINTS=" << jc.Size() << endl;

  if (jc.jcConstSpacing)        // Data is in constant spacing format
  {
    dataout << "##DELTAX=" << jc.jcdx << endl
            << "##XYDATA=(X++(Y..Y))";          // endl is inside the data loop

    // Saves the data
    double xx = jc.jcx0;
    for (unsigned int i = 0; i < jc.Size(); i++)
    {
      if (!(i % 10))            // Went through 10 values.
      {
        dataout << endl << long(xx);  // New line with integer version of old x
        xx += 10 * jc.jcdx;           // New x
      }

      long_y = long(yfactor * jc.jcDataY[i]);   // Converts the value to a large integer
      if (long_y >= 0)
        dataout << "+";           // Makes sure to save the positive sign.

      dataout << long_y;          // Saves the data.
    }
  }
  else                        // Data has arbitrary x's
  {
    dataout << "##XYDATA=(XY..XY)";
    for (unsigned int i = 0; i < jc.Size(); i++)
      dataout << endl << jc.jcDataX[i] << '\t' << jc.jcDataY[i];
  }

  dataout << endl << "##END=" << endl;

  return dataout;
}
//----------------------------------------------------------------------------------------
// Stream in - friend
//----------------------------------------------------------------------------------------
std::istream &operator>>(std::istream &datain, Tjcampdx &jc)
{
  bool EndFound = false;                        // "##END=" flag found
  std::vector <std::string> FileGuts;           // Vector to hold the file contents
  std::string OneLine;                          // Each line of the file

  while (!EndFound && !datain.eof())            // Loops until the stream is all read
  {
    getline(datain, OneLine);                   // Read each line...
    if (!OneLine.empty())                       // ...and if it is not empty...
    {
      FileGuts.push_back(OneLine);              // ...adds it to the vector.

      // Checks if the line is the end of data flag
      if (jc.GetKey(OneLine) == "##END")
          EndFound = true;
    }
  }

  if (!EndFound)                                // File does not contain END tag...
    throw Ejcampdx(Ejcampdx::END_NOT_FOUND);    // ...throws an exception and quits.

  try                               // Tries to...
    { jc.ParseJCampDX(FileGuts); }  // ...parses the file contents.
  catch (Ejcampdx e)                // If there is a problem...
    { throw; }                      // ...passes the buck.

  return datain;
}
//----------------------------------------------------------------------------------------
// Finds upper and lower indexes for sub-set in the range [x_0,x_1]
//----------------------------------------------------------------------------------------
void Tjcampdx::FindBounds(double x_0, double x_1, unsigned &LB, unsigned &UB) const
{
  double xi = x_1 > x_0 ? x_0 : x_1,   // Makes sure that desired low and...
         xf = x_0 < x_1 ? x_1 : x_0;   // ...high limits are in ascending order.

  if (xi > jcxf || xf < jcx0)                 // Checks if limits are fully outside data range
    throw Ejcampdx(Ejcampdx::OUT_OF_BOUNDS);  // Throws an exception on error

  xi = xi < jcx0 ? jcx0 : xi;       // Makes sure that limits are...
  xf = xf > jcxf ? jcxf : xf;       // ...within data range.

  LB = 0;           // Initial index for lower limit
  UB = Size() - 1;  // Initial index for high limit

  if (jcConstSpacing)       // Evenly spaced data
  {
    double x0 = jcx0,       // Starting point for x0...
           x1 = jcxf;       // ...and x1 search.  

    while (x0 < xi)           // Hunts for lower bound index... 
    {
      LB++;
      x0 += jcdx;
    }
    if (LB > 0) LB--;         // ...and tries to get an extra point. 

    while (x1 > xf)             // Hunts for upper bound index...
    {
      UB--;
      x1 -= jcdx;
    }
    if (UB < Size() - 1) UB++;  // ...and tries to get an extra point.
  }
  else                  // Arbitrarily spaced data
  {
    while (jcDataX[LB] < xi) 
      LB++;      // Finds lower bound index...
    if (LB > 0) LB--;                   // ...and tries to get an extra point.
    while (jcDataX[UB] > xf) 
      UB--;      // Finds upper bound index...
    if (UB < Size() - 1) UB++;          // ...and tries to get an extra point.
  }

  if (LB > 0 && fabs(x(LB-1) - xi) < jcamp::ZERO) LB--;                   // Checks if previous or...
  else if (LB != (Size() - 1) && fabs(x(LB+1) - xi) < jcamp::ZERO) LB++;  // ...next x's are a perfect match for LB.

  if (UB > 0 && fabs(x(UB-1) - xf) < jcamp::ZERO) UB--;                   // Checks if previous or...
  else if (UB != (Size() - 1) && fabs(x(UB+1) - xf) < jcamp::ZERO) UB++;  // ...next x's are a perfect match for UB.

  if (UB <= LB) UB = LB + 1;  // Makes sure that UB is not less than LB...
  if (UB == Size()) UB--;     // ...but not out of the range
}
//----------------------------------------------------------------------------------------
// Find extremes values for x and y
//----------------------------------------------------------------------------------------
void Tjcampdx::FindExtremes() const
{
  // Finds x extremes
  if (jcConstSpacing)         // For evenly spaced data
  {
    Extr->x_min = jcx0;
    Extr->x_max = jcx0 + (jcDataY.size() - 1) * jcdx;
  }
  else                        // For unevenly spaced data 
  {
    Extr->x_min = jcDataX.front();
    Extr->x_max = jcDataX.back();
  }
  if (Extr->x_min == Extr->x_max)   // Both extremes are the same
  {
    Extr->x_min = 0.995 * Extr->x_min;    // Sets them 1%
    Extr->x_max = 1.005 * Extr->x_max;    // ...appart
  }

  // Finds y extremes
  Extr->y_min = *std::min_element(jcDataY.begin(), jcDataY.end());
  Extr->y_max = *std::max_element(jcDataY.begin(), jcDataY.end());
  if (Extr->y_min == Extr->y_max)   // Both extremes are the same
  {
    Extr->y_min = 0.995 * Extr->y_min;    // Sets them 1%
    Extr->y_max = 1.005 * Extr->y_max;    // ...appart
  }
}
//----------------------------------------------------------------------------------------
// Loads x0 and x1 with the overlapping regions between jc and *this
//----------------------------------------------------------------------------------------
void Tjcampdx::FindOverlap(const Tjcampdx &jc, double &x0, double &x1) const
{
  switch (EquiRange(jc))
  {
    case EQUIVALENT:      // Evenly spaced and compatible
      return;         

    case SUB_SET:         // Data range is a sub set of this block
    case SAME_RANGE:      // Data are at the same range but do not have the same resolution or spacing type
      x0 = jc.jcx0;
      x1 = jc.jcxf;
      return;

    case SUPER_SET:       // Data range is a super set of this block
      x0 = jcx0;
      x1 = jcxf;
      return;

    case ABOVE_OVERLAP:   // Data are partially overlapping but goes above range
      x0 = jc.jcx0;
      x1 = jcxf;
      return;

    case BELOW_OVERLAP:   // Data are partially overlapping but goes below range
      x0 = jcx0;
      x1 = jc.jcxf;
      return;

    case NO_OVERLAP:      // Data are not compatible - no overlap
      throw Ejcampdx(Ejcampdx::DATA_NOT_COMPATIBLE);    // Throws an exception

    default:              // Should never get here
      throw Ejcampdx(Ejcampdx::UNKNOWN_COMPAT);         // Throws an exception
  }
}

//----------------------------------------------------------------------------------------
// Replace Line by the value attached to the key. Return the key in uppercase.
//----------------------------------------------------------------------------------------
std::string Tjcampdx::GetKey(std::string &Line)
{
  std::string Key;                     // Temporary std::string to be hold the key
  unsigned long Pos = Line.find("=");  // Looks for first occurence of '='

  if (Pos < std::string::npos)         // Found the sign '='
  {
    Key = Line.substr(0, Pos);                        // Put the key in the return value
    Line = Line.substr(Pos + 1, Line.length() - Pos); // Trim the key (and '=') from Line
  }
  trim(Line);           // Trims leading and trailing spaces and control characters
  std::transform(Key.begin(), Key.end(), Key.begin(), ::toupper);  // Converts the key to uppercase characters only
  return Key;                 // Returns the Key
}
//----------------------------------------------------------------------------------------
// Converts Contents buffer to Tjcampdx
//----------------------------------------------------------------------------------------
void Tjcampdx::ParseJCampDX(const std::vector <std::string> &Contents)
{
  // Looks for header location. The header starts at line 0 and goes to DataStart - 1.
  unsigned int DataStart = 0,                 // Assumes that the data starts at 1st line...
               DataEnd = Contents.size() - 1; // ...and uses all file.

  double XX, YY;                                // Aux double values for reading purposes
  std::string tmp;                                   // Aux temporary std::string
  std::istringstream Value;                     // Will hold each value to be read
  std::vector <std::string>::const_iterator it;            // Forward iterator
  std::vector <std::string>::const_reverse_iterator rit;   // Reverse iterator

  // Loops through the buffer (starting at the begining) counting the first n lines that
  // start with a '#' character.
  it = Contents.begin();      // Begining of the buffer
  while (((*it)[0] == '#') && (it != Contents.end()))
  {
    DataStart++;      // Increase the counter indicating the beginning of the data
    it++;             // Goes to the next std::string
  }

  // Loops through the buffer (starting at the end) looking for the line containing the
  // '##END=' flag. This is the line just past the data.
  rit = Contents.rbegin();          // End of the buffer
  tmp = *rit;                       // Gets the std::string in the reverse iterator (needed becauce GetKey changes its argument)
  while ((GetKey(tmp) != "##END") && (rit != Contents.rend()))
  {
    DataEnd--;                      // Decreases counter indicating the end of the data
    rit--;                          // Goes to previous std::string...
    tmp = *rit;                     // ...and updates the tmp variable.
  }

  // Parses the file header
  try { ParseJCampDXHeader(Contents, DataStart); }  // Tries to parse the header. On error...
  catch (Ejcampdx e) { throw; }                     // ...passes the buck.

  // Reads the values into appropriate vectors
  jcDataX.clear();
  jcDataY.clear();
  if (jcConstSpacing)             // Constant spacing -- only Y vector is used
  {
    unsigned long PosA, PosB;
    //char C;                       // Character to test for exponential
    // Loop through all data lines
    for (it = Contents.begin() + DataStart; it != Contents.begin() + DataEnd; it++)
    {
      PosA = it->find_first_of(" +-", 1);           // Looks for +, - or space ignoring first character
      if (PosA != std::string::npos && toupper((*it)[PosA - 1]) == 'E')  // Character preceeding sign was an exponential.
        PosA = it->find_first_of(" +-", PosA + 1);  // As it was not the beginning of a new number, looks for next spacer.

      while (PosA != std::string::npos)                  // Loopping inside each line
      {
        PosB = it->find_first_of(" +-", PosA + 1);  // Looks for next +, - or space
        
        if (PosB != std::string::npos && toupper((*it)[PosB - 1]) == 'E')  // Character preceeding sign was an exponential.
          PosB = it->find_first_of(" +-", PosB + 1);  // As it was not the beginning of a new number, looks for next spacer.
        
        Value.clear();                              // Clears the std::string stream
        Value.str(it->substr(PosA, PosB - PosA));   // Gets one number
        Value >> YY;                                // Reads value from this stream
        jcDataY.push_back(YY * jcYFactor);          // Loads data into container
        PosA = PosB;                                // Update PosA
      }
    }
  }
  else                            // Arbitrary spacing -- both X and Y vectors used
  {
    for (it = Contents.begin() + DataStart; it != Contents.begin() + DataEnd; it++)
    {
      Value.clear();            // Clears the std::string stream
      Value.str(*it);           // Loads a std::string from contents into the stream
      Value >> XX >> YY;        // Reads two values from this stream
      jcDataX.push_back(XX);    // Loads X and...
      jcDataY.push_back(YY * jcYFactor);  // ...(corrected) Y values into containers.
    }
    // Fills out values for completeness
    jcx0 = jcDataX.front();
    jcxf = jcDataX.back();
    jcdx = 0;                   // Set to zero as dx is undefined
  }
}
//----------------------------------------------------------------------------------------
// Analyses the header looking for recognizable keys -- Must have the following headers:
// ##JCAMP-DX which should be a version 4.nnn
// 3 out the following 4: ##FIRSTX, ##LASTX, ##DELTAX and ##NPOINTS -- if all are
// present then ##LASTX is only used for checking but not necessary if data is (XY..XY).
// ##XYDATA with the value (X++(Y..Y)) or (XY..XY)
//----------------------------------------------------------------------------------------
void Tjcampdx::ParseJCampDXHeader(const std::vector <std::string> &Contents, unsigned int DataStart)
{
  unsigned int MandFieldsRead = 0,        // Mandatory fields read
               npt;                       // Number of points
  std::string Key, Line;                       // The key and its value

  jcComment.clear();                      // Clears the object comment std::string

  // Goes through all the lines in the header
  for (std::vector <std::string>::const_iterator it = Contents.begin();
                                           it != Contents.begin() + DataStart; it++)
  {
    Line = *it;         // Copies the line from contents to a local variable
    Key = GetKey(Line); // Gets the key from Line and replaces Line by the key value
    
    if ((Key == "##JCAMP-DX") && (Line[0] == '4') && !(MandFieldsRead & mf_Version))  // Reads JCamp-DX version
      MandFieldsRead |= mf_Version;
    else if ((Key == "##FIRSTX") && !(MandFieldsRead & mf_FirstX))    // Reads First X value
    {
      MandFieldsRead |= mf_FirstX;  
      jcx0 = strtod(Line.c_str(),NULL);
    }
    else if ((Key == "##LASTX") && !(MandFieldsRead & mf_LastX))      // Reads Last X value
    {
      MandFieldsRead |= mf_LastX;  
      jcxf = strtod(Line.c_str(),NULL);
    }
    else if ((Key == "##DELTAX") && !(MandFieldsRead & mf_DeltaX))    // Reads Delta X value
    {
      MandFieldsRead |= mf_DeltaX;  
      jcdx = strtod(Line.c_str(),NULL);
    }
    else if ((Key == "##NPOINTS") && !(MandFieldsRead & mf_NPoints))  // Reads the number of points    
    {
      MandFieldsRead |= mf_NPoints;  
      npt = strtod(Line.c_str(),NULL);
    }
    else if ((Key == "##XYDATA") && !(MandFieldsRead & mf_XYData))    // Reads data type
    {
      if (Line == "(X++(Y..Y))")
      {
        MandFieldsRead |= mf_XYData;
        jcConstSpacing = true;
      }
      else if (Line == "(XY..XY)")
      {
        MandFieldsRead |= mf_XYData;
        jcConstSpacing = false;
      }
    }
    // All fields below are optional
    else if ((Key == "##TITLE") || (Key == "##SAMPLING PROCEDURE"))   // Title and sample procedure go into the comment variable
    {
      if (!jcComment.empty())               // Comment was already set...
        jcComment += " | " + Line;          // ...appends separator to the comment.

      jcComment += Line;                    // Appends the Key value to the comment
    }
    else if (Key == "##RESOLUTION")         // Resolution goes into the comment variable with extra "Res = " label
    {
      if (!jcComment.empty())
        jcComment += " | Res = ";
      jcComment += Line;
    }
    else if (Key == "##XUNITS") Units(Line);                        // X axis title
    else if (Key == "##YUNITS") jcYTitle = Line;                    // Y axis title
    else if (Key == "##YFACTOR") jcYFactor = strtod(Line.c_str(),NULL);   // Factor by which Y values should be multiplied
  }

  // Checks for the existence of mandatory JCamp fields and throws exception if needed
  if (!(MandFieldsRead & mf_Version) || !(MandFieldsRead & mf_XYData))
    throw Ejcampdx(Ejcampdx::NO_VERSION_OR_TYPE);

  // If data is in (XY..XY) format, then no info on x0, dx, etc... are needed
  if (!jcConstSpacing)
      return;

  // The data has a constant spacing. Checks for required variables.
  if ((MandFieldsRead & mf_FirstX) &&           // x0, dx and np found
      (MandFieldsRead & mf_DeltaX) &&
      (MandFieldsRead & mf_NPoints))
  {
    if (MandFieldsRead & mf_LastX)              // xf also present -- checks
    {
      double tmp = jcxf;                        // Backs up last x read from file
      jcxf = jcx0 + (npt - 1) * jcdx;           // Calculates last x from x0, dx and npt
      if (fabs(jcxf - tmp) > 0.01)              // Tolerates just a small difference
        throw Ejcampdx(Ejcampdx::INCOMPATIBLE_XF);
    }
    return;
  }
  else if ((MandFieldsRead & mf_FirstX) &&      // x0, dx and xf found
           (MandFieldsRead & mf_DeltaX) &&
           (MandFieldsRead & mf_LastX))
  {
    npt = 1 + int(0.1 + fabs(jcxf - jcx0) / jcdx);  // Calculates the number of points
    return;
  }
  else if ((MandFieldsRead & mf_DeltaX) &&      // dx, npt and xf found
           (MandFieldsRead & mf_NPoints) &&
           (MandFieldsRead & mf_LastX))
  {
    jcx0 = jcxf - (npt - 1) * jcx0;             // Calculates x0
    return;
  }
  else
    throw Ejcampdx(Ejcampdx::HEADER_FAILURE);   // Generic error
}
