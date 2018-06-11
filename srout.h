#ifndef CONVODIFF_H
#define CONVODIFF_H

#include "sopnamsp.h"
#include "machdefs.h"
#include <math.h>
#include <iostream>
#include <typeinfo>
#include <sstream>
#include <fstream>

#include "tvector.h"
#include "srandgen.h"
#include "fioarr.h"
#include "sopemtx.h"
#include "pexceptions.h"
#include "matharr.h"
#include "fabtwriter.h"
#include <hisprof.h>
#include <vector>
#include "histats.h"   // includes more/less all HiStats module include files
#include "fftwserver.h"

#include "myrg.h"

// Return the spectral poser for a given wave number wk 
double SpecShape(double wk, double Lx, double Ly, double rdiff, double alpha, long double falpha);

//----------- Function GenerateSpectra ---------
// Generate a random fourier spectra with power spectrum 
// according to SpecShape function
//void GenerateSpectra(TMatrix< complex<r_8> > fftarr, double Lx, double Ly, double rdiff);
void GenerateSpectra(TMatrix< complex<r_4> > fftarr, double Lx, double Ly, double rdiff, 
		     double alpha, long double falpha);
//void GenerateSpectra(TMatrix <r_8>  fftarr, double Lx, double Ly, double rdiff);

//----------- Function ComputePowerSpectra ---------
// Compute power spectrum as a function of wave number k 
// Input : FFT Array (Matrix)
// Output : power spectrum (profile histogram)
// if fgnfy==true, consider negtive frequencies along Y pour y > ny/2
//HProf ComputePowerSpectra(TMatrix< complex<r_8> >& fftarr, double Lx, double Ly, int nbin=0, bool fgnfy=true);
HProf ComputePowerSpectra(TMatrix< complex<r_4> >& fftarr, double Lx, double Ly, int nbin=0, bool fgnfy=true);

//----------- Function GiveLambda ------------
double GiveLambda0(string ppfname, char objname[]);

//----------- Function  ComputeScreen ---------
// Compute a screen = exp[ i ( A*Phi + Spherical Wave Term ) ]
//void ComputeScreen(TMatrix< complex<r_8> >& scr, string ppfname, char objname[], int mfac, double A, bool smooth);
void ComputeScreen(TMatrix< complex<r_4> >& scr, string ppfname, char objname[], int mfac, double A, 
                   bool smooth, bool vac=false);

//----------- Function  ComputeDiffractionPattern ---------
// Compute the diffraction pattern corresponding to a screen
// TMatrix<r_8> ComputeDiffractionPattern(TMatrix< complex<r_8> > & scr);
TMatrix<r_4> ComputeDiffractionPattern(TMatrix< complex<r_4> > & scr);

//-------------------- Limb Profile (ordinary space convoluton---------------------
// To be used by convolution subroutin
// This is a 2D matrix of the star limb profile which is projected in image plane.
// There is considerd two kinds of such a profile; one has a sharp edge, it is one
// inside the circle of projected star and zer outside. The other has a gossian
// distribution in both x and y directions with sigma=projected star radius. 
// the argument, rad, is the projected star radius in image pixel unit.
  
TMatrix<r_8> profile( int rad);

//--------------------- convolution ------------------------
// This subroutin considers the source finite size effect on the diffraction pattern by
// calculating the convolution product between he limb profile of the star and 
// the diffraction pattern of the point like source. The convolution integral is computed
// inside the region ofprojected limb profile (in the image plane), this region is
// effectivly considerd as a circle with radius of the projrcted star radius. So, a
// margin with the length of projected star radius is set around the image plane to avoid
// the mal effect of boardes.
// The argument 'pat' is the difraction pattern of a point like source and 'rsp' is the
// projected star radius in image pixel unit.  

TMatrix<r_8> convolution(TMatrix<r_8> pat, double rsp);


//---------------------- Limb Profile (for fourier Convolution) -----------------------
// 

// TMatrix< complex<r_8> > fftprofile( TMatrix <r_8> diffrac, double rad);  
TMatrix< complex<r_4> > fftprofile( TMatrix <r_4> diffrac, double rad);  

//---------------------- convolution (by fourier space) ------------------------

// TMatrix<r_8> fconvolution( TMatrix< complex<r_8> > & fftprof, double rad);
TMatrix<r_4> fconvolution( TMatrix< complex<r_4> > & fftprof, double rad);

//---------------------- swap --------------------
// To be used by MakeLightCurve
// This subroutin swap the initial and final coordinates valus

inline void swap(int &x1, int &y1, int &x0, int &y0)
{
  int a,b;

  a=y1; y1=y0; y0=a;
  b=x1; x1=x0; x0=b;
  // cout << "swaping"<<endl;
  // cout << "(xi , yi) = ( " << x0 << "," << y0 << " )" << endl;
  // cout << "(xf , yf) = ( " << x1 << "," << y1 << " )" << endl;
  // cout << "end of swaping"<<endl;

}


//--------------------- Slop ---------------------
// To be used by MakeLightCurve subroutin
// Calculation of the slope of a line between to given points

inline double slope(int x1, int y1, int x0, int y0)
{
  double s=(x1-x0)/(y1-y0);
  return(s);
}  


//-------------------- lcgen ---------------------
// To be used by MakeLightCurve subroutin
// This subroutin generates the light curve through a line of slope m.
// It makes the interpolation of intensity across the line inside the diffraction pattern

TVector<r_8> lcgen(int p1, int q1, int p0, int q0, double m, TMatrix<r_8> & convo);

//-------------------- MakeLightCurve -----------------
// This subroutin calculates the light curve through a given line which is characterized
// by initial and final points (x0,y0),(x1,y1) respectively

// DataTable MakeLightCurve(int x1, int y1, int x0, int y0, TMatrix<r_8> & conv, double del0k);
DataTable MakeLightCurve(int x1, int y1, int x0, int y0, TMatrix<r_4> & conv, double del0k);


//------------------------------ ContraRed --------------------------------

// This subroutin calculates the contrast( (Imax-Imin)/(Imax+Imin) ) as
// a function of the redused star radius. Initial and final points of light curve
// are (x0,y0),(x1,y1) respectively. pat is the intensity diffraction
// pattern. rfp is the fresnel radius in screen pixelsize unit times by
// squar root of pi.del is pixelsize in image plane in kilometer.

DataTable contrared(int x1, int y1, int x0, int y0, TMatrix<r_8> &pat, double rfp, double del0k);


// ----------------------------- Contrvstime ----------------------------
// to be used by contratime subroutin. dtlc is the light curve data table which has been computed
// by contratime. s is the numbers of members
// in each chosen block of the light curve 

TVector<r_8> contrvstime(DataTable dtlc, int s);

// ----------------------------- Contratime ------------------------------
// contrasr vs time scale. This subroutin calculates the light curve through the line between two
// given points of the image. rfp is the Fresnel radius in screen pixelsize unit times by squar root 
// of pi, del0k is the image pixel length in kilometer. rst, vel and tinit are reduced star radius,
// relative image velocity and initial(minimum) time scale input respectfully.
 
DataTable contratime(int x1, int y1, int x0, int y0, TMatrix<r_8> &convo, double rsp, double del0k, double vel, double tinit);

// ------------------------ DiffusionRadius ------------------------
// This subroutin calculates the diagram of screen phase (phi) root mean squar
// vs distance inside the screen.

// DataTable DiffusionRadius(TMatrix<r_8> &phi, double del1, double rdiff);
DataTable DiffusionRadius(TMatrix<r_4> &phi, double del1, double rdiff, double &reRdiff, double gambeta=1.);


// --------------------------- calrms -----------------------------
// To be used by DiffusionRadius subroutin. This subroutin calculates
// the root mean squar of phi: sqrt(<(phi(x,y)-phi(x',y'))^2>)

// double calrms(TMatrix<r_8> &phi, int dist);
double calrms(TMatrix<r_4> &phi, int dist, double gambeta=1.);

// ----------------------------- compouteStrucFunc -------------------------
// Returns the structure Function value at sepration Rdiff
// double ComputeStrucFunc(TMatrix<r_8> &phi, double del1, double rdiff);
double ComputeStrucFunc(TMatrix<r_4> &phi, double del1, double rdiff);


// ------------------------------- Contramindex2d ------------------------------
// calculation of contrast and standard deviation of the whole image.
double contramindex2d(TMatrix<r_4> pat);


// ----------------------------------- contradel1 -------------------------------
// This subrutin puts del1 tvectorand the related point like source contrast tvector 
// in a datatable
DataTable contradel1(TVector<r_8> d1, TVector<r_8> cont);


// --------------------------------- ComputeReSampling ----------------------------------
// This subroutin resample the ohase screen by dividing the N by two and fixing the 
// physical length of the screen, So, the screeen pixel size (del1)will be doubled.

// void ComputeReSampling(TMatrix<r_8> &phi, double &del1);
void ComputeReSampling(TMatrix<r_4> &phi, double &del1);


// -------------------------------- Compute Resized Matrix ------------------------
// Thid subroutin reproduces a new diffraction image from a larger one. The new image
// will have a smaller physical size by homothety coefficient of alpha. To reproduce a
// smaller image pixel values, the interpolation of the four nearest neighbor pixels of
// the initial image (diff) to the new pixel are considerd.

TMatrix <r_8> ComputeResizedMatrix(TMatrix <r_8> diff, double homoth);

// -------------------------------- Compute Lens Phi ----------------------------------------
// This subroutin is a test of computing the wave front after paasing through a simple lens.
// The lens is a biconvexe with the same curvature in both parts.
// f: focal length, D0: lens thikness. The lengths are in meter

void ComputeLensPhi(TMatrix <r_8> &phi, int N, double f, double D0, double lambda, double del1);
  
// ---------------- Modulation Index For a Light Curve -----------------------
double ModulIndex(const DataTable &lc);
  
// -------------------- Fourir Transform of the Light Curve -----------------------------
DataTable ComputeFFTlc (const DataTable &lc, double length);

// ---------------- ReSizePatt -----------------
TMatrix<r_8> ReSizePatt(const TMatrix<r_8> &currpatt, double A);

// ---------------------- ComputeCorrSpec ------------------
DataTable ComputeCorrSpec(double DTmin, int nbin,  double Tf, int Nf, 
			  const TVector<r_4> &fcurve, const TVector<r_4> &ftime, double sig);

// ---------------------- ComputeSigma2D ------------------
double ComputeSigma2D(TMatrix<r_4> phi);


// ---------------- Compute1dSpec -------------------
Vector Compute1dSpec(TMatrix< complex<r_4> > fftphi, double del, double xmim, 
		     double xmax, int nbin, double step);

// -------------- Finction intpol ------------------
double intpol(double c1, double c2, int n, double d, double rdiff);


// -------------------------------- Compute Resampled Matrix ------------------------
// This subroutin re-grid the images with \lambda > \lambda_0.
TMatrix <r_4> ComputeResampledMatrix(TMatrix <r_4> M0, TMatrix <r_4> M);

// --------------- ComputeNph ----------- 
double ComputeNph(double m, double Texp, int Nx, int Ny, double lambda, double qeff, double del0x, double del0y, double D);

// --------------------- ShiftPattern -----------------
TMatrix<r_4> ShiftPattern(TMatrix<r_4> diffrac, int dx, int dy, int Nx, int Ny, int d);
 
// ------------------ FillFitsHeader --------------------
void FillFitsHeader(FitsImg2DWriter &fitsimag, const char * headinfo);

// ------------------ ScalePattern --------------------
HProf ScalePattern(TMatrix<r_4> &diff, double &rmax, double threshold=1.e-3);

// ------------------ GiveDiffMax --------------------
double GiveDiffMax(TMatrix<r_4> diff, int &ii, int &jj);

// --------------- ComputeNph2 ----------- 
double ComputeNph2(double mag, double Texp, int Nx, int Ny, double lambda, double qeff, double del0x, double del0y, char *filter); 

#endif



