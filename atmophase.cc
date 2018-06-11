 
/*
      This program calculates the phase delay of the wave front 
              which passes through a turbulent cloud.            
*/

#include "sopnamsp.h"
#include "machdefs.h"
#include <math.h>
#include <iostream>
#include <typeinfo>

#include "tvector.h"
#include "srandgen.h"
#include "fioarr.h"
#include "sopemtx.h"
#include "pexceptions.h"
#include "matharr.h"

#include "resusage.h"

#include "histats.h"
#include "fftwserver.h"

#include "srout.h"

#include "sambainit.h"
#include <datatable.h>

// #include "tarrinit.h"

#include "timing.h"

#include "myrg.h"

int main(int narg, char* arg[])
{


SophyaInit();
InitTim();   // Initializing the CPU timer
Auto_Ini_Ranf(2);
int rc = 0;


// We define here the parameters which are decoded from the command line
 double Rdiff = 11.; //4.5;  // diffusion radius in centimeter 

 double alpha = 11./3;

 string ppfname = "atmophase.ppf";

if ((narg > 1) && (strcmp(arg[1],"-h")==0)) {
  cout << " phase/Usage: phase [OutPPFName [ RDiff ] ] \n" 
       << "  Rdiff in centimeters " << endl;
  return 0;
}

 if (narg > 1)  ppfname = arg[1];     // out put ppf name
 if (narg > 2)  Rdiff = atof(arg[2]); // add to float
 if (narg > 3)  alpha = atof(arg[3]); 

cout << " atmophase/Info : PPFName= " << ppfname << " Rdiff= " << Rdiff << endl;

// check resource usage for a given part of the program
  ResourceUsage res;
  res.Update();
  cout << " DEBUGGING:phase.cc: start of maine  " << endl;
  //  cout << " Memory size increase (KB):" << res.getDeltaMemorySize() << endl;
  cout << " Resource usage info : \n" << res << endl;

// Define the size of real arrays (N x N)

  int Nx = 4000;   
  int Ny = 4000; 

try {


/////////////////////////////////////////////////////
///////////feeding the physical parameters///////////
///////////////////// MKS System ////////////////////
/////////////////////////////////////////////////////
  double pi = 3.141592635;    // pi number  
  double nm2cm = 1.e-7;       // nanometer to centimeter
  double cm2micron = 1.e4;
  double lambda = 675*nm2cm;  // wave length
  double D = 154;             // Telescope diameter in cm
  double fr = 8.6;            // Focal ratio
  double f = D*fr;            // Telescoope focal length in cm 
  double Rf = sqrt(lambda*f/2/pi); //3.; // Fresnel radius in cm

  
  double del1 = .1; //.05; //sqrt(2*pi*Rf*Rf/Nx);  // phase pixelsize in cm  

  bool dosample = false;    // if true performs resampling

  double Lx = Nx * del1; // wave-front physical x-length in centimeter.
  double Ly = Ny * del1;
  double del0 = 2*pi*Rf*Rf/Lx; 
  
  POutPersist po(ppfname);            
  cout << endl;
  cout << " atmophase/Info: opening PPF file " << ppfname  << endl;
  cout << endl;
  
  cout << endl;
  cout << "-------------------------------------------------------------------" << endl;
  cout << " Number of pixels (Ny*Nx)= " << Ny << " * " << Nx << endl; 
  cout << " Lambda = " << lambda/nm2cm <<"(nm)"<< endl; 
  cout << " Focal length = " << f/100 << " m " << endl;
  cout << " Rdiff = " << Rdiff << " cm" << "  = " << Rdiff/del1 <<" pixels" << endl;
  cout << " Rf = " << Rf << " cm" << "  = " << Rf/del1 <<" pixels" << endl;
  cout << " pixel size = " << del1 <<" cm " << endl; 
  cout << " Image pixel size = " << del0*cm2micron << " micron " << endl;
  cout << " Wave-front Physical size (Ly*Lx) = " << Ly << " cm * " << Lx << " cm" << endl;
  cout << "-------------------------------------------------------------------" << endl;
  cout << endl;
  
  
  // Compute fourier transform and spectra
  FFTWServer ffts;                     
  ffts.setNormalize(); // false by default = No division by N*N through FFTForward             
  cout << " FFTServer info string = " << ffts.getInfo() << endl;
  cout << " FFTServer Normalization = " << ffts.getNormalize() << endl;
  cout << endl;
  
  
  // Create a complex matrix to generate a new image, using SpecShape function 
  // We use the fftimg sizes in order to recover NxN real image after FFTBackward
  cout << "1.Generating phase screen Fourier spectra GenerateSpectra() ... " << endl;
  
  //  TMatrix< complex<r_8> > fftphi(N, int(N/2+1));

  //  TMatrix< complex<r_8> > fftphi(Ny, int(Nx/2+1));
  TMatrix< complex<r_4> > fftphi(Ny, int(Nx/2+1));

  long double falpha = alpha * tgamma(-alpha/2) / tgamma(alpha/2) / pow(2,alpha);

//  long double falpha = pow(2.,alpha-2) * tgamma((alpha+2)/2.)*tgamma((alpha+2)/2.)
//                       * tgamma((alpha+4)/2.) * tgamma(alpha/2) * sin(pi*(alpha-2)/2.)
//                       / pow(pi,alpha) * tgamma(alpha+1);  
  //Abeta in NBD spectrum
    
  long double xx = 1./(2*pow(2*pi,alpha-1)* falpha); // compare us with nbd
  cout << " Turbulence exponent: " << alpha << endl;
  cout.setf(ios::fixed,ios::floatfield);
  cout.precision(8);
  cout << " falpha = " << falpha << endl;
  cout << " xx = " << xx << endl;
  cout << endl;  
  GenerateSpectra(fftphi, Lx, Ly, Rdiff, alpha, falpha); 


  //  GenerateSpectra2(fftphi, Lx/1000., Ly/1000., Rdiff/1000., alpha, falpha, del1/1000.); 
  // all lengths are input in kilometers 
    
/*
  cout << " Computing power spectra using ComputePowerSpectra(fftphi) ... " << endl;
  HProf cpsphi =  ComputePowerSpectra(fftphi, Lx/1000., Ly/1000., Nx/10);
  po << PPFNameTag("psphi") << cpsphi;
*/
//  po << PPFNameTag("fftphi") << fftphi;

  TMatrix< r_4 > phi(Ny, Nx);   
  //  TMatrix< r_4 > phi(fftphi.NRows(), fftphi.NRows());   


  // FFTBackward uses the input array as a workspace
  
  cout << "2.Computing phase Phi / calling ffts.FFTBackward(fftphi, phi)... " << endl; 
  ffts.FFTBackward(fftphi, phi,true); 

  phi /= (Lx*Ly); 

  // FFT calculates the serie of furier, by dividing the results by (1/(N*del1))^2
  // (for fourier backward transform) we obtain the approximated fourier integral 
  // with the correct dimension in result. (phase is dimensionless, radian) 
  // po << PPFNameTag("phi8kcut") << phi; cout << " CUT-OFF CHECK: objname: phi8kcut " << endl;


  
  /*
    // Computing phase screen of a simple convex lens
    
    TMatrix <r_8> phi;
    double au2meter = 149598.e6; // Astronomical unit to meter
    double f = z0;
    double D0 = 0.01;// 0.0001 * au2meter;
    cout << " Computing phase screen of a simple biconvex lens " << endl;
    cout << " f = " << f << "(m) = " << f / del1 << " screen pixelsize " << endl;
    // cout << " D0 = " << D0 << " (m) = " << D0 / del1 << " screen pixelsize " << endl; 
    ComputeLensPhi(phi, N, f, D0, lambda, del1);
    phi.Info()["L"] = Lx;
    po << PPFNameTag("phi8k") << phi;
    
  */
  
  
  cout << "3.Computing the structure function:  DataTable DiffusionRadius(phi,del1,Rdiff)" 
       << endl;
  double reRdiff = 0.; // recalculated Rdiff by computing the phase structure function
  DataTable diffusionrad;
  double gambeta = pow(2.,alpha-1) * tgamma((alpha+2)/2.)*tgamma((alpha+2)/2.)
                   * tgamma((alpha+4)/2.) / tgamma(alpha/2.) / tgamma(alpha+1.);
  gambeta = 1.;
  cout << " gammabeta=" << gambeta << endl;
  diffusionrad = DiffusionRadius(phi,del1,Rdiff,reRdiff,gambeta);

  if (fabs(Rdiff-reRdiff) < del1) reRdiff =  Rdiff;

  cout << " Best diffusion radius from computing the phase structure function: " 
       << reRdiff << " cm" << endl;
  po << PPFNameTag("rdiff") << diffusionrad;

  cout << "4.Rescaling the phase ... " << endl;

  double beta = reRdiff/Rdiff;
  beta =  pow(beta,(alpha-2.)/2.);
  cout << " phase scaling factor = " << beta << endl;
  phi *= beta; 
  reRdiff = 0.;
  diffusionrad.Clear();
  diffusionrad = DiffusionRadius(phi,del1,Rdiff,reRdiff,gambeta);
  cout << endl;
  cout << " Best diffusion radius after rescaling:" << reRdiff << " cm" << endl;

  cout << " Initial Fried radius = " << Rdiff*pow(6.88,.6) << " cm " << endl;  
  cout << " Fried radius = " << reRdiff*pow(6.88,.6) << " cm " << endl;  
  cout << " Deviation = " << fabs(Rdiff*pow(6.88,.6)-reRdiff*pow(6.88,.6))/(Rdiff*pow(6.88,.6))
       << endl;
//  cout << " Fried radius = " << reRdiff/.982 << " cm " << endl;  
  cout << endl;

//  cout << " Saving the phase structure function in data table 'Rdiff' ... " << endl;
//  cout << endl;
  po << PPFNameTag("rerdiff") << diffusionrad;


  phi.Info()["Ly"] = Ly;
  phi.Info()["Lx"] = Lx;
  phi.Info()["Rdiff"] = Rdiff;
  phi.Info()["Rf"] = Rf;
  phi.Info()["reRdiff"] = reRdiff;
  phi.Info()["lambda"] = lambda;
  phi.Info()["D"] = D;
  phi.Info()["fr"] = fr;
  //  po << PPFNameTag("phi8k") << phi; 
  po << PPFNameTag("phi") << phi; 



//////////// /////// Phi Resampling ///////////////////// 
  if(dosample){
     cout << "4.Phi Resamlpling for pixel numbers of N/2 and N/4. Physical size is fixed " 
	 << endl;
    ComputeReSampling(phi,del1);
    cout << " new del1 = " << del1/1000 << " (km) " << "    new N = " << phi.NCols() << endl;
    phi.Info()["Ly"] = Ly;
    phi.Info()["Lx"] = Lx;
    po << PPFNameTag("phi4k") << phi;
    cout << " Rdiff = " << Rdiff/del1 << "(screen pixelsize)" << endl;
    diffusionrad.Clear();
    cout << " Computing the structure function ... " << endl;
    cout << endl;
    diffusionrad = DiffusionRadius(phi,del1,Rdiff,reRdiff);
    cout << endl;
    {
      //      TMatrix< complex<r_8> > fftphi4k;
      //      TMatrix<r_8> cphi = phi;

      TMatrix< complex<r_4> > fftphi4k;
      TMatrix<r_4> cphi = phi;

      ffts.FFTForward(cphi,fftphi4k);
      cout << " Computing the phase power spectrum ... " << endl;
      cout << endl;
      HProf cpsphi =  ComputePowerSpectra(fftphi4k, Lx/1000., Ly/1000., Nx/20);
      po << PPFNameTag("psphi4k") << cpsphi;
    }
    
    ComputeReSampling(phi,del1);
    cout << " new del1 = " << del1/1000 << " (km) " << "    new N = " 
	 << phi.NCols() << endl;
    phi.Info()["Ly"] = Ly;
    phi.Info()["Lx"] = Lx;
    po << PPFNameTag("phi2k") << phi;
    cout << " Rdiff = " << Rdiff/del1 << "(screen pixelsize)" << endl;
    diffusionrad.Clear();
    cout << " Computing the structure function ... " << endl;
    cout << endl;
    diffusionrad = DiffusionRadius(phi,del1,Rdiff,reRdiff);
    cout << endl;
    {
      //      TMatrix< complex<r_8> > fftphi2k;
      //      TMatrix<r_8> cphi = phi;

      TMatrix< complex<r_4> > fftphi2k;
      TMatrix<r_4> cphi = phi;

      ffts.FFTForward(cphi,fftphi2k);
      cout << " Computing the phase power spectrum ... " << endl;
      cout << endl;
      HProf cpsphi =  ComputePowerSpectra(fftphi2k, Lx/1000., Ly/1000., Nx/40);
      po << PPFNameTag("psphi2k") << cpsphi;
    }
  } // end if(dosampling)
  
}// end of try
 
 
catch (PThrowable& exc) {
  cerr << " phase.cc catched Exception " << exc.Msg() << endl;
  rc = 77;
}  
catch (std::exception& sex) {
  cerr << "\n treccyl.cc std::exception :" 
       << (string)typeid(sex).name() << "\n msg= " 
       << sex.what() << endl;
}
catch (...) {
  cerr << " phase.cc catched unknown (...) exception  " << endl; 
  rc = 78; 
} 

cout << " ============ End of phase.cc ============== " << endl;
 cout << endl;

return rc;
}


