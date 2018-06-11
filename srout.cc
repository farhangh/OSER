#include "srout.h"
#include "resusage.h"
#include "myrg.h"


double pi = 3.141592635;     // pi number  

//----------- Function SpecShape ---------
// Return the spectral poser for a given wave number wk
// Lengths are in kilometers 
double SpecShape(double wk, double Lx, double Ly, double rdiff, double alpha, long double falpha) 
{

  //  double pi = 3.141592635;     // pi number  

  // rdiff: diffusion radius in kilometer

  // protectection for zero/negative wk
  if (wk < 1.e-10) return 0.;
  //  if (wk < 1.e-5) return 0.; cout << " Info/convodif.cc: SpecShape: CUT-OFF " 
  //  << endl;
  // spectra of type (w/scale)^alpha

 
  long double xx = 1./(2*pow(2*pi,alpha-1)* falpha); // us
//  long double xx = 1./(2*pow(2*pi,alpha-2)*falpha); // Lovelace
  double spectre = Lx * Ly * rdiff * rdiff * pow(rdiff*wk,-alpha) * xx;

//  double spectre = Lx*Ly * rdiff * rdiff * pow(rdiff*wk,-alpha) * falpha;

//  double spectre = .023 * Lx*Ly * rdiff * rdiff * pow(rdiff*wk,-alpha); //Wiener
  return (spectre);
}

//----------- Function GenerateSpectra -------
// Generate a random fourier spectra with power spectrum 
// according to SpecShape function
//void GenerateSpectra(TMatrix< complex<r_8> > fftarr, double Lx, double Ly, double rdiff)
void GenerateSpectra(TMatrix< complex<r_4> > fftarr, double Lx, double Ly, 
		     double rdiff, double alpha, long double falpha)
{
  // fftarr is supposed to represent 2-D fourier transform of 
  // a real input array. 
  // The second half of the array along Y (matrix rows) contain
  // negative frequencies

  cout << " convodif.cc:GenarateSpectra: alpha = " << alpha << endl;
  cout << endl;

  for(int ky=0; ky<fftarr.NRows(); ky++) {

    double kyy = ky;
    if (ky > fftarr.NRows()/2) kyy = fftarr.NRows()-ky;
    for(int kx=0; kx<fftarr.NCols(); kx++) {

      double wk = sqrt((double)(kx/Lx*kx/Lx+kyy/Ly*kyy/Ly));
      //      double wk = sqrt((double)(kx*kx+kyy*kyy))
      double amp = SpecShape(wk, Lx, Ly, rdiff, alpha, falpha);

      complex<r_8> za = complex<r_8> (GauRnd(0, sqrt(amp)), GauRnd(0., sqrt(amp))); 
      // <za*za> = 2*amp

      //      complex<r_4> za = complex<r_4> (GauRnd(0, sqrt(amp)), GauRnd(0., sqrt(amp)));
      //double za = GauRnd(0, sqrt(amp));

      fftarr(ky, kx) = (complex<r_4>)(za/sqrt(2.)); // /Lx/Ly ; 
    }
  }
}

//----------- Function ComputePowerSpectra ---------
// Compute power spectrum as a function of wave number k 
// Input : FFT Array (Matrix)
// Output : power spectrum (profile histogram)
// if fgnfy==true, consider negtive frequencies along Y pour y > ny/2

//HProf ComputePowerSpectra(TMatrix< complex<r_8> >& fftarr, double Lx, double Ly, int nbin, bool fgnfy)
HProf ComputePowerSpectra(TMatrix< complex<r_4> >& fftarr, double Lx, double Ly, int nbin, bool fgnfy)
{
  // fftarr is supposed to represent 2-D fourier transform of 
  // a real input array. 
  // The second half of the array along Y (matrix rows) contain
  // negative frequencies
   // Lx & Ly: screen sizes in kilometers

  int nbh = int(sqrt(2.)*fftarr.NCols()) ;

  // The profile histogram will contain the mean value of FFT amplitude
  // as a function of wave-number k = sqrt(kx*kx+ky*ky)
  if (nbin < 1) nbin = nbh;
  // HProf hp(-0.5, nbh-0.5, nbin);
  HProf hp(-7., -1., nbin,-1000,1000); // logarithmic scales

  // cout << " convodif.cc: ComoutePowerSpectras: fftarr = " << fftarr;

  for(int ky=0; ky<fftarr.NRows(); ky++) {
    double kyy = ky;
    if (fgnfy && (ky > fftarr.NRows()/2) ) kyy = fftarr.NRows()-ky;
    for(int kx=0; kx<fftarr.NCols(); kx++) {
      //      double wk = sqrt((double)(kx/Lx*kx/Lx+kyy/Ly*kyy/Ly));
      float wk = sqrt((float)(kx/Lx*kx/Lx+kyy/Ly*kyy/Ly));

      complex<r_8> za = (complex<r_8>)(fftarr(ky, kx)) * Lx * Ly; 
      // the FT of phi is dimensional because the power spectrum has dimension of length^4

      //      double amp = sqrt(za.real()*za.real()+za.imag()*za.imag());

      //     double amp = za.real()*za.real()+za.imag()*za.imag();
      float amp = (float)(za.real()*za.real()+za.imag()*za.imag());

      if((wk>=0.00001)&&(amp>0.00001))  hp.Add(log10(wk), log10(amp));
      //      hp.Add(wk, amp);
    }

  }
  return hp;
}


//----------- Function GiveLambda ------------
double GiveLambda0(string ppfname, char objname[])
{
  double kpc2m = 3.08e19;      // kilo parsec to meter
  double pi = 3.141592635;     // pi number  

//  TMatrix<r_8> phi;
  TMatrix<r_4> phi;

  PInPersist pin(ppfname);
  pin >> PPFNameTag(objname) >> phi; 

  double lambda0 = phi.Info()["lambda"];
  double Rdiff = phi.Info()["Rdiff"];
  double reRdiff = phi.Info()["reRdiff"];
  double Rf = phi.Info()["Rf"];
  double Rref = 2*pi*Rf*Rf / Rdiff ;
  double z1 = phi.Info()["z1"];
  double z0 = phi.Info()["z0"];

  cout << endl;
  cout << " info/convodif.cc:GiveLambda0: " << endl;
  cout << " lambda0 = " << lambda0 * 1.e9 << " nm " << endl;
  cout << " Rdiff = " << Rdiff / 1000. << " km " << endl;
  cout << " reRdiff = " << reRdiff / 1000. << " km " << endl;
  cout << " Rf = " << Rf / 1000. << " km " << endl;
  cout << " Rref = " << Rref / 1000. << " km " << endl;
  cout << " source distance = " << (z1+z0) / kpc2m << " kpc " << endl;
  cout << " cloud distance = " <<  z0 / kpc2m << " kpc " << endl;
  cout << endl;

  cout << " DEB:GiveLambda0 " << endl;
  
  return(lambda0);
}


//----------- Function  ComputeScreen ---------
// Compute a screen = exp[ i ( A*Phi + Spherical Wave Term ) ]
//void ComputeScreen(TMatrix< complex<r_8> > &scr, string ppfname, char objname[], 
//		   int mfac, double A, bool smooth)

void ComputeScreen(TMatrix< complex<r_4> > &scr, string ppfname, char objname[], 
		   int mfac, double A, bool smooth, bool vac)
{

  //  TMatrix<r_8> phi;
  TMatrix<r_4> phi;

  cout << " Info/srout.cc: ComputScreen: Opening ppf file: " << endl; 
  cout << "                    " << ppfname << ", ppf obj.name: " << objname << endl;


  PInPersist pin(ppfname);
  pin >> PPFNameTag(objname) >> phi; 

  scr.ReSize(phi.NRows(),phi.NCols());

  double rf = phi.Info()["Rf"] ;
  double Ly = phi.Info()["Ly"];
  double Lx = phi.Info()["Lx"];
  double D  = phi.Info()["D"];
  double fr = phi.Info()["fr"];
  double del1 = Ly / phi.NRows();
  double rfpix = rf / del1; 

  double vphi = 0.;
  double multij = 1.;
  double rfpix2 = rfpix*rfpix;
  double nny = phi.NRows()/2;  //*sycord;
  double nnx = phi.NCols()/2;  //*sxcord;
  
  int ND = int(D/del1/2);
//  cout << " DEB: srout.cc: ComputeScreen: D = Lx" << endl;

  double pi = 3.141592635;     // pi number


  //  double f; // dark obstacle
  //  cout << " Info/convodif.cc:computescreen: phi(y,x)=0 " << endl;
  //  cout << " Info/convodif.cc:computescreen: The Fiducial zone is multiplied by sin(pi/50*x) " 
  //       << endl;
  //  cout << " DEBUGGING:convodif.cc:computescreen: The screen is just: 1 + 0.1 * sin(pi/50*x1)" 
  //       << endl;
  //  cout << " Info/convodif.cc:computescreen: Light emission is restircted inside a squar hole " 
  //       << endl;
  //  cout << " DEBUGGING:convodif.cc:computescreen: no multipcliation by multij " << endl;


  if(smooth) cout << " Info/srout.cc:computescreen: smoothing included " << endl;
  else cout << " Info/srout.cc:computescreen: no smoothing " << endl;

  for(int y1 = 0; y1 < phi.NRows(); y1++) 
    {
      for(int x1 = 0; x1 < phi.NCols(); x1++) {
	
//	vphi = phi(y1, x1);

//      vphi += ((x1-nnx+.5)*(x1-nnx+.5)+(y1-nny+.5)*(y1-nny+.5))/ 2./rfpix2;


//	vphi *= A ; // different wave-lengths
	multij = pow(-1.,(double)(x1+y1));
	
       if((x1-nnx)*(x1-nnx)+(y1-nny)*(y1-nny) < ND*ND){
         if(vac) vphi = 0.;
         else vphi = phi(y1, x1);

         if((x1==0)&&(y1==0)) cout << " Info/srout.cc:ComputeScreen: plane wave " << endl;
	       scr(y1, x1) = complex<r_8>(cos(vphi), sin(vphi))*complex<r_8>(multij,0.);
//	 scr(y1, x1) = complex<r_8>(1.,0.); 
       }
       else scr(y1, x1) = complex<r_8>(0., 0.);

	
	if(smooth){
	  // Smoothing with each half of the function 0.5*sin(3pi/2-x)
	  double top = mfac*rfpix ;
	  if(y1<top) scr(y1,x1) *= 0.5*(sin(3*pi/2-pi*y1/top)+1);
	  else if(x1<top) scr(y1,x1) *= 0.5*(sin(3*pi/2-pi*x1/top)+1);

	  else if(y1>phi.NRows()-1-top) scr(y1,x1) *= 
	    0.5*(sin(3*pi/2-pi*(1+(y1-phi.NRows()+1+top)/top))+1);

	  else if(x1>phi.NCols()-1-top) scr(y1,x1) *= 
	    0.5*(sin(3*pi/2-pi*(1+(x1-phi.NCols()+1+top)/top))+1);
	  }

      }
    }

    scr.Info()["Lx"] =  Lx;
    scr.Info()["Ly"] =  Ly;
    scr.Info()["D"] =  D;
    scr.Info()["fr"] = fr;
    double Rdiff = phi.Info()["Rdiff"];
    double reRdiff = phi.Info()["reRdiff"];
    scr.Info()["Rdiff"] =  Rdiff / pow(A,1.2);
    scr.Info()["reRdiff"] =  reRdiff / pow(A,1.2);
    scr.Info()["Rf"] = rf / sqrt(A);
    double lambda0 = phi.Info()["lambda"];
    scr.Info()["lambda"] = lambda0 / A ;

}



//----------- Function  ComputeDiffractionPattern ---------
// Compute the diffraction pattern corresponding to a screen
// TMatrix<r_8> ComputeDiffractionPattern(TMatrix< complex<r_8> > & scr)
TMatrix<r_4> ComputeDiffractionPattern(TMatrix< complex<r_4> > & scr)
{

  ResourceUsage res;
/*
  res.Update();
  //  cout << " Memory size increase (KB):" << res.getDeltaMemorySize() << endl;
  cout << " convodif.cc:ComputeDiffractionPattern: " << endl;
  cout << " Resource usage info : \n" << res << endl;
*/

  // FFTFoward uses the input array as a workspace
  //  TMatrix< complex<r_8> > zdp; 
  TMatrix< complex<r_4> > zdp;   
  FFTWServer ffts;      
  ffts.setNormalize(false);               
  ffts.FFTForward(scr, zdp); 

  scr.ZeroSize(); // to gain more memory
  
  
/*
  res.Update();
  cout << " convodif.cc:ComputeDiffractionPattern: Resource usage info after fft: \n" << res << endl;
*/
  // We take the modulus of the computed FFT ...
  // and return the square of its elements
  //  TMatrix<r_8> ramp = module2(zdp);

  TMatrix<r_4> ramp(zdp.NRows(),zdp.NCols()); // = module2(zdp);
  for(int i=0; i<zdp.NRows();i++) 
    for(int j=0; j<zdp.NCols();j++) 
      ramp(i,j) = zdp(i,j).real()*zdp(i,j).real() + zdp(i,j).imag()*zdp(i,j).imag() ; 
  cout << " DEB: ramp(0,0) = " << ramp(0,0) << endl;
  cout << endl;
  //  ramp.Mul(ramp);
  res.Update();
  cout << " srout.cc:ComputeDiffractionPattern: ramp is computed: \n" << res << endl;
  return (ramp); 
}


//-------------------- Limb Profile (ordinary space convoluton) ---------------------
// To be used by convolution subroutin
// This is a 2D matrix of the star limb profile which is projected in image plane.
// There is considerd two kinds of such a profile; one has a sharp edge, it is one
// inside the circle of projected star and zer outside. The other has a gossian
// distribution in both x and y directions with sigma=projected star radius. 
// the argument, rad, is the projected star radius in image pixel unit.
  
TMatrix<r_8> profile( int rad)
{
  int dim=2*rad; // the dimention of the profile matrix
  TMatrix<r_8> pr(dim,dim); // the limb profile output matrix
  
  // sharp edge distribution calculation
  for(int j=0; j<dim; j++)
    for(int i=0; i<dim; i++){
      if ((i-rad)*(i-rad)+(j-rad)*(j-rad)<=rad*rad) pr(j,i)=1.;
      else pr(j,i)=0.;
    }

/*
  // Gaussian distribution calculation
  for(int j=0; j<dim; j++)
    for(int i=0; i<dim; i++)
      pr(j,i)=exp((j-rad)*(j-rad)/rad/rad/2)*exp((i-rad)*(i-rad)/rad/rad/2)
*/
  return(pr);
}

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

TMatrix<r_8> convolution(TMatrix<r_8> pat, double rsp)
{
 TMatrix<r_8> intens(pat.NRows(),pat.NCols()); 
 TMatrix<r_8> prof;

 int rad=int(rsp); // making the radius an integer number of pixels
 
 for(int j0=0; j0<pat.NRows(); j0++)
   for(int i0=0; i0<pat.NCols(); i0++){
     // seting the margins to zero
     if((j0<rsp)||(j0>pat.NRows()-rsp)||(i0<rsp)||(i0>pat.NCols()-rsp)) intens(j0,i0)=0.;
     else 
     {
       // convolution integral
       intens(j0,i0)=0.;
       prof = profile(rad); 
       for(int yy=0; yy<2*rad; yy++)
         for(int xx=0; xx<2*rad; xx++) 
           intens(j0,i0)+=prof(yy,xx)*pat(j0-(yy-rad),i0-(xx-rad));
      }// end of else
   } // end of for(i0)

 return(intens);

}



//---------------------- Limb Profile (for fourier Convolution) -----------------------
// To be used by fconvolution subroutin. In this case, the 2D limb profile matrix is at
// the same size of the image. In the case of a profile with sharp edge, the matrix elements
// of region at the center of the matrix within the projected star radis (in image pixelsize unit)
// are valued 1 and others are zeros. The production of fft of this profile and the fft of diffraction pattern 
// is calculated to be fed to fconvolution(). rad: reduced star radius in image pixel unit  

// TMatrix< complex<r_8> > fftprofile( TMatrix <r_8> diffrac, double rad) 
TMatrix< complex<r_4> > fftprofile( TMatrix <r_4> diffrac, double rad) 
{ 
  
  int ydim = diffrac.NRows();
  int xdim = diffrac.NCols();
  //  TMatrix<r_8> pr(ydim,xdim); 
  TMatrix<r_4> pr(ydim,xdim); 

  for(int j=0; j<ydim; j++) 
    for(int i=0; i<xdim; i++){ 
//       if ((i-xdim/2+1)*(i-xdim/2+1)+(j-ydim/2+1)*(j-ydim/2+1)<=rad*rad) pr(j,i)=1.; 
       if ((i-xdim/2)*(i-xdim/2)+(j-ydim/2)*(j-ydim/2)<=rad*rad) pr(j,i)=1.; 
      else pr(j,i)=0.; 
    }

  //  TMatrix < complex <r_8> > fftpr;
  //  TMatrix< complex<r_8> > fftdiffrac; 
  TMatrix < complex <r_4> > fftpr;
  TMatrix< complex<r_4> > fftdiffrac; 

  FFTWServer ffts;
  ffts.setNormalize(false); // false = No division by N*N through FFTForward  (true by def.)  
  ffts.FFTForward(pr, fftpr);
  ffts.FFTForward(diffrac, fftdiffrac);

 for(int j = 0; j < fftpr.NRows(); j++)
   for(int i = 0; i < fftpr.NCols(); i++){
      double multij = pow(-1.,(double)(j+i));
      //      fftpr(j,i) =  fftpr(j,i) * fftdiffrac(j,i) * complex<r_8>(multij,0.);
      fftpr(j,i) =  fftpr(j,i) * fftdiffrac(j,i) * complex<r_4>(multij,0.);
   }

  return(fftpr);   
}

//---------------------- convolution (by fourier space) ------------------------
// This subroutin perform the convolution between the star limb profile and the point like
// source intensity pattern.

// TMatrix<r_8> fconvolution( TMatrix< complex<r_8> > & fftprof, double rad)
TMatrix<r_4> fconvolution( TMatrix< complex<r_4> > & fftprof, double rad)
{

 FFTWServer ffts;  
 
 //  TMatrix<r_8> conv(fftprof.NRows(),2*fftprof.NCols()-2);
 TMatrix<r_4> conv(fftprof.NRows(),2*fftprof.NCols()-2);
 ffts.FFTBackward(fftprof,conv,true);

 conv /= conv.NRows() * conv.NCols(); // fftbackward normalization

 return(conv); 
}


// ------------------------------- Contramindex2d ------------------------------
// calculation of contrast and modlation index of the entire image either point like
// source diffraction pattern or convoluted one. marg is the pixel number to avoid
// the boarders.

// double contramindex2d(TMatrix<r_8> pat, int margx, int margy)
double contramindex2d(TMatrix<r_4> pat)
{
  double del1 = pat.Info()["del1"];
  double del0 = pat.Info()["del0x"];
 double imax = pat(pat.NRows()/2,pat.NCols()/2);
 double imin = pat(pat.NRows()/2,pat.NCols()/2);
 double totdiff = 0.;
 double totdiff2 = 0.;
 double mindex; // modulation index: intensity dispertion divided by average intensity
 double counter = 0.;
 for(int j=0; j<pat.NRows(); j++)
   for(int i=0; i<pat.NCols(); i++){
     if(imin>pat(j,i)) imin = pat(j,i);
     if(imax<pat(j,i)) imax = pat(j,i);
     totdiff += pat(j,i);
     totdiff2 += pat(j,i)*pat(j,i);
     counter++;
   }
 double meandiff = totdiff/counter;
 double meandiff2 = totdiff2/counter;
 mindex = sqrt(meandiff2-meandiff*meandiff);

 cout << "--------------------------------------------" << endl;
 // cout <<" imin = " << imin <<" imax = " << imax << endl; 
 cout <<" Average flux = " << meandiff << " w/micron^2 " <<endl;
 cout <<" Contrast = " << (imax-imin)/(imax+imin)*100 << " % " << endl;
 cout <<" Standard deviation = " << mindex << endl; 
 cout << "--------------------------------------------" << endl;
 //cout << " convodif.cc:contramindex2d:returning the contrast" << endl; return((imax-imin)/(imax+imin));
 cout << " Info/srout.cc:contramindex2d:returning the standard deviations" << endl; 
 return(mindex);
}


// --------------------------- calrms -----------------------------
// To be used by DiffusionRadius subroutin. This subroutin calculates
// the root mean squar of phi(x,y)-phi(x',y')

// double calrms(TMatrix<r_8> &phi, int dist)
double calrms(TMatrix<r_4> &phi, int dist, double gambeta)
{
  double res;
  //  double err;

  double res1 = 0.;
  double res2 = 0.;
  double res3 = 0.;
  double res4 = 0.;
  double res5 = 0.;
  double res6 = 0.;
  double res7 = 0.;
  double res8 = 0.;

  //  double err1 = 0.;
  //  double err2 = 0.;
  //  double err3 = 0.;
  //  double err4 = 0.;


  int counter = 500000;
  for(int i = 0; i < counter; i++){

    double xdrand = drand01(); // a random number in interval [0,1]
    double ydrand = drand01();

    xdrand *= phi.NCols()-2*dist-1.; 
    ydrand *= phi.NRows()-2*dist-1.;

    xdrand += dist+1;  // dist+1 <= xdrand <= N-dist
    ydrand += dist+1;

    int xrand = int(xdrand);
    int yrand = int(ydrand); 

    double f1 = phi(yrand+dist,xrand)-phi(yrand,xrand);
    double f2 = phi(yrand,xrand+dist)-phi(yrand,xrand);
    double f3 = phi(yrand,xrand-dist)-phi(yrand,xrand);
    double f4 = phi(yrand-dist,xrand)-phi(yrand,xrand);

    int dcos = int(dist/sqrt(2));
    if ((dist/sqrt(2)-dcos) > .5) dcos++;

    double f5 = phi(yrand+dcos,xrand+dcos)-phi(yrand,xrand);
    double f6 = phi(yrand-dcos,xrand+dcos)-phi(yrand,xrand);
    double f7 = phi(yrand+dcos,xrand-dcos)-phi(yrand,xrand);
    double f8 = phi(yrand-dcos,xrand-dcos)-phi(yrand,xrand);

    res1 += f1*f1;
    res2 += f2*f2;
    res3 += f3*f3;
    res4 += f4*f4;
    res5 += f5*f5;
    res6 += f6*f6;
    res7 += f7*f7;
    res8 += f8*f8;

    //    err1 += f1*f1*f1*f1;
    //    err2 += f2*f2*f2*f2;
    //    err3 += f3*f3*f3*f3;
    //    err4 += f4*f4*f4*f4;
  }

  res = res1+res2+res3+res4+res5+res6+res7+res8;
  //  res = sqrt(res/4*counter);
  res /= (8*counter);

  //  err = err1+err2+err3+err4;
  //  err /= (4*counter);
  //  err -= res;
  //  err = sqrt(err);
  //  cout << " Info/convodif.cc:calrms: error = " << err / sqrt(4*counter) << endl;
  //  cout << endl;
 
  return(res);
} 


// ----------------------- DiffusionRadius ------------------------
// This subroutin calculates the structure function of screen phase (phi)
// vs distance, in the form of sqrt(<(phi(x,y)-phi(x',y'))^2>) as a function
// of distancebetween (x,y) and (x',y'), inside the screen.


// DataTable DiffusionRadius(TMatrix<r_8> &phi, double del1, double rdiff)
DataTable DiffusionRadius(TMatrix<r_4> &phi, double del1, double rdiff, double &reRdiff, double gambeta)
{

  // rdiff: diffusion radius in meters
  // reRdiff: recalculated Rdiff by computing the phase structure function

  DataTable diffrad; //(2000);
  //  diffrad.AddDoubleColumn("r");   // x axis is the radius
  //  diffrad.AddDoubleColumn("rms"); // y axis is the root mean squar
  //  diffrad.AddDoubleColumn("rkcm"); // x axis considered in centimeter

  diffrad.AddFloatColumn("r");  
  diffrad.AddFloatColumn("rms");
  diffrad.AddFloatColumn("rkcm");

  double minSF = 10.;
  reRdiff = 0.; 
  int n =int(round(rdiff/del1)); // diffusion radius in the screen pixelsize
 
  double dat[3]; // the data table contains 3 columns
  

  //for(int dist = 2; dist<=100; dist++){ // distance in pixel
  for(int dist = 2; dist<=5*n; dist++){ // distance in pixel    
     double SF = calrms(phi, dist); // phase structure function
     dat[0] = dist;
     dat[1] = SF; 
     if( (SF>0.9*gambeta)&&(SF<1.1*gambeta) )
       cout <<" Structure Func. = " << SF << "  at " << dist << "th screen pixel" << endl;
     dat[2] = dist*del1; 
     diffrad.AddRow(dat);

     if ( fabs(SF - 1.) < fabs(minSF - 1.) ) {
       minSF = SF;
       reRdiff = dist * del1;
     }

  } // end of for

   //  if((rdiff/del1-n) > .5) n++;
  double RMSn1 = calrms(phi, n); 
  double RMSn2 = calrms(phi, n+1);
  double dphi = RMSn1 + (rdiff-n*del1)*(RMSn2-RMSn1)/del1; // interpolation
  cout << " Structure Func. = "  << dphi << " at Rdiff " << endl;
  cout << " Structure Func. = "  << calrms(phi,int(round(reRdiff/del1))) << " at reRdiff " << endl;
//  cout <<" Structure Func. if using Fried radius = "  << dphi/6.88 << " at Rdiff " <<endl;

  return(diffrad);
} // end of DiffractionRadius

// ----------------------------- compouteStrucFunc -------------------------
// Returns the structure Function value at sepration Rdiff
// double ComputeStrucFunc(TMatrix<r_8> &phi, double del1, double rdiff)
double ComputeStrucFunc(TMatrix<r_4> &phi, double del1, double rdiff)
{

  // rdiff: diffusion radius in meters

  int n =int(rdiff/del1); // diffusion radius in the screen pixelsize
  if((rdiff/del1-n) > .5) n++;
  double RMSn = calrms(phi, n); 
  double RMSn1 = calrms(phi, n+1); 

  cout <<" Structure Func. = "  << RMSn << " at " << n << " screen pixel (Rdiff) " 
       <<endl;
  cout <<" Structure Func. = "  << RMSn1 << " at " << n+1 << " screen pixel " 
       <<endl; 
  return(RMSn);
}


// ----------------------------------- contradel1 -------------------------------
// This subrutin puts del1 tvector and the related point like source contrast tvector 
// in a datatable

DataTable contradel1(TVector<r_8> d1, TVector<r_8> cont){

  DataTable conde(100);
  conde.AddDoubleColumn("del1");
  conde.AddDoubleColumn("contrast");
  double dat[3];
  for(int i=0; i<d1.Size();i++){
    dat[0] = d1(i)/1000;
    dat[1] = cont(i);
    conde.AddRow(dat);
  }
  return(conde);
} 


// --------------------------------- Compute ReSampling ----------------------------------
// This subroutin resample the ohase screen by dividing the N by two and fixing the 
// physical length of the screen, So, the screeen pixel size (del1)will be doubled.

// void ComputeReSampling(TMatrix<r_8> &phi, double &del1)
void ComputeReSampling(TMatrix<r_4> &phi, double &del1)
{

  int newNy = int(phi.NRows()/2);
  int newNx = int(phi.NCols()/2);

  //  TMatrix <r_8> newphi(newNy, newNx);
  TMatrix <r_4> newphi(newNy, newNx);

  //  cout << " DEBUGGING:convodif.cc:computerReSampling: " << endl;
  //  cout << " newphi is memory alocated: newphi.NCols() = " << newphi.NCols() 
  //       << " newphi.NRows() = " << newphi.NRows() << endl;

  for(int j=0; j<phi.NRows(); j+=2){
    int jp = int(j/2);
    for(int i=0; i<phi.NCols(); i+=2)
      {
	int ip = int(i/2);
	newphi(jp,ip) = ( phi(j,i) + phi(j+1,i) + phi(j,i+1) + phi(j+1,i+1) ) / 4.;
	//	cout << "newphi(" << jp << "," <<ip <<") = "<< newphi(jp,ip) << endl;
      }
  }

  del1 *= 2.;
  phi.ReSize(newNy,newNx);

  //  cout << " DEBUGGING:convodif.cc:computerReSampling: " << endl;
  //  cout << " phi is resized " << endl;

  phi = newphi;

  //  cout << " DEBUGGING:convodif.cc:computerReSampling: " << endl;
  //  cout << " phi is filled " << endl;

} 


// -------------------------------- Compute Resized Matrix ------------------------
// Thid subroutin reproduces a new diffraction image from a larger one. The new image
// will have a smaller physical size by homothety coefficient of alpha. To reproduce a
// smaller image pixel values, the interpolation of the four nearest neighbor pixels of
// the initial image (diff) to the new pixel are considerd.

TMatrix <r_8> ComputeResizedMatrix(TMatrix <r_8> diff, double homoth)
{
  int i,j; // pixel coordinate in initial image (diff)
  TMatrix <r_8> Ip(diff.NRows(),diff.NCols()); // resized image
  double l0 = diff.NCols()/2*(1-homoth); // physical length diffrence between initial
                                         // and new images

  for(int ip=0; ip<diff.NRows();ip++)
    for(int jp=0; jp<diff.NCols();jp++){

      i = int(l0+homoth*ip);
      j = int(l0+homoth*jp);
 
      Ip(ip,jp) = (i + 1 - l0 - homoth*ip)*(j + 1 - l0 - homoth*jp)*diff(i,j) 
	+ (l0 + homoth*ip - i)*(j + 1 - l0 - homoth*jp)*diff(i+1,j)
	+ (i + 1 - l0 - homoth*ip)*(l0 + homoth*jp - j)*diff(i,j+1)
	+ (l0 + homoth*ip -i)*(l0 + homoth*jp -j)*diff(i+1,j+1);

	}

  return(Ip);
}


// -------------------------------- Compute Resampled Matrix ------------------------
// This subroutin re-grid the images with \lambda > \lambda_0.
TMatrix <r_4> ComputeResampledMatrix(TMatrix <r_4> M0, TMatrix <r_4> M)
{

// Hereafter, the pattern of the smallest wavelength is called s-image and
// pattern with any other wavelength is called current image.

// Pixel size of the s-image  
double delx0 = M0.Info()["del0x"];
double dely0 = M0.Info()["del0y"];
//cout << " delx0 = " <<  delx0/1000 << " km " << endl;
//cout << " dely0 = " <<  dely0/1000 << " km " << endl;
//cout << endl;

//pixel size of the current image
double delx  = M.Info()["del0x"];
double dely  = M.Info()["del0y"];
//cout << " delx = " <<  delx/1000 << " km " << endl;
//cout << " dely = " <<  dely/1000 << " km " << endl;
//cout << endl;

int Nx0 = M0.NCols();
int Ny0 = M0.NRows();
//cout << " Nx0 = " <<  Nx0 << endl;
//cout << " Ny0 = " <<  Ny0 << endl;
//cout << endl;

int Nx = M.NCols();
int Ny = M.NRows();
//cout << " Nx = " <<  Nx << endl;
//cout << " Ny = " <<  Ny << endl;
//cout << endl;

// phisical size of the s-image 
double Lx0 =  Nx0 * delx0;
double Ly0 =  Ny0 * dely0;
//cout << " Lx0 = " <<  Lx0/1000 << " km " << endl;
//cout << " Ly0 = " <<  Ly0/1000 << " km " << endl;
//cout << endl;

// phisical size of the current image 
double Lx =  Nx * delx;
double Ly =  Ny * dely;
//cout << " Lx = " <<  Lx/1000 << " km " << endl;
//cout << " Ly = " <<  Ly/1000 << " km " << endl;
//cout << endl;

// Sampling the physical size of the s-image (Lx0) in current image step (delx)
int Mpx = int(Lx0/delx);
int Mpy = int(Ly0/dely);
//cout << " Mpx = " <<  Mpx << endl;
//cout << " Mpy = " <<  Mpy << endl;
//cout << endl;

// Putting the middle points (originjs) of s-image and current image on each other.  
int i0 = (Lx-Lx0)/2/delx;
int j0 = (Ly-Ly0)/2/dely;
//cout << " i0 = " <<  i0 << endl;
//cout << " j0 = " <<  j0 << endl;
//cout << endl;

// Mp: an extraction of the current image with the same physical size of the s-image and the same centre
TMatrix<r_4> Mp(Mpy,Mpx);
for(int jj = 0; jj < Mpy; jj++)
  for(int ii = 0; ii < Mpx; ii++){
    Mp(jj,ii) = M(j0+jj,i0+ii); 
  }

// Mf: the resampled current image with pixel size of the s-image
TMatrix<r_4> Mf(Ny0,Nx0);
double hx = delx0/delx;
double hy = dely0/dely;

// The equivalent pixel numbers in resampled (jp & ip) and current image (j & i)
// for the same physical length from the top left of the images.
for(int jp = 0; jp < Ny0; jp++)
  for(int ip = 0; ip < Nx0; ip++){
    int i = int(ip*hx);
    int j = int(jp*hy);
    
// To do bilinear interpolation, we need to find 4 nearest neighbour of the pixel Mf(jp,ip) from Mp(j,i).
// To search for nearest neighbours j+1 and j-1 members are neaded.   
    if ( (j>0)&&(j<Mpy-1) ){
      // six nearest pixels
      TMatrix<r_4> d(6,3);
      int k = 0;
      for(int nj=j-1;nj<j+2;nj++)
	for(int ni=i;ni<i+2;ni++){
	  d(k,0) = nj;
	  d(k,1) = ni;
	  d(k,2) = sqrt( (ni*delx-ip*delx0)*(ni*delx-ip*delx0) + (nj*dely-jp*dely0)*(nj*dely-jp*dely0) );
	  k++;
	}

//      if ( (jp==8000) && (ip==12000) ) cout << d << endl; 
      
      // Sorting the pixel coordinates and distance form (jp,ip) according to their smallest distance
      // from (jp,ip)
      for(k=0; k<6; k++){
	int c1 = d(k,0);
	int c2 = d(k,1);
	int c3 = d(k,2);
	for(int kk=k; kk<6; kk++)
	  if(c3>d(kk,2)) {
	    d(k,0) = d(kk,0); d(k,1) = d(kk,1); d(k,2) = d(kk,2);
	    d(kk,0) = c1; d(kk,1) = c2; d(kk,2) = c3;
	    c1 = d(k,0); c2 = d(k,1); c3 = d(k,2);
	  }
      }
//      if ( (jp==8000) && (ip==12000) ) cout << d << endl;   
//    Computing the bilinear interpolation of the fidr 4 nearest neghbours of the (jp,ip)
      if ( d(0,0) == d(1,0) ){
	double fj1 = (fabs(ip*delx0-d(0,1)*delx)*Mp(d(1,0),d(1,1))+fabs(ip*delx0-d(1,1)*delx)*Mp(d(0,0),d(0,1)))/delx;
	double fj2 = (fabs(ip*delx0-d(2,1)*delx)*Mp(d(3,0),d(3,1))+fabs(ip*delx0-d(3,1)*delx)*Mp(d(2,0),d(2,1)))/delx;
	Mf(jp,ip) = (fabs(jp*dely0-d(2,0)*dely)*fj1 + fabs(jp*dely0-d(0,0)*dely)*fj2)/dely;
      }
      else if ( d(0,1) == d(1,1) ){
	double fi1 = (fabs(jp*dely0-d(0,0)*dely)*Mp(d(1,0),d(1,1))+fabs(jp*dely0-d(1,0)*dely)*Mp(d(0,0),d(0,1)))/dely;
	double fi2 = (fabs(jp*dely0-d(2,0)*dely)*Mp(d(3,0),d(3,1))+fabs(jp*dely0-d(3,0)*dely)*Mp(d(2,0),d(2,1)))/dely;
	Mf(jp,ip) = (fabs(ip*delx0-d(2,1)*delx)*fi1 + fabs(ip*delx0-d(0,1)*delx)*fi2)/delx;
      }  
      else {cout << " DEB: Not a square! " << endl;}
      
    }


    else if(j==0){
      TMatrix<r_4> d(4,3);
      int k = 0;
      for(int nj=j;nj<j+2;nj++)
	for(int ni=i;ni<i+2;ni++){
	  d(k,0) = nj;
	  d(k,1) = ni;
	  d(k,2) = sqrt( (ni*delx-ip*delx0)*(ni*delx-ip*delx0) + (nj*dely-jp*dely0)*(nj*dely-jp*dely0) );
	  k++;
	}
      for(k=0; k<4; k++){
	int c1 = d(k,0);
	int c2 = d(k,1);
	int c3 = d(k,2);
	for(int kk=k; kk<4; kk++)
	  if(c3>d(kk,2)) {
	    d(k,0) = d(kk,0); d(k,1) = d(kk,1); d(k,2) = d(kk,2);
	    d(kk,0) = c1; d(kk,1) = c2; d(kk,2) = c3;
	    c1 = d(k,0); c2 = d(k,1); c3 = d(k,2);
	  }
      }
      if ( d(0,0) == d(1,0) ){
	double fj1 = (fabs(ip*delx0-d(0,1)*delx)*Mp(d(1,0),d(1,1))+fabs(ip*delx0-d(1,1)*delx)*Mp(d(0,0),d(0,1)))/delx;
	double fj2 = (fabs(ip*delx0-d(2,1)*delx)*Mp(d(3,0),d(3,1))+fabs(ip*delx0-d(3,1)*delx)*Mp(d(2,0),d(2,1)))/delx;
	Mf(jp,ip) = (fabs(jp*dely0-d(2,0)*dely)*fj1 + fabs(jp*dely0-d(0,0)*dely)*fj2)/dely;
      }
      else if ( d(0,1) == d(1,1) ){
	double fi1 = (fabs(jp*dely0-d(0,0)*dely)*Mp(d(1,0),d(1,1))+fabs(jp*dely0-d(1,0)*dely)*Mp(d(0,0),d(0,1)))/dely;
	double fi2 = (fabs(jp*dely0-d(2,0)*dely)*Mp(d(3,0),d(3,1))+fabs(jp*dely0-d(3,0)*dely)*Mp(d(2,0),d(2,1)))/dely;
	Mf(jp,ip) = (fabs(ip*delx0-d(2,1)*delx)*fi1 + fabs(ip*delx0-d(0,1)*delx)*fi2)/delx;
      }  
      else {cout << " DEB: Not a square! " << endl;}
  
    }
    
    
      
  }
  
return(Mf);

}






    
// -------------------------------- Compute Lens Phi ----------------------------------------
// This subroutin is a test of computing the wave front after paasing through a simple lens.
// The lens is a biconvex with the same curvature in both parts.
// f: focal length, D0: lens thikness. The lengths are in meter

void ComputeLensPhi(TMatrix <r_8> &phi, int N, double f, double D0, double lambda, double del1)
{

  f*=8;

  cout << " ((((((((((((((((((((((((((((((((((((((((((((((((( " << endl;
  cout << " DEBUGGING:CONVODIF.CC:COMPUTELENSPHI: " <<endl;
  cout << " f = " << f << endl;
  cout << " lambda = " << lambda << endl;
  cout << " del1 = " << del1 << endl; 
  cout << " ))))))))))))))))))))))))))))))))))))))))))))))))) " << endl;

  phi.ReSize(N,N);
  double pi = 3.141592635;

  for(int j = 0; j < phi.NRows(); j++)
    for(int i =0; i < phi.NCols(); i++){
      int r = 40; // lens radius
      int d = r; // lenses distance
      //      phi(j,i) = (2*pi/lambda) * ( D0 - 1 / (2*f/del1/del1) * ( (i-N/2)*(i-N/2) + (j-N/2)*(j-N/2) ) ); 
      //      phi(j,i) = - (2*pi/lambda) * ( (1 / (2*f/del1/del1)) * ( (i-N/2)*(i-N/2) + (j-N/2)*(j-N/2) ) ) ;
      //  phi(j,i) = - (2*pi/lambda) * ( (1 / (2*f/del1/del1)) * ( (i-N/4)*(i-N/4) + (j-N/2)*(j-N/2) ) ) ; // not centered lens

      // a plan with a samall lens inside, the lens radius = 40 pixels
      //      if( ((i-N/4)*(i-N/4) + (j-N/4)*(j-N/4)) > 1600 )  phi(j,i) =  (2*pi/lambda) * ( (1 / (2*f/del1/del1)) * ( r*r ) ) ;
      //      else phi(j,i) = - (2*pi/lambda) * ( (1 / (2*f/del1/del1)) * ( (i-N/4)*(i-N/4) + (j-N/4)*(j-N/4) ) ) ;

      //a plan with double lenses, biconvex and biconcave
      
      if( ((i-N/2+d)*(i-N/2+d) + (j-N/2)*(j-N/2)) > 1600 )  phi(j,i) =  0.; 
      if( ((i-N/2-d)*(i-N/2-d) + (j-N/2)*(j-N/2)) > 1600 )  phi(j,i) =  0.; 
      if( ((i-N/2+d)*(i-N/2+d) + (j-N/2)*(j-N/2)) <= 1600 ) phi(j,i) =  (2*pi/lambda) * (- (1 / (2*f/del1/del1)) * ( (i-N/2+d)*(i-N/2+d) + (j-N/2)*(j-N/2) - r*r) ); // biconvex
      if( ((i-N/2-d)*(i-N/2-d) + (j-N/2)*(j-N/2)) <= 1600 ) phi(j,i) =  (2*pi/lambda) * ( (1 / (2*f/del1/del1)) * ( (i-N/2-d)*(i-N/2-d) + (j-N/2)*(j-N/2) - r*r) ); // biconcave


    }

}
  
                       
// ---------------- ReSizePatt -----------------
TMatrix<r_8> ReSizePatt(const TMatrix<r_8> &currpatt, double A)
{

  // This code resizes the physical sizeof the input diffraction pattern (currpatt)
  // and computes the new intensity for each new pixel by interpolating among the first
  // four neighbours. both patterns are at the same center
  // A: lambda0 / lambda < 1

  int Ny = currpatt.NRows();
  int Nx = currpatt.NCols();

  ResourceUsage res;
/*
  res.Update();
  cout << " convodif.cc:ReSizePatt:Resource usage info before alocating patt: \n" << res << endl;
  cout << endl;
*/

  TMatrix<r_8> patt(Ny,Nx);

/*
  res.Update();
  cout << " convodif.cc:ReSizePatt:Resource usage info: \n" << res << endl;
  cout << endl;
*/
  double l0 = (1-A) * Ny / 2; 

  for(int ip = 0; ip < Ny; ip++)
    for(int jp = 0 ; jp < Nx; jp++){

      int i = int(l0 + ip*A);
      int j = int(l0 + jp*A);

      patt(ip,jp) = 
	(i+1-l0-A*ip)*(j+1-l0-A*jp)*currpatt(i,j) +
	(l0+A*ip-i)*(j+1-l0-A*jp)*currpatt(i+1,j) +
	(i+1-l0-A*ip)*(l0+A*jp-j)*currpatt(i,j+1) +
	(l0+A*ip-i)*(l0+A*jp-j)*currpatt(i+1,j+1);
    }

  return(patt);
}



// ------------ ComputeSigma2D ---------------------
double ComputeSigma2D(TMatrix<r_4> phi)
{
  double sig2 = 0.;
  double mean = 0.;
  int Ny = phi.NRows();
  int Nx = phi.NCols(); 

  for (int i = 0; i < Ny; i++){
    for (int j = 0; j < Nx; j++){

      mean += phi(i,j);
      sig2 += phi(i,j) * phi(i,j);

    }
    //    cout << "DEB:convodif.cc:CpmputeSigma2D:mean = " << mean << endl;
  }
  mean /= Nx*Ny;
  sig2 = sig2/Nx/Ny - mean*mean;
  //  cout << " convodif.cc:CpmputeSigma2D: phase mean = " << mean << endl;

  return(sig2);
}


// ---------------- Compute1dSpec -------------------
Vector Compute1dSpec(TMatrix< complex<r_4> > fftphi, double del, double xmin, double xmax, 
		     int nbin, double step)
{
  // Computing the 1d spectrum from the 2D generated fftphi
  // del in kilometere

  int Ny = fftphi.NRows() ;
  int Nx = fftphi.NCols() ;

  double Ly = Ny*del ;
  double Lx = 2*(Nx-1)*del ;

  Vector spec(nbin) ;

  cout << endl;
  cout << " convodif.cc:Compute1dSpec: Checking the parametres: " << endl;
  cout << " fftphi.NCols() = " << fftphi.NCols() << "  Lx = " << Lx << endl;
  cout << " fftphi.NRows() = " << fftphi.NRows() << "  Ly = " << Ly << endl;



  //  HProf sd(xmin,xmax,nbin,1.e-20,1.e24) ;
  HProf sd(xmin,xmax,nbin,1.e-30,1.e30) ;

  double sumfftphi = 0.;
  for(int ky = 0; ky < Ny; ky++){
    int kyy = ky ;
    if (ky > Ny/2) kyy = Ny - ky; // - (Ny/2+1);

    for(int kx = 0; kx < Nx; kx++){
      double wk = sqrt(kx*kx/Lx/Lx+kyy*kyy/Ly/Ly) ; 
      complex<r_4> za = fftphi(ky,kx) ;
      double Syx = (za.real()*za.real() + za.imag()*za.imag())/Lx/Ly ;
      double s1 = 2*pi*wk * Syx * step;
      sd.Add(wk,s1);
      
      sumfftphi += Syx;

    } // end for(kx)
  } // end for(ky)


  sumfftphi /= Lx*Ly;
  cout << " sumfftphi = " << sumfftphi << endl;
  cout << endl;

  sd.GetValue(spec) ;

  return(spec);
}


// -------------- Function intpol ------------------
double intpol(double c1, double c2, int n, double d, double rdiff)
{

  cout << "intpol: " << endl;
  cout << " c1 = " << c1 << "  c2 = " << c2 << "  n = " << n << endl;
  cout << " rdiff = " << rdiff << "  d= " << d << endl;
  cout << endl;
 
 
  double x = c1 + (rdiff-n*d)*(c2-c1)/(d);
  return x;

}

// --------------- ComputeNph ----------- 
double ComputeNph(double mag, double Texp, int Nx, int Ny, double lambda, double qeff, double del0x, double del0y, double D) 
// magnitude, exposure time, Nx, Ny, wanelength, quantum efficiency, mirro diameter
{
    double pi = 3.141592635;     // pi number  
    double nm2m = 1.e-9;         // nanometer to centimeter
    double m2micron = 1.e6;
    double micron2cm = 1.e-4;
    double h = 6.6261e-34;    // plank constant SI
    double kb = 1.3807e-23;   // Boltzmann constant SI
    double c = 2.9979e8;      // light speed SI
    double Tsolar = 5778;     // solar temperature in K
    double Rsolar = 6.9634e8;  // solar radius SI
    double pc2m = 3.0857e16;     // parsec to meter
    double dlambda = 3500*nm2m; //100*nm2m;   // passband width

      double lam = lambda * .01; // in meter 
      cout << " Lam = " << lam << " m" << endl;
      double Isolar = 2*h*c*c/lam/lam/lam/lam/lam/(exp(h*c/lam/kb/Tsolar)-1); // Planck radiation
                // flux per meter per stradian
      cout << " Isolar = " << Isolar << " w/m^2/lambda" << endl;
      
      double Lsolar = Isolar * 4*pi * Rsolar*Rsolar * dlambda; // solar luminosity in watt
      cout << " Lsolar = " << Lsolar << " w " << endl;
      
      double Lvega = 40.12 * Lsolar;
      double Dvega = 7.68 * pc2m;
      double Fvega = Lvega/4/pi/Dvega/Dvega;
      cout << " Vega Flux = " << Fvega << " w/m^2 " << endl; 
      Fvega /= (m2micron * m2micron);
      cout << " Vega flux = " << Fvega << " w/micron^2 " << endl;
      cout << endl;

//      double Flux = Fvega * pow(.3981,mag);
      double Flux = Fvega * pow(10.,-mag/2.5);    
      cout << " Source apparent magnitude = " << mag << endl;
      cout << " Source flux = " << Flux << " w/micron^2 " << endl;      
      double qeffFlux = Flux * qeff;
      cout << " Source effective flux = " << qeffFlux << " w/micron^2 " <<endl;
      
//      double intensity = effFlux * (del0x/micron2cm)*(del0x/micron2cm); //watt per pixel
//      double intensity = qeffFlux * pi*(D/micron2cm)*(D/micron2cm)/4.; //watt (crossing mirorr surface)
      double intensity = qeffFlux * Nx*Ny *del0x*del0y / (micron2cm*micron2cm); 
                      //watt (crossing mirorr surface)
      double energy = intensity * Texp; 
      double Nph = energy*lam/h/c; // total number of photons 


//      double totenergy = effFlux * (D/micron2cm)*(D/micron2cm)/4 * Texp;
      cout << " Exposure time  = " << Texp << " s" << endl;
      cout << " Crossing Intensity = " << intensity << " w" << endl;
      cout << " Energy = " << energy << " J" << endl;
//      cout << " Total Energy = " << totenergy << " J" << endl;
      cout << " Number of photons = " << Nph << endl;
      cout << endl;
      
      return Nph;
	
}


// --------------- ComputeNph2 ----------- 
double ComputeNph2(double mag, double Texp, int Nx, int Ny, double lambda, double qeff, double del0x, double del0y, char *filter) 
// magnitude, exposure time, Nx, Ny, wanelength, quantum efficiency, file name including filter wavelengths and transmissions
{
    double pi = 3.141592635;     // pi number  
    double nm2m = 1.e-9;         // nanometer to centimeter
    double m2micron = 1.e6;
    double micron2cm = 1.e-4;
    double a2m = 1.e-10;
    double h = 6.6261e-34;    // plank constant SI
    double kb = 1.3807e-23;   // Boltzmann constant SI
    double c = 2.9979e8;      // light speed SI
    double Tsolar = 5778;     // solar temperature in K
    double Rsolar = 6.9634e8;  // solar radius SI
    double pc2m = 3.0857e16;     // parsec to meter
    double dlambda = 100*nm2m;   // passband width

    double lam = lambda * .01; // in meter 
    cout << " Lam = " << lam << " m" << endl;

    ifstream inp(filter, ios::in | ios::binary);
    vector<double> x,y;
    double a,b;
    inp >> a; inp >> b;
    x.push_back(a);
    y.push_back(b);
    while(!inp.eof())
    {
     inp >> a; x.push_back(a);
     inp >> b; y.push_back(b);
//     cout << " ComputeNph2:" << j+1 << ". a=" << a << "  b=" << b << endl; 
//     cout << " ComputeNph2: x=" << x[j] << "  y=" << y[j] << endl; 
    } 
    x.pop_back();
    y.pop_back();
    int nv = x.size();
//    for(int i=0; i<nv; i++)
//      cout << " ComputeNph2: x=" << x[i] << "  y=" << y[i] << endl; 

    cout << " ComputeNph2: # of wavelength in " << filter << "=" << nv << endl;
    double lmin = x[0];
    double lmax = x[nv-1];
    double dpa = (lmax-lmin)/(nv-1); //angestrom
    cout << " ComputeNph2: lmin=" << lmin << ", lmax=" << lmax << ", dpa=" << dpa << " angestrom " << endl;    

    int N = 128; //3*nv;
    Vector xf(N), yf(N);
    for(int i=0; i<N; i++) 
    {
     if( i<nv) { yf(i) = 0.; xf(i) = 2*lmin - lmax + dpa*i;}
     else if (i>2*nv-1) { yf(i) = 0.; xf(i) = lmax + dpa*(i-2*nv+1);}
     else { xf(i) = x[i-nv]; yf(i) = y[i-nv];}
//     cout << i+1 << ", " << xf(i) << ", " << yf(i) << endl; 
    }

    x.clear();
    y.clear();
    
    Vector yp(N);
    for(int i=0; i<N; i++)
    { 
     double lambd = xf(i)*a2m;
     yp(i) =  2*h*c*c/lambd/lambd/lambd/lambd/lambd/(exp(h*c/lambd/kb/Tsolar)-1); // Planck radiation
                // flux per meter per stradian
//     cout << " DEB: ComputeNph2: i=" << i+1 << ", yp=" << yp(i) << endl; 
    }

//    cout << endl;
//    cout << "lamb trans   I " << endl;
//    for(int i=0; i<yp.NElts(); i++) cout << xf(i) << "  " << yf(i) << "  " << yp(i) << endl;

    FFTWServer fftws;
//    fftws.setNormalize(false); // false = No division by N*N through FFTForward  (true by def.)  
    TVector< complex< r_8 > > fftyf, fftyp, fftfp;
    fftws.FFTForward(yf,fftyf);
    fftws.FFTForward(yp,fftyp);

//    cout << fftyf << fftyp << endl;

    int nn = fftyp.NElts();
    fftfp.ReSize(nn);
    for(int i=0; i<nn; i++) fftfp(i) = fftyp(i)*fftyf(i);

//    cout << " fftyp:" << fftyp << endl;

//    TVector<complex<r_8> > fp; //(2*(nn-1)); 
    TVector<r_8>  fp(2*(nn-1)); 
    fftws.FFTBackward(fftfp,fp,true);
//    cout << " DEB: ComputeNph2: fp computed. " << endl;
//    cout << " fp:" << fp << endl;
    double Isolar = 0; // flux per meter per stradian
    for(int i=0; i<fp.Size(); i++) Isolar += fp(i); 
    cout << " Isolar = " << Isolar << " w/m^2/lambda" << endl;
      
    double Lsolar = Isolar * 4*pi * Rsolar*Rsolar * dpa*a2m; // * dlambda; // solar luminosity in watt
    cout << " Lsolar = " << Lsolar << " w " << endl;
      
      double Lvega = 40.12 * Lsolar;
      double Dvega = 7.68 * pc2m;
      double Fvega = Lvega/4/pi/Dvega/Dvega;
      cout << " Vega Flux = " << Fvega << " w/m^2 " << endl; 
      Fvega /= (m2micron * m2micron);
      cout << " Vega flux = " << Fvega << " w/micron^2 " << endl;
      cout << endl;

//      double Flux = Fvega * pow(.3981,mag);
      double Flux = Fvega * pow(10.,-mag/2.5);    
      cout << " Source apparent magnitude = " << mag << endl;
      cout << " Source flux = " << Flux << " w/micron^2 " << endl;      
      double qeffFlux = Flux * qeff;
      cout << " Source effective flux = " << qeffFlux << " w/micron^2 " <<endl;
      
//      double intensity = effFlux * (del0x/micron2cm)*(del0x/micron2cm); //watt per pixel
      double intensity = qeffFlux * pi*(154/micron2cm)*(154/micron2cm)/4.; //watt (crossing mirorr surface)
//      double intensity = qeffFlux * Nx*Ny *del0x*del0y / (micron2cm*micron2cm); 
                      //watt (crossing mirorr surface)
      double energy = intensity * Texp; 
      double Nph = energy*lam/h/c; // total number of photons 


//      double totenergy = effFlux * (D/micron2cm)*(D/micron2cm)/4 * Texp;
      cout << " Exposure time  = " << Texp << " s" << endl;
      cout << " Crossing Intensity = " << intensity << " w" << endl;
      cout << " Energy = " << energy << " J" << endl;
//      cout << " Total Energy = " << totenergy << " J" << endl;
      cout << " Number of photons = " << Nph << endl;
      cout << endl;
      
      return Nph;
	
}




// --------------------- ShiftPattern -----------------
TMatrix<r_4> ShiftPattern(TMatrix<r_4> const diffrac, int dx, int dy, int Nx, int Ny, int d)
// centered psf, shift in x, shift in y, mirror diameter in pixel
{
 TMatrix<r_4> diff;
 diff.ReSize(Ny,Nx);	
 if((dx-Nx/2)*(dx-Nx/2)+(dy-Ny/2)*(dy-Ny/2) < d*d/4) 
 	{
 		cout << " srout.cc: ShiftPattern: Error! Target out of field of view. " << endl;
 		exit;	
 	}
 	
// 	cout << "srout.cc DEB: Nx=" << diff.NCols() << ", Ny=" << diff.NRows() << endl;
 	
 	for(int j=0; j<Ny; j++)
  	for(int i=0; i<Nx; i++)
  	{
  		if((j+dy<Ny)&&(i+dx<Nx)) diff(j,i) = diffrac(j+dy,i+dx);
  		else if((j+dy<Ny)&&(i+dx>Nx)) diff(j,i) = diffrac(j+dy,i+dx-Nx);
  		else if((j+dy>Ny)&&(i+dx<Nx)) diff(j,i) = diffrac(j+dy-Ny,i+dx);
  		else if((j+dy>Ny)&&(i+dx>Nx)) diff(j,i) = diffrac(j+dy-Ny,i+dx-Nx);
  		
//  		cout << " j=" << j << ", i=" << i << ", diff=" << diff(i,j) << endl;		
  	}
  	
  	return diff;
	
}

// ------------------ FillFitsHeader --------------------
void FillFitsHeader(FitsImg2DWriter &fitsimag, const char * headinfo)
{

/*
 ifstream inp(headinfo, ios::in | ios::binary);
 if(!inp){cerr << " No Fits Header Info file " << endl; exit(0);}

 string st;
 char c1[50],c2[50],c3[50];
 while(!inp.eof()) 
 {
//  getline(inp,st);
//  istringstream istst(st);
//  istst >> c1 >> c2 >> c2;
  inp >> c1 >> c2 >> c3;
  cout << "DEB: srout.cc: FillFitsHeader: c1=" << c1 << ", c2=" << c2 << ", c3=" << c3 << endl;
  fitsimag.WriteKey(c1, c2, c3); 	
 }
*/
 fitsimag.WriteKey("FILTER", "i SDSS c6019");
 fitsimag.WriteKey("OBJECT", "cdfs_14_i");
 fitsimag.WriteKey("AIRMASS", 1.34900);
 fitsimag.WriteKey("EXPTIME", 360.000);
 fitsimag.WriteKey("PHOTFLAG", "Y");
 fitsimag.WriteKey("PHOT_K", 0.05);
 fitsimag.WriteKey("GAIN", 2.06);
 fitsimag.WriteKey("RDNOISE", 8.0);
 fitsimag.WriteKey("SATURATE", 44000.);
 fitsimag.WriteKey("MAGZPT",24.91 );
 fitsimag.WriteKey("EQUINOX", 2000.);
 fitsimag.WriteKey("CTYPE1","RA---TAN" );
 fitsimag.WriteKey("CTYPE2","DEC--TAN");
 fitsimag.WriteKey("CD1_1",-4.21315906076E-07 );
 fitsimag.WriteKey("CD2_1", 7.48998728615E-05);
 fitsimag.WriteKey("CD1_2", 7.48998728615E-05);
 fitsimag.WriteKey("CD2_2",4.21315906076E-07 );
 fitsimag.WriteKey("CRPIX1", 2271.30000000);
 fitsimag.WriteKey("CRPIX2",169.840000000);
 fitsimag.WriteKey("CRVAL1",52.2186620000);
 fitsimag.WriteKey("CRVAL2",-27.8610250000);
 fitsimag.WriteKey("LONPOLE", 180.000000000);
 fitsimag.WriteKey("LATPOLE",0.00000000000 );
 fitsimag.WriteKey("PV2_1",0.00000000000 );
}

// ------------------ ScalePattern --------------------
HProf ScalePattern(TMatrix<r_4> &diff, double &rmax, double threshold)
{
 int Nx = diff.NCols();
 int Ny = diff.NRows();

 // Computing maximum flux
 double max = -1.;
 int ii = 0;
 int jj = 0;
 for(int j=int(Ny/2); j<int(.75*Ny); j++)
  for(int i=int(Nx/2); i<int(.75*Nx); i++) 
   if (diff(j,i) > max ){ max=diff(j,i); ii=i; jj=i;} 
 
 cout << " srout.cc: ScalePattern: Threshold=" << threshold << endl; 
 cout << " srout.cc: ScalePattern: maxFlux=" << max << ", j= " << jj << ", i=" << ii << endl;

 float xmin = 0.;
 float xmax = sqrt(Nx*Nx+Ny*Ny)/2.;
 int nbin = int(xmax); 
 HProf hp(xmin, xmax, nbin);
 for(int j=0; j<Ny; j++)
  for(int i=0; i<Nx; i++)
  {
   float x = sqrt( (j-Ny/2)*(j-Ny/2)+(i-Nx/2)*(i-Nx/2) ) ;  
   hp.Add(x,diff(j,i)); 
  }

 TVector<r_8> val, absc;
 hp.GetAbsc(absc); 
 hp.GetValue(val);

 rmax = 0.;
 for(int k=1; k<nbin; k++)
 {
//  cout << " k=" << k << ", r=" << absc(k) << ", val=" << val(k) << endl;
  if(val(k)/max < threshold) {rmax = absc(k); break;}
//  if( fabs(val(k)-val(k-1)) < threshold) {rmax = absc(k); break;} 
 }
 double ave = 0.;
 int counter = 0; 
 for(int j=0; j<Ny; j++)
  for(int i=0; i<Nx; i++)
   if( ((j-Ny/2)*(j-Ny/2)+(i-Nx/2)*(i-Nx/2)) < rmax*rmax )
   {
    ++counter;
    ave += diff(j,i);
   }
 ave /= counter;
 cout << " srout.cc: ScalePattern: rmax=" << rmax << " pixel, average flux=" << ave << endl;

 diff /= ave;    

 return hp;
 
}


// ------------------ GiveDiffMax --------------------
double GiveDiffMax(TMatrix<r_4> diff, int &ii, int &jj)
{ 
 
 int Nx = diff.NCols();
 int Ny = diff.NRows();

 double max = -1.;
 ii = 0;
 jj = 0; 
 for(int j=int(Ny/2); j<int(.75*Ny); j++)
  for(int i=int(Nx/2); i<int(.75*Nx); i++) 
   if (diff(j,i) > max ){ max=diff(j,i); ii=i; jj=i;} 
 
 return max;
}


