// TPM for a smooth, spherical asteroid with output surface temperatures for
// various input parameters. The time-dependent (w/ optional depth dependence)
// heat diffusion equation is solved for an array of facets on the object.
//
// For a more general, parameterized version of this program, which doesn't include 
// depth-dependent parameters, see asteroidtpm.cpp.
//
// For use with IDL wrapper asteroidtpm_d.pro, which writes the input variables to a
// file which is to be read into this procedure.
//
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <math.h>
#include "therm.h"
#include "tpminput.h"
#include "tpmgeom.h"

#define PI   3.1415926
#define dtor 0.0174533 // degrees to radians

using namespace std;

int main(int argc, char* argv[]){

// input parameters

  float aorb = atof(argv[1]); // orbital semimajor axis
  cout <<"semimaj axis = "<< aorb <<" AU" << endl;

  float eorb = atof(argv[2]); // orbital eccentricity
  cout <<"eccentricity = "<< eorb << endl;

  float perh = atof(argv[3]); // rotation period (hours)
  cout <<"rot period   = "<< perh <<" hours" << endl;

  float A;      // bond albedo
  float emis;   // emissivity
  float obliq0; // sub-solar lat
  float solst0; // sub-solar lat

  int nxstep;   // # of depth steps per skindepth
  int ntrot;    // # of timesteps per rotation
  int norb;     // # of orbits (revolutions)
  int nlat;     // # of latitude bins

  float ti;     // thermal inertia
  float kappa;  // thermal diffusivity

  readintpm_orb(A,emis,obliq0,solst0,nxstep,ntrot,norb,nlat,ti,kappa);

// define constant values
  float sigma = 5.67e-8; // stefan-boltzmann constant
  float invalb = 1. - A;    // fraction of energy absorbed
  float obliq = obliq0*dtor; float solst = solst0*dtor; // convert degrees to radians
  int   bval = 0;

  float pers = perh*3600.; // period in seconds
  float dt = pers/ntrot;  // time step = rotation period/# of timesteps per rotation

  float skindepth = sqrt(kappa*pers/2./PI);
  float dx = skindepth/nxstep; // depth step = skin depth/# of steps per skindepth
  int   nx = 30*nxstep+1;   // # of depth steps

  //  float stepsize = 2.*dx*dx/dt; // used in cranknicol
  float fconv = dt*sqrt(kappa)/ti/dx; // convert flux (q) to dT at surface
  float kperx = sqrt(kappa)*ti/dx; // thermal conductivity/depth step at surface

// calculate orbital parameters
  float porb = sqrt(pow(aorb,3))*365.25; // period of revolution (in days)
  float nrot = porb*24./perh + 1.; // number of *solar* days per revolution (prograde)
  int   ntorb = ntrot*nrot+1; //  number of time steps over an orbit

// create vectors

  // over a rotation
  vector< vector< vector<float> > > // temp[lat][time][depth]
    temp(nlat,vector<vector<float> > (ntorb,vector<float> (nx,0.0)));
  vector<float> Temp0(nlat); // array of initial temperatures over latitude bins
  //  vector<float> kappa(nx);   // diffusivity[x]
  vector<float> lat(nlat);   // latitude
  vector<float> lon(ntrot);  // longitude
  vector<float> stepsize(nx); // stepsize in cranknicol
  //  vector<float> xarr(nx);    // depth array (not needed)

  // over an orbit
  vector<float> Fsun(ntorb); // insolation flux
  vector<double> sunang(ntorb); // solar angular diameter
  vector<float> nuang(ntorb); // true anomaly
  vector< vector<float> > fluxfrac(nlat,vector<float>(ntrot)); // fraction of incident solar energy

  vector<float> tss(ntorb); // sub-solar temperature
  vector<float> theta(ntorb); // thermal parameter

  vector<float> temp_t0(nx,0.0); vector<float> temp_t1(nx,0.0); // for cranknicol
  // used to find converged solution
  vector<float> denergy(nlat); vector<float> dflux(nlat); vector<float> dtemp(nlat);
  vector<float> Tsum(nlat);    vector<float> avgT(nlat);  vector<float> avgT0(nlat);

// populate vectors
for (int y=0;y<ntorb;y++){ // over an orbit

  double M = y*2.*PI/ntorb;
  double R = 2.*aorb - aorb*(1.-pow(eorb,2.))/(1.+eorb*cos(M));

  if (eorb == 0) nuang[y] = M; else
    nuang[y] = acos(((aorb*(1-pow(eorb,2.))/R)-1.)/eorb); // true anomaly
  if (M > PI) nuang[y] = 2.*PI - nuang[y];

  Fsun[y] = 1367./R/R; // flux at R with solar constant 1367 W/m^2
  sunang[y] = 4.*6.957e5/R/1.495978707e8; //apparent angular diameter of the Sun
  tss[y] = pow(Fsun[y]*invalb/emis/sigma,.25); // equilibrium temperature
  theta[y] = ti*sqrt(2.*PI/pers)/emis/sigma/pow(tss[y],3.); // thermal parameter
}

// set Temp0 to empirically derived Tdeep (Spencer et al. 1989; fig 2.)
 float ftss; if (theta[0] > 4.) ftss=.75;
 else ftss = .01*pow(log10(theta[0]),2)+.1*log10(theta[0])+.71;

for (int t=0;t<ntrot;t++){ // longitude loop
  lon[t]=2.*PI*t/ntrot; // get longitude array
}

for (int p=0;p<nlat;p++){ // latitude loop
  lat[p]=p*PI/(nlat-1)-PI/2.;
  Temp0[p]=ftss*tss[0]*pow(cos(lat[p]),.25); // initial temperature guess
  if (Temp0[p] != Temp0[p]) Temp0[p]=2.7; // if in permanent shadow, then T = 2.7K

  for (int j=0;j<nx;j++){ // depth loop
    if (p == 0) stepsize[j]=2.*dx*dx/dt*pow(1.02,j);// stepsize array
    temp[p][0][j]=Temp0[p];
  }
}

// begin rotations
 int t=1; int r=1; avgT0=Temp0; bool converged=false; bool lastrot=false; int latconv=0;

do { // time loop
      for (int p=0;p<nlat;p++){ // latitude loop

        // surface temperature gradient term = k(dt/dx)
	float tempgrad = kperx*(temp[p][(t-1) % ntorb][0]-temp[p][(t-1) % ntorb][1]);

	// energy radiated term = emis*sigma*T^4
	float radiated = emis*sigma*pow(temp[p][(t-1) % ntorb][0],4);

	double coszangle = coszang(obliq*cos(nuang[t % ntorb]+solst),nuang[t % ntorb]-PI,
				   lat[p],lon[t % ntrot],sunang[t % ntorb]);

	// energy balance to get surf temp from previous time step
	temp_t1[0] = temp[p][(t-1) % ntorb][0] + fconv*(Fsun[t % ntorb]*invalb*coszangle 
	              - tempgrad - radiated);

//	if (abs(coszangle) < cos(sunang[t % ntorb])) temp_t1[0] = temp[p][(t-1) % ntorb][0] + 
//		    fconv*(Fsun[t % ntorb]*invalb*(sunang[t % ntorb] - abs(coszangle))*0.2 -
//			   tempgrad - radiated);
//	else temp_t1[0] = temp[p][(t-1) % ntorb][0] + fconv*(Fsun[t % ntorb]*invalb*coszangle 
//						   - tempgrad - radiated);

	// solve for temperatures at depth
	temp_t0=temp[p][(t-1) % ntorb];
        cranknicol_orb(temp_t0, kappa, stepsize, temp_t1);
	temp[p][t % ntorb]=temp_t1;

	// get values for energy balance
	/* dE = -4*emis*sigma*T^3 dT or dT = dE/-4*emis*sigma*T^3 */
	denergy[p] += Fsun[t % ntorb]*coszangle - radiated; // dE
	dflux[p] += -4.*emis*sigma*pow(temp[p][t % ntorb][0],3.); // dF = -4*emis*sigma*T^3
	Tsum[p] += temp[p][t % ntorb][0]; // sum of surf T over a rotation
      }

      // print to screen the # of rotations			       
      cout << "\r# rotations = " << (t+1)/ntrot << flush;

	// adjust deep temperature
      if ((t+1) % ntorb == 0 && r < norb) {
	for (int p=0;p<nlat;p++) {
	  float dtemp = denergy[p]/dflux[p]; // dT = denergy/dF
	  float avgT = Tsum[p]/ntorb;

	  Temp0[p]=avgT;//-dtemp;
	  for (int j=0;j<nx;j++) temp[p][t % ntorb][j]=temp[p][t % ntorb][j] + 
				   (Temp0[p]-temp[p][t % ntorb][nx-1]);
	  denergy[p]=0.;}
	dflux=denergy; Tsum=denergy;}

      if ((t+1) % ntorb == 0) {cout << "\nrevolution " << r << " completed!" << endl;
	r++;} t++;
} while (t != norb*ntorb); // time continues

// write surface temperatures to file
 char buffer [100];
 int z = snprintf(buffer, 100, "surftemporb _%f.%f.%f.tbl", aorb, eorb, perh);

  ofstream surftemp;
  //  surftemp.open("surftemptpm_orb.tbl");
  surftemp.open(buffer);

  for (int phi=0;phi<nlat;phi++){
    for (int theta=0;theta<ntorb;theta++) surftemp << temp[phi][theta][0] << ' ';
    surftemp << endl;
  }
  surftemp.close();

// write temperature profiles to file
  ofstream deeptemp;
  deeptemp.open("deeptemptpm_orb.tbl");

  for (int theta=0;theta<ntorb;theta++){
    for (int z=0;z<nx;z++) deeptemp << temp[2][theta][z] <<' ';
    deeptemp << endl;
  }

  deeptemp.close();

  return 0;
}
