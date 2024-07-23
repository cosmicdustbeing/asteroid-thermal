// TPM for a smooth, spherical asteroid with output surface temperatures for
// various input parameters. The time-dependent heat diffusion equation is solved for
// an array of facets on the object.
//
// This version doesn't require specific information about the thermal inertia, 
// rotation period, Albedo, or Heliocentric distance. Instead, the heat conduction
// equation is parameterized in terms of the thermal parameter and the theoretical
// equilibrium sub-solar temperature in an effort to generalize the equations used
// and effectively reduce the number of input variables.
//
// For a version of this program with specific input parameters, see asteroidtpm_d.cpp
// and asteroidtpm_t.cpp.
//
// For use with IDL wrapper asteroidtpm.pro, which writes the input variables to a
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
#include "cnpy.h"

#define PI   3.1415926
#define dtor 0.0174533 // degrees to radians

using namespace std;

int main(){

// input parameters

  float theta;  // thermal parameter
  float sslat0; // sub-solar lat
  float sslon0; // sub-solar lon

  int nzstep;    // # of depth steps per skindepth
  int ntrot;     // # of timesteps per rotation
  int nrot;      // maximum # of rotations
  int nlat;      // # of latitude bins

  readintpm_g(theta,sslat0,sslon0,nzstep,ntrot,nrot,nlat);

// define constant values

  float sslat=sslat0*dtor; float sslon=sslon0*dtor; // convert to radians

  float dtp=2*PI/ntrot/theta;  // time step (a step around the circumference)
  int   ntp=ntrot*nrot+1; // total # of time steps

  float dz=(float) 1/nzstep; // depth step = skin depth/# of steps per skindepth
  int   nz=20*nzstep+1;   // # of depth steps (for 20 skin depths)

  float stepsize=2.*dz*dz/theta/dtp; // used in cranknicol
  float fconv=dtp/dz; // convert flux (q) to dT' at surface
  float perz=theta/dz; // thermal parameter/depth step at surface

// create vectors
  vector< vector< vector<float> > > // temp[lat][time][depth]
   temp(nlat,vector<vector<float> > (ntp,vector<float> (nz,0.0)));
  vector<float> Temp0(nlat); // array of initial temperatures over latitude bins
  vector<float> lat(nlat); // latitude
  vector<float> lon(ntrot);// longitude
  vector<float> zarr(nz);  // depth array
  vector< vector<float> > fluxfrac(nlat,vector<float>(ntrot)); // incident angle
  vector<float> temp_t0(nz,0.0); vector<float> temp_t1(nz,0.0); // for cranknicol
  // used to find converged solution
  vector<float> denergy(nlat); vector<float> dflux(nlat); vector<float> dtemp(nlat);
  vector<float> Tsum(nlat); vector<float> avgT(nlat); vector<float> avgT0(nlat);
  // set Temp0 to empirically derived Tdeep (Spencer et al. 1989; fig 2.)
  float ftss; if (theta > 2.) ftss=0.75;
  else ftss=.01*pow(log10(theta),2)+.1*log10(theta)+.71;

// populate vectors
for (int p=0;p<nlat;p++){ // latitude loop
  lat[p]=p*PI/(nlat-1)-PI/2.;
  Temp0[p]=ftss*pow(cos(lat[p]-sslat),.25); // initial temperature guess    
  if (Temp0[p] != Temp0[p]) Temp0[p]=0.01; // if in permanent shadow T' = 0.01

  for (int t=0;t<ntrot;t++){ // time(longitude) loop				   
    if (p == 0) lon[t]=2.*PI*t/ntrot; // get longitude array			   
    fluxfrac[p][t]=coszang(sslat,sslon,lat[p],lon[t],0.); // get solar incidence angle
  }
  for (int j=0;j<nz;j++){ // depth loop
    if (p == 0) zarr[j]=j*dz; // depth array
    temp[p][0][j]=Temp0[p];
  }
}
int t=1; avgT0=Temp0; bool converged=false; bool lastrot=false; int latconv=0;

// begin rotations
 cout <<"begin"<< endl;
do { // time loop
      for (int p=0;p<nlat;p++){ // latitude loop

        // surface temperature gradient term = theta*(dT'/dz)
	float tempgrad=perz*(temp[p][t-1][0]-temp[p][t-1][1]);

	// energy radiated term = T'^4
	float radiated=pow(temp[p][t-1][0],4.);

	// energy balance to get surf temp from previous time step
	temp_t1[0]=temp[p][t-1][0]+fconv*(fluxfrac[p][t % ntrot]-tempgrad-radiated);

	// solve for temperatures at depth
	temp_t0=temp[p][t-1];
        cranknicol_g(temp_t0, stepsize, temp_t1);
	temp[p][t]=temp_t1;

	// get values for energy balance
	denergy[p] += fluxfrac[p][t % ntrot]-radiated; // dE
	dflux[p] += -4.*pow(temp[p][t][0],3.); // dF = -4*T'^3
	Tsum[p] += temp[p][t][0]; // sum of surf T' over a rotation
      }

      if (converged) {lastrot=true; cout << "\ndone!" << endl;}

      // check to see if temperature converged/energy conserved after each rotation
      if ((t+1) % ntrot == 0 && !lastrot){
	for (int p=0;p<nlat;p++){
	  dtemp[p]=denergy[p]/dflux[p]; // dT' = dE/dF
	  avgT[p]=Tsum[p]/ntrot; // avg surf temp for each lat

	  // reset avgT & avgT0 for next time
	  if (abs(avgT[p]-avgT0[p]) < 1e-4 || abs(dtemp[p]) < 1e-5) {
	    avgT0[p]=avgT[p]; avgT[p]=0.0; latconv += 1;}
	  else {avgT0[p]=avgT[p]; avgT[p]=0.0; latconv = 0;}
	}

	// print to screen the # of rotations		
	cout << "\rrotation # " << (t+1)/ntrot << flush;

	// if converged solution then exit loop
	if (latconv >= nlat) converged=true;

	// if not, reset temperature guess
	else for (int p=0;p<nlat;p++){
	    Temp0[p]=Tsum[p]/ntrot-dtemp[p];
	    for (int j=0;j<nz;j++) temp[p][t][j]=temp[p][t][j]+
				     (Temp0[p]-temp[p][t][nz-1]);
	    denergy[p]=0.;
	  } dflux=denergy; Tsum=denergy;
      }

      if (t == ntp-1) cout << "\nerror: did not converge!" << endl; t++;
} while (t != ntp && !lastrot); // continue if not converged

 vector<float> surf(nlat*ntrot);
 for (int i=0;i<nlat;i++)
   for (int j=0;j<ntrot;j++) surf[i*ntrot+j] = temp[i][t-ntrot+j][0];
 
 cnpy::npz_save("surftemp_g.npz","smooth",&surf[0],{(unsigned long)nlat,(unsigned long)ntrot},"w");
// cnpy::npy_save("arr1.npy",&surf[0],{(unsigned long)nlat,(unsigned long)ntrot},"w");

// write surface temperatures to file
  ofstream surftemp;
  surftemp.open("surftemptpm_g.tbl");

  for (int phi=0;phi<nlat;phi++){
    for (int pi=0;pi<ntrot;pi++) surftemp << temp[phi][pi][0] << ' ';
    surftemp << endl;
  }
  surftemp.close();

//// write temperature profiles to file
//  ofstream deeptemp;
//  deeptemp.open("deeptemptpm_g.tbl");
//
//  for (int pi=t-ntrot;pi<t;pi++){
//    for (int x=0;x<nz;x++) deeptemp << temp[2][pi][x] <<' ';
//    deeptemp << endl;
//  }
//
//  deeptemp.close();

  return 0;
}
