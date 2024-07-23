// TPM for a rough, spherical asteroid with output surface temperatures for
// various input parameters. The time-dependent heat diffusion equation is solved for
// an array of facets on the object.
//
// Best used with IDL wrapper asteroidtpm.pro, which writes the input variables to a
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
#include "cratergeom.h"

#include "boost/multi_array.hpp"


#define PI   3.1415927
#define dtor 0.0174533 // degrees to radians

using namespace std;

int main(){

// input parameters

  float A;      // bond albedo
  float theta;  // thermal parameter
  float sslat0; // sub-solar lat
  float sslon0; // sub-solar lat
  int   nring;  // number of crater rings
  float cr_angle; // crater opening half-angle (degrees)

  int nzstep;    // # of depth steps per skindepth
  int ntrot;     // # of timesteps per rotation
  int nrot;      // maximum # of rotations
  int nlat;      // # of latitude bins

  readintpm_gr(A,theta,sslat0,sslon0,nring,cr_angle,nzstep,ntrot,nrot,nlat);

// define constant values

  float sslat=sslat0*dtor; float sslon=sslon0*dtor; // convert to radians

  float dtp=(float) 2*PI/ntrot/theta;  // time step (a step around the circumference)
  int   ntp=ntrot*nrot+1; // total # of time steps

  float cr_param=cr_angle*dtor; // crater parameter
  int   ncelem=0; for (int i=nring;i>0;i--){ncelem += i*4;} // number of crater elements
  float radtot; // for thermally scattered energy

  float fij=(1.-cos(cr_param))/2./ncelem;
  float albfac=A/(1.-A*cr_param/PI);

  float dz=(float) 1/nzstep; // depth step = skin depth/# of steps per skindepth
  int   nz=10*nzstep+1;   // # of depth steps (for 10 skin depths)

  float stepsize=2.*dz*dz/theta/dtp; // used in cranknicol
  float fconv=dtp/dz; // convert flux (q) to dT' at surface
  float perz=theta/dz; // thermal parameter/depth step at surface

// create vectors
  typedef boost::multi_array<float, 2> array_2d; typedef array_2d::index index;
  typedef boost::multi_array<float, 3> array_3d; typedef array_3d::index index;
  typedef boost::multi_array<float, 4> array_4d; typedef array_4d::index index;

  // temperatures
  array_4d temp(boost::extents[nlat][ntp][ncelem][nz]); // temp[lat][time][depth][crelem]
  array_3d temp_out(boost::extents[nlat][ntrot][ncelem]); // for output

  // array of initial temperatures over latitude bins
  array_2d Temp0(boost::extents[nlat][ncelem]);

  vector<float> lat(nlat); // latitude
  vector<float> lon(ntrot);// longitude
  vector<float> zarr(nz);  // depth array

  vector<double> nu(ncelem); // angle to element from bottom of crater
  vector<double> mu(ncelem); // azimuthal angle to element from north-line

  // cosine of normal vector with sun
  array_3d cr_cosincelem(boost::extents[nlat][ntrot][ncelem]);
  // incident energy flux [lat][time][cr_elem]
  array_3d cr_fluxfrac(boost::extents[nlat][ntrot][ncelem]);

  vector<float> temp_t0(nz); vector<float> temp_t1(nz); // for cranknicol

  // used to find converged solution
  array_2d denergy(boost::extents[nlat][ncelem]);
  array_2d dflux(boost::extents[nlat][ncelem]);
  array_2d dtemp(boost::extents[nlat][ncelem]);
  array_2d Tsum(boost::extents[nlat][ncelem]);
  array_2d avgT(boost::extents[nlat][ncelem]);
  array_2d avgT0(boost::extents[nlat][ncelem]);

  // set Temp0 to empirically derived Tdeep (Spencer et al. 1989; fig 2.)
  float ftss; if (theta > 4.) ftss=0.75; else ftss=log10(theta)*.08+.7;

// populate vectors
  make_crater(nring,ncelem,cr_param,nu,mu); // calculate crater geometry

for (int p=0;p<nlat;p++){ // latitude loop
  lat[p]=p*PI/(nlat-1)-PI/2.; // array of latitude bins

  for (int t=0;t<ntrot;t++){ // time (longitude) loop
    if (p == 0) lon[t]=2.*PI*t/ntrot; // make longitude array
    float sumcosinc=0.; // for relfected flux within craters

    for (int i=0;i<ncelem;i++){ // crater geometries for each facet
      cr_cosincelem[p][t][i]=crater_geom(cr_param,sslat,sslon,lat[p],lon[t],nu[i],mu[i]);
      cr_fluxfrac[p][t][i]=cr_cosincelem[p][t][i]; // direct solar energy
      sumcosinc += cr_cosincelem[p][t][i];
    }

    for (int i=0;i<ncelem;i++){ // crater element loop
      // reflected solar energy
      cr_fluxfrac[p][t][i] += fij*albfac*(sumcosinc-cr_cosincelem[p][t][i]);

      if (t == 0){ // populate temperature array with initial guess

	// for permanently shadowed elements
       	if (cr_cosincelem[p][t][i] == 0.){
	  if ((sslat > 0. && -lat[p] > PI/2-sslat) || (sslat < 0. && lat[p] > PI/2+sslat)
	   || (sslat == 0. && (p == 0 || p == nlat-1))) Temp0[p][i]=5e-4; // near poles
	  else // in partially-lit craters
	    Temp0[p][i]=ftss*pow(cos(lat[p]-sslat),.25)*cos(cr_param);
	} else Temp0[p][i]=((1.-ftss)*sin(lat[p])+ftss)*pow(cr_fluxfrac[p][t][i],.25);

	// ensure finite temperatures
       	if (Temp0[p][i] != Temp0[p][i]) Temp0[p][i]=1e-3;

	for (int k=0;k<nz;k++){ // depth loop
	  zarr[k]=k*dz; // depth array
	  temp[p][t][i][k]=Temp0[p][i]; // initial temperatures at depth (constant)
	}
      }
    }
  }
}
int t=1; avgT0=Temp0; bool converged=false;
 int elemconv=0; vector<int> latconv(nlat,0); int sumlatconv=0;

// begin rotations
 cout <<"begin"<< endl;
do { // time loop
    for (int p=0;p<nlat;p++){ // latitude loop
      if (latconv[p] == 0) for (int i=0;i<ncelem;i++){ // crater element loop

	// surface temperature gradient term = theta*(dT'/dz)
	float tempgrad=perz*(temp[p][t-1][i][0]-temp[p][t-1][i][1]);

	// energy radiated term = T'^4
	float radiated=pow(temp[p][t-1][i][0],4.);

	// thermally scattered energy within the crater
	if (i == 0) {radtot = 0.;
	  for (int j=0;j<ncelem;j++) radtot += pow(temp[p][t-1][j][0],4.);}
	float thermscat=fij*(radtot-radiated);

	// energy balance to get surf temp from previous time step
	temp_t1[0]=temp[p][t-1][i][0]+fconv*(cr_fluxfrac[p][t % ntrot][i]+thermscat
					     -tempgrad-radiated);

	// solve for temperatures at depth
	for (int k=0;k<nz;k++) temp_t0[k]=temp[p][t-1][i][k];
	cranknicol_g(temp_t0, stepsize, temp_t1);
	for (int k=0;k<nz;k++) temp[p][t][i][k]=temp_t1[k];

	// get values for energy balance
	denergy[p][i] += cr_fluxfrac[p][t % ntrot][i]+thermscat-radiated; // dE
	dflux[p][i] += -4.*pow(temp[p][t][i][0],3.); // dF = -4*T'^3
	Tsum[p][i] += temp[p][t][i][0]; // sum of surf T' over a rotation
      }
    }

    // check to see if temperature converged/energy conserved after each rotation
    if ((t+1) % ntrot == 0 && !converged){
      for (int p=0;p<nlat;p++){
	if (latconv[p] == 0){
	  for (int i=0;i<ncelem;i++){
	    dtemp[p][i]=denergy[p][i]/dflux[p][i]; // dT' = dE/dF
	    avgT[p][i]=Tsum[p][i]/ntrot; // avg surf temp for each lat

	    // reset avgT & avgT0 for next step
	    if (abs(avgT[p][i]-avgT0[p][i]) < 1e-4 || abs(dtemp[p][i]) < 1e-5
		|| avgT[p][i] < 5e-4) {
	      avgT0[p][i]=avgT[p][i]; avgT[p][i]=0.0; elemconv += 1;}
	    else {avgT0[p][i]=avgT[p][i]; avgT[p][i]=0.0; elemconv = 0;}
	  }
	  if (elemconv >= ncelem) {latconv[p]=1; sumlatconv += 1;
	    for (int n=0;n<ntrot;n++) for (int i=0;i<ncelem;i++) 
					temp_out[p][n][i]=temp[p][t+1-ntrot+n][i][0];}
	}
      }

      // print to screen the # of rotations
      //cout << "\rrotation # " << (t+1)/ntrot << flush;
      cout << "\rrotation #" << (t+1)/ntrot << " " << sumlatconv <<"/"<< nlat << flush;

      // if converged solution then exit loop
      if (sumlatconv >= nlat) {converged=true; cout << "\ndone!" << endl;}

      // if not, reset temperature guess
      else for (int p=0;p<nlat;p++){ if (latconv[p] == 0)
	  for (int i=0;i<ncelem;i++){
	    Temp0[p][i]=Tsum[p][i]/ntrot-dtemp[p][i];

	    for (int k=0;k<nz;k++) temp[p][t][i][k]=temp[p][t][i][k]*
				     (Temp0[p][i]/temp[p][t][i][nz-1]);
	    denergy[p][i]=0.;
	  }
	} dflux=denergy; Tsum=denergy;
    }

    if (t == ntp-1) cout << "\nerror: did not converge!" << endl; t++;
} while (t != ntp && !converged); // continue if not converged

// write surface temperatures of each crater element to separate files
for (int i=0;i<ncelem;i++){
  ofstream surftemp;
  surftemp.open("surftemptpm_gr_crel"+std::to_string(i+1)+".tbl");

  for (int phi=0;phi<nlat;phi++){
    for (int pi=0;pi<ntrot;pi++) surftemp << temp_out[phi][pi][i] << ' ';
    surftemp << endl;
  }
  surftemp.close();
}

  return 0;
}
