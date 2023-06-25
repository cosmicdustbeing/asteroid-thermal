// This header includes programs that increments a temperature-depth array by one 
// time-step using the Crank-Nicolson scheme for solving the heat diffusion equation 
//
// cranknicol: diffusivity in parameterized semi-infinite half-space
// cranknicol_d: depth dependent diffusivity in a semi-infinite half-space
// cranknicol_t: temperature & depth dependent diffusivity in a semi-infinite half-space
//
#include <vector>
#include "tdpndt.h"

using namespace std;

void tridiagsolv(const vector<float> &Aa, const vector<float> &Ab,
                 const vector<float> &Ac, const vector<float> &bvec,
		 vector<float> &x){

  int j,n=Aa.size();
  float bet;
  vector<float> gam(n);

  if (Aa[0] == 0.0) throw("Error 1 in trisolv");
  x[1]=bvec[0]/(bet=Aa[0]);

// forward substitution
  for (j=1; j<n; j++){
    gam[j]=Ab[j-1]/bet;
    bet=Aa[j]-Ac[j]*gam[j];
    if (bet == 0.0) throw("Error 2 in trisolv");
    x[j+1]=(bvec[j]-Ac[j]*x[j])/bet;
  }
// backsubstitution
  for (j=(n-2); j>=0; j--)
    x[j+1] -= gam[j+1]*x[j+2];
}

void cranknicol_orb(const vector<float> &temp_t0, const float &kappa,
		    const vector<float> &stepsize, vector<float> &temp_t1){
// input vectors: temp_t0 - nx size of temperature at time n0
//                stepsize - nx size 2*dx*dx/dt
//                temp_t1 - nx size of temperature at time t1 (w/top boundary)
//
// input param:   kappa - thermal diffusivity
//
  int nx=temp_t0.size(); // number of depth steps
  int i,k;

// set up A matrix as 3 vectors: Aa, Ab and Ac
  vector<float> Aa(nx); // Aa is diagonal
  vector<float> Ab(nx,0.0); // Ab is superdiagonal
  vector<float> Ac(nx,0.0); // Ac is subdiagonal

  for (i=1; i<(nx-1); i++){
    Aa[i] = stepsize[i-1]+2.*kappa;
    if (i == nx-2) Aa[i] = stepsize[i-1]+kappa; // for dt/dx = 0
    if (i < nx-2) Ab[i] = -1.*kappa; // Ab[nx-2] = 0
    if (i > 1) Ac[i] = -1.*kappa;    // Ac[0] = 0
  }
  Aa.erase (Aa.begin()); Aa.pop_back(); // trim vectors to size nx-2
  Ab.erase (Ab.begin()); Ab.pop_back();
  Ac.erase (Ac.begin()); Ac.pop_back();

// set up b vector for tridiag decomposition solver
  vector<float> bvec(nx);

  for (k=1; k<(nx-1); k++){
    bvec[k]=temp_t0[k-1]*kappa+
     temp_t0[k]*(stepsize[k-1]-2.*kappa)+
     temp_t0[k+1]*kappa;

  // upper boundary condition for next time step
    if (k == 1) {bvec[1]=bvec[1]+temp_t1[0]*kappa;}
  }
  bvec.erase (bvec.begin()); bvec.pop_back(); // trim bvec to size nx-2

// solve tridiagonal matrix
  tridiagsolv(Aa,Ab,Ac,bvec,temp_t1);

// impliment lower boundary condition dt/dx = 0
  temp_t1[nx-1]=temp_t1[nx-2];
}

void cranknicol_torb(const vector<float> &temp_t0, const float &kappa,
		     const float &kfrac,
		     const vector<float> &stepsize, vector<float> &temp_t1){
  // input vectors: temp_t0 - nx size of temperature at time n0
  //                stepsize - nx size 2*dx*dx/dt
  //                temp_t1 - nx size of temperature at time t1 (w/top boundary)
  //
  // input param:   kappa - thermal diffusivity
  //

  int nx=temp_t0.size(); // number of depth steps
  int i,k;

  // set up A matrix as 3 vectors: Aa, Ab and Ac
  vector<float> Aa(nx); // Aa is diagonal
  vector<float> Ab(nx,0.0); // Ab is superdiagonal
  vector<float> Ac(nx,0.0); // Ac is subdiagonal

   for (i=1; i<(nx-1); i++){
    Aa[i] = stepsize[i-1]+2.*kappa*(1+kfrac*pow(temp_t0[i]/350.,3.));
    if (i == nx-2) Aa[i] = stepsize[i-1]+
		     kappa*(1+kfrac*pow(temp_t0[i]/350.,3.));          // for dt/dx = 0
    if (i < nx-2) Ab[i] = -1.*kappa*(1+kfrac*pow(temp_t0[i]/350.,3.)); // Ab[nx-2] = 0
    if (i > 1) Ac[i] = -1.*kappa*(1+kfrac*pow(temp_t0[i]/350.,3.));    // Ac[0] = 0
  }

  Aa.erase (Aa.begin()); Aa.pop_back(); // trim vectors to size nx-2
  Ab.erase (Ab.begin()); Ab.pop_back();
  Ac.erase (Ac.begin()); Ac.pop_back();

  // set up b vector for tridiag decomposition solver
  vector<float> bvec(nx);

  for (k=1; k<(nx-1); k++){
    bvec[k]=temp_t0[k-1]*kappa*(1+kfrac*pow(temp_t0[k-1]/350.,3.))+
      temp_t0[k]*(stepsize[k-1]-2.*kappa*(1+kfrac*pow(temp_t0[k]/350.,3.)))+
      temp_t0[k+1]*kappa*(1+kfrac*pow(temp_t0[k+1]/350.,3.));

    // upper boundary condition for next time step
    if (k == 1) {bvec[1]=bvec[1]+temp_t1[0]*kappa*(1+kfrac*pow(temp_t0[1]/350.,3.));}
  }
  bvec.erase (bvec.begin()); bvec.pop_back(); // trim bvec to size nx-2

  // solve tridiagonal matrix
  tridiagsolv(Aa,Ab,Ac,bvec,temp_t1);

  // impliment lower boundary condition dt/dx = 0
  temp_t1[nx-1]=temp_t1[nx-2];
}
