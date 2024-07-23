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

void cranknicol(const vector<float> &temp_t0, const float &kappa,
		const float &stepsize, vector<float> &temp_t1){
  // input vectors: temp_t0 - nx size of temperature at time n0
  //                kappa   - nx size of thermal diffusivity
  //                temp_t1 - nx size of temperature at time t1 (w/top boundary)
  //
  // input param:   stepsize - 2*dx*dx/dt
  //

  int nx=temp_t0.size(); // number of depth steps
  int i,k;

  // set up A matrix as 3 vectors: Aa, Ab and Ac
  vector<float> Aa(nx); // Aa is diagonal
  vector<float> Ab(nx,0.0); // Ab is superdiagonal
  vector<float> Ac(nx,0.0); // Ac is subdiagonal

  for (i=1; i<(nx-1); i++){
    Aa[i] = stepsize+2.*kappa;
    if (i == nx-2) Aa[i] = stepsize+kappa; // for dt/dx = 0
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
      temp_t0[k]*(stepsize-2.*kappa)+
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

void cranknicol_g(const vector<float> &temp_t0, const float &stepsize, 		 
		  vector<float> &temp_t1){					 
  // input vectors: temp_t0 - nx size of temperature at time n0			 
  //                temp_t1 - nx size of temperature at time t1 (w/top boundary) 
  //										 
  // input param:   stepsize - 2*dx*dx/theta/dt					 
  //										 
  int nx=temp_t0.size(); // number of depth steps				 
  int i,k;									 
										 
  // set up A matrix as 3 vectors: Aa, Ab and Ac				 
  vector<float> Aa(nx); // Aa is diagonal					 
  vector<float> Ab(nx,0.0); // Ab is superdiagonal				 
  vector<float> Ac(nx,0.0); // Ac is subdiagonal				 
										 
  for (i=1; i<(nx-1); i++){							 
    Aa[i] = stepsize+2;								 
    if (i == nx-2) Aa[i] = stepsize+1; // for dt/dx = 0				 
    if (i < nx-2) Ab[i] = -1; // Ab[nx-2] = 0					 
    if (i > 1) Ac[i] = -1;    // Ac[0] = 0					 
  }										 
  Aa.erase (Aa.begin()); Aa.pop_back(); // trim vectors to size nx-2		 
  Ab.erase (Ab.begin()); Ab.pop_back();						 
  Ac.erase (Ac.begin()); Ac.pop_back();						 
										 
  // set up b vector for tridiag decomposition solver				 
  vector<float> bvec(nx);							 
										 
  for (k=1; k<(nx-1); k++){							 
    bvec[k]=temp_t0[k-1]+temp_t0[k]*(stepsize-2)+temp_t0[k+1];			 
										 
    // upper boundary condition for next time step				 
    if (k == 1) {bvec[1]=bvec[1]+temp_t1[0];}					 
  }										 
  bvec.erase (bvec.begin()); bvec.pop_back(); // trim bvec to size nx-2		 
										 
  // solve tridiagonal matrix							 
  tridiagsolv(Aa,Ab,Ac,bvec,temp_t1);						 
										 
  // impliment lower boundary condition dt/dx = 0				 
  temp_t1[nx-1]=temp_t1[nx-2];							 
}                                                                                

void cranknicol_d(const vector<float> &temp_t0, const vector<float> &kappa,
		  const float &stepsize, vector<float> &temp_t1){
// input vectors: temp_t0 - nx size of temperature at time n0
//                kappa   - nx size of thermal diffusivity
//                temp_t1 - nx size of temperature at time t1 (w/top boundary)
//
// input param:   stepsize - 2*dx*dx/dt
//
  int nx=temp_t0.size(); // number of depth steps
  int i,k;

// set up A matrix as 3 vectors: Aa, Ab and Ac
  vector<float> Aa(nx); // Aa is diagonal
  vector<float> Ab(nx,0.0); // Ab is superdiagonal
  vector<float> Ac(nx,0.0); // Ac is subdiagonal

  for (i=1; i<(nx-1); i++){
    Aa[i] = stepsize+kappa[i]+0.5*(kappa[i-1]+kappa[i+1]);
    if (i == nx-2) Aa[i] = stepsize+0.5*(kappa[i-1]+kappa[i]); // for dt/dx = 0
    if (i < nx-2) Ab[i] = -0.5*(kappa[i]+kappa[i+1]); // Ab[nx-2] = 0
    if (i > 1) Ac[i] = -0.5*(kappa[i-1]+kappa[i]);    // Ac[0] = 0
  }
  Aa.erase (Aa.begin()); Aa.pop_back(); // trim vectors to size nx-2
  Ab.erase (Ab.begin()); Ab.pop_back();
  Ac.erase (Ac.begin()); Ac.pop_back();

// set up b vector for tridiag decomposition solver
  vector<float> bvec(nx);

  for (k=1; k<(nx-1); k++){
    bvec[k]=temp_t0[k-1]*0.5*(kappa[k]+kappa[k-1])+
     temp_t0[k]*(stepsize-0.5*(kappa[k]+kappa[k-1])-0.5*(kappa[k]+kappa[k+1]))+
     temp_t0[k+1]*0.5*(kappa[k]+kappa[k+1]);

  // upper boundary condition for next time step
    if (k == 1) {bvec[1]=bvec[1]+temp_t1[0]*0.5*(kappa[1]+kappa[0]);}
  }
  bvec.erase (bvec.begin()); bvec.pop_back(); // trim bvec to size nx-2

// solve tridiagonal matrix
  tridiagsolv(Aa,Ab,Ac,bvec,temp_t1);

// impliment lower boundary condition dt/dx = 0
  temp_t1[nx-1]=temp_t1[nx-2];
}

//void cranknicol_t(const vector<float> &temp_t0, const float &stepsize,
//		  vector<float> &temp_t1){
//  // input vectors: temp_t0 - nx size of temperature at time t0
//  //                temp_t1 - nx size of temperature at time t1 (w/top boundary)
//  //
//  // input param:   stepsize - 2*dx*dx/dt
//  //
//  int nx=temp_t0.size(); // number of depth steps
//  int i,k;
//
//  // set up A matrix as 3 vectors: Aa, Ab and Ac
//  vector<float> Aa(nx); // Aa is diagonal
//  vector<float> Ab(nx,0.0); // Ab is superdiagonal
//  vector<float> Ac(nx,0.0); // Ac is subdiagonal
//
//  for (i=1; i<(nx-1); i++){
//    Aa[i] = stepsize-2*tdpndt(temp_t0[i]);
//    if (i == nx-2) Aa[i] = stepsize-tdpndt(temp_t0[i]); // for dt/dx = 0
//    if (i < nx-2) Ab[i] = tdpndt(temp_t0[i+1]); // Ab[nx-2] = 0
//    if (i > 1) Ac[i] = tdpndt(temp_t0[i-1]);    // Ac[0] = 0
//  }
//  Aa.erase (Aa.begin()); Aa.pop_back(); // trim vectors to size nx-2
//  Ab.erase (Ab.begin()); Ab.pop_back();
//  Ac.erase (Ac.begin()); Ac.pop_back();
//
//  // set up b vector for tridiag decomposition solver
//  vector<float> bvec(nx);
//
//  for (k=1; k<(nx-1); k++){
//    bvec[k]=temp_t0[k-1]*tdpndt(temp_t0[k-1])+
//      temp_t0[k]*(stepsize-2*(tdpndt(temp_t0[k]))+
//      temp_t0[k+1]*tdpndt(temp_t0[k-1]));
//
//    // upper boundary condition for next time step
//    if (k == 1) {bvec[1]=bvec[1]+temp_t1[0]*(tdpndt(temp_t0[0]));}
//  }
//  bvec.erase (bvec.begin()); bvec.pop_back(); // trim bvec to size nx-2
//
//  // solve tridiagonal matrix
//  tridiagsolv(Aa,Ab,Ac,bvec,temp_t1);
//
//  // impliment lower boundary condition dt/dx = 0
//  temp_t1[nx-1]=temp_t1[nx-2];
//}

void cranknicol_dt(vector<float> &temp_t0, vector<float> A,
                   vector<float> B, float &stepsize,
		   vector<float> &temp_t1, vector<float> por){
  // input vectors: temp_t0 - nx size of temperature at time n0
  //                kappa   - nx size of thermal diffusivity
  //                temp_t1 - nx size of temperature at time t1 (w/top boundary)
  //
  // input param:   stepsize - 2*dx*dx/dt
  //
  int nx=temp_t0.size(); // number of depth steps
  int i,k;

  // set up A matrix as 3 vectors: Aa, Ab and Ac
  vector<float> Aa(nx); // Aa is diagonal
  vector<float> Ab(nx,0.0); // Ab is superdiagonal
  vector<float> Ac(nx,0.0); // Ac is subdiagonal

  for (i=1; i<(nx-1); i++){
    float kappai_1 = tdpndt(temp_t0[i-1],A[i-1],B[i-1],por[i-1]);
    float kappai   = tdpndt(temp_t0[i],A[i],B[i],por[i]);
    float kappai1 = tdpndt(temp_t0[i+1],A[i+1],B[i+1],por[i+1]);

    Aa[i] = stepsize+kappai+0.5*(kappai_1+kappai1);
    if (i == nx-2) Aa[i] = stepsize+0.5*(kappai_1+kappai); // for dt/dx = 0
    if (i < nx-2) Ab[i] = -0.5*(kappai+kappai1); // Ab[nx-2] = 0
    if (i > 1) Ac[i] = -0.5*(kappai_1+kappai);    // Ac[0] = 0
  }
  Aa.erase (Aa.begin()); Aa.pop_back(); // trim vectors to size nx-2
  Ab.erase (Ab.begin()); Ab.pop_back();
  Ac.erase (Ac.begin()); Ac.pop_back();

  // set up b vector for tridiag decomposition solver
  vector<float> bvec(nx);

  for (k=1; k<(nx-1); k++){
    float kappak_1 = tdpndt(temp_t0[k-1],A[k-1],B[k-1],por[k-1]);
    float kappak = tdpndt(temp_t0[k],A[k],B[k],por[k]);
    float kappak1 = tdpndt(temp_t0[k+1],A[k+1],B[k+1],por[k+1]);

    bvec[k]=temp_t0[k-1]*0.5*(kappak+kappak_1)+
      temp_t0[k]*(stepsize-0.5*(kappak+kappak_1)-0.5*(kappak+kappak1))+
      temp_t0[k+1]*0.5*(kappak+kappak1);

    // upper boundary condition for next time step
    if (k == 1) {bvec[1]=bvec[1]+temp_t1[0]*0.5*
	(tdpndt(temp_t0[1],A[1],B[1],por[1])+
	 tdpndt(temp_t0[0],A[0],B[0],por[0]));}
  }
  bvec.erase (bvec.begin()); bvec.pop_back(); // trim bvec to size nx-2

  // solve tridiagonal matrix
  tridiagsolv(Aa,Ab,Ac,bvec,temp_t1);

  // impliment lower boundary condition dt/dx = 0
  temp_t1[nx-1]=temp_t1[nx-2];
}
