// This header file contains procedures that read input for different thermophysical 
// models and converts to model variables.
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

void readintpm_g(float &theta, float &sslat0, float &sslon0,
		 int &nzstep, int &ntrot, int &nrot, int &nlat){

  ifstream infile("inputtpm_g.txt");
  string line;

  while (getline(infile,line)){

    istringstream in(line);
    string var;
    in >> var;

    if (var == "theta:"){
      in >> theta; // thermal parameter
      cout <<"theta  = "<< theta << endl;
    }
    if (var == "sslat:"){
      in >> sslat0; // sub-solar latitude
      cout <<"ss lat = "<< sslat0 <<" degrees"<< endl;
    }
    if (var == "sslon:"){
      in >> sslon0; // sub-solar longitude
      cout <<"ss lon = "<< sslon0 <<" degrees"<< endl;
    }
    if (var == "nzstep:"){
      in >> nzstep; // # of depth steps per skindepth
    }
    if (var == "ntrot:"){
      in >> ntrot; // # of timesteps per rotation
    }
    if (var == "nrot:"){
      in >> nrot; // max # of rotations
    }
    if (var == "nlat:"){
      in >> nlat; // # of latitude bins
    }
  }
}

void readintpm_gr(float &A, float &theta, float &sslat0, float &sslon0,
		  int &nring, float &cr_angle, int &nzstep, int &ntrot,
		  int &nrot, int &nlat){

  ifstream infile("inputtpm_gr.txt");
  string line;

  while (getline(infile,line)){

    istringstream in(line);
    string var;
    in >> var;

    if (var == "A:"){
      in >> A; // bond albedo
      cout <<"Albedo = "<< A << endl;
    }
    if (var == "theta:"){
      in >> theta; // thermal parameter
      cout <<"theta  = "<< theta << endl;
    }
    if (var == "sslat:"){
      in >> sslat0; // sub-solar latitude
      cout <<"ss lat = "<< sslat0 <<" degrees"<< endl;
    }
    if (var == "sslon:"){
      in >> sslon0; // sub-solar longitude
      cout <<"ss lon = "<< sslon0 <<" degrees"<< endl;
    }
    if (var == "ntrot:"){
      in >> ntrot; // # of timesteps per rotation
    }
    if (var == "nring:"){
      in >> nring; // number of crater rings
    }
    if (var == "cr_angle:"){
      in >> cr_angle; // crater half-opening angle
      cout <<"crater angle = "<< cr_angle <<" degrees"<< endl;
    }
    if (var == "nzstep:"){
      in >> nzstep; // # of depth steps per skindepth
    }
    if (var == "nrot:"){
      in >> nrot; // max # of rotations
    }
    if (var == "nlat:"){
      in >> nlat; // # of latitude bins
    }
  }
}

void readintpm_geq(float &sslat0, float &sslon0, int &ntrot, int &nlat){

  ifstream infile("inputtpm_geq.txt");
  string line;

  while (getline(infile,line)){

    istringstream in(line);
    string var;
    in >> var;

    if (var == "sslat:"){
      in >> sslat0; // sub-solar latitude
      cout <<"ss lat  = "<< sslat0 <<" degrees"<< endl;
    }
    if (var == "sslon:"){
      in >> sslon0; // sub-solar longitude
      cout <<"ss lon  = "<< sslon0 <<" degrees"<< endl;
    }
    if (var == "ntrot:"){
      in >> ntrot; // # of timesteps per rotation
    }
    if (var == "nlat:"){
      in >> nlat; // # of latitude bins
    }
  }
}

void readintpm_greq(float &A, float &sslat0, float &sslon0,
		    int &nring, float &cr_angle, int &ntrot, int &nlat){

  ifstream infile("inputtpm_greq.txt");
  string line;

  while (getline(infile,line)){

    istringstream in(line);
    string var;
    in >> var;

    if (var == "A:"){
      in >> A; // bond albedo
      cout <<"Albedo  = "<< A << endl;
    }
    if (var == "sslat:"){
      in >> sslat0; // sub-solar latitude
      cout <<"ss lat  = "<< sslat0 <<" degrees"<< endl;
    }
    if (var == "sslon:"){
      in >> sslon0; // sub-solar longitude
      cout <<"ss lon  = "<< sslon0 <<" degrees"<< endl;
    }
    if (var == "ntrot:"){
      in >> ntrot; // # of timesteps per rotation
    }
    if (var == "nlat:"){
      in >> nlat; // # of latitude bins
    }
    if (var == "nring:"){
      in >> nring; // number of crater rings
    }
    if (var == "cr_angle:"){
      in >> cr_angle; // crater half-opening angle
      cout <<"crater angle = "<< cr_angle <<" degrees"<< endl;
    }
  }
}

void readintpm_ginf(float &sslat0, float &sslon0, int &ntrot, int &nlat){

  ifstream infile("inputtpm_ginf.txt");
  string line;

  while (getline(infile,line)){

    istringstream in(line);
    string var;
    in >> var;

    if (var == "sslat:"){
      in >> sslat0; // sub-solar latitude                                                 
      cout <<"ss lat  = "<< sslat0 <<" degrees"<< endl;
    }
    if (var == "sslon:"){
      in >> sslon0; // sub-solar longitude                                                
      cout <<"ss lon  = "<< sslon0 <<" degrees"<< endl;
    }
    if (var == "ntrot:"){
      in >> ntrot; // # of timesteps per rotation
    }
    if (var == "nlat:"){
      in >> nlat; // # of latitude bins
    }
  }
}

void readintpm_grinf(float &A, float &sslat0, float &sslon0, int &nring,
		     float &cr_angle, int &ntrot, int &nlat){

  ifstream infile("inputtpm_grinf.txt");
  string line;

  while (getline(infile,line)){

    istringstream in(line);
    string var;
    in >> var;

    if (var == "A:"){
      in >> A; // bond albedo
      cout <<"Albedo = "<< A << endl;
    }
    if (var == "sslat:"){
      in >> sslat0; // sub-solar latitude
      cout <<"ss lat = "<< sslat0 <<" degrees"<< endl;
    }
    if (var == "sslon:"){
      in >> sslon0; // sub-solar longitude
      cout <<"ss lon = "<< sslon0 <<" degrees"<< endl;
    }
    if (var == "ntrot:"){
      in >> ntrot; // # of timesteps per rotation
    }
    if (var == "nring:"){
      in >> nring; // number of crater rings
    }
    if (var == "cr_angle:"){
      in >> cr_angle; // crater half-opening angle
      cout <<"crater angle = "<< cr_angle <<" degrees"<< endl;
    }
    if (var == "nlat:"){
      in >> nlat; // # of latitude bins
    }
  }
}

void readintpm_d(float &R, float &A, float &emis,
	         float &sslat0, float &sslon0, float &perh,
	         int &nxstep, int &ntrot, int &nrot, int &nlat,
	         float &ti, vector<float> &dffsvty, vector<float> &boundary){

  ifstream infile("inputtpm_d.txt");
  string line;

  while (getline(infile,line)){

    istringstream in(line);
    string var;
    in >> var;

    if (var == "R:"){
      in >> R; // heliocentric distance
      cout <<"helio dist  = "<< R <<" AU"<< endl;
    }
    if (var == "A:"){
      in >> A; // Bond Albedo
      cout <<"bond albedo = "<< A << endl;
    }
    if (var == "emis:"){
      in >> emis; // emissivity
    }
    if (var == "sslat:"){
      in >> sslat0; // sub-solar latitude
      cout <<"ss lat      = "<< sslat0 <<" degrees"<< endl;
    }
    if (var == "sslon:"){
      in >> sslon0; // sub-solar longitude
      cout <<"ss lon      = "<< sslon0 <<" degrees"<< endl;
    }
    if (var == "period:"){
      in >> perh; // rotation period (hr)
      cout <<"rot. period = "<< perh <<" hr"<< endl;
    }
    if (var == "nxstep:"){
      in >> nxstep; // # of depth steps per skindepth
    }
    if (var == "ntrot:"){
      in >> ntrot; // # of timesteps per rotation
    }
    if (var == "nrot:"){
      in >> nrot; // max # of rotations
    }
    if (var == "nlat:"){
      in >> nlat; // # of latitude bins
    }
    if (var == "ti:"){
      in >> ti ; // thermal inertia
    }
    if (var == "kappa:"){
      float in_kap;
      while (in >> in_kap) dffsvty.push_back(in_kap); // diffusivity array
    }
    if (var == "boundary:"){
      float in_bound;
      while (in >> in_bound)boundary.push_back(in_bound); // boundary points
    }
  }
}

void readintpm_dt(float &R, float &A, float &emis,
		  float &sslat0, float &sslon0, float &perh,
		  int &nxstep, int &ntrot, int &nrot, int &nlat,
		  float &ti, vector<float> &A0, vector<float> &B0,
		  vector<float> &por,vector<float> &boundary){

  ifstream infile("inputtpm_dt.txt");
  string line;

  while (getline(infile,line)){

    istringstream in(line);
    string var;
    in >> var;

    if (var == "R:"){
      in >> R; // heliocentric distance
      cout <<"helio dist  = "<< R <<" AU"<< endl;
    }
    if (var == "A:"){
      in >> A; // Bond Albedo
      cout <<"Bond albedo = "<< A << endl;
    }
    if (var == "emis:"){
      in >> emis; // emissivity
    }
    if (var == "sslat:"){
      in >> sslat0; // sub-solar latitude
      cout <<"ss lat      = "<< sslat0 <<" degrees"<< endl;
    }
    if (var == "sslon:"){
      in >> sslon0; // sub-solar longitude
      cout <<"ss lon      = "<< sslon0 <<" degrees"<< endl;
    }
    if (var == "period:"){
      in >> perh; // rotation period (hr)
      cout <<"rot. period = "<< perh <<" hr"<< endl;
    }
    if (var == "nxstep:"){
      in >> nxstep; // # of depth steps per skindepth
    }
    if (var == "ntrot:"){
      in >> ntrot; // # of timesteps per rotation
    }
    if (var == "nrot:"){
      in >> nrot; // max # of rotations
    }
    if (var == "nlat:"){
      in >> nlat; // # of latitude bins
    }
    if (var == "ti:"){
      in >> ti ; // thermal inertia
      cout <<"T.I.        = "<< ti <<" J/m^2/K/s^0.5"<< endl;
    }
    if (var == "Ax:"){
      float in_kap;
      while (in >> in_kap) A0.push_back(in_kap); // diffusivity array         
    }
    if (var == "Bx:"){
      float in_kap;
      while (in >> in_kap) B0.push_back(in_kap); // diffusivity array
    }
    if (var == "por:"){
      float inpor;
      while (in >> inpor) por.push_back(inpor); // porosity array
    }
//    if (var == "kappa:"){
//      float in_kap;
//      while (in >> in_kap)
//        vector<float> dffsvty; dffsvty.push_back(in_kap); // diffusivity array
//    }
    if (var == "boundary:"){
      float inbou;
      while (in >> inbou) boundary.push_back(inbou);
    }
  }
}


void readintpm_shape(float &R, float &A, float &emis, float &sslat0, float &sslon0, int &nzstep,
		     int &ntrot, int &nrot, float &ti, float &kappa){

  ifstream infile("inputtpm_shape.txt");
  string line;

  while (getline(infile,line)){

    istringstream in(line);
    string var;
    in >> var;

    if (var == "R:"){
      in >> R; // Heliocentric distance (au)
      cout <<"R           = "<< R << endl;
    }
    if (var == "A:"){
      in >> A; // Bond albedo
      cout <<"Bond albedo = "<< A << endl;
    }
    if (var == "emis:"){
      in >> emis; // emissivity
      cout <<"emissivity  = "<< emis << endl;
    }
    if (var == "sslat:"){
      in >> sslat0; // sub-solar latitude
      cout <<"SS Lat    = "<< sslat0 <<" degrees"<< endl;
    }
    if (var == "sslon:"){
      in >> sslon0; // sub-solar longitude
      cout <<"SS Lon    = "<< sslon0 <<" degrees"<< endl;
    }
    if (var == "nzstep:"){
      in >> nzstep; // # of depth steps per skindepth
    }
    if (var == "ntrot:"){
      in >> ntrot; // # of timesteps per rotation
    }
    if (var == "nrot:"){
      in >> nrot; // max # of rotations
    }
    if (var == "ti:"){
      in >> ti; // thermal inertia
      cout <<"therm inert  = "<< ti <<" (SI)"<< endl;
    }
    if (var == "kappa:"){
      in >> kappa; // thermal diffusivity
    }
  }
}
