// This header file contains procedures that read input for different thermophysical 
// models and converts to model variables.
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

void readintpm_orb(float &A, float &emis,
		   float &obliq, float &solst,
		   int &nxstep, int &ntrot, int &norb, int &nlat,
		   float &ti, float &dffsvty){

  ifstream infile("inputtpm_orb.txt");
  string line;

  while (getline(infile,line)){

    istringstream in(line);
    string var;
    in >> var;

//    if (var == "aorb:"){
//      in >> aorb; // heliocentric distance
//      cout <<"semimajor axis = "<< aorb <<" AU"<< endl;
//    }
//    if (var == "eorb:"){
//      in >> eorb; // eccentricity
//      cout <<"eccentricity   = "<< eorb << endl;
//    }
    if (var == "A:"){
      in >> A; // Bond Albedo
      cout <<"bond albedo = "<< A << endl;
    }
    if (var == "emis:"){
      in >> emis; // emissivity
    }
    if (var == "obliquity:"){
      in >> obliq; // spin obliquity
      cout <<"obliquity   = "<< obliq <<" degrees"<< endl;
    }
    if (var == "solstice:"){
      in >> solst; // position of northern summer solstice (rel. to perihelion)
      if (solst < 0) cout <<"solstice    = "<< -1*solst <<" degrees before q"<< endl;
      if (solst >= 0) cout <<"solstice    = "<< solst <<" degrees after q"<< endl;
    }
//    if (var == "period:"){
//      in >> perh; // rotation period (hr)
//      cout <<"rot. period = "<< perh <<" hr"<< endl;
//    }
    if (var == "nxstep:"){
      in >> nxstep; // # of depth steps per skindepth
    }
    if (var == "ntrot:"){
      in >> ntrot; // # of timesteps per rotation
    }
    if (var == "orbits:"){
      in >> norb; // # of orbits (revoltions)
    }
    if (var == "nlat:"){
      in >> nlat; // # of latitude bins
    }
    if (var == "ti:"){
      in >> ti ; // thermal inertia
    }
    if (var == "kappa:"){
      in >> dffsvty ; // thermal diffusivity
    }
  }
}

void readintpm_torb(float &aorb, float &eorb, float &A, float &emis,
		    float &obliq, float &solst, float &perh,
		    int &nxstep, int &ntrot, int &nlat,
		    float &ti, float &dffsvty, float &kfrac){

  ifstream infile("inputtpm_torb.txt");
  string line;

  while (getline(infile,line)){

    istringstream in(line);
    string var;
    in >> var;

    if (var == "aorb:"){
      in >> aorb; // heliocentric distance
      cout <<"semimajor axis = "<< aorb <<" AU"<< endl;
    }
    if (var == "eorb:"){
      in >> eorb; // eccentricity
      cout <<"eccentricity   = "<< eorb << endl;
    }
    if (var == "A:"){
      in >> A; // Bond Albedo
      cout <<"bond albedo = "<< A << endl;
    }
    if (var == "emis:"){
      in >> emis; // emissivity
    }
    if (var == "obliquity:"){
      in >> obliq; // spin obliquity
      cout <<"obliquity   = "<< obliq <<" degrees"<< endl;
    }
    if (var == "solstice:"){
      in >> solst; // position of northern summer solstice (rel. to perihelion)
      if (solst < 0) cout <<"solstice    = "<< -1*solst <<" degrees before q"<< endl;
      if (solst >= 0) cout <<"solstice    = "<< solst <<" degrees after q"<< endl;
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
    if (var == "nlat:"){
      in >> nlat; // # of latitude bins
    }
    if (var == "ti:"){
      in >> ti ; // thermal inertia
    }
    if (var == "kappa:"){
      in >> dffsvty ; // thermal diffusivity
    }
    if (var == "chi:"){
      in >> kfrac ; //
    }
  }
}
