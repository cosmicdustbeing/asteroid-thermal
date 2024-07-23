// This header contains procedures for calculating various angles used in various TPMs
#include <cmath>

using namespace std;

// transpose spherical coordinates (r theta phi) to cartesian (x y z)
void sph2cart(const float r, const float theta, const float phi, vector<float> &outvec){

  outvec[0] = r*cos(phi)*cos(theta);
  outvec[1] = r*cos(phi)*sin(theta);
  outvec[2] = r*sin(phi);
}

// transpose cartesian coordinates (x y z) to spherical (r theta phi)
void cart2sph(const float x, const float y, const float z, vector<float> &outvec){

  outvec[0] = sqrt(x*x + y*y + z*z);
  outvec[1] = atan(y/x);
  outvec[2] = asin(z/outvec[0]);
}

// calculate the cosine of the zenith angle
float coszang(const float &sublat, const float &sublon,
	      const float &lat, const float &lon, const float &sunang){

  float coszang;

  coszang=sin(sublat)*sin(lat)+
    cos(sublat)*cos(lat)*cos(abs(lon-sublon));

  if (abs(coszang) < sin(sunang/2.))
    coszang -= (coszang - sin(sunang/2.)*sqrt(1.-pow(coszang,2.)))/2.;

  if (coszang < 0.) coszang = 0.;

  return coszang;
}

void rotate_x(vector<float> &in, const float &angle, vector<float> &out){

  out[0] = cos(angle)*in[0] + sin(angle)*in[1];
  out[1] = 
  out[2] = in[2];

}

void rotate_z(vector<float> &in, const float &angle, vector<float> &out){

  out[0] = in[0];
  out[1] = cos(angle)*in[1] + sin(angle)*in[2];
  out[2] = cos(angle)*in[2] - sin(angle)*in[1];

}
