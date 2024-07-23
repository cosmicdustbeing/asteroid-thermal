// This header contains procedures for calculating geometry for spherical craters in
// the rough surface TPM
#include <cmath>
#include <vector>

#define PI   3.1415927
#define dtor 0.0174533 // degrees to radians

using namespace std;

float cr_coszang(const float &sublat, const float &sublon,
		  const float &lat, const float &lon){

  float coszang;

  coszang=sin(sublat)*sin(lat)+
    cos(sublat)*cos(lat)*cos(abs(lon-sublon));
//  if (abs(lat) == 1.5708) coszang = -coszang;

  return coszang;
}

float cr_sinzang(const float &coszang){

  return sqrt(1.-pow(coszang,2.));
}

float cr_cosphi(const float &sublat, const float &lat,
		 const float &coszang, const float &sinzang){

  float numer; float denom;

  numer=sin(sublat)-coszang*sin(lat);
  denom=sinzang*cos(lat);

  if (sinzang*cos(lat) == 0.) denom += 1;

  return numer/denom;
}

float cr_sinphi(const float &sublat, const float &sublon,
		 const float &lon, const float &sinzang){

  float numer; float denom;

  numer=sin(lon-sublon)*cos(sublat);
  if (sinzang == 0) denom=1.; else denom=sinzang;

  return numer/denom;
}

void cr_ang(const float &sublat, const float &sublon,
	    const float &lat, const float &lon,
	    float &coszang, float &sinzang, float &cosphi, float &sinphi){

  coszang=cr_coszang(sublat,sublon,lat,lon);
  sinzang=cr_sinzang(coszang);
  cosphi=cr_cosphi(sublat,lat,coszang,sinzang);
  sinphi=cr_sinphi(sublat,sublon,lon,sinzang);
}

float crater_geom(const float &crparam, const float &sublat, const float &sublon,
		   const float &lat, const float &lon, const double &nui,
		   const double &mui){

  float cosinc,sininc,cosphiz,sinphiz; // for angle calculation

  cr_ang(sublat,sublon,lat,lon,cosinc,sininc,cosphiz,sinphiz);

  double sinnui=sin(nui); double cosnui=cos(nui);
  double sinmunorth=sin(mui); double cosmunorth=cos(mui);

  double tau=sinnui;
  float z=sin(crparam); float z2=pow(z,2.);
  float x=cos(crparam); float x2=pow(x,2.);

  //  float taninc=sininc/cosinc;
  float taninc=tan(acos(cosinc));

  // convert munorth to mu, the azimuthal angle from element to sun-line
  double sinmu=sinmunorth*cosphiz-cosmunorth*sinphiz;
  double cosmu=cosmunorth*cosphiz+sinmunorth*sinphiz;

  // calculate y1 (length of shadow from center-line of crater)
  // and other variables involved in shadow calculation
  double y1=tau*cosmu;
  double zprime2=z2-pow((tau*sinmu),2.);
  double rprime2=x2+zprime2;

  // test if shadowed
  double a=sqrt(zprime2)+x*taninc;
  float cosinc2=pow(cosinc,2.);
  double y2s=rprime2/cosinc2-pow(a,2.); if (y2s < 0.) y2s=0.;
  double y2=cosinc2*(a-taninc*sqrt(y2s));

  float  shadow; if (y1 >= y2) shadow=1; else shadow=0;

  float coszelem=(cosinc*cosnui-sininc*sinnui*cosmu)*(1.-shadow);
  if (coszelem < 0.) coszelem=0.;
  if (coszelem == -0.) coszelem=0;
  if (cosinc < 1e-5) coszelem=0.;

  return coszelem;

}

void make_crater(const int &nring, const int &nelem, float &cr_param,
		 vector<double> &nui, vector<double> &mui){

  vector<int>    ninring(nring); // array w/# of elelemts in each ring
  vector<double> cosomega(nring); // cosines of the angles between the crater
                                  // normals and the dividing lines between rings
  vector<double> nu(nring); // angle to ring from bottom of crater

  double revtot=nelem; // for omega calculation

  for (int i=nring-1;i>=0;i--){ // loop over the rings (backward)
    ninring[i]=(i+1)*4.; // 4*i elements in ith ring

    if (i == nring-1){cosomega[i]=abs(cos(cr_param));
      if (cosomega[i] < 1e-4) cosomega[i] = 1e-4;}
    else {
      double num=ninring[i+1]+cosomega[i+1]*revtot;
      double num2=revtot+ninring[i+1];
      cosomega[i]=num/num2; // Eq. 8 in Emery et al. (1998)
    }
    revtot -= ninring[i]; // subtract # of elements in larger ring
  }

  int beg=0; int end=ninring[0]; // start/finish index for ring elements
  for (int i=0;i<nring;i++){ // loop over the rings (forward)
    if (i == 0){nu[i]=0.5*acos(cosomega[i]);}
    if (i != 0){nu[i]=0.5*(acos(cosomega[i-1])+acos(cosomega[i]));}
    for (int k=beg;k<end;k++){ // loop over elements in each ring
      nui[k]=nu[i]; // set element nu to ring nu
      mui[k]=(k-beg+0.5)*(2.*PI/ninring[i]); // mu angle for crater element
      if (mui[k] > PI) mui[k] -= 2.*PI;
    }
    if (i < (nring-1)){beg += ninring[i]; end += ninring[i+1];}
    if (i == (nring-1)){beg = nelem-ninring[i] ;end = nelem;}
  }
}
