// This header contains procedures for calculating various angles used in various TPMs
#include <cmath>

using namespace std;

double coszang(const double &sublat, const double &sublon,
	       const double &lat, const double &lon, double &sunang){

  double coszang;

  coszang=sin(sublat)*sin(lat)+
    cos(sublat)*cos(lat)*cos(abs(lon-sublon));

  if (abs(coszang) < sin(sunang/2.)) coszang += cos(sunang/2.)*coszang + 
				   sin(sunang/2.)*sqrt(1-pow(coszang,2.));

  if (coszang < 0.) coszang = 0.0;

  return coszang;
}
