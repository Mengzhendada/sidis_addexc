#ifndef Exc_structure_function_H
#define Exc_structure_function_H
 
#include <cmath>
#include <string>
#include <vector>
#include "sidis/kinematics.hpp"
#include "sidis/vector.hpp"
using namespace sidis::kin;
using namespace sidis;
namespace sidis{namespace exc{
const double m_n = 0.9395612928;
const double ALPHA = 7.2973525664e-3L;
double Get_thetacm(double W, double Q2, double t);
struct DataPoint {
	double Q2,W, theta;
	std::vector<double> values;
};
void printDataPoint(DataPoint dp);
double linearInterpolated(double x,double x0,double x1,double y0,double y1);
//std::vector<DataPoint> readDataFromFile(const std::string& filename);
std::vector<double> Get_exc_sf(double W, double Q2, double t);
struct EXC_A{
  double A1r;
  double A1i;
  double A2r;
  double A2i;
  double A3r;
  double A3i;
  double A4r;
  double A4i;
  double A5r;
  double A5i;
  double A6r;
  double A6i;
  EXC_A(double W,double Q2,double t);
};
struct EXC_SF_F{
  double r1;
  double r2;
  double f1r;
  double f1i;
  double f2r;
  double f2i;
  double f3r;
  double f3i;
  double f4r;
  double f4i;
  double f5r;
  double f5i;
  double f6r;
  double f6i;
  EXC_SF_F(KinematicsRad kin);
};
struct EXC_SF_combine{
  double f11;
  double f12r;
  double f12i;
  double f13r;
  double f13i;
  double f14r;
  double f14i;
  double f15r;
  double f15i;
  double f16r;
  double f16i;
  double f22;
  double f23r;
  double f23i;
  double f24r;
  double f24i;
  double f25r;
  double f25i;
  double f26r;
  double f26i;
  double f33;
  double f34r;
  double f34i;
  double f35r;
  double f35i;
  double f36r;
  double f36i;
  double f44;
  double f45r;
  double f45i;
  double f46r;
  double f46i;
  double f55;
  double f56r;
  double f56i;
  double f66;
  EXC_SF_combine(EXC_SF_F exc_sf_f);

};


struct EXC_SF{

  //C.23 exclusive structure functions in orthogonal basis
  //p is plus, m is minus, o is over, 'ex' in C.23 is ignored in variable name
  //_000 means no target polarization _123 corresponds to eta1, eta2, eta3
  double r1;
  double r2;
  double r3;
  double r4;
  double r5;
  double rex;
  double H00pH22_000;
  double H11mH22opt2_000;
  double H22_000;
  double H01r_000;
  double H01i_000;
  double H00pH22_010;
  double H11mH22opt2_010;
  double H22_010;
  double H01r_010;
  double H01i_010;
  double H02r_100;
  double H02i_100;
  double H12r_100;
  double H12i_100;
  double H02r_001;
  double H02i_001;
  double H12r_001;
  double H12i_001;
  explicit EXC_SF(EXC_SF_combine exc_sf_com, KinematicsRad kin);
};
struct EXCUU{
  //equations in 42, the generalized exclusive structure functions of exclusive processes that contribute to the cross section of exclusive radiative tail
  double H1_000;
  double H2_000;
  double H3_000;
  double H4_000;
  explicit EXCUU(KinematicsRad kin);
};
struct EXCLU{
  double H5_000;
  explicit EXCLU(KinematicsRad kin);
};
struct EXCUT{
	math::Vec3 H1_010;
	math::Vec3 H2_010;
	math::Vec3 H3_010;
	math::Vec3 H4_010;
	math::Vec3 H6_100;
	math::Vec3 H8_100;
	explicit EXCUT(KinematicsRad kin);
};
struct EXCLT{
	math::Vec3 H5_010;
	math::Vec3 H7_100;
	math::Vec3 H9_100;
	explicit EXCLT(KinematicsRad kin);
};
struct EXCUL{
	math::Vec3 H6_001;
	math::Vec3 H8_001;
	explicit EXCUL(KinematicsRad kin);
};
struct EXCLL{
	math::Vec3 H7_001;	
	math::Vec3 H9_001;
	explicit EXCLL(KinematicsRad kin);
};

}}//namespace 
#endif
