#include "sidis/Exc_structure_function.hpp"
#include "sidis/constant.hpp"
#include "sidis/extra/math.hpp"
//#include "sidis/phinom.hpp"
#include <iostream>
#include <sstream>
#include <fstream>


using namespace std;
using namespace sidis::math;

///Add updates from exclusive contribution arXiv:2310.17961v1 (2023)

///Fine structure constant. Not sure about the Q2 dependence for exclusive structure function
double ALPHA = 7.2973525664e-3L;

double Get_thetacm(double W, double Q2, double t){
  //The exclusive structure function grid depends on thetacm
  double m_pi = 0.1349770;//define values here is too dirty
  //double m_n = 0.9395612928;
  double m_p = 0.9382720813;

  double mpi2=m_pi*m_pi;
  double mn2=m_n*m_n;
  double mp2=m_p*m_p;
  double W2 = W*W;
  double numerator = (W2 + mpi2 - mn2) * (W2 - Q2 - mp2) + 2.0 * W2 * (t + Q2 - mpi2);
  double term1 = std::pow(W2 + mpi2 - mn2, 2) - 4.0 * mpi2 * W2;
  double term2 = std::pow(W2 - Q2 - mp2, 2) + 4.0 * Q2 * W2;
  double denominator = std::sqrt(term1 * term2);

  double thcm = std::acos(numerator / denominator)*180/3.14159;
  return thcm;
}

double Get_Interpolated(double y0,double y1,double x0,double x1,double x){
  return y0+(y1-y0)*(x-x0)/(x1-x0);//Why would I write this as a function :(
}

//Just interpolate the average between two theta
std::vector<double> Get_exc_sf(double W, double Q2, double t){
  double thetacm = Get_thetacm(W,Q2,t);
  //std::cout<<"thetacm "<<thetacm<<std::endl;
  std::ifstream file("src/exc_sf_set/exclu.grid");
  std::vector<double> A_exc;
  //std::array<6> A_real;
  //std::array<6> A_ima;
  if(!file.is_open()){
    std::cout<<"Error opening file."<<std::endl;
    return A_exc;
  }

  std::string line;
  std::getline(file,line);//Ignore the first row
  while (std::getline(file,line)){
    std::istringstream iss(line);
    //std::cout<<line<<std::endl;
    double Key_Q2,Key_W,Key_theta;
    iss>>Key_Q2>>Key_W>>Key_theta;
    //std::cout<<"check keys"<<Key_Q2<<" "<<Key_W<<" "<<Key_theta<<std::endl;

    if(std::abs(Key_Q2-Q2)<0.3&&std::abs(Key_W/1000.0-W)<0.02&&std::abs(Key_theta-thetacm)<3){
      //std::cout<<"check keys"<<Key_Q2<<" "<<Key_W<<" "<<Key_theta<<std::endl;
      std::string nextline;
      std::getline(file,nextline);
      //std::cout<<line<<std::endl;
      //std::cout<<nextline<<std::endl;
      std::istringstream iss_next(nextline);
      double Keynext_Q2,Keynext_W,Keynext_theta;
      iss_next>>Keynext_Q2>>Keynext_W>>Keynext_theta;
      double number,number_next;
      while (iss>>number && iss_next>>number_next){
        double ave = Get_Interpolated(number,number_next,Key_theta,Keynext_theta,thetacm);
        //std::cout<<number<<" "<<number_next<<"ave "<<ave<<std::endl;
        A_exc.push_back(ave);
      }
      break;//found the row, exit the loop
    }
  }

  file.close();
  if(A_exc.empty()) {std::cout<<"Warning! Exclusive structure function! Couldn't find grids."<<std::endl;}
  return A_exc;

}

//Grab the value of the current row
//std::vector<double> Get_exc_sf(double W, double Q2, double t){
//  double thetacm = Get_thetacm(W,Q2,t);
//  std::cout<<"thetacm "<<thetacm<<std::endl;
//  std::ifstream file("exclu.grid");
//  std::vector<double> A_exc;
//  //std::array<6> A_real;
//  //std::array<6> A_ima;
//  if(!file.is_open()){
//    std::cout<<"Error opening file."<<std::endl;
//    return A_exc;
//  }
//  
//  std::string line;
//  std::getline(file,line);//Ignore the first row
//  while (std::getline(file,line)){
//    std::istringstream iss(line);
//    //std::cout<<line<<std::endl;
//    double Key_Q2,Key_W,Key_theta;
//    iss>>Key_Q2>>Key_W>>Key_theta;
//    //std::cout<<"check keys"<<Key_Q2<<" "<<Key_W<<" "<<Key_theta<<std::endl;
//
//    if(std::abs(Key_Q2-Q2)<0.3&&std::abs(Key_W/1000.0-W)<0.02&&std::abs(Key_theta-thetacm)<3){
//      std::cout<<"check keys"<<Key_Q2<<" "<<Key_W<<" "<<Key_theta<<std::endl;
//      double number;
//      while (iss>>number){
//        A_exc.push_back(number);
//      }
//      break;//found the row, exit the loop
//    }
//    else{std::cout<<"Couldn't find the corresponding Exclusive Values"<<std::endl;} 
//  }
//
//  file.close();
//  if(A_exc.empty()) {std::cout<<"Warning! Exclusive structure function! Couldn't find grids."<<std::endl;}
//  return A_exc;
//
//}

EXC_A::EXC_A(double W, double Q2, double t){
  A1r = Get_exc_sf(W,Q2,t)[0];
  A1i = Get_exc_sf(W,Q2,t)[1];
  A2r = Get_exc_sf(W,Q2,t)[2];
  A2i = Get_exc_sf(W,Q2,t)[3];
  A3r = Get_exc_sf(W,Q2,t)[4];
  A3i = Get_exc_sf(W,Q2,t)[5];
  A4r = Get_exc_sf(W,Q2,t)[6];
  A4i = Get_exc_sf(W,Q2,t)[7];
  A5r = Get_exc_sf(W,Q2,t)[8];
  A5i = Get_exc_sf(W,Q2,t)[9];
  A6r = Get_exc_sf(W,Q2,t)[10];
  A6i = Get_exc_sf(W,Q2,t)[11];

}

EXC_SF_F::EXC_SF_F(KinematicsRad kin){
  //C.20
  EXC_A excA(sqrt(kin.shift_W_sq),kin.shift_Q_sq,kin.shift_t);
  r1=kin.shift_r1;//sq(sqrt(kin.shift_W_sq)+kin.M)+kin.shift_Q_sq;
  r2=kin.shift_r2;//sq(sqrt(kin.shift_W_sq)-kin.M)+kin.shift_Q_sq;
  f1r= 1.0/(2*sqrt(2.0*PI*ALPHA))*(excA.A1r+(sqrt(kin.shift_W_sq)-kin.M)*excA.A4r+(kin.shift_Q_sq*excA.A6r+kin.shift_V_m*(excA.A3r-excA.A4r))/(sqrt(kin.shift_W_sq)-kin.M));
  f1i= 1.0/(2*sqrt(2.0*PI*ALPHA))*(excA.A1i+(sqrt(kin.shift_W_sq)-kin.M)*excA.A4i+(kin.shift_Q_sq*excA.A6i+kin.shift_V_m*(excA.A3i-excA.A4i))/(sqrt(kin.shift_W_sq)-kin.M));
  f2r=1.0/(2*sqrt(2*PI*ALPHA))*(-excA.A1r+(sqrt(kin.shift_W_sq)+kin.M)*excA.A4r+(kin.shift_Q_sq*excA.A6r+kin.shift_V_m*(excA.A3r-excA.A4r))/(sqrt(kin.shift_W_sq)+kin.M));
  f2i=1.0/(2*sqrt(2*PI*ALPHA))*(-excA.A1i+(sqrt(kin.shift_W_sq)+kin.M)*excA.A4i+(kin.shift_Q_sq*excA.A6i+kin.shift_V_m*(excA.A3i-excA.A4i))/(sqrt(kin.shift_W_sq)+kin.M));
  f3r=1.0/(2*sqrt(2*PI*ALPHA))*(excA.A3r-excA.A4r+(sqrt(kin.shift_W_sq)-kin.M)*excA.A2r+(kin.shift_Q_sq*(excA.A2r-2*excA.A5r))/(2*(sqrt(kin.shift_W_sq)+kin.M)));
  f3i=1.0/(2*sqrt(2*PI*ALPHA))*(excA.A3i-excA.A4i+(sqrt(kin.shift_W_sq)-kin.M)*excA.A2i+(kin.shift_Q_sq*(excA.A2i-2*excA.A5i))/(2*(sqrt(kin.shift_W_sq)+kin.M)));
  f4r=1.0/(2*sqrt(2*PI*ALPHA))*(excA.A3r-excA.A4r+(sqrt(kin.shift_W_sq)+kin.M)*excA.A2r+(kin.shift_Q_sq*(excA.A2r-2*excA.A5r))/(2*(sqrt(kin.shift_W_sq)-kin.M)));
  f4i=1.0/(2*sqrt(2*PI*ALPHA))*(excA.A3i-excA.A4i+(sqrt(kin.shift_W_sq)+kin.M)*excA.A2i+(kin.shift_Q_sq*(excA.A2i-2*excA.A5i))/(2*(sqrt(kin.shift_W_sq)-kin.M)));
  f5r=1.0/(8*sqrt(kin.shift_W_sq)*sqrt(2*PI*ALPHA))*((kin.shift_W_sq+sq(kin.mh)-sq(m_n))*(kin.shift_Q_sq*(excA.A2r-2*excA.A5r)+2*(kin.M+sqrt(kin.shift_W_sq))*((sqrt(kin.shift_W_sq-kin.M)*excA.A2r+excA.A3r-excA.A4r))+4*sqrt(kin.shift_W_sq)*(excA.A4r-2*sqrt(kin.shift_W_sq)*excA.A2r-excA.A3r))-(sq(sq(kin.M)+kin.shift_Q_sq-kin.shift_W_sq)+4*kin.shift_Q_sq*kin.shift_W_sq)*excA.A2r+2*r1*((sqrt(kin.shift_W_sq)-kin.M)*(excA.A4r-excA.A6r)+excA.A1r));
  f5i=1.0/(8*sqrt(kin.shift_W_sq)*sqrt(2*PI*ALPHA))*((kin.shift_W_sq+sq(kin.mh)-sq(m_n))*(kin.shift_Q_sq*(excA.A2i-2*excA.A5i)+2*(kin.M+sqrt(kin.shift_W_sq))*((sqrt(kin.shift_W_sq-kin.M)*excA.A2i+excA.A3i-excA.A4i))+4*sqrt(kin.shift_W_sq)*(excA.A4i-2*sqrt(kin.shift_W_sq)*excA.A2i-excA.A3i))-(sq(sq(kin.M)+kin.shift_Q_sq-kin.shift_W_sq)+4*kin.shift_Q_sq*kin.shift_W_sq)*excA.A2i+2*r1*((sqrt(kin.shift_W_sq)-kin.M)*(excA.A4i-excA.A6i)+excA.A1i));
  f6r=1.0/(8*sqrt(kin.shift_W_sq)*sqrt(2*PI*ALPHA))*(kin.shift_lambda_Y_sqrt*excA.A2r-2*r2*(excA.A1r+(sqrt(kin.shift_W_sq)+kin.M)*excA.A6r)+(kin.shift_W_sq+sq(kin.mh)-m_n*m_n)*(2*kin.shift_Q_sq*excA.A5r+2*(sqrt(kin.shift_W_sq)-kin.M)*excA.A3r+2*kin.M*excA.A4r+(2*sq(kin.M)-2*kin.shift_W_sq-kin.shift_Q_sq)*excA.A2r)+2*(kin.M*(sq(kin.M)+kin.shift_Q_sq-kin.shift_W_sq)-sqrt(kin.shift_W_sq)*(sq(kin.M)-sq(m_n)+kin.shift_t))*excA.A4r+kin.shift_V_m*(2*(kin.shift_W_sq-sq(kin.M)-kin.shift_Q_sq)*excA.A5r-4*sqrt(kin.shift_W_sq)*excA.A3r+(3*sq(kin.M)+5*kin.shift_W_sq+3*kin.shift_Q_sq)*excA.A2r));
  f6i=1.0/(8*sqrt(kin.shift_W_sq)*sqrt(2*PI*ALPHA))*(kin.shift_lambda_Y_sqrt*excA.A2i-2*r2*(excA.A1i+(sqrt(kin.shift_W_sq)+kin.M)*excA.A6i)+(kin.shift_W_sq+sq(kin.mh)-m_n*m_n)*(2*kin.shift_Q_sq*excA.A5i+2*(sqrt(kin.shift_W_sq)-kin.M)*excA.A3i+2*kin.M*excA.A4i+(2*sq(kin.M)-2*kin.shift_W_sq-kin.shift_Q_sq)*excA.A2i)+2*(kin.M*(sq(kin.M)+kin.shift_Q_sq-kin.shift_W_sq)-sqrt(kin.shift_W_sq)*(sq(kin.M)-sq(m_n)+kin.shift_t))*excA.A4i+kin.shift_V_m*(2*(kin.shift_W_sq-sq(kin.M)-kin.shift_Q_sq)*excA.A5i-4*sqrt(kin.shift_W_sq)*excA.A3i+(3*sq(kin.M)+5*kin.shift_W_sq+3*kin.shift_Q_sq)*excA.A2i));
}

EXC_SF_combine::EXC_SF_combine(EXC_SF_F exc_sf_f){
  //C.22
  f11=sq(exc_sf_f.f1r)+sq(exc_sf_f.f1i);
  f12r=exc_sf_f.f1r*exc_sf_f.f2r+exc_sf_f.f1i*exc_sf_f.f2i;
  f12i=exc_sf_f.f1r*exc_sf_f.f2r-exc_sf_f.f1i*exc_sf_f.f2i;
  f13r=exc_sf_f.f1r*exc_sf_f.f3r+exc_sf_f.f1i*exc_sf_f.f3i;
  f13i=exc_sf_f.f1r*exc_sf_f.f3r-exc_sf_f.f1i*exc_sf_f.f3i;
  f14r=exc_sf_f.f1r*exc_sf_f.f4r+exc_sf_f.f1i*exc_sf_f.f4i;
  f14i=exc_sf_f.f1r*exc_sf_f.f4r-exc_sf_f.f1i*exc_sf_f.f4i;
  f15r=exc_sf_f.f1r*exc_sf_f.f5r+exc_sf_f.f1i*exc_sf_f.f5i;
  f15i=exc_sf_f.f1r*exc_sf_f.f5r-exc_sf_f.f1i*exc_sf_f.f5i;
  f16r=exc_sf_f.f1r*exc_sf_f.f6r+exc_sf_f.f1i*exc_sf_f.f6i;
  f16i=exc_sf_f.f1r*exc_sf_f.f6r-exc_sf_f.f1i*exc_sf_f.f6i;
  f22=sq(exc_sf_f.f2r)+sq(exc_sf_f.f2i);
  f23r=exc_sf_f.f2r*exc_sf_f.f3r+exc_sf_f.f2i*exc_sf_f.f3i;
  f23i=exc_sf_f.f2r*exc_sf_f.f3r-exc_sf_f.f2i*exc_sf_f.f3i;
  f24r=exc_sf_f.f2r*exc_sf_f.f4r+exc_sf_f.f2i*exc_sf_f.f4i;
  f24i=exc_sf_f.f2r*exc_sf_f.f4r-exc_sf_f.f2i*exc_sf_f.f4i;
  f25r=exc_sf_f.f2r*exc_sf_f.f5r+exc_sf_f.f2i*exc_sf_f.f5i;
  f25i=exc_sf_f.f2r*exc_sf_f.f5r-exc_sf_f.f2i*exc_sf_f.f5i;
  f26r=exc_sf_f.f2r*exc_sf_f.f6r+exc_sf_f.f2i*exc_sf_f.f6i;
  f26i=exc_sf_f.f2r*exc_sf_f.f6r-exc_sf_f.f2i*exc_sf_f.f6i;
  f33=sq(exc_sf_f.f3r)+sq(exc_sf_f.f3i);
  f34r=exc_sf_f.f3r*exc_sf_f.f4r+exc_sf_f.f3i*exc_sf_f.f4i;
  f34i=exc_sf_f.f3r*exc_sf_f.f4r-exc_sf_f.f3i*exc_sf_f.f4i;
  f35r=exc_sf_f.f3r*exc_sf_f.f5r+exc_sf_f.f3i*exc_sf_f.f5i;
  f35i=exc_sf_f.f3r*exc_sf_f.f5r-exc_sf_f.f3i*exc_sf_f.f5i;
  f36r=exc_sf_f.f3r*exc_sf_f.f6r+exc_sf_f.f3i*exc_sf_f.f6i;
  f36i=exc_sf_f.f3r*exc_sf_f.f6r-exc_sf_f.f3i*exc_sf_f.f6i;
  f44=sq(exc_sf_f.f4r)+sq(exc_sf_f.f4i);
  f45r=exc_sf_f.f4r*exc_sf_f.f5r+exc_sf_f.f4i*exc_sf_f.f5i;
  f45i=exc_sf_f.f4r*exc_sf_f.f5r-exc_sf_f.f4i*exc_sf_f.f5i;
  f46r=exc_sf_f.f4r*exc_sf_f.f6r+exc_sf_f.f4i*exc_sf_f.f6i;
  f46i=exc_sf_f.f4r*exc_sf_f.f6r-exc_sf_f.f4i*exc_sf_f.f6i;
  f55=sq(exc_sf_f.f5r)+sq(exc_sf_f.f5i);
  f56r=exc_sf_f.f5r*exc_sf_f.f6r+exc_sf_f.f5i*exc_sf_f.f6i;
  f56i=exc_sf_f.f5r*exc_sf_f.f6r-exc_sf_f.f5i*exc_sf_f.f6i;
  f66=sq(exc_sf_f.f6r)+sq(exc_sf_f.f6i);
}

EXC_SF::EXC_SF(EXC_SF_combine exc_SF_com,KinematicsRad kin){
  r1=kin.shift_r1;//sq(sqrt(kin.shift_W_sq)+kin.M)+kin.shift_Q_sq;
  r2=kin.shift_r2;//sq(sqrt(kin.shift_W_sq)-kin.M)+kin.shift_Q_sq;
  r3=kin.shift_r3;//sq(sqrt(kin.shift_W_sq)+m_n)-sq(kin.mh);
  r4=kin.shift_r4;//sq(sqrt(kin.shift_W_sq)-m_n)-sq(kin.mh);
  r5=kin.shift_r5;//kin.shift_W_sq*(sq(kin.M)+sq(kin.mh)+sq(m_n)-kin.shift_Q_sq-2*kin.shift_t)+(sq(kin.M)+kin.shift_Q_sq)*(sq(kin.mh)-sq(m_n))-sq(kin.shift_W_sq);
  H00pH22_000= 
  (1 / kin.shift_W_sq) * (
    2 * kin.shift_Q_sq * kin.shift_W_sq * ( (r3 / r1) * exc_SF_com.f55 + (r4 / r2) * exc_SF_com.f66 ) +
    (r5 / (r1 * r2)) * ( (kin.shift_W_sq - sq(kin.M)) * kin.shift_lambda_Y_sqrt * exc_SF_com.f12r - 4 * kin.shift_Q_sq * kin.shift_W_sq * exc_SF_com.f56r ) +0.5 * (
      sq(sqrt(kin.shift_W_sq) - kin.M) * r1 * r3 * exc_SF_com.f11 +
      sq(sqrt(kin.shift_W_sq) + kin.M) * r2 * r4 * exc_SF_com.f22)
    );
  H11mH22opt2_000=
  (1 / (2 * kin.shift_W_sq)) * (sq(sqrt(kin.shift_W_sq) + kin.M) * r2 * (4 * sqrt(kin.shift_W_sq) * exc_SF_com.f23r + r3 * exc_SF_com.f33) + 2 * (sq(kin.M) - kin.shift_W_sq) * r5 * exc_SF_com.f34r +sq(sqrt(kin.shift_W_sq) - kin.M) * r1 * (4 * sqrt(kin.shift_W_sq) * exc_SF_com.f14r + r4 * exc_SF_com.f44));
  // 1/(2*kin.shift_W_sq)*(sq(sqrt(kin.shift_W_sq)+kin.M)*r2*(4*sqrt(kin.shift_W_sq)*exc_SF_com.f23r+r3*exc_SF_com.f33)+2*(sq(kin.M)-kin.shift_W_sq)*r5*exc_SF_com.f34r+sq(sqrt(kin.shift_W_sq)-kin.M)*r1*(4*sqrt(kin.shift_W_sq)*exc_SF_com.f14r+r4*exc_SF_com.f44));
  H22_000=(1 / (2 * kin.shift_W_sq)) * (sq(sqrt(kin.shift_W_sq) - kin.M) * r1 * r3 * exc_SF_com.f11 + sq(sqrt(kin.shift_W_sq) + kin.M) * r2 * r4 * exc_SF_com.f22 + 2 * (kin.shift_W_sq - sq(kin.M)) * r5 * exc_SF_com.f12r);
  H01r_000=(sqrt(kin.shift_Q_sq) * kin.shift_ph_t / (sqrt(kin.shift_W_sq) * kin.shift_lambda_Y_sqrt)) * (
    (sqrt(kin.shift_W_sq) - kin.M) * (r1 * (2 * sqrt(kin.shift_W_sq) * exc_SF_com.f16r + r4 * exc_SF_com.f46r) - r5 * exc_SF_com.f45r) +
    (sqrt(kin.shift_W_sq) + kin.M) * (r2 * (2 * sqrt(kin.shift_W_sq) * exc_SF_com.f25r + r3 * exc_SF_com.f35r) - r5 * exc_SF_com.f36r)
    );
  H01i_000=(sqrt(kin.shift_Q_sq) * kin.shift_ph_t / (sqrt(kin.shift_W_sq) * kin.shift_lambda_Y_sqrt)) * (
    (sqrt(kin.shift_W_sq) - kin.M) * (r5 * exc_SF_com.f45i - r1 * (2 * sqrt(kin.shift_W_sq) * exc_SF_com.f16i) + r4 * exc_SF_com.f46i) +
    (sqrt(kin.shift_W_sq) + kin.M) * (r5 * exc_SF_com.f36i - r2 * (2 * sqrt(kin.shift_W_sq) * exc_SF_com.f25i + r3 * exc_SF_com.f35i)));
  H00pH22_010 = (2  * kin.shift_ph_t / (sqrt(kin.shift_W_sq) * kin.shift_lambda_Y_sqrt)) * ((kin.shift_W_sq - sq(kin.M)) * kin.shift_lambda_Y_sqrt * exc_SF_com.f12i + 4 * kin.shift_Q_sq * kin.shift_W_sq * exc_SF_com.f56i);
  H11mH22opt2_010= (1. / (kin.shift_W_sq * sq(kin.shift_ph_t) * kin.shift_lambda_Y_sqrt)) * (r5 * (sq(sqrt(kin.shift_W_sq) - kin.M) * r1 * exc_SF_com.f14i - sq(sqrt(kin.shift_W_sq) + kin.M) * r2 * exc_SF_com.f23i) + (kin.shift_W_sq - sq(kin.M)) * kin.shift_lambda_Y_sqrt * (r4 * exc_SF_com.f24i + 2 * sqrt(kin.shift_W_sq) * (sq(kin.shift_ph_t) * exc_SF_com.f34i - 2 * exc_SF_com.f12i)) - r3 * exc_SF_com.f13i);
  H22_010=(2 * kin.shift_ph_t * kin.shift_lambda_Y_sqrt / sqrt(kin.shift_W_sq)) *(kin.shift_W_sq - sq(kin.M)) * exc_SF_com.f12i;  
  H01r_010=(sqrt(kin.shift_Q_sq) / (r1 * r2 * sqrt(kin.shift_W_sq))) * (
    (sqrt(kin.shift_W_sq) - kin.M) * (r1 * r5 * exc_SF_com.f16i - kin.shift_lambda_Y * (r3 * exc_SF_com.f15i + 2 * sqrt(kin.shift_W_sq) * sq(kin.shift_ph_t) * exc_SF_com.f45i)) +
    (sqrt(kin.shift_W_sq) + kin.M) * (kin.shift_lambda_Y * (r4 * exc_SF_com.f26i + 2 * sqrt(kin.shift_W_sq) * sq(kin.shift_ph_t) * exc_SF_com.f36i) - r2 * r5 * exc_SF_com.f25i)
    );
  H01i_010=(sqrt(kin.shift_Q_sq) / (r1 * r2 * sqrt(kin.shift_W_sq))) * (
    (sqrt(kin.shift_W_sq) - kin.M) * (r1 * r5 * exc_SF_com.f16r - kin.shift_lambda_Y * (r3 * exc_SF_com.f15r + 2 * sqrt(kin.shift_W_sq) * sq(kin.shift_ph_t) * exc_SF_com.f45r)) +
    (sqrt(kin.shift_W_sq) + kin.M) * (kin.shift_lambda_Y * (r4 * exc_SF_com.f26r + 2 * sqrt(kin.shift_W_sq) * sq(kin.shift_ph_t) * exc_SF_com.f36r) - r2 * r5 * exc_SF_com.f25r));
  H02r_100=(sqrt(kin.shift_Q_sq) / (r1 * r2 * sqrt(kin.shift_W_sq))) * (
    (sqrt(kin.shift_W_sq) - kin.M) * (r3 * kin.shift_lambda_Y * exc_SF_com.f15i - r1 * r5 * exc_SF_com.f16i) +   (sqrt(kin.shift_W_sq) + kin.M) * (r2 * (r5 * exc_SF_com.f25i - r1 * r4 * exc_SF_com.f26i))
    );
  H02i_100=(sqrt(kin.shift_Q_sq) / (r1 * r2 * sqrt(kin.shift_W_sq))) * (
    (sqrt(kin.shift_W_sq) - kin.M) * (r3 * kin.shift_lambda_Y * exc_SF_com.f15r - r1 * r5 * exc_SF_com.f16r) +   (sqrt(kin.shift_W_sq) + kin.M) * (r2 * (r5 * exc_SF_com.f25r - r1 * r4 * exc_SF_com.f26r)));
  H12r_100=(kin.shift_ph_t  / (2 * sq(kin.shift_W_sq) * kin.shift_lambda_Y_sqrt)) * (
    (kin.shift_W_sq - sq(kin.M)) * (r1 * r2 * (4 * sqrt(kin.shift_W_sq) * exc_SF_com.f12i - r4 * exc_SF_com.f24i) + r3 * kin.shift_lambda_Y * exc_SF_com.f13i) +
    r5 * (sq(sqrt(kin.shift_W_sq) + kin.M) * r2 * exc_SF_com.f23i - sq(sqrt(kin.shift_W_sq) - kin.M) * r1 * exc_SF_com.f14i)
    );
  H12i_100= (kin.shift_ph_t  / (2 * kin.shift_W_sq * kin.shift_lambda_Y_sqrt)) * (
    (kin.shift_W_sq - sq(kin.M)) * (r3 * kin.shift_lambda_Y * exc_SF_com.f13r - r1 * r2 * r4 * exc_SF_com.f24r) +
    r5 * (sq(sqrt(kin.shift_W_sq) + kin.M) * r2 * exc_SF_com.f23r - sq(sqrt(kin.shift_W_sq) - kin.M) * r1 * exc_SF_com.f14r));
  H02r_001=(-2* sqrt(kin.shift_Q_sq) * kin.shift_ph_t / (kin.shift_lambda_Y_sqrt)) * (
    (sqrt(kin.shift_W_sq) + kin.M) * r2 * exc_SF_com.f25i + (sqrt(kin.shift_W_sq) - kin.M) * r1 * exc_SF_com.f16i);
  H02i_001=(-2 * sqrt(kin.shift_Q_sq) * kin.shift_ph_t / (kin.shift_lambda_Y_sqrt)) * (
    (sqrt(kin.shift_W_sq) + kin.M) * r2 * exc_SF_com.f25r + (sqrt(kin.shift_W_sq) - kin.M) * r1 * exc_SF_com.f16r);
  H12r_001= (-sq(kin.shift_ph_t) / sqrt(kin.shift_W_sq)) * (
    sq(sqrt(kin.shift_W_sq) - kin.M) * r1 * exc_SF_com.f14i + sq(sqrt(kin.shift_W_sq) + kin.M) * r2 * exc_SF_com.f23i
    );
  H12i_001= (-1.0 / (2 * kin.shift_W_sq * r1 * r2)) * (
    sq(sqrt(kin.shift_W_sq) - kin.M) * r1 * (2 *sqrt(kin.shift_W_sq)* kin.shift_lambda_Y * sq(kin.shift_ph_t) * exc_SF_com.f14r + r1 * r2 * r3 * exc_SF_com.f11) +
    sq(sqrt(kin.shift_W_sq) + kin.M) * r2 * (2 *sqrt(kin.shift_W_sq)* kin.shift_lambda_Y * sq(kin.shift_ph_t) * exc_SF_com.f23r + r1 * r2 * r4 * exc_SF_com.f22) +
    2 * (sq(kin.shift_W_sq) - sq(kin.M)) * r5 * kin.shift_lambda_Y * exc_SF_com.f12r);

}

//Generalized exclusive structure functions of exclusive process equation 42
EXCUU::EXCUU(EXC_SF exc_sf, KinematicsRad kin){
  H1_000=exc_sf.H22_000;
  H2_000=1/kin.shift_lambda_Y_sqrt*(4*kin.shift_Q_sq*(exc_sf.H00pH22_000)-4*sqrt(kin.shift_Q_sq)/kin.shift_ph_t*exc_sf.H01r_000+sq(kin.shift_rex)*(exc_sf.H11mH22opt2_000));
  H3_000=exc_sf.H11mH22opt2_000;
  H4_000=1/kin.shift_lambda_Y_sqrt*(-kin.shift_rex*(exc_sf.H11mH22opt2_000)+2*sqrt(kin.shift_Q_sq)/kin.shift_ph_t*exc_sf.H01r_000);
}
EXCLU::EXCLU(EXC_SF exc_sf, KinematicsRad kin){
  H5_000=2*sqrt(kin.shift_Q_sq)/(kin.shift_ph_t*kin.shift_lambda_Y_sqrt)*exc_sf.H01i_000;
}
EXCUT::EXCUT(EXC_SF exc_sf, KinematicsRad kin){
  H1_010=exc_sf.H22_010;
  H2_010=1/kin.shift_lambda_Y_sqrt*(4*kin.shift_Q_sq*(exc_sf.H00pH22_010)-4*sqrt(kin.shift_Q_sq)/kin.shift_ph_t*exc_sf.H01r_010+sq(kin.shift_rex)*(exc_sf.H11mH22opt2_010));
  H3_010=exc_sf.H11mH22opt2_010;
  H4_010=1/kin.shift_lambda_Y_sqrt*(-kin.shift_rex*(exc_sf.H11mH22opt2_010)+2*sqrt(kin.shift_Q_sq)/kin.shift_ph_t*exc_sf.H01r_010);
  H6_100=2/kin.shift_ph_t/kin.shift_lambda_Y*(2*sqrt(kin.shift_Q_sq)*exc_sf.H02r_100-kin.shift_rex/kin.shift_ph_t*exc_sf.H12r_100);
  H8_100=2/sq(kin.shift_ph_t)/kin.shift_lambda_Y_sqrt*exc_sf.H12r_100;
}
EXCLT::EXCLT(EXC_SF exc_sf,KinematicsRad kin){
  H5_010=2*sqrt(kin.shift_Q_sq)/(kin.shift_ph_t*kin.shift_lambda_Y_sqrt)*exc_sf.H01i_010;
  H7_100=2/kin.shift_ph_t/kin.shift_lambda_Y*(2*sqrt(kin.shift_Q_sq)*exc_sf.H02i_100-kin.shift_rex/kin.shift_ph_t*exc_sf.H12i_100);
  H9_100=2/sq(kin.shift_ph_t)/kin.shift_lambda_Y_sqrt*exc_sf.H12i_100;
}
EXCUL::EXCUL(EXC_SF exc_sf, KinematicsRad kin){
  kin.shift_rex=2*(kin.shift_Q_sq*(sq(kin.M)-sq(m_n)+kin.shift_S_x+kin.shift_t)+kin.shift_S_x*kin.shift_V_m)/kin.shift_lambda_Y_sqrt;
  H6_001=2/kin.shift_ph_t/kin.shift_lambda_Y*(2*sqrt(kin.shift_Q_sq)*exc_sf.H02r_001-kin.shift_rex/kin.shift_ph_t*exc_sf.H12r_001);
  H8_001=2/sq(kin.shift_ph_t)/kin.shift_lambda_Y_sqrt*exc_sf.H12r_001;
}
EXCLL::EXCLL(EXC_SF exc_sf, KinematicsRad kin){
  kin.shift_rex=2*(kin.shift_Q_sq*(sq(kin.M)-sq(m_n)+kin.shift_S_x+kin.shift_t)+kin.shift_S_x*kin.shift_V_m)/kin.shift_lambda_Y_sqrt;
  H7_001=2/kin.shift_ph_t/kin.shift_lambda_Y*(2*sqrt(kin.shift_Q_sq)*exc_sf.H02i_001-kin.shift_rex/kin.shift_ph_t*exc_sf.H12i_001);
  H9_001=2/sq(kin.shift_ph_t)/kin.shift_lambda_Y_sqrt*exc_sf.H12i_001;

}
