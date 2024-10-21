#include "sidis/Exc_structure_function.hpp"
#include "sidis/constant.hpp"
#include "sidis/extra/math.hpp"
//#include "sidis/phenom.hpp"
#include "sidis/frame.hpp"
#include "sidis/transform.hpp"
#include <iostream>
#include <sstream>
#include <fstream>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>


using namespace std;
using namespace sidis::math;

///Add updates from exclusive contribution arXiv:2310.17961v1 (2023)

///Fine structure constant. Not sure about the Q2 dependence for exclusive structure function
namespace sidis{namespace exc{ 
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

//double Get_Interpolated(double y0,double y1,double x0,double x1,double x){
//  return y0+(y1-y0)*(x-x0)/(x1-x0);//Why would I write this as a function :(
//}
double linearInterpolate(double x0, double x1, double f0, double f1, double x) {
        if(x1-x0==0){
		return f0;
	}
	return f0 + (f1 - f0) * (x - x0) / (x1 - x0);
}
void printDataPoint(DataPoint dp){
	std::cout<<" Q2: "<<dp.Q2<<" W: "<<dp.W<<" theta "<<dp.theta<<std::endl;
	std::cout<<"Values: ";
	for (size_t i = 0; i < dp.values.size(); ++i) {
        std::cout << dp.values[i] << " ";
    }
    std::cout << std::endl;
}

// Function to perform trilinear interpolation for a given Q2, W, and theta
std::vector<double> Get_exc_sf(double W, double Q2, double t) {
  double theta = Get_thetacm(W,Q2,t);
  //std::cout<<" check in Get_exc_sf theta "<<theta<<std::endl;
  //std::cout<<" check in Get_exc_sf W "<<W<<" Q2 "<<Q2<<" t "<<t<<std::endl;
	double a2=1.15,a30=-1.23,a31=0.16,a3;
	double Q2corr,Wcorr,Q2arg,W2arg;
	if(Q2>5){
		Q2corr=std::pow(5.0,a2)/std::pow(Q2,a2);
		Q2arg = 5.0;
	}
	else{
		Q2corr=1.0;
		Q2arg=Q2;
	}
	if(W>2){
		a3=a30+a31*theta;
		if(theta<50) a3 = a30+a31*50.0;
		else if(theta>180) a3 = a30+a31*100.0;
		Wcorr = std::pow(2.0,a3)/std::pow(W*W,(a3/2.0));
		W2arg = 4.0;
	}
	else{
		W2arg = W*W;
		Wcorr = 1.0;
	}
	

	//I think I should read the file somewhere else just once
	std::vector<DataPoint> data;
	std::ifstream file("src/exc_sf_set/exclu.grid");
	if(!file.is_open()){
		std::cout<<"Error opening file."<<std::endl;
	}
	std::string line;
	std::getline(file,line);//Ignore the first row
	while (std::getline(file,line)){
		std::istringstream iss(line);
		DataPoint dp;
		iss>>dp.Q2>>dp.W>>dp.theta;
		double val;
		while (iss >>val){
			dp.values.push_back(val);//push the rest of the values
		}
		data.push_back(dp);
	}
	// This is a simplified trilinear interpolation between the 8 surrounding points

	// Find the 8 surrounding points
	DataPoint lowerQ2_lowerW_lowerTheta, lowerQ2_lowerW_upperTheta;
	DataPoint lowerQ2_upperW_lowerTheta, lowerQ2_upperW_upperTheta;
	DataPoint upperQ2_lowerW_lowerTheta, upperQ2_lowerW_upperTheta;
	DataPoint upperQ2_upperW_lowerTheta, upperQ2_upperW_upperTheta;

	bool foundtheta=false;
	bool foundQ2=false;
	bool foundW = false;

	// Find the 8 surrounding points for Q2, W, and theta (you can optimize this search)
	for (size_t i = 0; i < data.size() - 1; ++i) {

		// Check for exact match on Q2 and W
		if (Q2 == data[i].Q2 && W == data[i].W / 1000) {

			// We only need to interpolate over theta
			if ((theta - data[i].theta)>0&&(theta - data[i].theta) < 3 ) {
				//std::cout << "Exact match for Q2 and W at index " << i << std::endl;
				lowerQ2_lowerW_lowerTheta = data[i];
				lowerQ2_lowerW_upperTheta = (std::abs(theta - data[i].theta) == 0) ? data[i] : data[i + 1];
				lowerQ2_upperW_lowerTheta = data[i];
				lowerQ2_upperW_upperTheta = (std::abs(theta - data[i].theta) == 0) ? data[i] : data[i + 1];
				upperQ2_lowerW_lowerTheta = data[i];
				upperQ2_lowerW_upperTheta = (std::abs(theta - data[i].theta) == 0) ? data[i] : data[i + 1];
				upperQ2_upperW_lowerTheta = data[i];
				upperQ2_upperW_upperTheta = (std::abs(theta - data[i].theta) == 0) ? data[i] : data[i + 1];
				foundtheta=true;
				foundQ2=true;
				foundW=true;
			}
		}
		else if((Q2-data[i].Q2)==0 && abs(W-data[i].W/1000)>0){
			if ((W - data[i].W / 1000) > 0 && (1000*W - data[i].W) < 20 &&
			    (theta - data[i].theta) > 3&&(theta - data[i].theta) < 3) {
				//std::cout<<" why W lower " <<abs(1000*W-data[i].W)<<std::endl;
				//std::cout << "Exact match for Q2 at index " << i << " lower W "<<std::endl;

				lowerQ2_lowerW_lowerTheta = data[i];
				lowerQ2_lowerW_upperTheta = (std::abs(theta - data[i].theta) == 0) ? data[i] : data[i + 1];
				upperQ2_lowerW_lowerTheta = data[i];
				upperQ2_lowerW_upperTheta = (std::abs(theta - data[i].theta) == 0) ? data[i] : data[i + 1];
				foundtheta=true;
				foundQ2=true;
				foundW=true;
			}
			else if ((data[i].W / 1000 - W) > 0 && (data[i].W - 1000*W) < 20 &&
				 (theta - data[i].theta) > 3&&(theta - data[i].theta) < 3) {
				//std::cout<<" why W upper " <<abs(1000*W-data[i].W)<<std::endl;
				//std::cout << "Exact match for Q2 at index " << i << " upper W "<< std::endl;
				lowerQ2_upperW_lowerTheta = data[i];
				lowerQ2_upperW_upperTheta = (std::abs(theta - data[i].theta) == 0) ? data[i] : data[i + 1];
				upperQ2_upperW_lowerTheta = data[i];
				upperQ2_upperW_upperTheta = (std::abs(theta - data[i].theta) == 0) ? data[i] : data[i + 1];
				foundtheta=true;
				foundQ2=true;
				foundW=true;
			}
		}//Q2 exactly match
		// Check for exact match on W (but not Q2)
		else if (W == data[i].W / 1000 && abs(Q2-data[i].Q2)>0) {

			// We only need to interpolate over Q2 and theta
			if ((Q2 - data[i].Q2) > 0 && (Q2 - data[i].Q2) < 0.3 && (theta - data[i].theta)>0&&(theta - data[i].theta) < 3) {
				//std::cout << "Exact match for W at index " << i << std::endl;
				lowerQ2_lowerW_lowerTheta = data[i];
				lowerQ2_lowerW_upperTheta = (std::abs(theta - data[i].theta) == 0) ? data[i] : data[i + 1];
				lowerQ2_upperW_lowerTheta = data[i];
				lowerQ2_upperW_upperTheta = (std::abs(theta - data[i].theta) == 0) ? data[i] : data[i + 1];
				foundtheta=true;
				foundQ2=true;
				foundW=true;
			}
			else if ((data[i].Q2 - Q2) > 0 && (data[i].Q2 - Q2) < 0.3 &&  (theta - data[i].theta)>0&&(theta - data[i].theta) < 3) {
				//std::cout << "Exact match for W at index " << i << std::endl;
				upperQ2_lowerW_lowerTheta = data[i];
				upperQ2_lowerW_upperTheta = (std::abs(theta - data[i].theta) == 0) ? data[i] : data[i + 1];
				upperQ2_upperW_lowerTheta = data[i];
				upperQ2_upperW_upperTheta = (std::abs(theta - data[i].theta) == 0) ? data[i] : data[i + 1];
				foundtheta=true;
				foundQ2=true;
				foundW=true;
			}
		}
		//regular
		else if (
			(Q2-data[i].Q2)>0 && (Q2-data[i].Q2)<0.3 && 
			(W-data[i].W/1000)>0 && (W-data[i].W/1000)<0.02 &&
			(theta-data[i].theta)>0&&(theta-data[i].theta) < 3
			) {
			//std::cout<< " check in Get_exc_sf for loop lowerQ2 lowerW "<<i<<std::endl;

			lowerQ2_lowerW_lowerTheta = data[i];
			lowerQ2_lowerW_upperTheta = (std::abs(theta - data[i].theta) == 0) ? data[i] : data[i + 1];
			foundtheta = true;
		}
		else if(
			(Q2-data[i].Q2)>0 && (Q2-data[i].Q2)<0.3 && 
			(data[i].W/1000-W)>0 &&(data[i].W/1000-W)<0.02 &&
			(theta-data[i].theta)>0&&(theta-data[i].theta) < 3
		       ){
			//std::cout<< " check in Get_exc_sf for loop lowerQ2 upperW "<<i<<std::endl;
			lowerQ2_upperW_lowerTheta = data[i];
			lowerQ2_upperW_upperTheta = (std::abs(theta - data[i].theta) == 0) ? data[i] : data[i + 1];
			foundW = true;
		}
		else if(
			(data[i].Q2-Q2)>0 && (data[i].Q2-Q2)<0.3 && 
			(W-data[i].W/1000)>0 && (W-data[i].W/1000)<0.02 &&
			(theta-data[i].theta)>0&&(theta-data[i].theta) < 3
		       ){
			//std::cout<< " check in Get_exc_sf for loop upperQ2 lowerW "<<i<<std::endl;
			upperQ2_lowerW_lowerTheta = data[i];
			upperQ2_lowerW_upperTheta = (std::abs(theta - data[i].theta) == 0) ? data[i] : data[i + 1];
			foundQ2=true;
		}
		else if(
			(data[i].Q2-Q2)>0 && (data[i].Q2-Q2)<0.3 && 
			(data[i].W/1000-W)>0 && (data[i].W/1000-W)<0.02 &&
			(theta-data[i].theta)>0&&(theta-data[i].theta) < 3
		       ){
			//std::cout<< " check in Get_exc_sf for loop upperQ2 upperW "<<i<<std::endl;
			upperQ2_upperW_lowerTheta = data[i];
			upperQ2_upperW_upperTheta = (std::abs(theta - data[i].theta) == 0) ? data[i] : data[i + 1];
		}
	}

	if (!(foundtheta&&foundQ2&&foundW)) {
		std::cerr << "Could not find appropriate data points for interpolation." << std::endl;
		return std::vector<double>(12, 0.0); // Return zeros if not found
	}

	//print out the eight edges to check
	//std::cout<<"lowerQ2 lowerW lowerTheta"<<std::endl;
	//printDataPoint(lowerQ2_lowerW_lowerTheta); 
	//std::cout<<"lowerQ2 lowerW upperTheta"<<std::endl;
	//printDataPoint(lowerQ2_lowerW_upperTheta); 
	//std::cout<<"lowerQ2 upperW lowerTheta"<<std::endl;
	//printDataPoint(lowerQ2_upperW_lowerTheta); 
	//std::cout<<"lowerQ2 upperW upperTheta"<<std::endl;
	//printDataPoint(lowerQ2_upperW_upperTheta); 
	//std::cout<<"upperQ2 lowerW lowerTheta"<<std::endl;
	//printDataPoint(upperQ2_lowerW_lowerTheta); 
	//std::cout<<"upperQ2 lowerW upperTheta"<<std::endl;
	//printDataPoint(upperQ2_lowerW_upperTheta); 
	//std::cout<<"upperQ2 upperW lowerTheta"<<std::endl;
	//printDataPoint(upperQ2_upperW_lowerTheta); 
	//std::cout<<"upperQ2 upperW upperTheta"<<std::endl;
	//printDataPoint(upperQ2_upperW_upperTheta); 

	// Interpolate between the 8 corner points
	std::vector<double> interpolatedValues(12, 0.0);

	for (size_t i = 0; i < 12; ++i) {
		// Interpolate along the theta axis
		double lowerQ2_lowerW = linearInterpolate(lowerQ2_lowerW_lowerTheta.theta, lowerQ2_lowerW_upperTheta.theta, lowerQ2_lowerW_lowerTheta.values[i], lowerQ2_lowerW_upperTheta.values[i], theta);
		//std::cout<<"check each 1D"<< lowerQ2_lowerW<<" x0 "<<lowerQ2_lowerW_lowerTheta.theta<<" x1 "<< lowerQ2_lowerW_upperTheta.theta<<" y0 "<< lowerQ2_lowerW_lowerTheta.values[i]<<" y1 "<< lowerQ2_lowerW_upperTheta.values[i]<<" theta "<< theta<<std::endl;
		double lowerQ2_upperW = linearInterpolate(lowerQ2_upperW_lowerTheta.theta, lowerQ2_upperW_upperTheta.theta, lowerQ2_upperW_lowerTheta.values[i], lowerQ2_upperW_upperTheta.values[i], theta);
		double upperQ2_lowerW = linearInterpolate(upperQ2_lowerW_lowerTheta.theta, upperQ2_lowerW_upperTheta.theta, upperQ2_lowerW_lowerTheta.values[i], upperQ2_lowerW_upperTheta.values[i], theta);
		double upperQ2_upperW = linearInterpolate(upperQ2_upperW_lowerTheta.theta, upperQ2_upperW_upperTheta.theta, upperQ2_upperW_lowerTheta.values[i], upperQ2_upperW_upperTheta.values[i], theta);

		// Interpolate along the W axis
		double lowerQ2 = linearInterpolate(lowerQ2_lowerW_lowerTheta.W/1000, lowerQ2_upperW_lowerTheta.W/1000, lowerQ2_lowerW, lowerQ2_upperW, W);
		double upperQ2 = linearInterpolate(upperQ2_lowerW_lowerTheta.W/1000, upperQ2_upperW_lowerTheta.W/1000, upperQ2_lowerW, upperQ2_upperW, W);

		// Interpolate along the Q2 axis
		interpolatedValues[i] = linearInterpolate(lowerQ2_lowerW_lowerTheta.Q2, upperQ2_upperW_lowerTheta.Q2, lowerQ2, upperQ2, Q2);
	}

	return interpolatedValues;
}


EXC_A::EXC_A(double W, double Q2, double t){
	//std::cout<<" check in EXC_A, W: "<<W<<" Q2: "<<Q2<<" t: "<<t<<std::endl;
        std::vector<double> As = Get_exc_sf(W,Q2,t);
	A1r = As[0];
	A1i = As[1];
	A2r = As[2];
	A2i = As[3];
	A3r = As[4];
	A3i = As[5];
	A4r = As[6];
	A4i = As[7];
        A5r = As[8];
        A5i = As[9];
        A6r = As[10];
        A6i = As[11];


	//A1r = Get_exc_sf(W,Q2,t)[0];
	//A1i = Get_exc_sf(W,Q2,t)[1];
	//A2r = Get_exc_sf(W,Q2,t)[2];
	//A2i = Get_exc_sf(W,Q2,t)[3];
	//A3r = Get_exc_sf(W,Q2,t)[4];
	//A3i = Get_exc_sf(W,Q2,t)[5];
	//A4r = Get_exc_sf(W,Q2,t)[6];
	//A4i = Get_exc_sf(W,Q2,t)[7];
	//A5r = Get_exc_sf(W,Q2,t)[8];
	//A5i = Get_exc_sf(W,Q2,t)[9];
	//A6r = Get_exc_sf(W,Q2,t)[10];
	//A6i = Get_exc_sf(W,Q2,t)[11];

}

EXC_SF_F::EXC_SF_F(KinematicsRad kin){
  //C.20
  EXC_A excA(sqrt(kin.shiftexc_W_sq),kin.shiftexc_Q_sq,kin.shiftexc_t);
  //std::cout<<"in EXC_SF_F W: "<<sqrt(kin.W_sq)<<" Q2: "<<kin.Q_sq<<" t: "<<kin.t<<std::endl;
  //std::cout<<"in EXC_SF_F shift W: "<<sqrt(kin.shiftexc_W_sq)<<" Q2: "<<kin.shiftexc_Q_sq<<" t: "<<kin.shiftexc_t<<std::endl;
  r1=kin.shiftexc_r1;//sq(sqrt(kin.shiftexc_W_sq)+kin.M)+kin.shiftexc_Q_sq;
  r2=kin.shiftexc_r2;//sq(sqrt(kin.shiftexc_W_sq)-kin.M)+kin.shiftexc_Q_sq;
  f1r= 1.0/(2*sqrt(2.0*PI*ALPHA))*(excA.A1r+(sqrt(kin.shiftexc_W_sq)-kin.M)*excA.A4r+(kin.shiftexc_Q_sq*excA.A6r+kin.shiftexc_V_m*(excA.A3r-excA.A4r))/(sqrt(kin.shiftexc_W_sq)-kin.M));
  f1i= 1.0/(2*sqrt(2.0*PI*ALPHA))*(excA.A1i+(sqrt(kin.shiftexc_W_sq)-kin.M)*excA.A4i+(kin.shiftexc_Q_sq*excA.A6i+kin.shiftexc_V_m*(excA.A3i-excA.A4i))/(sqrt(kin.shiftexc_W_sq)-kin.M));
  f2r=1.0/(2*sqrt(2*PI*ALPHA))*(-excA.A1r+(sqrt(kin.shiftexc_W_sq)+kin.M)*excA.A4r+(kin.shiftexc_Q_sq*excA.A6r+kin.shiftexc_V_m*(excA.A3r-excA.A4r))/(sqrt(kin.shiftexc_W_sq)+kin.M));
  f2i=1.0/(2*sqrt(2*PI*ALPHA))*(-excA.A1i+(sqrt(kin.shiftexc_W_sq)+kin.M)*excA.A4i+(kin.shiftexc_Q_sq*excA.A6i+kin.shiftexc_V_m*(excA.A3i-excA.A4i))/(sqrt(kin.shiftexc_W_sq)+kin.M));
  f3r=1.0/(2*sqrt(2*PI*ALPHA))*(excA.A3r-excA.A4r+(sqrt(kin.shiftexc_W_sq)-kin.M)*excA.A2r+(kin.shiftexc_Q_sq*(excA.A2r-2*excA.A5r))/(2*(sqrt(kin.shiftexc_W_sq)+kin.M)));
  f3i=1.0/(2*sqrt(2*PI*ALPHA))*(excA.A3i-excA.A4i+(sqrt(kin.shiftexc_W_sq)-kin.M)*excA.A2i+(kin.shiftexc_Q_sq*(excA.A2i-2*excA.A5i))/(2*(sqrt(kin.shiftexc_W_sq)+kin.M)));
  f4r=1.0/(2*sqrt(2*PI*ALPHA))*(excA.A3r-excA.A4r+(sqrt(kin.shiftexc_W_sq)+kin.M)*excA.A2r+(kin.shiftexc_Q_sq*(excA.A2r-2*excA.A5r))/(2*(sqrt(kin.shiftexc_W_sq)-kin.M)));
  f4i=1.0/(2*sqrt(2*PI*ALPHA))*(excA.A3i-excA.A4i+(sqrt(kin.shiftexc_W_sq)+kin.M)*excA.A2i+(kin.shiftexc_Q_sq*(excA.A2i-2*excA.A5i))/(2*(sqrt(kin.shiftexc_W_sq)-kin.M)));
  f5r=1.0/(8*sqrt(kin.shiftexc_W_sq)*sqrt(2*PI*ALPHA))*((kin.shiftexc_W_sq+sq(kin.mh)-sq(m_n))*(kin.shiftexc_Q_sq*(excA.A2r-2*excA.A5r)+2*(kin.M+sqrt(kin.shiftexc_W_sq))*((sqrt(kin.shiftexc_W_sq-kin.M)*excA.A2r+excA.A3r-excA.A4r))+4*sqrt(kin.shiftexc_W_sq)*(excA.A4r-2*sqrt(kin.shiftexc_W_sq)*excA.A2r-excA.A3r))-(sq(sq(kin.M)+kin.shiftexc_Q_sq-kin.shiftexc_W_sq)+4*kin.shiftexc_Q_sq*kin.shiftexc_W_sq)*excA.A2r+2*r1*((sqrt(kin.shiftexc_W_sq)-kin.M)*(excA.A4r-excA.A6r)+excA.A1r));
  f5i=1.0/(8*sqrt(kin.shiftexc_W_sq)*sqrt(2*PI*ALPHA))*((kin.shiftexc_W_sq+sq(kin.mh)-sq(m_n))*(kin.shiftexc_Q_sq*(excA.A2i-2*excA.A5i)+2*(kin.M+sqrt(kin.shiftexc_W_sq))*((sqrt(kin.shiftexc_W_sq-kin.M)*excA.A2i+excA.A3i-excA.A4i))+4*sqrt(kin.shiftexc_W_sq)*(excA.A4i-2*sqrt(kin.shiftexc_W_sq)*excA.A2i-excA.A3i))-(sq(sq(kin.M)+kin.shiftexc_Q_sq-kin.shiftexc_W_sq)+4*kin.shiftexc_Q_sq*kin.shiftexc_W_sq)*excA.A2i+2*r1*((sqrt(kin.shiftexc_W_sq)-kin.M)*(excA.A4i-excA.A6i)+excA.A1i));
  f6r=1.0/(8*sqrt(kin.shiftexc_W_sq)*sqrt(2*PI*ALPHA))*(kin.shiftexc_lambda_Y_sqrt*excA.A2r-2*r2*(excA.A1r+(sqrt(kin.shiftexc_W_sq)+kin.M)*excA.A6r)+(kin.shiftexc_W_sq+sq(kin.mh)-m_n*m_n)*(2*kin.shiftexc_Q_sq*excA.A5r+2*(sqrt(kin.shiftexc_W_sq)-kin.M)*excA.A3r+2*kin.M*excA.A4r+(2*sq(kin.M)-2*kin.shiftexc_W_sq-kin.shiftexc_Q_sq)*excA.A2r)+2*(kin.M*(sq(kin.M)+kin.shiftexc_Q_sq-kin.shiftexc_W_sq)-sqrt(kin.shiftexc_W_sq)*(sq(kin.M)-sq(m_n)+kin.shiftexc_t))*excA.A4r+kin.shiftexc_V_m*(2*(kin.shiftexc_W_sq-sq(kin.M)-kin.shiftexc_Q_sq)*excA.A5r-4*sqrt(kin.shiftexc_W_sq)*excA.A3r+(3*sq(kin.M)+5*kin.shiftexc_W_sq+3*kin.shiftexc_Q_sq)*excA.A2r));
  f6i=1.0/(8*sqrt(kin.shiftexc_W_sq)*sqrt(2*PI*ALPHA))*(kin.shiftexc_lambda_Y_sqrt*excA.A2i-2*r2*(excA.A1i+(sqrt(kin.shiftexc_W_sq)+kin.M)*excA.A6i)+(kin.shiftexc_W_sq+sq(kin.mh)-m_n*m_n)*(2*kin.shiftexc_Q_sq*excA.A5i+2*(sqrt(kin.shiftexc_W_sq)-kin.M)*excA.A3i+2*kin.M*excA.A4i+(2*sq(kin.M)-2*kin.shiftexc_W_sq-kin.shiftexc_Q_sq)*excA.A2i)+2*(kin.M*(sq(kin.M)+kin.shiftexc_Q_sq-kin.shiftexc_W_sq)-sqrt(kin.shiftexc_W_sq)*(sq(kin.M)-sq(m_n)+kin.shiftexc_t))*excA.A4i+kin.shiftexc_V_m*(2*(kin.shiftexc_W_sq-sq(kin.M)-kin.shiftexc_Q_sq)*excA.A5i-4*sqrt(kin.shiftexc_W_sq)*excA.A3i+(3*sq(kin.M)+5*kin.shiftexc_W_sq+3*kin.shiftexc_Q_sq)*excA.A2i));
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
  r1=kin.shiftexc_r1;//sq(sqrt(kin.shiftexc_W_sq)+kin.M)+kin.shiftexc_Q_sq;
  r2=kin.shiftexc_r2;//sq(sqrt(kin.shiftexc_W_sq)-kin.M)+kin.shiftexc_Q_sq;
  r3=kin.shiftexc_r3;//sq(sqrt(kin.shiftexc_W_sq)+m_n)-sq(kin.mh);
  r4=kin.shiftexc_r4;//sq(sqrt(kin.shiftexc_W_sq)-m_n)-sq(kin.mh);
  r5=kin.shiftexc_r5;//kin.shiftexc_W_sq*(sq(kin.M)+sq(kin.mh)+sq(m_n)-kin.shiftexc_Q_sq-2*kin.shiftexc_t)+(sq(kin.M)+kin.shiftexc_Q_sq)*(sq(kin.mh)-sq(m_n))-sq(kin.shiftexc_W_sq);
  H00pH22_000= 
  (1 / kin.shiftexc_W_sq) * (
    2 * kin.shiftexc_Q_sq * kin.shiftexc_W_sq * ( (r3 / r1) * exc_SF_com.f55 + (r4 / r2) * exc_SF_com.f66 ) +
    (r5 / (r1 * r2)) * ( (kin.shiftexc_W_sq - sq(kin.M)) * kin.shiftexc_lambda_Y_sqrt * exc_SF_com.f12r - 4 * kin.shiftexc_Q_sq * kin.shiftexc_W_sq * exc_SF_com.f56r ) +0.5 * (
      sq(sqrt(kin.shiftexc_W_sq) - kin.M) * r1 * r3 * exc_SF_com.f11 +
      sq(sqrt(kin.shiftexc_W_sq) + kin.M) * r2 * r4 * exc_SF_com.f22)
    );
  H11mH22opt2_000=
  (1 / (2 * kin.shiftexc_W_sq)) * (sq(sqrt(kin.shiftexc_W_sq) + kin.M) * r2 * (4 * sqrt(kin.shiftexc_W_sq) * exc_SF_com.f23r + r3 * exc_SF_com.f33) + 2 * (sq(kin.M) - kin.shiftexc_W_sq) * r5 * exc_SF_com.f34r +sq(sqrt(kin.shiftexc_W_sq) - kin.M) * r1 * (4 * sqrt(kin.shiftexc_W_sq) * exc_SF_com.f14r + r4 * exc_SF_com.f44));
  // 1/(2*kin.shiftexc_W_sq)*(sq(sqrt(kin.shiftexc_W_sq)+kin.M)*r2*(4*sqrt(kin.shiftexc_W_sq)*exc_SF_com.f23r+r3*exc_SF_com.f33)+2*(sq(kin.M)-kin.shiftexc_W_sq)*r5*exc_SF_com.f34r+sq(sqrt(kin.shiftexc_W_sq)-kin.M)*r1*(4*sqrt(kin.shiftexc_W_sq)*exc_SF_com.f14r+r4*exc_SF_com.f44));
  H22_000=(1 / (2 * kin.shiftexc_W_sq)) * (sq(sqrt(kin.shiftexc_W_sq) - kin.M) * r1 * r3 * exc_SF_com.f11 + sq(sqrt(kin.shiftexc_W_sq) + kin.M) * r2 * r4 * exc_SF_com.f22 + 2 * (kin.shiftexc_W_sq - sq(kin.M)) * r5 * exc_SF_com.f12r);
  H01r_000=(sqrt(kin.shiftexc_Q_sq) * kin.shiftexc_ph_t / (sqrt(kin.shiftexc_W_sq) * kin.shiftexc_lambda_Y_sqrt)) * (
    (sqrt(kin.shiftexc_W_sq) - kin.M) * (r1 * (2 * sqrt(kin.shiftexc_W_sq) * exc_SF_com.f16r + r4 * exc_SF_com.f46r) - r5 * exc_SF_com.f45r) +
    (sqrt(kin.shiftexc_W_sq) + kin.M) * (r2 * (2 * sqrt(kin.shiftexc_W_sq) * exc_SF_com.f25r + r3 * exc_SF_com.f35r) - r5 * exc_SF_com.f36r)
    );
  H01i_000=(sqrt(kin.shiftexc_Q_sq) * kin.shiftexc_ph_t / (sqrt(kin.shiftexc_W_sq) * kin.shiftexc_lambda_Y_sqrt)) * (
    (sqrt(kin.shiftexc_W_sq) - kin.M) * (r5 * exc_SF_com.f45i - r1 * (2 * sqrt(kin.shiftexc_W_sq) * exc_SF_com.f16i) + r4 * exc_SF_com.f46i) +
    (sqrt(kin.shiftexc_W_sq) + kin.M) * (r5 * exc_SF_com.f36i - r2 * (2 * sqrt(kin.shiftexc_W_sq) * exc_SF_com.f25i + r3 * exc_SF_com.f35i)));
  H00pH22_010 = (2  * kin.shiftexc_ph_t / (sqrt(kin.shiftexc_W_sq) * kin.shiftexc_lambda_Y_sqrt)) * ((kin.shiftexc_W_sq - sq(kin.M)) * kin.shiftexc_lambda_Y_sqrt * exc_SF_com.f12i + 4 * kin.shiftexc_Q_sq * kin.shiftexc_W_sq * exc_SF_com.f56i);
  H11mH22opt2_010= (1. / (kin.shiftexc_W_sq * sq(kin.shiftexc_ph_t) * kin.shiftexc_lambda_Y_sqrt)) * (r5 * (sq(sqrt(kin.shiftexc_W_sq) - kin.M) * r1 * exc_SF_com.f14i - sq(sqrt(kin.shiftexc_W_sq) + kin.M) * r2 * exc_SF_com.f23i) + (kin.shiftexc_W_sq - sq(kin.M)) * kin.shiftexc_lambda_Y_sqrt * (r4 * exc_SF_com.f24i + 2 * sqrt(kin.shiftexc_W_sq) * (sq(kin.shiftexc_ph_t) * exc_SF_com.f34i - 2 * exc_SF_com.f12i)) - r3 * exc_SF_com.f13i);
  H22_010=(2 * kin.shiftexc_ph_t * kin.shiftexc_lambda_Y_sqrt / sqrt(kin.shiftexc_W_sq)) *(kin.shiftexc_W_sq - sq(kin.M)) * exc_SF_com.f12i;  
  H01r_010=(sqrt(kin.shiftexc_Q_sq) / (r1 * r2 * sqrt(kin.shiftexc_W_sq))) * (
    (sqrt(kin.shiftexc_W_sq) - kin.M) * (r1 * r5 * exc_SF_com.f16i - kin.shiftexc_lambda_Y * (r3 * exc_SF_com.f15i + 2 * sqrt(kin.shiftexc_W_sq) * sq(kin.shiftexc_ph_t) * exc_SF_com.f45i)) +
    (sqrt(kin.shiftexc_W_sq) + kin.M) * (kin.shiftexc_lambda_Y * (r4 * exc_SF_com.f26i + 2 * sqrt(kin.shiftexc_W_sq) * sq(kin.shiftexc_ph_t) * exc_SF_com.f36i) - r2 * r5 * exc_SF_com.f25i)
    );
  H01i_010=(sqrt(kin.shiftexc_Q_sq) / (r1 * r2 * sqrt(kin.shiftexc_W_sq))) * (
    (sqrt(kin.shiftexc_W_sq) - kin.M) * (r1 * r5 * exc_SF_com.f16r - kin.shiftexc_lambda_Y * (r3 * exc_SF_com.f15r + 2 * sqrt(kin.shiftexc_W_sq) * sq(kin.shiftexc_ph_t) * exc_SF_com.f45r)) +
    (sqrt(kin.shiftexc_W_sq) + kin.M) * (kin.shiftexc_lambda_Y * (r4 * exc_SF_com.f26r + 2 * sqrt(kin.shiftexc_W_sq) * sq(kin.shiftexc_ph_t) * exc_SF_com.f36r) - r2 * r5 * exc_SF_com.f25r));
  H02r_100=(sqrt(kin.shiftexc_Q_sq) / (r1 * r2 * sqrt(kin.shiftexc_W_sq))) * (
    (sqrt(kin.shiftexc_W_sq) - kin.M) * (r3 * kin.shiftexc_lambda_Y * exc_SF_com.f15i - r1 * r5 * exc_SF_com.f16i) +   (sqrt(kin.shiftexc_W_sq) + kin.M) * (r2 * (r5 * exc_SF_com.f25i - r1 * r4 * exc_SF_com.f26i))
    );
  H02i_100=(sqrt(kin.shiftexc_Q_sq) / (r1 * r2 * sqrt(kin.shiftexc_W_sq))) * (
    (sqrt(kin.shiftexc_W_sq) - kin.M) * (r3 * kin.shiftexc_lambda_Y * exc_SF_com.f15r - r1 * r5 * exc_SF_com.f16r) +   (sqrt(kin.shiftexc_W_sq) + kin.M) * (r2 * (r5 * exc_SF_com.f25r - r1 * r4 * exc_SF_com.f26r)));
  H12r_100=(kin.shiftexc_ph_t  / (2 * sq(kin.shiftexc_W_sq) * kin.shiftexc_lambda_Y_sqrt)) * (
    (kin.shiftexc_W_sq - sq(kin.M)) * (r1 * r2 * (4 * sqrt(kin.shiftexc_W_sq) * exc_SF_com.f12i - r4 * exc_SF_com.f24i) + r3 * kin.shiftexc_lambda_Y * exc_SF_com.f13i) +
    r5 * (sq(sqrt(kin.shiftexc_W_sq) + kin.M) * r2 * exc_SF_com.f23i - sq(sqrt(kin.shiftexc_W_sq) - kin.M) * r1 * exc_SF_com.f14i)
    );
  H12i_100= (kin.shiftexc_ph_t  / (2 * kin.shiftexc_W_sq * kin.shiftexc_lambda_Y_sqrt)) * (
    (kin.shiftexc_W_sq - sq(kin.M)) * (r3 * kin.shiftexc_lambda_Y * exc_SF_com.f13r - r1 * r2 * r4 * exc_SF_com.f24r) +
    r5 * (sq(sqrt(kin.shiftexc_W_sq) + kin.M) * r2 * exc_SF_com.f23r - sq(sqrt(kin.shiftexc_W_sq) - kin.M) * r1 * exc_SF_com.f14r));
  H02r_001=(-2* sqrt(kin.shiftexc_Q_sq) * kin.shiftexc_ph_t / (kin.shiftexc_lambda_Y_sqrt)) * (
    (sqrt(kin.shiftexc_W_sq) + kin.M) * r2 * exc_SF_com.f25i + (sqrt(kin.shiftexc_W_sq) - kin.M) * r1 * exc_SF_com.f16i);
  H02i_001=(-2 * sqrt(kin.shiftexc_Q_sq) * kin.shiftexc_ph_t / (kin.shiftexc_lambda_Y_sqrt)) * (
    (sqrt(kin.shiftexc_W_sq) + kin.M) * r2 * exc_SF_com.f25r + (sqrt(kin.shiftexc_W_sq) - kin.M) * r1 * exc_SF_com.f16r);
  H12r_001= (-sq(kin.shiftexc_ph_t) / sqrt(kin.shiftexc_W_sq)) * (
    sq(sqrt(kin.shiftexc_W_sq) - kin.M) * r1 * exc_SF_com.f14i + sq(sqrt(kin.shiftexc_W_sq) + kin.M) * r2 * exc_SF_com.f23i
    );
  H12i_001= (-1.0 / (2 * kin.shiftexc_W_sq * r1 * r2)) * (
    sq(sqrt(kin.shiftexc_W_sq) - kin.M) * r1 * (2 *sqrt(kin.shiftexc_W_sq)* kin.shiftexc_lambda_Y * sq(kin.shiftexc_ph_t) * exc_SF_com.f14r + r1 * r2 * r3 * exc_SF_com.f11) +
    sq(sqrt(kin.shiftexc_W_sq) + kin.M) * r2 * (2 *sqrt(kin.shiftexc_W_sq)* kin.shiftexc_lambda_Y * sq(kin.shiftexc_ph_t) * exc_SF_com.f23r + r1 * r2 * r4 * exc_SF_com.f22) +
    2 * (sq(kin.shiftexc_W_sq) - sq(kin.M)) * r5 * kin.shiftexc_lambda_Y * exc_SF_com.f12r);

}

//Generalized exclusive structure functions of exclusive process equation 42
EXCUU::EXCUU(KinematicsRad kin){
	EXC_SF_F exc_sf_f(kin);
	EXC_SF_combine exc_sf_combine(exc_sf_f);
	EXC_SF exc_sf(exc_sf_combine,kin);
	H1_000=exc_sf.H22_000;
	H2_000=1/kin.shiftexc_lambda_Y_sqrt*(4*kin.shiftexc_Q_sq*(exc_sf.H00pH22_000)-4*sqrt(kin.shiftexc_Q_sq)/kin.shiftexc_ph_t*exc_sf.H01r_000+sq(kin.shiftexc_rex)*(exc_sf.H11mH22opt2_000));
	H3_000=exc_sf.H11mH22opt2_000;
	H4_000=1/kin.shiftexc_lambda_Y_sqrt*(-kin.shiftexc_rex*(exc_sf.H11mH22opt2_000)+2*sqrt(kin.shiftexc_Q_sq)/kin.shiftexc_ph_t*exc_sf.H01r_000);
}
EXCLU::EXCLU(KinematicsRad kin){
	EXC_SF_F exc_sf_f(kin);
	EXC_SF_combine exc_sf_combine(exc_sf_f);
	EXC_SF exc_sf(exc_sf_combine,kin);
	H5_000=2*sqrt(kin.shiftexc_Q_sq)/(kin.shiftexc_ph_t*kin.shiftexc_lambda_Y_sqrt)*exc_sf.H01i_000;
}
EXCUT::EXCUT(KinematicsRad kin){
	EXC_SF_F exc_sf_f(kin);
	EXC_SF_combine exc_sf_combine(exc_sf_f);
	EXC_SF exc_sf(exc_sf_combine,kin);
	Transform3 shiftexc_rot = frame::hadron_from_shiftexc(kin);
	H1_010=dot(shiftexc_rot,Vec3(
			0.,
			exc_sf.H22_010,
			0.));
	H2_010=dot(shiftexc_rot,Vec3(
			0.,
			1/kin.shiftexc_lambda_Y_sqrt*(4*kin.shiftexc_Q_sq*(exc_sf.H00pH22_010)-4*sqrt(kin.shiftexc_Q_sq)/kin.shiftexc_ph_t*exc_sf.H01r_010+sq(kin.shiftexc_rex)*(exc_sf.H11mH22opt2_010)),
			0.));
	H3_010=dot(shiftexc_rot,Vec3(
		0.,
		exc_sf.H11mH22opt2_010,
		0.));
	H4_010=dot(shiftexc_rot,Vec3(
			0.,
			1/kin.shiftexc_lambda_Y_sqrt*(-kin.shiftexc_rex*(exc_sf.H11mH22opt2_010)+2*sqrt(kin.shiftexc_Q_sq)/kin.shiftexc_ph_t*exc_sf.H01r_010),
			0.));
	H6_100=dot(shiftexc_rot,Vec3(
			2/kin.shiftexc_ph_t/kin.shiftexc_lambda_Y*(2*sqrt(kin.shiftexc_Q_sq)*exc_sf.H02r_100-kin.shiftexc_rex/kin.shiftexc_ph_t*exc_sf.H12r_100),
			0.,
			0.));
	H8_100=dot(shiftexc_rot,Vec3(
			2/sq(kin.shiftexc_ph_t)/kin.shiftexc_lambda_Y_sqrt*exc_sf.H12r_100,
			0.,
			0.));
}
EXCLT::EXCLT(KinematicsRad kin){
	EXC_SF_F exc_sf_f(kin);
	EXC_SF_combine exc_sf_combine(exc_sf_f);
	EXC_SF exc_sf(exc_sf_combine,kin);
	Transform3 shiftexc_rot = frame::hadron_from_shiftexc(kin);
	H5_010=dot(shiftexc_rot,Vec3(
			0.,
			2*sqrt(kin.shiftexc_Q_sq)/(kin.shiftexc_ph_t*kin.shiftexc_lambda_Y_sqrt)*exc_sf.H01i_010,
			0.));
	H7_100=dot(shiftexc_rot,Vec3(
			2/kin.shiftexc_ph_t/kin.shiftexc_lambda_Y*(2*sqrt(kin.shiftexc_Q_sq)*exc_sf.H02i_100-kin.shiftexc_rex/kin.shiftexc_ph_t*exc_sf.H12i_100),
			0.,
			0.));
	H9_100=dot(shiftexc_rot,Vec3(
			2/sq(kin.shiftexc_ph_t)/kin.shiftexc_lambda_Y_sqrt*exc_sf.H12i_100,
			0.,
			0.));
}
EXCUL::EXCUL(KinematicsRad kin){
	EXC_SF_F exc_sf_f(kin);
	EXC_SF_combine exc_sf_combine(exc_sf_f);
	EXC_SF exc_sf(exc_sf_combine,kin);
	Transform3 shiftexc_rot = frame::hadron_from_shiftexc(kin);
	H6_001=dot(shiftexc_rot,Vec3(
			0.,
			0.,
			2/kin.shiftexc_ph_t/kin.shiftexc_lambda_Y*(2*sqrt(kin.shiftexc_Q_sq)*exc_sf.H02r_001-kin.shiftexc_rex/kin.shiftexc_ph_t*exc_sf.H12r_001)));
	H8_001=dot(shiftexc_rot,Vec3(
			0.,
			0.,
			2/sq(kin.shiftexc_ph_t)/kin.shiftexc_lambda_Y_sqrt*exc_sf.H12r_001));
}
EXCLL::EXCLL(KinematicsRad kin){
	EXC_SF_F exc_sf_f(kin);
	EXC_SF_combine exc_sf_combine(exc_sf_f);
	EXC_SF exc_sf(exc_sf_combine,kin);
	Transform3 shiftexc_rot = frame::hadron_from_shiftexc(kin);
	H7_001=dot(shiftexc_rot,Vec3(
			0.,
			0.,
			2/kin.shiftexc_ph_t/kin.shiftexc_lambda_Y*(2*sqrt(kin.shiftexc_Q_sq)*exc_sf.H02i_001-kin.shiftexc_rex/kin.shiftexc_ph_t*exc_sf.H12i_001)));
	H9_001=dot(shiftexc_rot,Vec3(
			0.,
			0.,
			2/sq(kin.shiftexc_ph_t)/kin.shiftexc_lambda_Y_sqrt*exc_sf.H12i_001));

}
}}//namespace
