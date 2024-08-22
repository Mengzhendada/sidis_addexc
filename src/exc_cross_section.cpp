#include "sidis/cross_section.hpp"

#include <cmath>
#include <limits>
#include <string>

#include "sidis/bound.hpp"
#include "sidis/constant.hpp"
#include "sidis/cut.hpp"
#include "sidis/frame.hpp"
#include "sidis/hadronic_coeff.hpp"
#include "sidis/kinematics.hpp"
#include "sidis/leptonic_coeff.hpp"
#include "sidis/phenom.hpp"
#include "sidis/structure_function.hpp"
#include "sidis/vector.hpp"
#include "sidis/extra/integrate.hpp"
#include "sidis/extra/math.hpp"

using namespace sidis;
using namespace sidis::cut;
using namespace sidis::had;
using namespace sidis::kin;
using namespace sidis::lep;
using namespace sidis::math;
using namespace sidis::ph;
using namespace sidis::sf;
using namespace sidis::xs;

// Macro that computes the cross-section from the base cross-sections in an
// optimized way. For example, if the polarization is zero, then the
// cross-section can be computed just using the UU base cross-section.

// This one works with base cross-
// sections where the XL, XT1, and XT2 cases are all grouped together into an
// XP case.
#define SIDIS_MACRO_XS_FROM_BASE_P(name, Lep, Had, kin, sf, b, lambda_e, eta) ([&]() { \
	/* Create a mask describing the polarization state. */ \
	unsigned pol_mask = (((lambda_e) != 0.) << 1) \
		| (((eta).x != 0. || (eta).y != 0. || (eta).z != 0.) << 0); \
	Real uu = 0.; \
	Vec3 up = VEC3_ZERO; \
	Real lu = 0.; \
	Vec3 lp = VEC3_ZERO; \
	switch (pol_mask) { \
	case 0:  /* 00 */ \
		{ \
			Lep##UU lep((kin)); \
			Had##UU had((kin), (sf)); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
		} \
		break; \
	case 1:  /* 01 */ \
		{ \
			Lep##UP lep((kin)); \
			Had##UP had((kin), (sf)); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
			up = name##_base_up((b), lep.uu, lep.up, had.up); \
		} \
		break; \
	case 2:  /* 10 */ \
		{ \
			Lep##LU lep((kin)); \
			Had##LU had((kin), (sf)); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
			lu = name##_base_lu((b), lep.lu, had.lu); \
		} \
		break; \
	case 3:  /* 11 */ \
		{ \
			Lep##LP lep((kin)); \
			Had##LP had((kin), (sf)); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
			up = name##_base_up((b), lep.uu, lep.up, had.up); \
			lu = name##_base_lu((b), lep.lu, had.lu); \
			lp = name##_base_lp((b), lep.lu, lep.lp, had.lp); \
		} \
		break; \
	} \
	return uu + dot(up, (eta)) + (lambda_e)*(lu + dot(lp, (eta))); \
}())



Real xs::exc(KinematicsRad const& kin, Phenom const& phenom, SfSet const& sf, Real lambda_e, Vec3 eta) {
	Rad b(kin, phenom);
	return SIDIS_MACRO_XS_FROM_BASE_P(exc, LepRad, HadRad, kin, sf, b, lambda_e, eta);
}



EstErr xs::exc_integ(Kinematics const& kin, Phenom const& phenom, SfSet const& sf, Real lambda_e, Vec3 eta, Real k_0_bar, IntegParams params) {
	CutRad cut;
	cut.k_0_bar = Bound(k_0_bar, INF);
	EstErr xs_integ = integrate<3>(
		[&](std::array<Real, 3> x) {
			KinematicsRad kin_exc;
			Real jac;
			if (!take(cut, kin, x.data(), &kin_exc, &jac)) {
				return 0.;
			}
			Rad b(kin_exc, phenom);
			Real xs = jac * SIDIS_MACRO_XS_FROM_BASE_P(exc, LepRad, HadRad, kin_exc, sf, b, lambda_e, eta);
			if (std::isnan(xs)) {
				return 0.;
			} else {
				return xs;
			}
		},
		std::array<Real, 3>{ 0., 0., 0. },
		std::array<Real, 3>{ 1., 1., 1. },
		params);
	return xs_integ;
}

Real xs::exc(KinematicsRad const& kin, SfSet const& sf, Real lambda_e, Vec3 eta) {
	return exc(kin, Phenom(kin.project()), sf, lambda_e, eta);
}
EstErr xs::exc_integ(Kinematics const& kin, SfSet const& sf, Real lambda_e, Vec3 eta, Real k_0_bar, IntegParams params) {
	return exc_integ(kin, Phenom(kin), sf, lambda_e, eta, k_0_bar, params);
}




// Exclusive radiative base functions.
EXC::EXC(KinematicsRad const& kin, Phenom const& phenom) {
	// Equation [36].
	coeff = -(std::pow(phenom.alpha_qed, 3)*kin.S*sq(kin.S_x))
		/(512.*std::pow(PI,5)*kin.M*kin.ph_l*kin.lambda_S*kin.lambda_Y_sqrt);
	Rex = kin.Rex;//Need to define this Rex
}

Real xs::exc_base_uu(EXC const& b, LepRadBaseUU const& lep, EXCUU const& exc) {
	return b.coeff*(
		1./b.Rex*(
			lep.theta_011*exc.H1_000
			+ lep.theta_021*exc.H2_000
			+ lep.theta_031*exc.H3_000
			+ lep.theta_041*exc.H4_000)
		+ (
			lep.theta_012*exc.H1_000
			+ lep.theta_022*exc.H2_000
			+ lep.theta_032*exc.H3_000
			+ lep.theta_042*exc.H4_000)
		+ b.Rex*(
			lep.theta_013*exc.H1_000
			+ lep.theta_023*exc.H2_000
			+ lep.theta_033*exc.H3_000
			+ lep.theta_043*exc.H4_000));
}
//For Exclusive case, I don't think that we have to combine UL and UT, but it's combined in SIDIS case, so here it is
Vec3 xs::exc_base_up(EXC const& b, LepRadBaseUU const& lep_uu, LepRadBaseUP const & lep_up, EXCUT const& excut, EXCUL const& excul){
        //in SIDIS case, had.H1 is Vec3. Here in EXC hadron structure function, I wrote everything in double, so I have to build Vec3 here
	return b.coeff*(
		1/b.Rex*(
			lep_uu.theta_011*Vec3(0.,excut.H1_010,0.)
			+lep_uu.theta_021*Vec3(0.,excut.H2_010,0.)
			+lep_uu.theta_031*Vec3(0.,excut.H3_010,0.)
			+lep_uu.theta_041*Vec3(0.,excut.H4_010,0.)
			+lep_up.theta_061*Vec3(excut.H6_100,0.,excul.H6_001)
			+lep_up.theta_081*Vec3(excut.H8_100,0.,excul.H8_001)
			)
		+(
			lep_uu.theta_012*Vec3(0.,excut.H1_010,0.)
			+lep_uu.theta_022*Vec3(0.,excut.H2_010,0.)
			+lep_uu.theta_032*Vec3(0.,excut.H3_010,0.)
			+lep_uu.theta_042*Vec3(0.,excut.H4_010,0.)
			+lep_up.theta_062*Vec3(excut.H6_100,0.,excul.H6_001)
			+lep_up.theta_082*Vec3(excut.H8_100,0.,excul.H8_001)
			)
		+ b.Rex*(
			lep_uu.theta_013*Vec3(0.,excut.H1_010,0.)
			+lep_uu.theta_023*Vec3(0.,excut.H2_010,0.)
			+lep_uu.theta_033*Vec3(0.,excut.H3_010,0.)
			+lep_uu.theta_043*Vec3(0.,excut.H4_010,0.)
			+lep_up.theta_063*Vec3(excut.H6_100,0.,excul.H6_001)
			+lep_up.theta_083*Vec3(excut.H8_100,0.,excul.H8_001)
			)
		+ b.Rex*b.Rex*(
                        lep_up.theta_064*Vec3(excut.H6_100,0.,excul.H6_001)
                        lep_up.theta_084*Vec3(excut.H8_100,0.,excul.H8_001)
			)
		);
}
//Vec3 xs::exc_base_up(EXC const& b, LepRadBaseUU const& lep_uu, LepRadBaseUP const& lep_up, HadRadBaseUP const& had) {
//	return b.coeff*(
//		1./b.Rex*(
//			lep_uu.theta_011*had.H_1
//			+ lep_uu.theta_021*had.H_2
//			+ lep_uu.theta_031*had.H_3
//			+ lep_uu.theta_041*had.H_4
//			+ lep_up.theta_061*had.H_6
//			+ lep_up.theta_081*had.H_8)
//		+ (
//			lep_uu.theta_012*had.H_1
//			+ lep_uu.theta_022*had.H_2
//			+ lep_uu.theta_032*had.H_3
//			+ lep_uu.theta_042*had.H_4
//			+ lep_up.theta_062*had.H_6
//			+ lep_up.theta_082*had.H_8)
//		+ b.R*(
//			lep_uu.theta_013*had.H_1
//			+ lep_uu.theta_023*had.H_2
//			+ lep_uu.theta_033*had.H_3
//			+ lep_uu.theta_043*had.H_4
//			+ lep_up.theta_063*had.H_6
//			+ lep_up.theta_083*had.H_8)
//		+ b.R*b.R*(
//			lep_up.theta_064*had.H_6
//			+ lep_up.theta_084*had.H_8));
//}
Real xs::exc_base_lu(EXC const& b, LepRadBaseLU const& lep, EXCLU const& exc) {
	return b.coeff*(
		1./b.Rex*(lep.theta_051 + lep.theta_151)*exc.H5_000
		+ (lep.theta_052 + lep.theta_152)*exc.H5_000
		+ b.Rex*(lep.theta_053 + lep.theta_153)*exc.H5_000);
}
Vec3 xs::exc_base_lp(EXC const& b, LepRadBaseLU const& lep_lu, LepRadBaseLP const& lep_lp, EXCLT const& exclt, EXCLL const& excll){
	return b.coeff*(
		1./b.Rex*(
			(lep_lu.theta_051 + lep_lu.theta_151)*Vec3(0., exclt.H5_010, 0.)
			+ (lep_lp.theta_071 + lep_lp.theta_171)*Vec3(exclt.H7_100,0.,excll.H7_001)
			+ (lep_lp.theta_091 + lep_lp.theta_191)*Vec3(exclt.H9_100,0.,exclt.H9_001)
			)
		+ (
		 	(lep_lu.theta_052 + lep_lu.theta_152)*Vec3(0.,exclt.H5_010,0.)
			+ (lep_lp.theta_072 + lep_lp.theta_172)*Vec3(exclt.H7_100,0.,excll.H7_001)
			+ (lep_lp.theta_092 + lep_lp.theta_192)*Vec3(exclt.H9_100,0.,excll.H9_001)
			)
		+ b.Rex*(
		 	(lep_lu.theta_053 + lep_lu.theta_153)*Vec3(0.,exclt.H5_010,0.)
			+ (lep_lp.theta_073 + lep_lp.theta_173)*Vec3(exclt.H7_100,0.,excll.H7_001)
			+ (lep_lp.theta_093 + lep_lp.theta_193)*Vec3(exclt.H9_100,0.,excll.H9_001)
			)
		+ b.Rex*b.Rex*(
			(lep_lp.theta_074 + lep_lp.theta_174)*Vec3(exclt.H7_100,0.,excll.H7_001)
			+ (lep_lp.theta_094 + lep_lp.theta_194)*Vec3(exclt.H9_100,0.,excll.H9_001)
			)
		);
		

}

