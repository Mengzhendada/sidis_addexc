#include "sidis/cross_section.hpp"

#include <cmath>
#include <limits>
#include <string>
#include <iostream>
#include <fstream>
#include <array>

#include "sidis/bound.hpp"
#include "sidis/constant.hpp"
#include "sidis/cut.hpp"
#include "sidis/frame.hpp"
#include "sidis/hadronic_coeff.hpp"
#include "sidis/kinematics.hpp"
#include "sidis/leptonic_coeff.hpp"
#include "sidis/phenom.hpp"
#include "sidis/structure_function.hpp"
#include "sidis/Exc_structure_function.hpp"
#include "sidis/vector.hpp"
#include "sidis/extra/integrate.hpp"
#include "sidis/extra/math.hpp"

using namespace sidis;
using namespace sidis::cut;
using namespace sidis::had;
using namespace sidis::exc;
using namespace sidis::kin;
using namespace sidis::lep;
using namespace sidis::math;
using namespace sidis::ph;
using namespace sidis::sf;
using namespace sidis::xs;

// Macro that computes the cross-section from the base cross-sections in an
// optimized way. For example, if the polarization is zero, then the
// cross-section can be computed just using the UU base cross-section.
#define SIDIS_MACRO_XS_FROM_BASE(name, Lep, Had, kin, sf, b, lambda_e, eta) ([&]() { \
	/* Create a mask describing the polarization state. */ \
	unsigned pol_mask = (((lambda_e) != 0.) << 3) \
		| (((eta).x != 0.) << 2) \
		| (((eta).y != 0.) << 1) \
		| (((eta).z != 0.) << 0); \
	Real uu = 0.; \
	Vec3 up = VEC3_ZERO; \
	Real lu = 0.; \
	Vec3 lp = VEC3_ZERO; \
	switch (pol_mask) { \
	case 0:  /* 0000 */ \
		{ \
			Lep##UU lep((kin)); \
			Had##UU had((kin), (sf)); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
		} \
		break; \
	case 1:  /* 0001 */ \
		{ \
			Lep##UP lep((kin)); \
			Had##UL had((kin), (sf)); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
			up.z = name##_base_ul((b), lep.up, had.ul); \
		} \
		break; \
	case 2:  /* 0010 */ \
		{ \
			Lep##UU lep((kin)); \
			Had##UT had((kin), (sf)); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
			up.y = name##_base_ut2((b), lep.uu, had.ut); \
		} \
		break; \
	case 3:  /* 0011 */ \
		{ \
			Lep##UP lep((kin)); \
			Had##UP had((kin), (sf)); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
			up.y = name##_base_ut2((b), lep.uu, had.ut); \
			up.z = name##_base_ul((b), lep.up, had.ul); \
		} \
		break; \
	case 4:  /* 0100 */ \
		{ \
			Lep##UP lep((kin)); \
			Had##UT had((kin), (sf)); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
			up.x = name##_base_ut1((b), lep.up, had.ut); \
		} \
		break; \
	case 5:  /* 0101 */ \
		{ \
			Lep##UP lep((kin)); \
			Had##UP had((kin), (sf)); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
			up.x = name##_base_ut1((b), lep.up, had.ut); \
			up.z = name##_base_ul((b), lep.up, had.ul); \
		} \
		break; \
	case 6:  /* 0110 */ \
		{ \
			Lep##UP lep((kin)); \
			Had##UT had((kin), (sf)); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
			up.x = name##_base_ut1((b), lep.up, had.ut); \
			up.y = name##_base_ut2((b), lep.uu, had.ut); \
		} \
		break; \
	case 7:  /* 0111 */ \
		{ \
			Lep##UP lep((kin)); \
			Had##UP had((kin), (sf)); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
			up.x = name##_base_ut1((b), lep.up, had.ut); \
			up.y = name##_base_ut2((b), lep.uu, had.ut); \
			up.z = name##_base_ul((b), lep.up, had.ul); \
		} \
		break; \
	case 8:  /* 1000 */ \
		{ \
			Lep##LU lep((kin)); \
			Had##LU had((kin), (sf)); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
			lu = name##_base_lu((b), lep.lu, had.lu); \
		} \
		break; \
	case 9:  /* 1001 */ \
		{ \
			Lep##LP lep((kin)); \
			Had##LL had((kin), (sf)); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
			up.z = name##_base_ul((b), lep.up, had.ul); \
			lu = name##_base_lu((b), lep.lu, had.lu); \
			lp.z = name##_base_ll((b), lep.lp, had.ll); \
		} \
		break; \
	case 10: /* 1010 */ \
		{ \
			Lep##LU lep((kin)); \
			Had##LT had((kin), (sf)); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
			up.y = name##_base_ut2((b), lep.uu, had.ut); \
			lu = name##_base_lu((b), lep.lu, had.lu); \
			lp.y = name##_base_lt2((b), lep.lu, had.lt); \
		} \
		break; \
	case 11: /* 1011 */ \
		{ \
			Lep##LP lep((kin)); \
			Had##LP had((kin), (sf)); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
			up.y = name##_base_ut2((b), lep.uu, had.ut); \
			up.z = name##_base_ul((b), lep.up, had.ul); \
			lu = name##_base_lu((b), lep.lu, had.lu); \
			lp.y = name##_base_lt2((b), lep.lu, had.lt); \
			lp.z = name##_base_ll((b), lep.lp, had.ll); \
		} \
		break; \
	case 12: /* 1100 */ \
		{ \
			Lep##LP lep((kin)); \
			Had##LT had((kin), (sf)); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
			up.x = name##_base_ut1((b), lep.up, had.ut); \
			lu = name##_base_lu((b), lep.lu, had.lu); \
			lp.x = name##_base_lt1((b), lep.lp, had.lt); \
		} \
		break; \
	case 13: /* 1101 */ \
		{ \
			Lep##LP lep((kin)); \
			Had##LP had((kin), (sf)); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
			up.x = name##_base_ut1((b), lep.up, had.ut); \
			up.z = name##_base_ul((b), lep.up, had.ul); \
			lu = name##_base_lu((b), lep.lu, had.lu); \
			lp.x = name##_base_lt1((b), lep.lp, had.lt); \
			lp.z = name##_base_ll((b), lep.lp, had.ll); \
		} \
		break; \
	case 14: /* 1110 */ \
		{ \
			Lep##LP lep((kin)); \
			Had##LT had((kin), (sf)); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
			up.x = name##_base_ut1((b), lep.up, had.ut); \
			up.y = name##_base_ut2((b), lep.uu, had.ut); \
			lu = name##_base_lu((b), lep.lu, had.lu); \
			lp.x = name##_base_lt1((b), lep.lp, had.lt); \
			lp.y = name##_base_lt2((b), lep.lu, had.lt); \
		} \
		break; \
	case 15: /* 1111 */ \
		{ \
			Lep##LP lep((kin)); \
			Had##LP had((kin), (sf)); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
			up.x = name##_base_ut1((b), lep.up, had.ut); \
			up.y = name##_base_ut2((b), lep.uu, had.ut); \
			up.z = name##_base_ul((b), lep.up, had.ul); \
			lu = name##_base_lu((b), lep.lu, had.lu); \
			lp.x = name##_base_lt1((b), lep.lp, had.lt); \
			lp.y = name##_base_lt2((b), lep.lu, had.lt); \
			lp.z = name##_base_ll((b), lep.lp, had.ll); \
		} \
		break; \
	} \
	return uu + dot(up, (eta)) + (lambda_e)*(lu + dot(lp, (eta))); \
}())

// Similar to `SIDIS_MACRO_XS_FROM_BASE`, except this one works with base cross-
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

// This variant of `SIDIS_MACRO_XS_FROM_BASE_P` allows for an "endpoint" set of
// structure functions to be provided, for endpoint-subtraction-related
// calculations.
#define SIDIS_MACRO_XS_FROM_BASE_P_0(name, Lep, Had, kin, sf, had_0, b, lambda_e, eta) ([&]() { \
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
			HadUU had_0_uu; \
			had_0_uu.uu = (had_0).uu; \
			Had##UU had((kin), (sf), had_0_uu); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
		} \
		break; \
	case 1:  /* 01 */ \
		{ \
			Lep##UP lep((kin)); \
			HadUP had_0_up; \
			had_0_up.uu = (had_0).uu; \
			had_0_up.ul = (had_0).ul; \
			had_0_up.ut = (had_0).ut; \
			Had##UP had((kin), (sf), had_0_up); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
			up = name##_base_up((b), lep.uu, lep.up, had.up); \
		} \
		break; \
	case 2:  /* 10 */ \
		{ \
			Lep##LU lep((kin)); \
			HadLU had_0_lu; \
			had_0_lu.uu = (had_0).uu; \
			had_0_lu.lu = (had_0).lu; \
			Had##LU had((kin), (sf), had_0_lu); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
			lu = name##_base_lu((b), lep.lu, had.lu); \
		} \
		break; \
	case 3:  /* 11 */ \
		{ \
			Lep##LP lep((kin)); \
			Had##LP had((kin), (sf), (had_0)); \
			uu = name##_base_uu((b), lep.uu, had.uu); \
			up = name##_base_up((b), lep.uu, lep.up, had.up); \
			lu = name##_base_lu((b), lep.lu, had.lu); \
			lp = name##_base_lp((b), lep.lu, lep.lp, had.lp); \
		} \
		break; \
	} \
	return uu + dot(up, (eta)) + (lambda_e)*(lu + dot(lp, (eta))); \
}())

namespace {

Real delta_vert_rad_0(Kinematics const& kin) {
	// Equation [1.3].
	Real Q_m_sq = kin.Q_sq + 2.*sq(kin.m);
	Real S_prime = kin.S - kin.Q_sq - kin.V_1;
	Real X_prime = kin.X + kin.Q_sq - kin.V_2;
	Real lambda_m = kin.Q_sq*(kin.Q_sq + 4.*sq(kin.m));
	Real lambda_S_prime = sq(S_prime) - 4.*sq(kin.m)*kin.mx_sq;
	Real lambda_X_prime = sq(X_prime) - 4.*sq(kin.m)*kin.mx_sq;
	Real lambda_m_sqrt = std::sqrt(lambda_m);
	Real lambda_S_prime_sqrt = std::sqrt(lambda_S_prime);
	Real lambda_X_prime_sqrt = std::sqrt(lambda_X_prime);

	// Differences of the form `√λ/|S| - 1`.
	Real diff_m = sqrt1p_1m((4.*sq(kin.m))/kin.Q_sq);
	Real diff_S_prime = sqrt1p_1m(-(4.*sq(kin.m)*kin.mx_sq)/sq(S_prime));
	Real diff_X_prime = sqrt1p_1m(-(4.*sq(kin.m)*kin.mx_sq)/sq(X_prime));
	Real sum_m = 2. + diff_m;
	Real sum_S_prime = 2. + diff_S_prime;
	Real sum_X_prime = 2. + diff_X_prime;

	// Equation [1.C10].
	Real L_m = 1./lambda_m_sqrt*std::log(sum_m/diff_m);
	Real L_S_prime = 1./lambda_S_prime_sqrt*std::log(-sum_S_prime/diff_S_prime);
	Real L_X_prime = 1./lambda_X_prime_sqrt*std::log(-sum_X_prime/diff_X_prime);

	// Equation [1.40].
	Real rho = 1./lambda_m_sqrt*(
		(Q_m_sq + lambda_m_sqrt)*S_prime
		- 2.*sq(kin.m)*X_prime);
	Real S_phi = Q_m_sq/lambda_m_sqrt*(
		lambda_S_prime*sq(L_S_prime)/4. - lambda_X_prime*sq(L_X_prime)/4.
		+ dilog(1. - 1./(sum_S_prime*S_prime)*rho)
		+ dilog(1. - (sum_S_prime*S_prime)/(4.*sq(kin.m)*kin.mx_sq)*rho)
		- dilog(1. - (sum_X_prime*X_prime)/(kin.mx_sq*sq(sum_m)*kin.Q_sq)*rho)
		- dilog(1. - (4.*sq(kin.m))/(sum_X_prime*X_prime*sq(sum_m)*kin.Q_sq)*rho));

	// Equation [1.52].
	Real delta = 0.5*S_prime*L_S_prime + 0.5*X_prime*L_X_prime + S_phi - 2.
		+ (1.5*kin.Q_sq + 4.*sq(kin.m))*L_m
		- Q_m_sq/lambda_m_sqrt*(
			0.5*lambda_m*sq(L_m)
			+ 2.*dilog((2.*lambda_m_sqrt)/(kin.Q_sq + lambda_m_sqrt))
			- 0.5*sq(PI));
	return delta;
}

}

Real xs::born(Kinematics const& kin, Phenom const& phenom, SfSet const& sf, Real lambda_e, Vec3 eta) {
	Born b(kin, phenom);
	return SIDIS_MACRO_XS_FROM_BASE(born, LepBorn, Had, kin, sf, b, lambda_e, eta);
}

Real xs::amm(Kinematics const& kin, Phenom const& phenom, SfSet const& sf, Real lambda_e, Vec3 eta) {
	Amm b(kin, phenom);
	return SIDIS_MACRO_XS_FROM_BASE(amm, LepAmm, Had, kin, sf, b, lambda_e, eta);
}

Real xs::nrad_ir(Kinematics const& kin, Phenom const& phenom, SfSet const& sf, Real lambda_e, Vec3 eta, Real k_0_bar) {
	Nrad b(kin, phenom, k_0_bar);
	return SIDIS_MACRO_XS_FROM_BASE(nrad_ir, LepNrad, Had, kin, sf, b, lambda_e, eta);
}

Real xs::rad(KinematicsRad const& kin, Phenom const& phenom, SfSet const& sf, Real lambda_e, Vec3 eta) {
	Rad b(kin, phenom);
	return SIDIS_MACRO_XS_FROM_BASE_P(rad, LepRad, HadRad, kin, sf, b, lambda_e, eta);
}

Real xs::rad_f(KinematicsRad const& kin, Phenom const& phenom, SfSet const& sf, Real lambda_e, Vec3 eta) {
	Rad b(kin, phenom);
	return SIDIS_MACRO_XS_FROM_BASE_P(rad_f, LepRad, HadRadF, kin, sf, b, lambda_e, eta);
}

EstErr xs::nrad_integ(Kinematics const& kin, Phenom const& phenom, SfSet const& sf, Real lambda_e, Vec3 eta, Real k_0_bar, IntegParams params) {
	// The soft part of the radiative cross-section (below `k_0_bar`) is bundled
	// into the return value here.
	Real xs_nrad_ir = nrad_ir(kin, phenom, sf, lambda_e, eta, k_0_bar);
	// TODO: The integration parameters should be modified here to account for
	// the `xs_nrad_ir` contribution.
	EstErr xs_rad_f = rad_f_integ(kin, phenom, sf, lambda_e, eta, k_0_bar, params);
	return { xs_nrad_ir + xs_rad_f.val, xs_rad_f.err };
}

EstErr xs::rad_f_integ(Kinematics const& kin, Phenom const& phenom, SfSet const& sf, Real lambda_e, Vec3 eta, Real k_0_bar, IntegParams params) {
	HadLP had_0(kin, sf);
	CutRad cut;
	cut.k_0_bar = Bound(0., k_0_bar);
	EstErr xs_integ = integrate<3>(
		[&](std::array<Real, 3> x) {
			KinematicsRad kin_rad;
			Real jac;
			if (!take(cut, kin, x.data(), &kin_rad, &jac)) {
				return 0.;
			}
			Rad b(kin_rad, phenom);
			Real xs = jac * SIDIS_MACRO_XS_FROM_BASE_P_0(rad_f, LepRad, HadRadF, kin_rad, sf, had_0, b, lambda_e, eta);
			if (std::isnan(xs)) {
				// If the result is `NaN`, it most likely means we went out of
				// the allowed region for the structure function grids (or we
				// are in a kinematically disallowed region). In that case, just
				// return zero.
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

EstErr xs::rad_integ(Kinematics const& kin, Phenom const& phenom, SfSet const& sf, Real lambda_e, Vec3 eta, Real k_0_bar, IntegParams params) {
	CutRad cut;
	cut.k_0_bar = Bound(k_0_bar, INF);
	EstErr xs_integ = integrate<3>(
		[&](std::array<Real, 3> x) {
			KinematicsRad kin_rad;
			Real jac;
			if (!take(cut, kin, x.data(), &kin_rad, &jac)) {
				return 0.;
			}
			Rad b(kin_rad, phenom);
			Real xs = jac * SIDIS_MACRO_XS_FROM_BASE_P(rad, LepRad, HadRad, kin_rad, sf, b, lambda_e, eta);
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

Real xs::born(Kinematics const& kin, SfSet const& sf, Real lambda_e, Vec3 eta) {
	return born(kin, Phenom(kin), sf, lambda_e, eta);
}
Real xs::amm(Kinematics const& kin, SfSet const& sf, Real lambda_e, Vec3 eta) {
	return amm(kin, Phenom(kin), sf, lambda_e, eta);
}
Real xs::nrad_ir(Kinematics const& kin, SfSet const& sf, Real lambda_e, Vec3 eta, Real k_0_bar) {
	return nrad_ir(kin, Phenom(kin), sf, lambda_e, eta, k_0_bar);
}
Real xs::rad(KinematicsRad const& kin, SfSet const& sf, Real lambda_e, Vec3 eta) {
	return rad(kin, Phenom(kin.project()), sf, lambda_e, eta);
}
Real xs::rad_f(KinematicsRad const& kin, SfSet const& sf, Real lambda_e, Vec3 eta) {
	return rad_f(kin, Phenom(kin.project()), sf, lambda_e, eta);
}
EstErr xs::nrad_integ(Kinematics const& kin, SfSet const& sf, Real lambda_e, Vec3 eta, Real k_0_bar, IntegParams params) {
	return nrad_integ(kin, Phenom(kin), sf, lambda_e, eta, k_0_bar, params);
}
EstErr xs::rad_f_integ(Kinematics const& kin, SfSet const& sf, Real lambda_e, Vec3 eta, Real k_0_bar, IntegParams params) {
	return rad_f_integ(kin, Phenom(kin), sf, lambda_e, eta, k_0_bar, params);
}
EstErr xs::rad_integ(Kinematics const& kin, SfSet const& sf, Real lambda_e, Vec3 eta, Real k_0_bar, IntegParams params) {
	return rad_integ(kin, Phenom(kin), sf, lambda_e, eta, k_0_bar, params);
}

// Radiative corrections to Born cross-section.
Real xs::delta_vert_rad_ir(Kinematics const& kin, Real k_0_bar) {
	// Paragraph following equation [1.C17].
	Real k0_max = (kin.mx_sq - sq(kin.Mth))/(2.*kin.mx);
	if (!(k_0_bar > 0.)) {
		return -INF;
	}
	Real Q_m_sq = kin.Q_sq + 2.*sq(kin.m);
	Real diff_m = sqrt1p_1m((4.*sq(kin.m))/kin.Q_sq);
	Real sum_m = 2. + diff_m;
	Real lambda_m = kin.Q_sq*(kin.Q_sq + 4.*sq(kin.m));
	Real lambda_m_sqrt = std::sqrt(lambda_m);
	Real L_m = 1./lambda_m_sqrt*std::log(sum_m/diff_m);
	Real delta_0 = delta_vert_rad_0(kin);
	// This comes from subtracting `delta_H` (equation [1.38]) from `delta_VR`
	// (equation [1.52]).
	Real delta_shift = 2.*(Q_m_sq*L_m - 1.)*std::log(
		k_0_bar < k0_max ?
		(2.*k_0_bar)/kin.m :
		(kin.mx_sq - sq(kin.Mth))/(kin.m*kin.mx));
	return delta_0 + delta_shift;
}
Real xs::delta_rad_ir_hard(Kinematics const& kin, Real k_0_bar) {
	Real k0_max = (kin.mx_sq - sq(kin.Mth))/(2.*kin.mx);
	if (!(k_0_bar > 0.)) {
		return INF;
	} else if (!(k_0_bar < k0_max)) {
		return 0.;
	}
	Real Q_m_sq = kin.Q_sq + 2.*sq(kin.m);
	Real diff_m = sqrt1p_1m((4.*sq(kin.m))/kin.Q_sq);
	Real sum_m = 2. + diff_m;
	Real lambda_m = kin.Q_sq*(kin.Q_sq + 4.*sq(kin.m));
	Real lambda_m_sqrt = std::sqrt(lambda_m);
	Real L_m = 1./lambda_m_sqrt*std::log(sum_m/diff_m);
	// Equation [1.38].
	Real delta = 2.*(Q_m_sq*L_m - 1.)*std::log(
		(kin.mx_sq - sq(kin.Mth))/(2.*k_0_bar*kin.mx));
	return delta;
}

Real xs::delta_vac_lep(Kinematics const& kin) {
	// Equation [1.50].
	Real ms[3] = { MASS_E, MASS_MU, MASS_TAU };
	Real delta = 0.;
	for (unsigned idx = 0; idx < 3; ++idx) {
		Real m = ms[idx];
		Real lambda_sqrt = std::sqrt(kin.Q_sq*(kin.Q_sq + 4.*sq(m)));
		Real diff_m = sqrt1p_1m((4.*sq(m))/kin.Q_sq);
		Real sum_m = 2. + diff_m;
		Real L_m = 1./lambda_sqrt*std::log(sum_m/diff_m);
		delta += 2./3.L*(kin.Q_sq
			+ 2.*sq(m))*L_m
			- 10./9.L
			+ (8.*sq(m))/(3.*kin.Q_sq)*(1. - 2.*sq(m)*L_m);
	}
	return delta;
}

// Born base functions.
Born::Born(Kinematics const& kin, Phenom const& phenom) :
	// Equation [1.15]. The `Q^4` factor has been absorbed into `C_1`.
	coeff((sq(phenom.alpha_qed)*kin.S*sq(kin.S_x))/(8.*kin.M*kin.ph_l*kin.lambda_S)) { }

Real xs::born_base_uu(Born const& b, LepBornBaseUU const& lep, HadBaseUU const& had) {
	return b.coeff*(
		lep.theta_1*had.H_10
		+ lep.theta_2*had.H_20
		+ lep.theta_3*had.H_30
		+ lep.theta_4*had.H_40);
}
Real xs::born_base_ul(Born const& b, LepBornBaseUP const& lep, HadBaseUL const& had) {
	return b.coeff*(lep.theta_6*had.H_63 + lep.theta_8*had.H_83);
}
Real xs::born_base_ut1(Born const& b, LepBornBaseUP const& lep, HadBaseUT const& had) {
	return b.coeff*(lep.theta_6*had.H_61 + lep.theta_8*had.H_81);
}
Real xs::born_base_ut2(Born const& b, LepBornBaseUU const& lep, HadBaseUT const& had) {
	return b.coeff*(
		lep.theta_1*had.H_12
		+ lep.theta_2*had.H_22
		+ lep.theta_3*had.H_32
		+ lep.theta_4*had.H_42);
}
Real xs::born_base_lu(Born const& b, LepBornBaseLU const& lep, HadBaseLU const& had) {
	return b.coeff*lep.theta_5*had.H_50;
}
Real xs::born_base_ll(Born const& b, LepBornBaseLP const& lep, HadBaseLL const& had) {
	return b.coeff*(lep.theta_7*had.H_73 + lep.theta_9*had.H_93);
}
Real xs::born_base_lt1(Born const& b, LepBornBaseLP const& lep, HadBaseLT const& had) {
	return b.coeff*(lep.theta_7*had.H_71 + lep.theta_9*had.H_91);
}
Real xs::born_base_lt2(Born const& b, LepBornBaseLU const& lep, HadBaseLT const& had) {
	return b.coeff*lep.theta_5*had.H_52;
}

// AMM base functions.
Amm::Amm(Kinematics const& kin, Phenom const& phenom) {
	// Equation [1.53]. The `Q^4` factor has been absorbed into `C_1`.
	Real lambda_m = kin.Q_sq*(kin.Q_sq + 4.*sq(kin.m));
	Real lambda_m_sqrt = std::sqrt(lambda_m);
	Real diff_m = sqrt1p_1m((4.*sq(kin.m))/kin.Q_sq);
	Real sum_m = 2. + diff_m;
	Real L_m = 1./lambda_m_sqrt*std::log(sum_m/diff_m);
	coeff = L_m*kin.Q_sq*(std::pow(phenom.alpha_qed, 3)*sq(kin.m)*kin.S*sq(kin.S_x))
		/(16.*PI*kin.M*kin.ph_l*kin.lambda_S);
}

Real xs::amm_base_uu(Amm const& b, LepAmmBaseUU const& lep, HadBaseUU const& had) {
	return b.coeff*(
		lep.theta_1*had.H_10
		+ lep.theta_2*had.H_20
		+ lep.theta_3*had.H_30
		+ lep.theta_4*had.H_40);
}
Real xs::amm_base_ul(Amm const& b, LepAmmBaseUP const& lep, HadBaseUL const& had) {
	return b.coeff*(lep.theta_6*had.H_63 + lep.theta_8*had.H_83);
}
Real xs::amm_base_ut1(Amm const& b, LepAmmBaseUP const& lep, HadBaseUT const& had) {
	return b.coeff*(lep.theta_6*had.H_61 + lep.theta_8*had.H_81);
}
Real xs::amm_base_ut2(Amm const& b, LepAmmBaseUU const& lep, HadBaseUT const& had) {
	return b.coeff*(
		lep.theta_1*had.H_12
		+ lep.theta_2*had.H_22
		+ lep.theta_3*had.H_32
		+ lep.theta_4*had.H_42);
}
Real xs::amm_base_lu(Amm const& b, LepAmmBaseLU const& lep, HadBaseLU const& had) {
	return b.coeff*lep.theta_5*had.H_50;
}
Real xs::amm_base_ll(Amm const& b, LepAmmBaseLP const& lep, HadBaseLL const& had) {
	return b.coeff*(lep.theta_7*had.H_73 + lep.theta_9*had.H_93);
}
Real xs::amm_base_lt1(Amm const& b, LepAmmBaseLP const& lep, HadBaseLT const& had) {
	return b.coeff*(lep.theta_7*had.H_71 + lep.theta_9*had.H_91);
}
Real xs::amm_base_lt2(Amm const& b, LepAmmBaseLU const& lep, HadBaseLT const& had) {
	return b.coeff*lep.theta_5*had.H_52;
}

// Non-radiative infrared-divergence-free base functions.
Nrad::Nrad(Kinematics const& kin, Phenom const& phenom, Real k_0_bar) {
	Born born(kin, phenom);
	Amm amm(kin, phenom);
	Real born_factor = 1. + phenom.alpha_qed/PI*(
		delta_vert_rad_ir(kin, k_0_bar)
		+ delta_vac_lep(kin)
		+ phenom.delta_vac_had);
	coeff_born = born_factor*born.coeff;
	coeff_amm = amm.coeff;
}

Real xs::nrad_ir_base_uu(Nrad const& b, LepNradBaseUU const& lep, HadBaseUU const& had) {
	return
		(b.coeff_born*lep.born.theta_1 + b.coeff_amm*lep.amm.theta_1)*had.H_10
		+ (b.coeff_born*lep.born.theta_2 + b.coeff_amm*lep.amm.theta_2)*had.H_20
		+ (b.coeff_born*lep.born.theta_3 + b.coeff_amm*lep.amm.theta_3)*had.H_30
		+ (b.coeff_born*lep.born.theta_4 + b.coeff_amm*lep.amm.theta_4)*had.H_40;
}
Real xs::nrad_ir_base_ul(Nrad const& b, LepNradBaseUP const& lep, HadBaseUL const& had) {
	return
		(b.coeff_born*lep.born.theta_6 + b.coeff_amm*lep.amm.theta_6)*had.H_63
		+ (b.coeff_born*lep.born.theta_8 + b.coeff_amm*lep.amm.theta_8)*had.H_83;
}
Real xs::nrad_ir_base_ut1(Nrad const& b, LepNradBaseUP const& lep, HadBaseUT const& had) {
	return
		(b.coeff_born*lep.born.theta_6 + b.coeff_amm*lep.amm.theta_6)*had.H_61
		+ (b.coeff_born*lep.born.theta_8 + b.coeff_amm*lep.amm.theta_8)*had.H_81;
}
Real xs::nrad_ir_base_ut2(Nrad const& b, LepNradBaseUU const& lep, HadBaseUT const& had) {
	return
		(b.coeff_born*lep.born.theta_1 + b.coeff_amm*lep.amm.theta_1)*had.H_12
		+ (b.coeff_born*lep.born.theta_2 + b.coeff_amm*lep.amm.theta_2)*had.H_22
		+ (b.coeff_born*lep.born.theta_3 + b.coeff_amm*lep.amm.theta_3)*had.H_32
		+ (b.coeff_born*lep.born.theta_4 + b.coeff_amm*lep.amm.theta_4)*had.H_42;
}
Real xs::nrad_ir_base_lu(Nrad const& b, LepNradBaseLU const& lep, HadBaseLU const& had) {
	return (b.coeff_born*lep.born.theta_5 + b.coeff_amm*lep.amm.theta_5)*had.H_50;
}
Real xs::nrad_ir_base_ll(Nrad const& b, LepNradBaseLP const& lep, HadBaseLL const& had) {
	return
		(b.coeff_born*lep.born.theta_7 + b.coeff_amm*lep.amm.theta_7)*had.H_73
		+ (b.coeff_born*lep.born.theta_9 + b.coeff_amm*lep.amm.theta_9)*had.H_93;
}
Real xs::nrad_ir_base_lt1(Nrad const& b, LepNradBaseLP const& lep, HadBaseLT const& had) {
	return
		(b.coeff_born*lep.born.theta_7 + b.coeff_amm*lep.amm.theta_7)*had.H_71
		+ (b.coeff_born*lep.born.theta_9 + b.coeff_amm*lep.amm.theta_9)*had.H_91;
}
Real xs::nrad_ir_base_lt2(Nrad const& b, LepNradBaseLU const& lep, HadBaseLT const& had) {
	return (b.coeff_born*lep.born.theta_5 + b.coeff_amm*lep.amm.theta_5)*had.H_52;
}

// Radiative base functions.
Rad::Rad(KinematicsRad const& kin, Phenom const& phenom) {
	// Equation [1.43].
	coeff = -(std::pow(phenom.alpha_qed, 3)*kin.S*sq(kin.S_x))
		/(64.*sq(PI)*kin.M*kin.ph_l*kin.lambda_S*kin.lambda_Y_sqrt);
	R = kin.R;
}

Real xs::rad_base_uu(Rad const& b, LepRadBaseUU const& lep, HadRadBaseUU const& had) {
	return b.coeff*(
		1./b.R*(
			lep.theta_011*had.H_10
			+ lep.theta_021*had.H_20
			+ lep.theta_031*had.H_30
			+ lep.theta_041*had.H_40)
		+ (
			lep.theta_012*had.H_10
			+ lep.theta_022*had.H_20
			+ lep.theta_032*had.H_30
			+ lep.theta_042*had.H_40)
		+ b.R*(
			lep.theta_013*had.H_10
			+ lep.theta_023*had.H_20
			+ lep.theta_033*had.H_30
			+ lep.theta_043*had.H_40));
}
Vec3 xs::rad_base_up(Rad const& b, LepRadBaseUU const& lep_uu, LepRadBaseUP const& lep_up, HadRadBaseUP const& had) {
	return b.coeff*(
		1./b.R*(
			lep_uu.theta_011*had.H_1
			+ lep_uu.theta_021*had.H_2
			+ lep_uu.theta_031*had.H_3
			+ lep_uu.theta_041*had.H_4
			+ lep_up.theta_061*had.H_6
			+ lep_up.theta_081*had.H_8)
		+ (
			lep_uu.theta_012*had.H_1
			+ lep_uu.theta_022*had.H_2
			+ lep_uu.theta_032*had.H_3
			+ lep_uu.theta_042*had.H_4
			+ lep_up.theta_062*had.H_6
			+ lep_up.theta_082*had.H_8)
		+ b.R*(
			lep_uu.theta_013*had.H_1
			+ lep_uu.theta_023*had.H_2
			+ lep_uu.theta_033*had.H_3
			+ lep_uu.theta_043*had.H_4
			+ lep_up.theta_063*had.H_6
			+ lep_up.theta_083*had.H_8)
		+ b.R*b.R*(
			lep_up.theta_064*had.H_6
			+ lep_up.theta_084*had.H_8));
}
Real xs::rad_base_lu(Rad const& b, LepRadBaseLU const& lep, HadRadBaseLU const& had) {
	return b.coeff*(
		1./b.R*(lep.theta_051 + lep.theta_151)*had.H_50
		+ (lep.theta_052 + lep.theta_152)*had.H_50
		+ b.R*(lep.theta_053 + lep.theta_153)*had.H_50);
}
Vec3 xs::rad_base_lp(Rad const& b, LepRadBaseLU const& lep_lu, LepRadBaseLP const& lep_lp, HadRadBaseLP const& had) {
	return b.coeff*(
		1./b.R*(
			(lep_lu.theta_051 + lep_lu.theta_151)*had.H_5
			+ (lep_lp.theta_071 + lep_lp.theta_171)*had.H_7
			+ (lep_lp.theta_091 + lep_lp.theta_191)*had.H_9)
		+ (
		 	(lep_lu.theta_052 + lep_lu.theta_152)*had.H_5
			+ (lep_lp.theta_072 + lep_lp.theta_172)*had.H_7
			+ (lep_lp.theta_092 + lep_lp.theta_192)*had.H_9)
		+ b.R*(
		 	(lep_lu.theta_053 + lep_lu.theta_153)*had.H_5
			+ (lep_lp.theta_073 + lep_lp.theta_173)*had.H_7
			+ (lep_lp.theta_093 + lep_lp.theta_193)*had.H_9)
		+ b.R*b.R*(
			(lep_lp.theta_074 + lep_lp.theta_174)*had.H_7
			+ (lep_lp.theta_094 + lep_lp.theta_194)*had.H_9));
}

Real xs::rad_f_base_uu(Rad const& b, LepRadBaseUU const& lep, HadRadFBaseUU const& had) {
	return b.coeff*(
		(
			lep.theta_011*had.H_10_diff
			+ lep.theta_021*had.H_20_diff
			+ lep.theta_031*had.H_30_diff
			+ lep.theta_041*had.H_40_diff)
		+ (
			lep.theta_012*had.H_10
			+ lep.theta_022*had.H_20
			+ lep.theta_032*had.H_30
			+ lep.theta_042*had.H_40)
		+ b.R*(
			lep.theta_013*had.H_10
			+ lep.theta_023*had.H_20
			+ lep.theta_033*had.H_30
			+ lep.theta_043*had.H_40));
}
Vec3 xs::rad_f_base_up(Rad const& b, LepRadBaseUU const& lep_uu, LepRadBaseUP const& lep_up, HadRadFBaseUP const& had) {
	return b.coeff*(
		(
			lep_uu.theta_011*had.H_1_diff
			+ lep_uu.theta_021*had.H_2_diff
			+ lep_uu.theta_031*had.H_3_diff
			+ lep_uu.theta_041*had.H_4_diff
			+ lep_up.theta_061*had.H_6_diff
			+ lep_up.theta_081*had.H_8_diff)
		+ (
			lep_uu.theta_012*had.H_1
			+ lep_uu.theta_022*had.H_2
			+ lep_uu.theta_032*had.H_3
			+ lep_uu.theta_042*had.H_4
			+ lep_up.theta_062*had.H_6
			+ lep_up.theta_082*had.H_8)
		+ b.R*(
			lep_uu.theta_013*had.H_1
			+ lep_uu.theta_023*had.H_2
			+ lep_uu.theta_033*had.H_3
			+ lep_uu.theta_043*had.H_4
			+ lep_up.theta_063*had.H_6
			+ lep_up.theta_083*had.H_8)
		+ b.R*b.R*(
			lep_up.theta_064*had.H_6
			+ lep_up.theta_084*had.H_8));
}
Real xs::rad_f_base_lu(Rad const& b, LepRadBaseLU const& lep, HadRadFBaseLU const& had) {
	return b.coeff*(
		(lep.theta_051 + lep.theta_151)*had.H_50_diff
		+ (lep.theta_052 + lep.theta_152)*had.H_50
		+ b.R*(lep.theta_053 + lep.theta_153)*had.H_50);
}
Vec3 xs::rad_f_base_lp(Rad const& b, LepRadBaseLU const& lep_lu, LepRadBaseLP const& lep_lp, HadRadFBaseLP const& had) {
	return b.coeff*(
		(
			(lep_lu.theta_051 + lep_lu.theta_151)*had.H_5_diff
			+ (lep_lp.theta_071 + lep_lp.theta_171)*had.H_7_diff
			+ (lep_lp.theta_091 + lep_lp.theta_191)*had.H_9_diff)
		+ (
		 	(lep_lu.theta_052 + lep_lu.theta_152)*had.H_5
			+ (lep_lp.theta_072 + lep_lp.theta_172)*had.H_7
			+ (lep_lp.theta_092 + lep_lp.theta_192)*had.H_9)
		+ b.R*(
		 	(lep_lu.theta_053 + lep_lu.theta_153)*had.H_5
			+ (lep_lp.theta_073 + lep_lp.theta_173)*had.H_7
			+ (lep_lp.theta_093 + lep_lp.theta_193)*had.H_9)
		+ b.R*b.R*(
			(lep_lp.theta_074 + lep_lp.theta_174)*had.H_7
			+ (lep_lp.theta_094 + lep_lp.theta_194)*had.H_9));
}

//EstErr xs::exc_integ(Kinematics const& kin, SfSet const& sf, Real lambda_e, Vec3 eta, Real k_0_bar, IntegParams params) {
//	return exc_integ(kin, Phenom(kin), sf, lambda_e, eta);
//}
Real xs::EXC_XS_FROM_BASE_P(KinematicsRad kin,Real lambda_e,Vec3& eta) 
{
    EXC b(kin);
    // Create a mask describing the polarization state.
    unsigned pol_mask = (((lambda_e) != 0.) << 1) 
                      | (((eta.x != 0. || eta.y != 0. || eta.z != 0.) << 0));

    Real uu = 0.;
    Vec3 up = VEC3_ZERO;
    Real lu = 0.;
    Vec3 lp = VEC3_ZERO;

    switch (pol_mask) {
        case 0:  // 00
            {
                LepRadUU leprad(kin);
                EXCUU excuu(kin);
                uu = exc_base_uu(b, leprad.uu, excuu);
            }
            break;
        case 1:  // 01
            {
                LepRadUP leprad(kin);
                EXCUU excuu(kin);
                EXCUT excut(kin);
		EXCUL excul(kin);
                uu = exc_base_uu(b, leprad.uu, excuu);
                up = exc_base_up(b, leprad.uu, leprad.up, excut,excul);
            }
            break;
        case 2:  // 10
            {
                LepRadLU leprad(kin);
                EXCUU excuu(kin);
                EXCLU exclu(kin);
                uu = exc_base_uu(b, leprad.uu, excuu);
                lu = exc_base_lu(b, leprad.lu, exclu);
            }
            break;
        case 3:  // 11
            {
                LepRadLP leprad(kin);
                EXCUU excuu(kin);
                EXCUT excut(kin);
                EXCUL excul(kin);
                EXCLU exclu(kin);
                EXCLT exclt(kin);
                EXCLL excll(kin);
                uu = exc_base_uu(b, leprad.uu, excuu);
                up = exc_base_up(b, leprad.uu, leprad.up, excut,excul);
                lu = exc_base_lu(b, leprad.lu, exclu);
                lp = exc_base_lp(b, leprad.lu, leprad.lp, exclt,excll);
            }
            break;
    }

    return uu + dot(up, eta) + (lambda_e)*(lu + dot(lp, eta));
}

//Real integrate2D(xs::EXC_XS_FROM_BASE_P func, 
//                   double x1min, double x1max, double x2min, double x2max, 
//                   int nx1, int nx2) {
//    double hx1 = (x1max - x1min) / nx1;
//    double hx2 = (x2max - x2min) / nx2;
//    Real integral = 0.0;
//
//    for (int i = 0; i < nx1; ++i) {
//        
//        for (int j = 0; j < nx2; ++j) {
//            integral += func * hx1 * hx2;
//        }
//    }
//
//    return integral;
//}

//Define the exc integrand function
//int integrand(unsigned ndim, const double *x, void *data, unsigned fdim, double *fval){
//	KinematicsRad *kinrad = static_cast<KinematicsRad*>(data);
//	KinematicsRad kin_temp(*kinrad,x[0],x[1],1.0);
//	fval[0] = EXC_XS_FROM_BASE_P(kin_temp,lambda_e,eta);
//	return 0;
//}

Real xs::exc_integ(Kinematics const& kin, Real lambda_e, Vec3 eta) {
	/*
	std::cout<<"check taus and taup Q2: "<<kin.Q_sq<<" S: "<<kin.S<<" X: "<<kin.X<<std::endl;
	double taus = -kin.Q_sq/kin.S;
	double taup = kin.Q_sq/kin.X;
	std::cout<<" check taus and taup taus: "<<taus<<" taup "<<taup<<std::endl;
	std::cout<<" check taumin:: "<<kin.tau_min<<" taumax "<<kin.tau_max<<std::endl;
	constexpr std::size_t D = 2;	
	auto func = [&](std::array<Real, 2> x) -> Real {
		KinematicsRad kinrad(kin, x[0], x[1], 1.0);
		Real xs = EXC_XS_FROM_BASE_P(kinrad, lambda_e, eta);
		return std::isnan(xs) ? 0.0 : xs;
	};
	// Define the GSL-compatible integrand function
	auto gsl_integrand = [](double* xs, size_t dim, void* data) -> double {
		auto* func_ptr = static_cast<decltype(func)*>(data);  // Cast void* to lambda pointer
		std::array<double, D> point;
		std::move(xs, xs + D, point.begin());  // Convert raw input to std::array
		return (*func_ptr)(point);  // Call the lambda function
	};
	const gsl_rng_type* T;
	gsl_rng* r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_monte_vegas_state* s = gsl_monte_vegas_alloc(D);

	// Prepare GSL function structure
	gsl_monte_function func_gsl = {gsl_integrand, D, &func};

	double xl1[D] = {kin.tau_min,0};
	double xu1[D] = {taup+(taup-taus)/2,0.2};
	double xl2[D] = {kin.tau_min,6};
	double xu2[D] = {taup+(taup-taus)/2,2*PI};
        	
	double result1, error1, result2, error2;
	size_t calls = 10000;

	// Integrate over the first region
	gsl_monte_vegas_integrate(&func_gsl, xl1, xu1, D, 1000, r, s, &result1, &error1);  // Warm-up
	gsl_monte_vegas_integrate(&func_gsl, xl1, xu1, D, calls, r, s, &result1, &error1);
	std::cout << "Integral near (x1, x2): " << result1 << " ± " << error1 << std::endl;

	// Reset the state and integrate over the second region
	// Reset the state by freeing and reallocating
	gsl_monte_vegas_free(s);
	s = gsl_monte_vegas_alloc(D);

	gsl_monte_vegas_integrate(&func_gsl, xl2, xu2, D, 1000, r, s, &result2, &error2);  // Warm-up
	gsl_monte_vegas_integrate(&func_gsl, xl2, xu2, D, calls, r, s, &result2, &error2);
	std::cout << "Integral near (x3, x4): " << result2 << " ± " << error2 << std::endl;

	// Combine results from both regions
	double total_result = result1 + result2;
	double total_error = std::sqrt(error1 * error1 + error2 * error2);

	std::cout << "Total estimated integral: " << total_result << std::endl;
	std::cout << "Total estimated error: " << total_error << std::endl;

	// Clean up
	gsl_monte_vegas_free(s);
	gsl_rng_free(r);

	return total_result;
	*/
	/*
	constexpr std::size_t D = 2;	
	auto func = [&](std::array<Real, 2> x) -> Real {
		KinematicsRad kinrad(kin, x[0], x[1], 1.0);
		Real xs = EXC_XS_FROM_BASE_P(kinrad, lambda_e, eta);
		return std::isnan(xs) ? 0.0 : xs;
	};
	// Define the GSL-compatible integrand function
	auto gsl_integrand = [](double* xs, size_t dim, void* data) -> double {
		auto* func_ptr = static_cast<decltype(func)*>(data);  // Cast void* to lambda pointer
		std::array<double, D> point;
		std::move(xs, xs + D, point.begin());  // Convert raw input to std::array
		return (*func_ptr)(point);  // Call the lambda function
	};

	// Define the integration bounds
	double xl[D] = {kin.tau_min, 0.0};  // Lower bounds
	double xu[D] = {kin.tau_max, 2 * M_PI};  // Upper bounds

	// Initialize the VEGAS state
	gsl_monte_vegas_state* s = gsl_monte_vegas_alloc(D);

	// Set up the GSL Monte Carlo function structure
	gsl_monte_function func_gsl = {gsl_integrand, D, &func};

	// Initialize the random number generator
	const gsl_rng_type* T;
	gsl_rng* r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);

	// Perform the VEGAS integration
	double result, error;
	size_t initial_calls = 1000;  // Initial sampling points for grid setup
	size_t refinement_calls = 50000;  // Refined sampling points for accurate result

	// First burn-in iteration to establish the sampling grid
	gsl_monte_vegas_integrate(&func_gsl, xl, xu, D, initial_calls, r, s, &result, &error);
	std::cout << "First burn-in: " << result << " +/- " << error << std::endl;
	std::cout << "Chi-squared: " << gsl_monte_vegas_chisq(s) << std::endl;

	// Optional: Refine the grid with multiple iterations
	for (int i = 0; i < 3; ++i) {
		std::cout << "Refining grid, iteration " << i + 1 << "..." << std::endl;
		gsl_monte_vegas_integrate(&func_gsl, xl, xu, D, refinement_calls, r, s, &result, &error);
		std::cout << "Iteration " << i + 1 << ": " << result << " +/- " << error << std::endl;
		std::cout << "Chi-squared: " << gsl_monte_vegas_chisq(s) << std::endl;

		// Check if chi-squared is close to 1.0, indicating convergence
		if (std::abs(gsl_monte_vegas_chisq(s) - 1.0) < 0.1) {
			std::cout << "Converged!" << std::endl;
			break;
		}
	}

	// Output final result and error
	std::cout << "Final result: " << result << " +/- " << error << std::endl;

	// Clean up
	gsl_monte_vegas_free(s);
	gsl_rng_free(r);

	return result;
	*/
	//Vegas
	constexpr std::size_t D = 2;	
	auto func = [&](std::array<Real, 2> x) -> Real {
		KinematicsRad kinrad(kin, x[0], x[1], 1.0);
		Real xs = EXC_XS_FROM_BASE_P(kinrad, lambda_e, eta);
		return std::isnan(xs) ? 0.0 : xs;
	};

	// Define the GSL-compatible integrand function
	auto gsl_integrand = [](double* xs, size_t dim, void* data) -> double {
		auto* func_ptr = static_cast<decltype(func)*>(data);  // Cast void* to lambda pointer
		std::array<double, D> point;
		std::move(xs, xs + D, point.begin());  // Convert raw input to std::array
		return (*func_ptr)(point);  // Call the lambda function
	};

	// Define the integration bounds
	double xl[D] = {kin.tau_min, -PI};  // Lower bounds
	double xu[D] = {kin.tau_max, PI};  // Upper bounds

	// Initialize the VEGAS state
	gsl_monte_vegas_state* s = gsl_monte_vegas_alloc(D);

	// Set up the GSL Monte Carlo function structure
	gsl_monte_function func_gsl = {gsl_integrand, D, &func};

	// Initialize the random number generator
	const gsl_rng_type* T;
	gsl_rng* r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);

	// Perform the VEGAS integration
	double result, error;
	size_t calls = 10000;  // Number of sampling points

	// First iteration to set up the sampling grid
	gsl_monte_vegas_integrate(&func_gsl, xl, xu, D, 1000, r, s, &result, &error);
        std::cout<<" First burn "<<result<<std::endl;
	// Optional: Refine the grid with a second iteration
	gsl_monte_vegas_integrate(&func_gsl, xl, xu, D, calls, r, s, &result, &error);

	// Output the result and error
	std::cout << "Estimated integral: " << result << std::endl;
	std::cout << "Estimated error: " << error << std::endl;

	// Clean up
	gsl_monte_vegas_free(s);
	gsl_rng_free(r);

	return result;
	/*
	// Define bounds using std::array, but correctly pass them to cubature::Point
	cubature::Point<2, double> lower = {kin.tau_min, -PI};
	cubature::Point<2, double> upper = {kin.tau_max,  PI};

	// Output variables for the integral and error
	double integral[1] = {0.0};
	double error[1] = {0.0};

	// Perform the cubature integration with the correct number of arguments
	cubature::EstErr<Real> result = cubature::cubature<2>(
		[&](std::array<Real, 2> x) -> Real {
			KinematicsRad kinrad(kin, x[0], x[1], 1.0);
			Real xs = EXC_XS_FROM_BASE_P(kinrad, lambda_e, eta);
			return std::isnan(xs) ? 0.0 : xs;
		},
		lower, upper,  // Bounds
		10000, 0.0, 0.0  // Max evaluations and tolerances
		);

	// Check if the integration was successful
	if (result.err==0) {
		std::cerr << "Cubature failed with error: " << result.err << std::endl;
		return 0.0;  // Return 0 on failure
	}

	std::cout<<"Check Cubature results "<<integral[0]<<std::endl;
	// Return the computed integral value
	return integral[0];
	*/
	/*	
	//Get xs for each tau phi bins
	double taumax = kin.tau_max;
	double taumin = kin.tau_min;

	int ntau = 100;
	int nphi = 100;
	double htau = (taumax-taumin)/ntau;
	double hphi = (2*PI)/nphi;
	Real integral = 0.0;
	std::ofstream file("check_dtau_dphi.txt");
	for (int i = 0; i< ntau; ++i){
		double dtau = taumin+(i+0.5)*htau;
		for (int j=0; j<nphi;++j){
			double dphi=0+(j+0.5)*hphi;
			KinematicsRad kinrad(kin,dtau,dphi,1.0);//Kinrad(kin,dtau,dphi,1.0_ don't worry, the last number is R, which for exc case, I hard coded Rex, not using this 1.0.
			Real xs =  EXC_XS_FROM_BASE_P(kinrad, lambda_e, eta);
			
			file<<dtau<<" "<<dphi<<" "<<xs<<std::endl;;
			integral+=xs*dtau*dphi;
			std::cout<<" check integ j"<<j<<std::endl;
			std::cout<<" check integ"<<xs<<std::endl;
			std::cout<<" check W2 "<<(kinrad.shift_W_sq)<<" Q2 "<<kinrad.shift_Q_sq<<std::endl;
			std::cout<<dtau<<" "<<dphi<<" "<<xs<<std::endl;
		}
		std::cout<<"check integ i "<<i<<std::endl;
	}
	file.close();
	return integral;
	*/
}

// Exclusive radiative base functions.
EXC::EXC(KinematicsRad const& kin) {
	// Equation [36].
	coeff = -(std::pow(ALPHA, 3)*kin.S*sq(kin.S_x))
		/(512.*std::pow(PI,5)*kin.M*kin.ph_l*kin.lambda_S*kin.lambda_Y_sqrt*(1-kin.tau-kin.mu)*sq(kin.shiftexc_Q_sq));
	Rex = kin.Rex;
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
			lep_uu.theta_011*excut.H1_010
			+lep_uu.theta_021*excut.H2_010
			+lep_uu.theta_031*excut.H3_010
			+lep_uu.theta_041*excut.H4_010
			+lep_up.theta_061*excut.H6_100
			+lep_up.theta_061*excul.H6_001
			+lep_up.theta_081*excut.H8_100
			+lep_up.theta_081*excul.H8_001
			)
		+(
			lep_uu.theta_012*excut.H1_010
			+lep_uu.theta_022*excut.H2_010
			+lep_uu.theta_032*excut.H3_010
			+lep_uu.theta_042*excut.H4_010
			+lep_up.theta_062*excut.H6_100
			+lep_up.theta_062*excul.H6_001
			+lep_up.theta_082*excut.H8_100
			+lep_up.theta_082*excul.H8_001
			)
		+ b.Rex*(
			lep_uu.theta_013*excut.H1_010
			+lep_uu.theta_023*excut.H2_010
			+lep_uu.theta_033*excut.H3_010
			+lep_uu.theta_043*excut.H4_010
			+lep_up.theta_063*excut.H6_100
			+lep_up.theta_063*excul.H6_001
			+lep_up.theta_083*excut.H8_100
			+lep_up.theta_083*excul.H8_001
			)
		+ b.Rex*b.Rex*(
                        lep_up.theta_064*excut.H6_100
                        +lep_up.theta_064*excul.H6_001
                        +lep_up.theta_084*excut.H8_100
                        +lep_up.theta_084*excul.H8_001
			)
		);
}
Real xs::exc_base_lu(EXC const& b, LepRadBaseLU const& lep, EXCLU const& exc) {
	return b.coeff*(
		1./b.Rex*(lep.theta_051 + lep.theta_151)*exc.H5_000
		+ (lep.theta_052 + lep.theta_152)*exc.H5_000
		+ b.Rex*(lep.theta_053 + lep.theta_153)*exc.H5_000);
}
Vec3 xs::exc_base_lp(EXC const& b, LepRadBaseLU const& lep_lu, LepRadBaseLP const& lep_lp, EXCLT const& exclt, EXCLL const& excll){
	return b.coeff*(
		1./b.Rex*(
			(lep_lu.theta_051 + lep_lu.theta_151)*exclt.H5_010
			+ (lep_lp.theta_071 + lep_lp.theta_171)*exclt.H7_100
			+ (lep_lp.theta_071 + lep_lp.theta_171)*excll.H7_001
			+ (lep_lp.theta_091 + lep_lp.theta_191)*exclt.H9_100
			+ (lep_lp.theta_091 + lep_lp.theta_191)*excll.H9_001
			)
		+ (
		 	(lep_lu.theta_052 + lep_lu.theta_152)*exclt.H5_010
			+ (lep_lp.theta_072 + lep_lp.theta_172)*exclt.H7_100
			+ (lep_lp.theta_072 + lep_lp.theta_172)*excll.H7_001
			+ (lep_lp.theta_092 + lep_lp.theta_192)*exclt.H9_100
			+ (lep_lp.theta_092 + lep_lp.theta_192)*excll.H9_001
			)
		+ b.Rex*(
		 	(lep_lu.theta_053 + lep_lu.theta_153)*exclt.H5_010
			+ (lep_lp.theta_073 + lep_lp.theta_173)*exclt.H7_100
			+ (lep_lp.theta_073 + lep_lp.theta_173)*excll.H7_001
			+ (lep_lp.theta_093 + lep_lp.theta_193)*exclt.H9_100
			+ (lep_lp.theta_093 + lep_lp.theta_193)*excll.H9_001
			)
		+ b.Rex*b.Rex*(
			(lep_lp.theta_074 + lep_lp.theta_174)*exclt.H7_100
			+(lep_lp.theta_074 + lep_lp.theta_174)*excll.H7_001
			+ (lep_lp.theta_094 + lep_lp.theta_194)*exclt.H9_100
			+ (lep_lp.theta_094 + lep_lp.theta_194)*excll.H9_001
			)
		);
		

}


