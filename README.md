# SIDIS-RC EvGen

This is an event generator for polarized semi-inclusive deep inelastic
scattering (SIDIS) with QED radiative corrections, following closely work done
in [1]. The project is divided into a library component `sidis` for calculating
SIDIS cross-sections, and a binary component `sidisgen` for generating events.

A detailed description of the generator can be found on
[the arxiv](https://arxiv.org/abs/2210.03785).

## Quick start

### Generator

The `sidisgen` generator uses a parameter file. The allowed options are listed
by `sidisgen --help`. As an example:

```csv
mc.num_events        10000
file.gen             gen.root
file.event           event.root
phys.sf_set          prokudin
phys.rc_method       approx
phys.soft_threshold  0.01
phys.mass_threshold  1.073249081
setup.beam           e
setup.target         p
setup.hadron         pi+
setup.beam_energy    10.6
setup.beam_pol       0
setup.target_pol     0 0 0
```

Then, initialize a FOAM to approximate the differential cross-section in the
specified kinematic region.

```bash
sidisgen --initialize <param-file>
```

Finally, use the FOAM for Monte-Carlo event generation.

```bash
sidisgen --generate <param-file>
```

To process the resulting ROOT file, see `examples/process_events.cpp`.

### Library

Many demonstrations of the different features of the `sidis` library can be
found in the `examples` folder. To get started quickly:

```cpp
#include <iostream>
#include <sidis/sidis.hpp>
#include <sidis/sf_set/prokudin.hpp>

sidis::Real const PI = sidis::PI;
sidis::Real const M_TH = sidis::MASS_P + sidis::MASS_PI_0;

int main() {
	sidis::part::Particles particles(
		sidis::part::Nucleus::P,   // Target nucleus.
		sidis::part::Lepton::E,    // Beam lepton.
		sidis::part::Hadron::PI_P, // Leading hadron.
		M_TH                       // Threshold mass of undetected part.
	);
	sidis::Real S = 2. * 10.6 * particles.M; // Kinematic variable `S = 2 p k1`.
	sidis::kin::PhaseSpace phase_space {
		0.2,      // Bjorken x.
		0.9,      // Bjorken y.
		0.3,      // Bjorken z.
		2.,       // Transverse momentum of hadron, squared.
		0.5 * PI, // Azimuthal angle of hadron.
		0.,       // Azimuthal angle of transverse target polarization.
	};
	sidis::kin::Kinematics kin(particles, S, phase_space);
	sidis::Real beam_pol = 0.;
	sidis::math::Vec3 target_pol(0., 0., 0.);
	// Compute structure functions with WW-type approximation.
	sidis::sf::set::ProkudinSfSet sf;
	sidis::Real born_xs = sidis::xs::born(kin, sf, beam_pol, target_pol);
	std::cout << "Born unpolarized cross-section is " << born_xs << std::endl;
	return 0;
}
```

## Build

For building the generator, ROOT must be installed on your system. For building
the documentation, Doxygen must be installed. For building the tests, Catch2
must be installed.

```bash
# Clone the project including submodules.
git clone https://github.com/duanebyer/sidis.git
cd sidis
git submodule init
git submodule update
# Make a build directory.
mkdir build
cd build
# Configure the build.
cmake -DCMAKE_BUILD_TYPE=Release ..
# Build.
make
# Optionally install the library and binary files to your system. Use
# `CMAKE_INSTALL_PREFIX` during the configure to choose the install location.
make install
```

The following CMake configuration options may be of use:
* `Sidis_REAL_TYPE`: The floating point type to use for all cross-section
  calculations. Only the standard C++ floating point types `float`, `double`,
  and `long double` are supported.
* `Sidis_BUILD_TESTS`: Whether to build the tests.
* `Sidis_BUILD_EXAMPLES`: Whether to build the examples.
* `Sidis_BUILD_DOXYGEN`: Whether to build the Doxygen documentation.
* `Sidis_BUILD_APPS`: Whether to build the `sidisgen` binary.
* `Sidis_ENABLE_IPO`: Whether to build with interprocedural optimization.

## Acknowledgements

This software makes use of the following libraries:
* [WW-SIDIS](https://github.com/prokudin/WW-SIDIS): Structure functions for
  proton using WW-type approximation [2].
* [mstwpdf](https://mstwpdf.hepforge.org/): Parton distribution functions for
  proton [3].
* [ROOT](https://root.cern/): Plotting and data analysis.
* [FOAM](http://jadach.web.cern.ch/jadach/Foam/Index.html): Monte-Carlo event
  generation using spatial partitioning.
* [cubature](https://github.com/stevengj/cubature): Multi-dimensional numerical
  integration.
* [GSL](https://www.gnu.org/software/gsl/): Multi-dimensional Monte-Carlo
  integration.

## References

[1] I. Akushevich and A. Ilyichev, 2019. Lowest order QED radiative effects in
    polarized SIDIS. Phys. Rev. D 100(3), 033005.

[2] S. Bastami, H. Avakian, A. V. Efremov, A. Kotzinian, B. U. Musch, B.
    Parsamyan, A. Prokudin, M. Schlegel, P. Schweitzer. Semi-inclusive deep-
	inelastic scattering in Wandzura-Wilczek-type approximation. J. High Energy
	Phys. 1906(2019), 007.

[3] A. D. Martin, W. J. Stirling, R. S. Thorne, G. Watt. Parton distributions
    for the LHC. Eur. Phys. J. C. 63(2009), 189-285.

