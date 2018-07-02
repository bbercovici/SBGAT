# Small Bodies Geophysical Analysis Tool

The Small Bodies Geophysical Analysis Tool (SBGAT) implementation.

* SBGAT's classes inherit from VTK's filters so as to facilitate manipulation and visualization tasks. SBGAT comes in the form of a backend library *SbgatCore* and a frontend component *SbgatGui* exposing *SbgatCore*'s classes through a GUI.

* The latest stable release of SBGAT can be retrieved from the `master` branch of the present repository, or by using Homebrew as described below.

* The `develop` branch features code that is undergoing active development or debugging.

## Wiki

The SBGAT User's Wiki can be found [here](https://github.com/bbercovici/SBGAT/wiki).

## Installation: 

### Mac users

[Homebrew](https://brew.sh/) can be used to install SBGAT's components and dependencies. Homebrew install options differ whether you want to use *SbgatGui* or just the core classes of *SbgatCore*. This is justified by *SbgatCore*'s dependency on *VTK*, as *VTK* needs to know whether it must be linked against *Qt* at build time.

#### With SbgatGui

    brew tap bbercovici/self
    brew update
    brew install vtk --with-qt
    brew install sbgat-core
    brew install sbgat-gui

The *SbgatGui* executable will be simlinked to `/usr/local/bin` .

#### Without SbgatGui

    brew tap bbercovici/self
    brew update
    brew install vtk
    brew install sbgat-core

If you change your mind and decide you want the *SbgatGui* component, run:

    brew uninstall sbgat-core
    brew uninstall vtk
    
before reinstalling.

### Linux & Mac users

[Refer to the detailed installation instructions](https://github.com/bbercovici/SBGAT/wiki/2:-Compile-and-install-SBGAT-dependencies).

## Getting updates

### Mac users

Assuming that *SbgatCore* was installed with Homebrew:

    brew update
    brew upgrade sbgat-core

If installed, after updating Homebrew, *SbgatGui* can be also upgraded:

    brew upgrade sbgat-gui


### Linux & Mac users

Check each of SBGAT's dependencies repository and SBGAT's repository itself for updates. Assuming that the current directory is the original Git local repository for the component you wish to update, run :

    git pull
    cd build
    cmake ..
    make
    make install

to apply the update (if any).

## Changelog

### [SBGAT 1.07.1](https://github.com/bbercovici/SBGAT/releases/tag/1.07.1)

#### New:
- Created module `SBGATKeplerianTraj`, a wrapper around the `OrbitConversions` library. This means that SBGAT can generate Keplerian trajectories

#### Improvements
-  `SBGATObsLightcurve` now handles multi-body configurations, allowing generation of lightcurves of binary asteroid systems under the assumption that the secondary is undergoing a keplerian orbit about the primary 

#### Bug fixes: 


### [SBGAT 1.06.1](https://github.com/bbercovici/SBGAT/releases/tag/1.06.1)

#### New:
- Added `SBGATObsLightcurve` to `SbgatCore` , a module enabling the generation of instantaneous-exposure lightcurves in a fixed-spin scenario. This module assumes constant small-body spin and phase angle between the sun, the small body and the observer.
- `SBGATObsRadar` now throws an instance of `std::runtime_error` if the specified bin sizes are incompatible with the collected data that may yield an empty histogram dimension
- Observations from `SBGATObsRadar` and `SBGATObsLightcurve` can be penalized by incidence so as to diminish the weight of a given measurement. `SBGATObsRadar` weighs by the `cos` of the angle between the observer and the surface normal, while `SBGATObsLightcurve` weighs by the product of the `cos` of the angle between the observer and the surface normal and the `cos` of the angle between the sun and the surface normal

#### Improvements
- Simulated Range/Range-rate images and lightcurves rely on area-weighted surface sampling : `N * surface_area/max_surface_area` points are sampled for each facet, where `max_surface_area` is the surface area of the largest facet in the shape and `surface_area` that of the considered facet
- Removed more deprecated functionalities

#### Bug fixes: 
- Fixed bug in `SbgatGui` that was allowing users to bin radar observations before effectively collecting them.
- Saved radar images now have correct color levels


### [SBGAT 1.05.2](https://github.com/bbercovici/SBGAT/releases/tag/1.05.2)

- Adds `SBGATObsRadar` to `SbgatCore`, a class emulating range/range-rate radar measurements. The corresponding menu and action are also available in `SbgatGui`
- If `gcc` exists in Homebrew's Cellar, SBGAT and its dependencies will be compiled using this OpenMP compliant compiler, giving better performance on multithreaded platforms. [This functionality had to be postponed due to Qt 5.10 incompability with recent gcc versions on MacOS](https://bugreports.qt.io/browse/QTBUG-66585). 


### [SBGAT 1.05.1](https://github.com/bbercovici/SBGAT/releases/tag/1.05.1)

This new release of SBGAT allows import/export of gravity spherical harmonics from/into SBGAT by means of Niels Lohmann's Modern C++ JSON library. This functionality is available from SbgatCore's classes and SbgatGui as well.

### SBGAT 1.04.3

* SBGAT 1.04.3 marks the shift to the Homebrew package manager as a way to greatly facilitate SBGAT's distribution and update. It is of course still possible to download each of SBGAT's dependencies separatly and manually build from source.

### SBGAT 1.04.2

* SBGAT 1.04.2 enables the computation of the spherical harmonics expansion directly from SbgatGUI

	* Added an action *Compute Gravity Spherical Harmonics* under *Analyses*
	* Added an action *Align Shape* under *Small Body*. This action aligns the barycenter of the selected shape with (0,0,0) and its principal axes with the rendering window axes. This is a prerequisite for meaningful YORP or gravity spherical harmonics computations.
	* Added an action *Save Shape Model* under *Small Body*. This action exports the selected shape model in its current state to an .obj file of choice.

 [**Users must update their versions of RigidBodyKinematics to reflect the latest changes**](https://github.com/bbercovici/RigidBodyKinematics) 

### SBGAT 1.04.1

* [SBGAT and its dependencies are now distributed under the MIT license](https://choosealicense.com/licenses/mit/). 
* No new functionalities besides updated license information. 

### SBGAT 1.04.0

* SBGAT 1.04.0 can now be used to compute the spherical harmonics expansion of the exterior gravity field about a constant-density polyhedron

	* Added *SBGATSphericalHarmo*, a SBGAT filter enabling the computation and evaluation of the spherical harmonics coefficients of the exterior gravity field caused by a constant density shape represented by a `vtkPolydata`. This class effectively provides a wrapper around *SHARMLib*, a library developed by Benjamin Bercovici from the original works of Yu Takahashi and Siamak Hesar at the University of Colorado Boulder. For more details, see [Spherical harmonic coefficients for the potential of a constant-density polyhedron](https://www.sciencedirect.com/science/article/pii/S0098300497001106) by Werner, R. a. (1997).
	* Added a test where the spherical harmonics expansion is computed and evaluated around KW4. The test succeeds if the acceleration error relative to the polyhedron gravity model is less that 0.0001 %

 **Note that *SHARMLib* is now a dependency of SBGAT and should be installed prior to compiling the newest *SBGATCore* and *SBGATGui*. [Instructions are provided on the corresponding wiki page.](https://github.com/bbercovici/SBGAT/wiki/2:-Compile-and-install-SBGAT-dependencies#sharmlib)** 


### SBGAT 1.03.0

* SBGAT 1.03.0 sees the introduction of YORP coefficients computation

	* Added *SBGATSrpYorp*, a SBGAT filter enabling the computation of the YORP force and torque Fourier coefficients from a VTK Polydata. This class effectively provides a wrapper around *YORPLib*, a library developed by Jay W. McMahon at the University of Colorado Boulder that implements the analytical results derived in [The dynamical evolution of uniformly rotating asteroids subject to YORP](https://doi.org/10.1016/j.icarus.2006.12.015) by Scheeres, D. J. (2007).
	* YORP coefficients computation can be performed from within *SBGATGui* through the *Analyses* drop-down menu.  

 **Note that *YORPLib* is now a dependency of SBGAT and should be installed prior to compiling the latest *SBGATCore* and *SBGATGui*. [Instructions are provided on the corresponding wiki page.](https://github.com/bbercovici/SBGAT/wiki/2:-Compile-and-install-SBGAT-dependencies#yorplib)** 

### SBGAT 1.02.1

* SBGAT 1.02.1 marks the transition to VTK as SBGAT's backbone. 

	* Added *SBGATMassProperties*, a SBGAT filter computing the surface area, volume, inertia and center of mass of a constant density polyhedron (see [Dobrovolskis, A. R. (1996). Inertia of Any Polyhedron. Icarus, 124, 698–704. ](https://doi.org/10.1006/icar.1996.0243]))
	* Added *SBGATPolyhedronGravityModel*, a SBGAT filter computing the acceleration and potential of a constant density, topologically-closed polyhedron polyhedron (see [Werner, R.A. & Scheeres, D.J. Celestial Mech Dyn Astr (1996) 65: 313.](https://doi.org/10.1007/BF00053511]))
	* Added validation tests 

### SBGAT 1.02.0

* SBGAT 1.02 is a first take at fully leveraging *VTK* data structures for visual props representation and operation. Current features of *SBGATGui* include: 
	* Small body shape model import from `.obj` files
	* Trajectory loaded from time-XYZ ascii files. This capability may eventually be replaced by SPICE kernels
	* Spacecraft shape model import from `.obj` files
	* Spacecraft displacement along previously loaded trajectory
	* Addition/removal of light sources at arbitrary positions
	* Computation of geometric measures such as surface area, volume, bounding boxes, center-of-mass and inertia tensor, the last two assuming a constant density distribution


![Visualization of gravity slopes on KW4 Alpha](http://i.imgur.com/fEvACWu.png)
![SBGATGui example](https://i.imgur.com/x0tb7hL.jpg)
![Visualization of a trajectory in Itokawa's body-fixed frame](https://i.imgur.com/xXRy1DY.png)


## Documentation

The SBGAT code documentation can be found [here](https://bbercovici.github.io/sbgat-doc/index.html). [It was generated with Doxygen and hosted on GitHub using the method described here](https://visualstudiomagazine.com/articles/2015/03/01/github-pages.aspx) 

## License

[This software is distributed under the MIT License](https://choosealicense.com/licenses/mit/)

Created by Benjamin Bercovici
