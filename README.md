# Small Bodies Geophysical Analysis Tool

The Small Bodies Geophysical Analysis Tool (SBGAT) implementation.

* SBGAT's classes inherit from VTK's filters so as to facilitate manipulation and visualization tasks. SBGAT comes in the form of a backend library *SbgatCore* and a frontend component *SbgatGui* exposing *SbgatCore*'s classes through a GUI.

* The latest stable release of SBGAT can be retrieved from the `master` branch of the present repository, or by using Homebrew as described below.

* The `develop` branch features code that is undergoing active development or debugging.

## Wiki

The SBGAT User's Wiki can be found [here](https://github.com/bbercovici/SBGAT/wiki), and the technical documentation [here](https://bbercovici.github.io/sbgat-doc/index.html)

## Installation: 

### Mac users

[Homebrew](https://brew.sh/) can be used to install SBGAT's components and dependencies. 

    brew tap bbercovici/self
    brew update
    brew install sbgat-core
    brew install sbgat-gui

The *SbgatGui* executable will be simlinked to `/usr/local/bin` .

### Linux & Mac users

[Refer to the detailed installation instructions](https://github.com/bbercovici/SBGAT/wiki/2:-Compiling-and-installing-SBGAT-dependencies).

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

### [SBGAT 2.05.1](https://github.com/bbercovici/SBGAT/releases/tag/2.05.1)

#### New
- Uncertainty quantification in the surface polyhedron gravity model can now be computed from within `SBGATGui`

#### Improvements
- Users now must load one shape model in `SBGATGui` before being able to open up any of the analyses/observations windows.

#### Bug fixes:
- `SBGATMassPropertiesUQ`'s master script has been corrected with the proper SIM_PREFIX


### [SBGAT 2.04.1](https://github.com/bbercovici/SBGAT/releases/tag/2.04.1)

#### New
- Introduced two new base classes, `SBGATFilter` and `SBGATFilterUQ`, the latter being tailored for uncertainty quantification in the considered geophysical properties
- A Python script producing shape "slices" has been added

<<<<<<< HEAD
### Improvement
=======
### Improvements
>>>>>>> develop
- `SBGATMassProperties` now inherits from `SBGATFilter`
- `SBGATPolyhedronGravityModel` now inherits from `SBGATMassProperties`
- `SBGATMassPropertiesUQ` now inherits from `SBGATFilterUQ`
- `SBGATPolyhedronGravityModelUQ` now inherits from `SBGATMassPropertiesUQ`
- Clicking on a facet in `SbgatGui` now results in displaying the facet center coordinates.
- Populated `Examples` folder with illustrative snippets for  `SBGATPolyhedronGravityModelUQ` and `SBGATMassPropertiesUQ`. These tests should be run from a Python 3 process (`run ../master_script.py`) called from their `build` directory.

## Regression
- As of April 7th 2019, OpenMP is no longer found by CMake on Mac.



### [SBGAT 2.02.4](https://github.com/bbercovici/SBGAT/releases/tag/2.02.4)

#### New
- SBGAT 2.02.4 introduces `SBGATPolyhedronGravityModelUQ`, a class dedicated to uncertainty quantification in Polyhedron Gravity potentials and accelerations arising from a stochastic shape. This class enables the evaluation of the variance in the potential and covariance in the acceleration at any point in space (excluding shape edges)
- Expanded and consolidated the Tests suite.

#### Improvements
- **Length unit consistency has been completely overhauled**. 
  * All of SBGAT's filters (`SBGATMassproperties`, `SBGATSphericalHarmo`, `SBGATPolyhedronGravityModel`,...) will return results using meters as length unit. For instance, calling `GetAcceleration` from the `SBGATPolyhedronGravityModel` class will always return an acceleration in `m/s^2`. 
  * Similarly, any method from the aforementioned filters expects an input position to be expressed in meters.
  * Classes documentation has been updated to reflect this change
  * Overall consistency is enforced by manually specifying the unit in which a given shape model is specified through the use of `SetScaleMeters` or `SetScaleKiloMeters`. This way, a shape whose coordinates are expressed in kilometers can be connected to an instance of `SBGATMassproperties` or any other filter and used in a completely transparent manner as long as `SetScaleKiloMeters` is called on the filter before `Update()`
  * `SBGATMassproperties`, `SBGATSphericalHarmo` and `SBGATPolyhedronGravityModel` are initialized by default for a shape whose length units are in meters (i.e the `scaleFactor` member is set to `1` by default, and set to `1000` if `SetScaleKiloMeters()` is called)
- NO MORE "UPDATE()" FROM GETTERS
- Work is in progress to revamp the main page of the doxygen documentation

### [SBGAT 2.02.3](https://github.com/bbercovici/SBGAT/releases/tag/2.02.3)

#### Bug fixes
- Fixed bug in `SBGATSphericalHarmo` that could have caused the evaluation of the spherical harmonics over a non-barycentered shape to be incorrect. 
- Pushed fix to latest version of `SHARMLib` dependency to address the same issue
- Modified `CMakeLists.txt` in Tests to fix issue caused by a conflicting header being sometimes included by one of VTK's dependencies


### [SBGAT 2.02.2](https://github.com/bbercovici/SBGAT/releases/tag/2.02.2)

#### New
- Shape models can now be modified from within SBGATGui, by selecting a vertex and applying a Gaussian interpolation of the vertex displacement to a k-neighbor neighborhood.

#### Improvements
- Camera is now positioned at the correct distance from the targeted shape body upon loading
- Improved visual aspect of selected facet
- Improved visual aspect of loaded shapes

#### Bug fixes:
- Fixed bug in SBGATGui that prevented proper alignment of the shape model with its principal axes

### [SBGAT 2.02.1](https://github.com/bbercovici/SBGAT/releases/tag/2.02.1)

#### New
- Second-moments about the mean in the volume, center-of-mass, unit-densiy inertia tensor, principal moments, principal dimensions and orientation of the principal axes can now be evaluated through the static methods of `SBGATShapeUncertainty`
- `SBGATShapeUncertainty` provides two methods to evaluate the statistical moments: a monte-carlo method `ComputeInertiaStatisticsMC`, and another method `ComputeInertiaStatistics` leveraging a linearized uncertainty model as proposed by Bercovici and McMahon (`Inertia statistics of an uncertain small body shape, ICARUS 2019 (In Review`). Both methods take the same argument (correleation length `l` and the standard deviation governing the error in the control point coordinates, directed along the average normal at these points). Both these methods assume that the shape error can be described as normally distributed and zero-mean.
- Added a new function in `Tests` to illustrate the convergence of the Monte-Carlo shape uncertainty sampling to the analytical prediction

#### Improvements
- `SBGATMassProperties::SaveMassProperties` is now saving the average radius (`r_avg = cbrt(3/4 *Volume/pi)`) to the JSON file
- `SBGATMassproperties::GetPrincipalAxes` could return 4 dcms, all representative of the same inertia ellipsoid. To enforce stability in the principal axes extractions, `SBGATMassproperties::GetPrincipalAxes` now returns the dcm that has the smallest-norm corresponding MRP.  
- `SBGAT` will try to link against `OpenMP` by default. This behaviour can now be disabled by passing the `--without-libomp` flag to Homebrew or by passing the `-DNO_OMP:BOOL=TRUE` flag to CMake

### Bug fixes:
- Fixed potential bug in `SBGATPolyhedronGravityModel` involving a parallel computing block where a variable with no viable reduction clause was being operated on
- `SBGATMassproperties::GetPrincipalAxes` is now returning the DCM `[PB]` converting from the body coordinate frame to the body principal frame. It was previously returning `[BP]`

**Note that *ShapeUQLib* is now a dependency of SBGAT and should be installed prior to compiling the newest *SBGATCore* and *SBGATGui*. [Instructions are provided on the corresponding wiki page.](https://github.com/bbercovici/SBGAT/wiki/2:-Compiling-and-installing-SBGAT-dependencies#shapeuqlib)** 

### [SBGAT 2.01.2](https://github.com/bbercovici/SBGAT/releases/tag/2.01.2)

#### New

- `SBGATMassProperties` now offers a method to save the computed mass properties to a JSON file (`SBGATMassProperties::SaveMassProperties`)
- A static method evaluating and saving the mass properties of the provided shape is now provided (`SBGATMassProperties::ComputeAndSaveMassProperties`)
- The `Measures` menu in `SbgatGUI` has been augmented with a `Save geometric measures` action

#### Improvements
- **The inertia tensor normalization has changed.** When computing the mass properties of a given small body, the following normalization is now applied to the inertia tensor: `I_norm = I / (mass * r_avg ^ 2)`  where `r_avg = cbrt(3/4 *Volume/pi)`. `r_avg` is now computed along with the other properties within `SBGATMassProperties`.
- The parallel axis theorem is no-longer applied to the small body. That is, the inertia will always be expressed about (0,0,0).
- Several GUI minors bug fixes

### [SBGAT 2.01.1](https://github.com/bbercovici/SBGAT/releases/tag/2.01.1)

#### New:
- Facets can now be individually selected in `SbgatGui` by clicking on them. If the surface PGM of the selected shape is available, the results for the selected facet will be shown in the console

#### Improvements
- The `Set Shape Mapper` option was renamed to `Set Results Overlay`
- Several GUI minors bug fixes

### Bug fixes:
- Fixed bug in `SbgatCore` that was due to an unecessary rescaling of the computed potentials and accelerations within `SbgatCore::ComputeSurfacePGM`. This bug was manifesting itself when calling `SbgatCore::ComputeSurfacePGM` with a shape model whose coordinates were expressed in kilometers. This bug was not affecting `SbgatGui` since it automatically rescales input shapes to meters upon loading.

### [SBGAT 1.12.2](https://github.com/bbercovici/SBGAT/releases/tag/1.12.2)

#### New:
- A previously evaluated surface Polyhedron Gravity Model can now be loaded from a JSON file via the static method `SBGATPolyhedronGravityModel::LoadSurfacePGM`
- `SbgatGui` can now overlay previously computed surface PGM results over a corresponding shape model (aka featuring the same number of facets as the one used to generate the surface PGM)

#### Improvements
- The surface polyhedron gravity model now computes: 
    * gravitational slopes
    * inertial gravity potentials
    * body-fixed gravity potentials
    * inertial gravity acceleration magnitudes
    * body-fixed gravity acceleration magnitudes

### [SBGAT 1.12.1](https://github.com/bbercovici/SBGAT/releases/tag/1.12.1)

#### New:
- The gravity-gradient matrix deriving from the polyhedron gravity model can now be evaluated
- The evaluation of the surface polyhedron gravity model can now be saved to a JSON file through `SBGATPolyhedronGravityModel::SaveSurfacePGM`
- This static method is available in `SbgatGUI` as the evaluation of the PGM will now require the specification of an output file

#### Improvements
- Stronger typing of inputs in `SBGATPolyhedronGravityModel`
- `SbgatGui` will now ask users whethers a loaded shape should be barycentered/principal-axis aligned, and apply the corresponding transform if answered `yes`


### [SBGAT 1.11.1](https://github.com/bbercovici/SBGAT/releases/tag/1.11.1)

#### New:
- The gravity-gradient matrix deriving from the spherical harmonics gravity model can now be evaluated.
- The partial derivative of the spherical harmonics gravity model with respect to the gravity spherical harmonics coefficients can now be evaluated.
- A static method `SBGATPolyhedronGravityModel::ComputeSurfacePGM` has been added to `SbgatCore` to facilitate the evaluation of the polyhedron gravity model at the surface of a small body shape.

#### Improvements
- The CMake configuration of SBGAT will no longer failed if OpenMP cannot be found.
- A warning will now be issued in `SbgatGui` if the surface polyhedron gravity model is evaluated with a zero rotation period.
- `SBGATPolyhedronGravityModel` and `SBGATSphericalHarmo` are now returning potentials and accelerations evaluated with the same unit of length as the shape model they were associated with.

### Bug fixes:
- Fixed bug in `SbgatGui` that was due to a `vtkPolydataMapper` not being properly assigned to the correct `ModelDataWrapper` after aligning the shape

### [SBGAT 1.10.1](https://github.com/bbercovici/SBGAT/releases/tag/1.10.1)

#### New:
- The Polyhedron Gravity Model can now be evaluated at the surface of loaded shape models from `SbgatGui`s
- The surface PGM can then be overlaid in the form of surface slopes, gravitional potential, inertial and body-fixed accelerations

#### Improvements:

- The SBGAT documentation should be back online at [the following location](https://bbercovici.github.io/sbgat-doc/index.html)

#### Bug fixes:
- Fixed a major regression that was causing *SbgatGui* to crash when computing lightcurves


### [SBGAT 1.09.1](https://github.com/bbercovici/SBGAT/releases/tag/1.09.1)

#### New:

- Built-from source GCC no longer required for OpenMP support on Mac since CLANG 10 now provides the required definitions
- Added a method simultaneously computing gravity potential and acceleration in `SBGATPolyhedronGravityModel` and updated corresponding tests 

#### Improvements
-  Added constant qualifiers to `SBGATPolyhedronGravityModel` inputs

#### Bug fixes: 
- Captions should now properly show in *SbgatGui*

### [SBGAT 1.08.1](https://github.com/bbercovici/SBGAT/releases/tag/1.08.1)

#### New:

- Created `SBGATObs` base-class from which `SBGATObsLightcurve` and `SBGATObsRadar` derive
- Added new functionalities to *SbgatGui*
- Created the dependency `OrbitConversions` to generate Keplerian trajectories or convert 6 dof states between Cartesians and Keplerian representations.

#### Improvements
-  `SBGATObsLightcurve` and `SBGATObsRadar` now handle primary/secondary asteroid systems, allowing generation of lightcurves and radar observations of binary asteroid systems. Relative trajectories can be computed under the assumption that the secondary is undergoing a keplerian orbit about the primary, or be loaded from an external file.
- Removed dependency of `SBGATSphericalHarmo` to density and mass of the considered shape model

#### Bug fixes: 
- Fixed a bug in `SBGATPolyhedronGravityModel.cpp` where the edge extraction would sometimes fail. Fix consisted in filtering the input through a `vtkPolyDataCleaner` before handing it to `vtkExtractEdges`.


### [SBGAT 1.06.1](https://github.com/bbercovici/SBGAT/releases/tag/1.06.1)

#### New:
- Added `SBGATObsLightcurve` to *SbgatCore* , a module enabling the generation of instantaneous-exposure lightcurves in a fixed-spin scenario. This module assumes constant small-body spin and phase angle between the sun, the small body and the observer.
- `SBGATObsRadar` now throws an instance of `std::runtime_error` if the specified bin sizes are incompatible with the collected data that may yield an empty histogram dimension
- Observations from `SBGATObsRadar` and `SBGATObsLightcurve` can be penalized by incidence so as to diminish the weight of a given measurement. `SBGATObsRadar` weighs by the `cos` of the angle between the observer and the surface normal, while `SBGATObsLightcurve` weighs by the product of the `cos` of the angle between the observer and the surface normal and the `cos` of the angle between the sun and the surface normal

#### Improvements
- Simulated Range/Range-rate images and lightcurves rely on area-weighted surface sampling : `N * surface_area/max_surface_area` points are sampled for each facet, where `max_surface_area` is the surface area of the largest facet in the shape and `surface_area` that of the considered facet
- Removed more deprecated functionalities

#### Bug fixes: 
- Fixed bug in *SbgatGui* that was allowing users to bin radar observations before effectively collecting them.
- Saved radar images now have correct color levels


### [SBGAT 1.05.2](https://github.com/bbercovici/SBGAT/releases/tag/1.05.2)

- Adds `SBGATObsRadar` to *SbgatCore*, a class emulating range/range-rate radar measurements. The corresponding menu and action are also available in *SbgatGui*
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
