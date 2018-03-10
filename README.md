# SBGAT
The Small Bodies Geophysical Analysis Tool implementation

Public repository storing the source code of the implementation of the Small Body Geophysical Analysis Tool (SBGAT). 

* SBGAT's classes inherit from VTK's filters so as to facilitate manipulation and visualization tasks. SBGAT comes into a backend SBGATCore and frontend SBGATGui libraries. 

* All of SBGAT's functionality will be eventually available from within and without the GUI.

* The latest stable release of SBGAT can be retrieved from the `master` branch of the present repository. 

* The `develop` branch features code is undergoing active development or debugging.

## Wiki
The SBGAT User's Wiki can be found [here](https://github.com/bbercovici/SBGAT/wiki)

## Documentation
The SBGAT code documentation can be found [here](https://bbercovici.github.io/SBGAT-doc/index.html) 

It was generated with Doxygen and hosted on GitHub using the method described [here](https://visualstudiomagazine.com/articles/2015/03/01/github-pages.aspx) 

## Changelog


### SBGAT 1.03.0
SBGAT 1.03.0 sees the introduction of YORP coefficients computation as a new capability

* Added *SBGATSrpYorp*, a SBGAT filter enabling the computation of the force and torque YORP Fourier coefficients from a VTK Polydata. This class effectively provides a wrapper around *YORPLib*, a library developed by Jay W. McMahon at the University of Colorado Boulder that implements the analytical results derived in *The dynamical evolution of uniformly rotating asteroids subject to YORP* by Scheeres, D. J. (2007), Icarus, 188(2), 430–450. https://doi.org/10.1016/j.icarus.2006.12.015 
* YORP coefficients computation can be performed from within SBGATGui through the *Analyses* drop-down menu.  

 **Note that *YORPLib* is now a dependency of SBGAT and should be installed prior to compiling the newest SBGATCore and SBGATGui  [Instructions are provided on the corresponding wiki page](https://github.com/bbercovici/SBGAT/wiki/2)-Compile-and-install-SBGAT-dependencies#yorplib)** 

### SBGAT 1.02.1
SBGAT 1.02.1 marks the transition to VTK as SBGAT's backbone. 

* Added *SBGATMassProperties*, a SBGAT filter computing the surface area, volume, inertia and center of mass of a constant density polyhedron (see [Dobrovolskis, A. R. (1996). Inertia of Any Polyhedron. Icarus, 124, 698–704. ](https://doi.org/10.1006/icar.1996.0243]))
* Added *SBGATPolyhedronGravityModel*, a SBGAT filter computing the acceleration and potential of a constant density, topologically-closed polyhedron polyhedron (see [Werner, R.A. & Scheeres, D.J. Celestial Mech Dyn Astr (1996) 65: 313.](https://doi.org/10.1007/BF00053511]))
* Added validation tests 

### SBGAT 1.02.0
SBGAT 1.02 is a first take at fully leveraging VTK data structures for visual props representation and operation. 

Current features of SBGATGui include: 
* Small body shape model import from `.obj` files
* Trajectory loaded from time-XYZ ascii files. This capability will eventually replaces by SPICE kernels
* Spacecraft shape model import from `.obj` files
* Spacecraft displacement along previously loaded trajectory
* Addition/removal of light sources at arbitrary positions
* Computation of geometric measures such as surface area, volume, bounding boxes, center-of-mass and inertia tensor, the last two assuming a constant density distribution


![Visualization of gravity slopes on KW4 Alpha](http://i.imgur.com/fEvACWu.png)
![SBGATGui example](https://i.imgur.com/x0tb7hL.jpg)
![Visualization of a trajectory in Itokawa's body-fixed frame](https://i.imgur.com/xXRy1DY.png)



Created by Benjamin Bercovici
