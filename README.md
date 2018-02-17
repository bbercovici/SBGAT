# SBGAT
The Small Bodies Geophysical Analysis Tool implementation

Public repository storing the source code of the implementation of the Small Body Geophysical Analysis Tool (SBGAT). 

The latest stable release of SBGAT can be retrieved from the `master` branch of the present repository. 

The `develop` branch features code is undergoing active development or debugging.

## Wiki
The SBGAT User's Wiki can be found [here](https://github.com/bbercovici/SBGAT/wiki)

## Documentation
The SBGAT code documentation can be found [here](https://bbercovici.github.io/sbgat-doc/index.html) 

It was generated with Doxygen and hosted on GitHub using the method described [here](https://visualstudiomagazine.com/articles/2015/03/01/github-pages.aspx) 

## Changelog

### Sbgat 1.02.1
Sbgat 1.02.1 marks the transition to VTK as SBGAT's backbone. 

* Polyhedron Gravity Model potential and accelerations can be generated from an arbitrary closed vtkPolydata representing a polyhedron
* SBGATMassProperties, a SBGAT filter computing the surface area, volume, inertia and center of mass of a constant density polyhedron
* SBGATPolyhedronGravityModel, a SBGAT filter computing the acceleration and potential of a constant density polyhedron. 
* More comprehensive validation tests 

### Sbgat 1.02.0
Sbgat 1.02 is a first take at fully leveraging VTK data structures for visual props representation and operation. 

Current features of SbgatGUI include: 
* Small body shape model import from `.obj` files
* Trajectory loaded from time-XYZ ascii files. This capability will eventually replaces by SPICE kernels
* Spacecraft shape model import from `.obj` files
* Spacecraft displacement along previously loaded trajectory
* Addition/removal of light sources at arbitrary positions
* Computation of geometric measures such as surface area, volume, bounding boxes, center-of-mass and inertia tensor, the last two assuming a constant density distribution


![Visualization of gravity slopes on KW4 Alpha](http://i.imgur.com/fEvACWu.png)
![SbgatGui example](https://i.imgur.com/x0tb7hL.jpg)
![Visualization of a trajectory in Itokawa's body-fixed frame](https://i.imgur.com/xXRy1DY.png)



Created by Benjamin Bercovici
