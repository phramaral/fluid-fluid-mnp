# fluid-fluid-mnp
Files containing structure of systems with interaction of water, heptane and functionalized magnetic nanoparticle.

In each directory, we have xyz files of each molecule belonging to the systems. In addition to these, there is code written by us (Python) that creates uniformly distributed points on the surface of a sphere of radius R (fibonacci_sphere.py). These points are the anchor points on the NPM to which the functional group is linked. In addition to this code, we have a code (FORTRAN) that links the functional groups to the anchor points (np_graphitization.f90).
