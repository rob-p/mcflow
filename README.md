About this project:
========================================

This project provides a C++ implementation of the methods described in the
Siggraph 1999 paper by Desbrun et al. entitled "Implicit fairing of irregular
meshes using diffusion and curvature flow."[1] It aims to provide a readable
implementation, but also focuses on performance.  The algorithm hinges on the
solution of a linear system of equations, and this project provides two
different implementations to compute this solution.  One implementation
(apply_mc_flow) solves the system iteratively using a PBCG solver, while the
other (apply_mc_flow_cholesky) performs a Cholesky factorization and then solves
the system directly.  In general, the direct method is much faster.  In the
future, I may add an implementation which solves the system using PBCG on the
GPU; it will be interesting to see how this approach performs.

Dependencies:
========================================

The code depends on a number of external libraries.  The dependencies are as
follows:

OpenMesh (>= v2.0-RC4): http://www.openmesh.org/index.php
Boost (>= v1.42): http://www.boost.org
GMM (>= v4.1): http://download.gna.org/getfem/html/homepage/index.html
TAUCS (>= v2.2): http://www.tau.ac.il/~stoledo/taucs/

The last dependency, TAUCS, is only required if you wish to compile and use the
Cholesky factorization based solver.  Additionally, this project (and most of my
others) use the CMake build system (http://www.cmake.org/).

Building:
========================================

If you have the dependencies, building the software should be fairly
straightforward.  Assuming that '>' represents the command prompt below, and
you're currently in the top-level directory, the project can be made as follows:

> mkdir build
> cd build
> cmake ..
> make

This will produce two executables; apply_mc_flow and apply_mc_flow_cholesky.  As
mentioned above, the former solves the necessary linear system using a PBCG
method while the second performs a Cholesky factorization followed by a direct
solve.  If you don't have TAUCS installed (it can be a pain to build), then you
can install just the PBCG based variant of the software by replacing

> make

with

> make apply_mc_flow

Conversely, if you just want to make the Cholesky based version, you can just
"make apply_mc_flow_cholesky".

Usage:
========================================

Using the program is fairly straightforward as well.  However for a good
description of the different parameters, I refer you to the paper by Desbrun et
al.[1].  The program takes as input a mesh (in any of the formats that can be
read by OpenMesh), and outputs a processed (smoothed) version of that mesh to
the desired file.  It additionally takes a parameter denoting the amount of
smoothing to be done each iteration, as well as the number of smoothing
iterations to perform.  Finally, you can choose between applying mean curvature
flow or Gaussian curvature flow.


References:
========================================

[1] @inproceedings{Desbrun:2003:Implicit,
     author = {Desbrun, Mathieu and Meyer, Mark and Schr\"{o}der, Peter
              and Barr, Alan H.},
     title = {Implicit fairing of irregular meshes using diffusion and
             curvature flow},
     booktitle = {SIGGRAPH '99: Proceedings of the 26th annual conference
                 on Computer graphics and interactive techniques},
     year = {1999},
     isbn = {0-201-48560-5},
     pages = {317--324},
     doi = {http://doi.acm.org/10.1145/311535.311576},
     publisher = {ACM Press/Addison-Wesley Publishing Co.},
     address = {New York, NY, USA},
    }
