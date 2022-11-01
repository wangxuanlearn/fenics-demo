#!/bin/bash
# Execute this file to recompile locally
c++ -Wall -shared -fPIC -std=c++11 -O3 -fno-math-errno -fno-trapping-math -ffinite-math-only -I/usr/lib/petscdir/petsc3.12/x86_64-linux-gnu-real/include -I/usr/lib/slepcdir/slepc3.12/x86_64-linux-gnu-real/include -I/usr/lib/x86_64-linux-gnu/openmpi/include -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi -I/usr/include/hdf5/openmpi -I/usr/include/eigen3 -I/home/wangxuan/.cache/dijitso/include dolfin_subdomain_b52c1a2e7ea15353a421996e34af254a.cpp -L/usr/lib/x86_64-linux-gnu/openmpi/lib -L/usr/lib/petscdir/petsc3.12/x86_64-linux-gnu-real/lib -L/usr/lib/slepcdir/slepc3.12/x86_64-linux-gnu-real/lib -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi -L/home/wangxuan/.cache/dijitso/lib -Wl,-rpath,/home/wangxuan/.cache/dijitso/lib -lmpi -lmpi_cxx -lpetsc_real -lslepc_real -lm -ldl -lz -lsz -lhdf5 -lboost_timer -ldolfin -olibdijitso-dolfin_subdomain_b52c1a2e7ea15353a421996e34af254a.so