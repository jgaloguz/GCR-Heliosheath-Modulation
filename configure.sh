#!/bin/bash

autoreconf
automake --add-missing
./configure CXXFLAGS="-Ofast" --with-mpi=$1 --with-execution=PARALLEL --with-trajectory=PARKER --with-time_flow=BACKWARD --with-rkmethod=0 --with-server=SELF