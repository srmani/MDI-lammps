module load intel/2017u4
module load intel-mpi/2017u4

#replace the standard cori2 makefile with one for Edison
cp Makefile.cori2 ../src/MAKE/MACHINES

#also replace the default MPI makefile, since the QMMM package uses it
cp Makefile.cori2 ../src/MAKE/Makefile.mpi

#replace the QM/MM makefile
cp Makefile.qmmm ../lib/qmmm/Makefile.ifort

#add a makefile for espresso
cp Makefile.espresso ../lib/qmmm/Makefile.espresso

cd ../src
make yes-standard
make no-gpu
make no-kim
make no-kokkos
make no-latte
make no-meam
make no-mscg
make no-poems
make no-python
make no-voronoi
make yes-user-qmmm
make yes-user-misc

cd ../lib/qmmm
make -j 32 -f Makefile.ifort

cd ../reax
make -j 32 -f Makefile.ifort

cd ../../src
make -j 32 cori2

cd ../lib/qmmm
make -j 32 -f Makefile.espresso all
make -j 32 -f Makefile.ifort all
