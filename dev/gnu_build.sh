#replace the standard cori2 makefile with one for Edison
cp Makefile.cori2_gnu ../src/MAKE/MACHINES/Makefile.cori2

#also replace the default MPI makefile, since the QMMM package uses it
cp Makefile.cori2_gnu ../src/MAKE/Makefile.mpi

#replace the QM/MM makefile
cp Makefile.qmmm ../lib/qmmm/Makefile.gfortran

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
make no-user-qmmm
make yes-user-mdi

#build mdi
cd ../lib/mdi
python Install.py -m gcc

cd ../reax
#make -j 32 -f Makefile.gfortran
python Install.py -m gfortran

cd ../../src
make -j 32 mpi
cp lmp_mpi lmp_cori2

#cd ../lib/qmmm
#make -j 32 -f Makefile.espresso all
#make -j 32 -f Makefile.ifort all
