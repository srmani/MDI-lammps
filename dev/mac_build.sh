cp Makefile.mac ../src/MAKE/MACHINES/Makefile.mac

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
make no-mpiio
make no-reax
make yes-user-mdi

make mpi-stubs

#build mdi
cd ../lib/mdi
python Install.py -m gcc

#cd ../reax
#make -j 32 -f Makefile.gfortran

cd ../../src
make -j 8 mac
