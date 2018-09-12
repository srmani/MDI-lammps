module load intel/2017u4
module load intel-mpi/2017u4

cd ../src
make no-user-qmmm
make no-user-mdi
make yes-user-mdi

#build mdi
cd ../lib/mdi
rm -r build
mkdir build
cd build
cmake ..
make
cd ../
cp build/molssi_driver_interface/libmdi.a .

cd ../../src
make -j 32 cori2
