module load intel/2017u4
module load intel-mpi/2017u4

cd ../src
make no-user-qmmm
make no-user-mdi
make yes-user-mdi

#build mdi
cd ../lib/mdi
python Install.py -m icc

cd ../../src
make -j 32 cori2
