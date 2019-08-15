cd ../src
make no-user-qmmm
make no-user-mdi
make yes-user-mdi

#build mdi
cd ../lib/mdi
python Install.py -m gcc

cd ../../src
make -j 32 mpi
cp lmp_mpi lmp_cori2
