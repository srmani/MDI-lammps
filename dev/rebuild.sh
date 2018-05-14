module load intel/2017u4
module load intel-mpi/2017u4

#replace the QM/MM makefile
cp Makefile.qmmm ../lib/qmmm/Makefile.ifort

cd ../src
make no-user-qmmm
make no-user-misc
make yes-user-qmmm
make yes-user-misc
make -j 32 cori2
