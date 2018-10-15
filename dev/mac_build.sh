cp Makefile.mac ../src/MAKE/MACHINES/Makefile.mac

cd ../src
make mpi-stubs
make no-user-mdi

#build mdi
cd ../lib/mdi
python Install.py -m gcc

cd ../../src
make -j 8 mac
