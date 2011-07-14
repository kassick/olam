#!/bin/bash

function saverastro {
	mkdir -p rastros/$1
	mv *.rst rastros/$1
}

np=${NP:-8}
echo $np
rm -rf /tmp/sfcfiles 
mkdir /tmp/sfcfiles 
rm -rf /tmp/hist 
mkdir /tmp/hist




echo SFC
mpirun -np $np ./olam-3.3-alterado-mpi -f OLAMIN.makesfc.local

saverastro makesfc

echo GRID

mpirun -np $np ./olam-3.3-alterado-mpi -f OLAMIN.makegrid.local

saverastro makegrid


echo RUN
mpirun -np $np ./olam-3.3-alterado-mpi -f OLAMIN.86400.1800.local

saverastro run.86400.1800

