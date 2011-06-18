#!/bin/bash

function saverastro {
	mkdir -p rastros/$1
	mv *.rst $1
}

nodefile=$1
np=`cat $nodefile|wc -l`
echo $nodefile $np
for i in `cat $nodefile`; do
	ssh $i "rm -rf /tmp/sfcfiles ;  mkdir /tmp/sfcfiles ; rm -rf /tmp/hist ; mkdir /tmp/hist"
done


echo GRID

mpirun -np $np -machinefile $nodefile ./olam-3.3-alterado-mpi -f OLAMIN.makegrid.local

saverastro makegrid

#copy gridfile
for i in `cat $nodefile`; do
	if [  "$i" == "`hostname`" ]; then
		echo no
		continue
	fi
	scp /tmp/sfcfiles/* $i:/tmp/sfcfiles;
done


echo SFC
mpirun -np $np -machinefile $nodefile ./olam-3.3-alterado-mpi -f OLAMIN.makesfc.local

saverastro makesfc


echo RUN
mpirun -np $np -machinefile $nodefile ./olam-3.3-alterado-mpi -f OLAMIN.86400.1800.local

saverastro run.86400.1800

