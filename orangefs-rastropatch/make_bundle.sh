#!/bin/sh

make
make kmod

rm -rf dst
mkdir dst
sudo make DESTDIR=`pwd`/dst install
sudo make DESTDIR=`pwd`/dst kmod_install
cd dst
tar czf $1 .
cd ..


echo create bundle $1

