#!/bin/sh

make
make kmod

echo "Clean dst?"
read yn
if [ "$yn" = "yes" ]; then
  rm -rf dst
fi
#rm -rf dst
mkdir dst
sudo make DESTDIR=`pwd`/dst install
sudo make DESTDIR=`pwd`/dst kmod_install
cd dst
tar czf $1 .
cd ..


echo create bundle $1

