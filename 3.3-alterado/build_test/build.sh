#!/bin/sh

export OMPI_MPIFC=ifort
export OMPI_MPIC=icc
export OMPI_CC=icc
export OMPI_FC=ifort
export CC=icc
exprot FC=ifort


make

