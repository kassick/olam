#!/bin/sh
# File: "/home/kassick/Work/olam/trace2paje/generate.sh"
# Created: "Ter, 31 Mai 2011 12:26:51 -0300 (kassick)"
# Updated: "Seg, 17 Out 2011 19:58:04 -0200 (kassick)"
# $Id$
# Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 

# Path to compiled librastro
export RASTRO_PREFIX=$HOME/Work/librastro
export LANG=C
rm -rf CMakeCache.txt
rm -rf CMakeFiles
cmake CMakeLists.txt

touch CMakeCache.txt
