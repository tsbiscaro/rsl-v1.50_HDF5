#!/bin/bash
rm -f ltmain.sh
rm -f aclocal.m4
rm -f Makefile
rm -f configure
libtoolize
aclocal
autoconf
automake
