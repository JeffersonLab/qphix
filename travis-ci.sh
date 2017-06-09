#!/bin/bash
# Copyright © 2017 Martin Ueding <dev@martin-ueding.de>

set -e
set -u
set -x

#sudo apt-get update
#sudo apt-get upgrade -y python3
sudo apt-get install -y python3-jinja2

python3 --version
python3 -c 'import sys; print("\n".join(sys.path))'
dpkg-query -L python3-jinja2

export PYTHONPATH=/usr/lib/python3/dist-packages

export SRCDIR=${PWD}
export SANDBOX_PATH=${PWD}/sandbox
echo Sandbox Path is ${SANDBOX_PATH}

# AVX build
export BUILDDIR=${SANDBOX_PATH}/build/avx
export INSTALLDIR=${SANDBOX_PATH}/install/avx
mkdir -p ${BUILDDIR}
mkdir -p ${INSTALLDIR}

pushd ${BUILDDIR}
CXX=g++ CXXFLAGS="-g -O2" cmake -Disa=avx -Dtarget_cxx=g++ -Dtarget_cxxflags="-march=sandybridge -O3" -DCMAKE_INSTALL_PREFIX=${INSTALLDIR} ${SRCDIR}
make VERBOSE=1 -j 10 
