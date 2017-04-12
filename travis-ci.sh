#!/bin/bash
# Copyright © 2017 Martin Ueding <dev@martin-ueding.de>

set -e
set -u
set -x



# The submodules are downloaded via SSH. This means that there has to be some
# SSH key registered with GitHub. The Travis CI virtual machine does not have
# any key. Therefore this will download a SSH key which is registered as a
# read-only deploy key with that repository. The virtual machine will then use
# this key to authenticate against GitHub. Only then it can download the public
# readable submodules. Another approach would be to switch to HTTPS for the
# submodules.
#wget -O ~/.ssh/id_rsa https://raw.githubusercontent.com/martin-ueding/ssh-access-dummy/master/dummy
#wget -O ~/.ssh/id_rsa.pub https://raw.githubusercontent.com/martin-ueding/ssh-access-dummy/master/dummy.pub
#chmod 0600 ~/.ssh/id_rsa
#chmod 0600 ~/.ssh/id_rsa.pub

# There is `#pragma omp simd` in the code. That needs OpenMP 4.0. That is
# supported from GCC 4.9. In Ubuntu Trusty, which is used by Travis CI, there
# is only 4.8. Therefore the newer version of GCC needs to be installed.
#sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test
sudo apt-get update
#sudo apt-get install -y gcc-4.9 g++-4.9

sudo apt-get install -y python3-jinja2

mkdir -p qphix-output/.git

./generate-and-copy-to-qphix scalar qphix-output
