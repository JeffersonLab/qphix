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

mkdir -p qphix-output/.git

./generate-and-copy-to-qphix scalar qphix-output
