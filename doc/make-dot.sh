#!/bin/bash
# Copyright © 2017 Martin Ueding <dev@martin-ueding.de>

set -e
set -u

for i in ./*.dot; do
    dot "$i" -Tpdf -o "${i%.*}.pdf"
done
