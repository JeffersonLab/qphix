#!/bin/bash
# Copyright © 2017 Martin Ueding <dev@martin-ueding.de>

# Generates the code for a given ISA and then copies the generated code to the
# QPhiX repository.

if (( $# == 0 )); then
    echo "Call it like this: $0 ISA"
    exit 1
fi

set -e
set -u
set -x

isa="$1"
qphix="../qphix"

if ! [[ -d "$qphix/.git" ]]; then
    echo "The QPhiX Git repository could not be found at $qphix. Please make sure it is there or change the variable in this script."
    exit 1
fi

# Generate the kernel code.
make -j "$(nproc)" "$isa"

# Generate the specialization code.
pushd jinja
./generate_files.py "$isa"
popd

# Copy the kernel files to QPhiX.
mkdir -p "$qphix/include/qphix/$isa/generated"
cp "$isa/"* "$qphix/include/qphix/$isa/generated/"

# Copy the specialization files to QPhiX.
cp "jinja/$isa/"* "$qphix/include/qphix/$isa/"
cp "jinja/"*.h "$qphix/include/qphix/"
