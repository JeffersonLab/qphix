#!/bin/bash
# Copyright © 2017 Martin Ueding <dev@martin-ueding.de>

# Compiles and tests QPhiX on the Travis CI infrastructure.

set -e
set -u
set -x

fold_start() {
echo -en 'travis_fold:start:'$1'\\r'
}

fold_end() {
echo -en 'travis_fold:end:'$1'\\r'
}

fold_start more_cpu_info
lscpu
cat /proc/cpuinfo
fold_end more_cpu_info

# The submodules are downloaded via SSH. This means that there has to be some
# SSH key registered with GitHub. The Travis CI virtual machine does not have
# any key. Therefore this will download a SSH key which is registered as a
# read-only deploy key with that repository. The virtual machine will then use
# this key to authenticate against GitHub. Only then it can download the public
# readable submodules. Another approach would be to switch to HTTPS for the
# submodules.
fold_start get_ssh_key
wget -O ~/.ssh/id_rsa https://raw.githubusercontent.com/martin-ueding/ssh-access-dummy/master/dummy
wget -O ~/.ssh/id_rsa.pub https://raw.githubusercontent.com/martin-ueding/ssh-access-dummy/master/dummy.pub
chmod 0600 ~/.ssh/id_rsa
chmod 0600 ~/.ssh/id_rsa.pub
fold_end get_ssh_key

cd ..

# There is `#pragma omp simd` in the code. That needs OpenMP 4.0. That is
# supported from GCC 4.9. In Ubuntu Trusty, which is used by Travis CI, there
# is only 4.8. Therefore the newer version of GCC needs to be installed.
fold_start update_gcc
sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test
sudo apt-get update
sudo apt-get install -y gcc-6 g++-6 ccache libopenmpi-dev

ls -l /usr/lib/ccache
fold_end update_gcc

basedir=$PWD

cc_name=mpicc
cxx_name=mpic++

export OMPI_CC=gcc-6
export OMPI_CXX=g++-6

color_flags=""
openmp_flags="-fopenmp"
base_flags="-O2 -finline-limit=50000 $color_flags"
cxx11_flags="--std=c++11"
disable_warnings_flags="-Wno-all -Wno-pedantic"
qphix_flags="-Drestrict=__restrict__"
qphix_configure=""

sourcedir="$basedir/sources"
mkdir -p "$sourcedir"

# Directory for the installed files (headers, libraries, executables).
prefix="$basedir/local"
mkdir -p "$prefix"

# Directory for building. The GNU Autotools support out-of-tree builds which
# allow to use different compilers on the same codebase.
build="$basedir/build"
mkdir -p "$build"

mv qphix $sourcedir/

PATH=$prefix/bin:$PATH

base_cxxflags="$base_flags"
base_cflags="$base_flags"
base_configure="--prefix=$prefix --disable-shared --enable-static CC=$(which $cc_name) CXX=$(which $cxx_name)"

# Clones a git repository if the directory does not exist. It does not call
# `git pull`. After cloning, it deletes the `configure` and `Makefile` that are
# shipped by default such that they get regenerated in the next step.
clone-if-needed() {
    local url="$1"
    local dir="$2"
    local branch="$3"

    if ! [[ -d "$dir" ]]
    then
        git clone "$url" --recursive -b "$branch"

        pushd "$dir"
        rm -f configure Makefile
        popd
    fi
}

# If the user has not given a variable `SMP` in the environment, use as many
# processes to compile as there are cores in the system.
make_smp_template="-j $(nproc)"
make_smp_flags="${SMP-$make_smp_template}"

# Runs `make && make install` with appropriate flags that make compilation
# parallel on multiple cores. A sentinel file is created such that `make` is
# not invoked once it has correctly built.
make-make-install() {
    if ! [[ -f build-succeeded ]]; then
        fold_start $repo.make
        make $make_smp_flags
        fold_end $repo.make

        fold_start $repo.make_install
        make install
        fold_end $repo.make_install

        touch build-succeeded
        pushd $prefix/lib
        rm -f *.so *.so.*
        popd
    fi
}

# Prints a large heading such that it is clear where one is in the compilation
# process. This is not needed but occasionally helpful.
print-fancy-heading() {
    set +x
    echo "######################################################################"
    echo "# $*"
    echo "######################################################################"
    set -x

    if [[ -d "$sourcedir/$repo" ]]; then
        pushd "$sourcedir/$repo"
        git branch
        popd
    fi
}

# I have not fully understood this here. I *feel* that there is some cyclic
# dependency between `automake --add-missing` and the `autoreconf`. It does not
# make much sense. Perhaps one has to split up the `autoreconf` call into the
# parts that make it up. Using this weird dance, it works somewhat reliably.
autotools-dance() {
    #automake --add-missing --copy || autoreconf -f || automake --add-missing --copy
    autoreconf -vif
}

# Invokes the various commands that are needed to update the GNU Autotools
# build system. Since the submodules are also Autotools projects, these
# commands need to be invoked from the bottom up, recursively. The regular `git
# submodule foreach` will do a traversal from the top. Due to the nested nature
# of the GNU Autotools, we need to have depth-first traversal. Assuming that
# the directory names do not have anything funny in them, the parsing of the
# output can work.
autoreconf-if-needed() {
    if ! [[ -f configure ]]; then
        if [[ -f .gitmodules ]]; then
            for module in $(git submodule foreach --quiet --recursive pwd | tac); do
                pushd "$module"
                aclocal
                autotools-dance
                popd
            done
        fi

        aclocal
        autotools-dance
    fi
}

cd "$sourcedir"

###############################################################################
#                                   libxml2                                   #
###############################################################################


repo=libxml2
fold_start $repo.download
print-fancy-heading $repo
clone-if-needed https://git.gnome.org/browse/libxml2 $repo v2.9.4
fold_end $repo.download

fold_start $repo.autoreconf
pushd $repo
cflags="$base_cflags"
cxxflags="$base_cxxflags"
if ! [[ -f configure ]]; then
    mkdir -p m4
    pushd m4
    ln -fs /usr/share/aclocal/pkg.m4 .
    popd
    NOCONFIGURE=yes ./autogen.sh
fi
popd
fold_end $repo.autoreconf

fold_start $repo.configure
mkdir -p "$build/$repo"
pushd "$build/$repo"
if ! [[ -f Makefile ]]; then
    if ! $sourcedir/$repo/configure $base_configure \
            --without-zlib \
            --without-python \
            --without-readline \
            --without-threads \
            --without-history \
            --without-reader \
            --without-writer \
            --with-output \
            --without-ftp \
            --without-http \
            --without-pattern \
            --without-catalog \
            --without-docbook \
            --without-iconv \
            --without-schemas \
            --without-schematron \
            --without-modules \
            --without-xptr \
            --without-xinclude \
            CFLAGS="$cflags" CXXFLAGS="$cxxflags"; then
        cat config.log
        exit 1
    fi
fi
fold_end $repo.configure
make-make-install
popd

###############################################################################
#                                     QMP                                     #
###############################################################################

repo=qmp
fold_start $repo.download
print-fancy-heading $repo
clone-if-needed https://github.com/usqcd-software/qmp.git $repo master
fold_end $repo.download

fold_start $repo.autoreconf
pushd $repo
cflags="$base_cflags $openmp_flags --std=c99"
cxxflags="$base_cxxflags $openmp_flags $cxx11_flags"
autoreconf-if-needed
popd
fold_end $repo.autoreconf

fold_start $repo.configure
mkdir -p "$build/$repo"
pushd "$build/$repo"
if ! [[ -f Makefile ]]; then
    if ! $sourcedir/$repo/configure $base_configure \
            --with-qmp-comms-type=MPI \
            CFLAGS="$cflags" CXXFLAGS="$cxxflags"; then
        cat config.log
        exit 1
    fi
fi
fold_end $repo.configure
make-make-install
popd

###############################################################################
#                                    QDP++                                    #
###############################################################################

repo=qdpxx
fold_start $repo.download
print-fancy-heading $repo
clone-if-needed https://github.com/usqcd-software/qdpxx.git $repo devel
fold_end $repo.download

fold_start $repo.autoreconf
pushd $repo
cflags="$base_cflags $openmp_flags --std=c99"
cxxflags="$base_cxxflags $openmp_flags $cxx11_flags"
autoreconf-if-needed
popd
fold_end $repo.autoreconf

fold_start $repo.configure
mkdir -p "$build/$repo"
pushd "$build/$repo"
if ! [[ -f Makefile ]]; then
    if ! $sourcedir/$repo/configure $base_configure \
            --enable-openmp \
            --enable-sse --enable-sse2 \
            --enable-parallel-arch=parscalar \
            --enable-precision=double \
            --with-qmp="$prefix" \
            --with-libxml2="$prefix/bin/xml2-config" \
            CFLAGS="$cflags" CXXFLAGS="$cxxflags"; then
        cat config.log
        exit 1
    fi
fi
fold_end $repo.configure
make-make-install
popd

###############################################################################
#                                    QPhiX                                    #
###############################################################################

repo=qphix
print-fancy-heading $repo

fold_start $repo.autoreconf
pushd $repo
cflags="$base_cflags $openmp_flags $qphix_flags"
cxxflags="$base_cxxflags $openmp_flags $cxx11_flags $qphix_flags"
autoreconf-if-needed
popd
fold_end $repo.download

fold_start $repo.configure
case "$QPHIX_ARCH" in
    SCALAR)
        archflag=
        soalen=1
        ;;
    AVX)
        archflag=-march=sandybridge
        soalen=2
        ;;
    AVX2)
        archflag=-march=haswell
        soalen=2
        ;;
    AVX512)
        archflag=-march=knl
        soalen=4
        ;;
    *)
        echo "Unsupported QPHIX_ARCH"
        exit 1;
esac

mkdir -p "$build/$repo"
pushd "$build/$repo"
if ! [[ -f Makefile ]]; then
    if ! $sourcedir/$repo/configure $base_configure \
            $qphix_configure \
            --enable-testing \
            --enable-proc=$QPHIX_ARCH \
            --enable-soalen=$soalen \
            --enable-clover \
            --enable-twisted-mass \
            --enable-tm-clover \
            --enable-openmp \
            --enable-mm-malloc \
            --enable-parallel-arch=parscalar \
            --with-qdp="$prefix" \
            CFLAGS="$cflags $archflag" CXXFLAGS="$cxxflags $archflag"; then
        cat config.log
        exit 1
    fi
fi
fold_end $repo.configure
make-make-install
popd

###############################################################################

# Only run the tests on architectures that are supported by Travis CI.
case "$QPHIX_ARCH" in
    SCALAR)
        ;;
    AVX)
        ;;
    AVX2)
        exit 0
        ;;
    AVX512)
        exit 0
        ;;
esac

export OMP_NUM_THREADS=2

pushd $build/qphix/tests

l=16
args="-by 8 -bz 8 -c 2 -sy 1 -sz 1 -pxy 1 -pxyz 0 -minct 1 -x $l -y $l -z $l -t $l -prec f -geom 1 1 1 2"

tests=(
t_clov_dslash
t_dslash
t_twm_dslash
t_twm_clover
)

for runner in "${tests[@]}"
do
    fold_start testing.$runner
    mpirun -n 2 ./$runner $args
    fold_end testing.$runner
done

