#!/bin/bash
# This script is used to build manylinux wheels. It should be run in a docker container:
# docker run --rm -v `pwd`:/io dmeliza/libtfr-manylinux /io/dev/build-wheels.sh
set -e -x
export STATIC_LAPACK=1

# Compile wheels
for PYBIN in /opt/python/*/bin; do
    "${PYBIN}/pip" install -r /io/requirements.txt
    "${PYBIN}/pip" wheel /io/ -w wheelhouse/
done

# Bundle external shared libraries into the wheels
for whl in wheelhouse/*.whl; do
    auditwheel repair "$whl" -w /io/wheelhouse/
done

# Install packages and test
for PYBIN in /opt/python/*/bin; do
    "${PYBIN}/pip" install libtfr --no-index -f /io/wheelhouse
    (cd "$HOME"; "${PYBIN}/nosetests" -w /io/test/)
done
