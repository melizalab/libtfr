FROM quay.io/pypa/manylinux1_x86_64

RUN yum install -y fftw3-devel gfortran

RUN wget http://www.netlib.org/lapack/lapack-3.7.1.tgz \
 && echo "84c4f7163b52b1bf1f6ca2193f6f48ed3dec0fab  lapack-3.7.1.tgz" | sha1sum -c - \
 && mkdir -p /usr/src/lapack \
 && tar zxvf lapack-3.7.1.tgz -C /usr/src/lapack --strip-components=1 \
 && rm lapack-3.7.1.tgz

WORKDIR /usr/src/lapack
RUN cp make.inc.example make.inc \
 && sed -i -e "s/frecursive$/frecursive -fPIC/" make.inc \
 && sed -i -e "s/SRC$/SRC double/" Makefile \
 && make -j $(nproc) lapack_install lib blaslib

WORKDIR /root
CMD /io/dev/build-wheels.sh
