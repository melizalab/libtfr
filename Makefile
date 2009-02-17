# Makefile 
#
#

CC=/usr/bin/gcc
LD=/usr/bin/ld
ARCHFLAGS = -arch i386
CFLAGS = -fPIC -Wall -O2 -g2 ${ARCHFLAGS}
#LDFLAGS = -lm -llapack -lfftw3_threads -lfftw3 -lpthread
LDFLAGS = -lm -llapack -lfftw3

LIBNAME = sono2
DYNLIB_FLAGS = -shared -Wl,-soname,lib${LIBNAME}.so.1 -o lib${LIBNAME}.so.1.0.0

SRCS = mtm.c tfr.c sonogram.c
HDR = mtm.h sonogram.h
OBJS = $(SRCS:.c=.o)

PREFIX=/usr/local
MODE=644
OWNER=root
GROUP=wheel
INSTALL=install -o ${OWNER} -g ${GROUP}

%.o : %.c
	$(CC) $(CFLAGS) -c $<

lib: $(OBJS)
	$(LD) -dylib $(ARCHFLAGS) -install_name lib${LIBNAME} -o lib${LIBNAME}.dylib $(LDFLAGS) $(OBJS)

test: lib test.o
	$(CC) $(CFLAGS) $(LDFLAGS) -ldataio $(OBJS) test.o -o sono_test