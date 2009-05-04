# Makefile 
#
#

CC=/usr/bin/gcc
LD=/usr/bin/ld
CFLAGS = -fPIC -Wall -O2 -g2 -std=c99
# Default is to try to link with ATLAS lapack library; this won't work on OS X (try SCONS)
# If you don't have LAPACK enable the NO_LAPACK compile flag
LDFLAGS = -lm -llapack -lcblas -latlas -lfftw3
#CFLAGS = $(CFLAGS) -DNO_LAPACK

LIBNAME = tfrspec

SRCS = mtm.c tfr.c
HDR = tfr.h
OBJS = $(SRCS:.c=.o)

PREFIX=/usr/local
MODE=644
OWNER=root
GROUP=wheel
INSTALL=install -o ${OWNER} -g ${GROUP}

%.o : %.c
	$(CC) $(CFLAGS) -c $<

lib: $(OBJS)
	$(CC) -shared -Wl,-soname,lib${LIBNAME}.so -o lib${LIBNAME}.so $(LDFLAGS) $(OBJS)
	ar rc lib${LIBNAME}.a $(OBJS)
	ranlib lib${LIBNAME}.a

install: lib$(LIBNAME).so
	$(INSTALL) -m $(MODE) lib$(LIBNAME).so $(PREFIX)/lib
	$(INSTALL) -m $(MODE) lib$(LIBNAME).a $(PREFIX)/lib	
	$(INSTALL) -m $(MODE) $(HDR) $(PREFIX)/include

clean:
	rm -f $(OBJS) lib$(LIBNAME).a lib$(LIBNAME).so