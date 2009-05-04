# Makefile 
#
#

CC=/usr/bin/gcc
LD=/usr/bin/ld
CFLAGS = -fPIC -Wall -O2 -g2
LDFLAGS = -lm -llapack -lfftw3

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
	$(LD) -dylib -install_name lib${LIBNAME} -o lib${LIBNAME}.so $(LDFLAGS) $(OBJS)
	ar rc lib${LIBNAME}.a $(OBJS)
	ranlib lib${LIBNAME}.a

install: lib$(LIBNAME).so
	$(INSTALL) -m $(MODE) lib$(LIBNAME).so $(PREFIX)/lib
	$(INSTALL) -m $(MODE) lib$(LIBNAME).a $(PREFIX)/lib	
	$(INSTALL) -m $(MODE) $(HDR) $(PREFIX)/include
	