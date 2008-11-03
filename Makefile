# Makefile 
#
#


###########################################################
### Comment out the following if you don't have FFTW
### (version 2.1.3 or so) installed.
###########################################################
FFTW_DEFINES = -DHAVE_FFTW
FFTW_LIBS = -lfftw3
###########################################################

MTM_DEFINES = 
MTM_TARGETS = multitaper.o
MTM_OBJECTS = multitaper.o

CC=/usr/bin/gcc
LIB_CFLAGS = -O2 -g2 -arch i386 -arch ppc

LIBRARY = sono

SRCS = sonogram.c mem.c multitaper.c fftmaster.c \
		b512.c bitrev.c dint.c dintime.c idint.c idintime.c tab.c 
HDR = sonogram.h
OBJS = $(SRCS:.c=.o)

TVER := $(shell grep BUILD version.h | cut -d' ' -f3)
BUILDVER := $(shell expr 1 + ${TVER})

PREFIX=/usr/local
MODE=644
OWNER=root
GROUP=wheel
INSTALL=install -o ${OWNER} -g ${GROUP}

CFLAGS = -fPIC -Wall ${LIB_CFLAGS} ${FFTW_DEFINES} ${MTM_DEFINES}

all: lib$(LIBRARY).a

lib$(LIBRARY).a: $(SRCS) $(HDR) ${OBJS} $(MTM_TARGETS)
	rm -f lib$(LIBRARY).a
	ar clq lib$(LIBRARY).a ${OBJS} $(MTM_OBJECTS)
	ranlib lib$(LIBRARY).a

mtm/libmtm.a:
	make -C mtm OPT_CFLAGS='${OPT_CFLAGS}'

clean:
	make -C mtm clean OPT_CFLAGS='${OPT_CFLAGS}'
	rm -f core $(OBJS) lib$(LIBRARY).a .tmp.h

spotless: clean
	make -C mtm spotless
	rm -f .depend

release:
	-rm -f version.h
	-echo "#define LIBSONO_BUILD "$(BUILDVER) > version.h
	-echo "#define LIBSONO_VERSION "`date +%Y%m%d` >> version.h

tar: 
	make spotless
	make release
	tar -C ../ --exclude libsono/misc --one-file-system -cpf ../libsono_`date +%Y%m%d`.tar libsono 
	bzip2 -9 ../libsono_`date +%Y%m%d`.tar
	mv -f ../libsono_`date +%Y%m%d`.tar.bz2 /var/www/html/software/

install: lib$(LIBRARY).a
	${INSTALL} -m ${MODE} lib$(LIBRARY).a ${PREFIX}/lib/lib$(LIBRARY).a
	cat version.h $(HDR) > .tmp.h
	${INSTALL} -m ${MODE} .tmp.h ${PREFIX}/include/$(HDR)
	rm -f .tmp.h

dep depend:
	$(CPP) -M $(CFLAGS) $(SRCS) > .depend

#
# include a dependency file if one exists
#
ifeq (.depend,$(wildcard .depend))
include .depend
endif

