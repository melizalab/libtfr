#
# Makefile for mtm
#
OPT_CFLAGS = -g2 -Wall
CFLAGS = ${OPT_CFLAGS}
LIBRARY = mtm

CSRCS = amtm.c
COBJS = $(SRCS:.c=.o)
FSRCS = findpeaks.f spec.f miscsubs.f 
FOBJS = $(FSRCS:.f=.o)
FPPSRCS = fieldio.F
FPPOBJS = $(FPPSRCS:.F=.o)
FSRCS2 = cosort.f random.f gamma.f beta.f tred2.f imtql2.f
FOBJS2 = $(FSRCS2:.f=.o)
OBJS = $(COBJS) $(FOBJS) $(FPPOBJS) $(FOBJS2)
SRCS = $(CSRCS) $(FSRCS) $(FPPSRCS) $(FSRCS2)

all : lib${LIBRARY}.a 

lib${LIBRARY}.a: $(OBJS) 
	rm -f lib${LIBRARY}.a
	ar clq lib${LIBRARY}.a $(OBJS) 

clean:
	-rm -f *.o libmtm.a

spotless: clean
	-rm -f .depend

dep depend:
	$(CPP) -M $(CFLAGS) $(SRCS) > .depend

#
# include a dependency file if one exists
#
ifeq (.depend,$(wildcard .depend))
include .depend
endif

