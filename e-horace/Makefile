## change the command to reflect your LHAPDF installation
## the report numbers have been obtained with LHAPDF 5.9.1
eMELALIB = -L/usr/local/lib/lib/ -leMELA
LHAPDFLIB= -L/usr/local/lib/ -lLHAPDF

## USE FIRST LINE FOR REPORT NUMBERS, OTHERWISE THE 2nd ONE
REPORT=-DREPORT
#REPORT=

## the compilers must be specifiend in the file compilers.in
include compilers.in

#DEBUG = 1

FFLAGS = -O3 -funroll-loops

FFLAGSdebug = -g -fbounds-check

HORVER       = 3.2
EXE          = horace
LIBSVDIR     = libvirtual
LIBZSVDIR    = libzvirtual

LIBSV        = $(LIBSVDIR)/squared_me.a
LIBSVUU      = $(LIBZSVDIR)/uubar/squared_me.a
LIBSVDD      = $(LIBZSVDIR)/ddbar/squared_me.a

#LOOPTOOLSDIR = LoopTools2-withD0cm
LOOPTOOLSDIR = LoopTools
LOOPTOOLS    = $(LOOPTOOLSDIR)/lib/libooptools.a

#CUBADIR      = Cuba-1.2
CUBADIR      = Cuba

CUBA         = -L$(CUBADIR) -lcuba
CUBALIB      = $(CUBADIR)/libcuba.a

EXTRADEPS = shared.inc compilers.in

OBJECTS = matching.o routines.o PDF.o event.o distributions.o\
          subtraction.o cuts.o QEDPS.o storage.o pdfsampling.o\
          hadr5n16.o sampling.o sq.o sqz.o phind.o\
          virtual.o softphotonintegral.o init.o ranlux.o

SOURCES = matching.f routines.f PDF.f event.f distributions.f\
          subtraction.f cuts.F QEDPS.f storage.f pdfsampling.f\
          hadr5n16.f sampling.f sq.f sqz.f phind.f\
          virtual.F softphotonintegral.f init.F ranlux.f\
          form/formme.f  form/bornme.f form/bornz.f form/ggborn.f form/zrad.f\
          main.f *.inc compilers.in

EXTRASOURCES = Makefile input-report-* libvirtual/ libzvirtual/ $(LOOPTOOLSDIR)/\
               $(CUBADIR)/ README-horace-*

ifeq ($(DEBUG),1)
 FFLAGS=$(FFLAGSdebug)
endif

F77 += $(FFLAGS)

LDMSG1 = 'Before running, make sure the environment variable LD_LIBRARY_PATH includes'
LDMSG2 = 'the LHAPDF library path. '

default: $(EXE)
debug:   $(EXE)

clean:
	rm -f *.o

deepclean: clean
	rm -f *.so $(EXE) *~ *H*.a

cleancuba:
	cd $(CUBADIR) && make distclean

verydeepclean: deepclean
	cd $(LOOPTOOLSDIR) && make clean && rm -rf lib/ bin/
	cd $(LIBZSVDIR)/ddbar/ && make clean
	cd $(LIBZSVDIR)/uubar/ && make clean
	cd $(LIBSVDIR)/ && make clean
	make cleancuba

pack:
	mkdir -p release/horace-$(HORVER)-report/ &&\
        cp -ra --parents $(SOURCES) $(EXTRASOURCES) release/horace-$(HORVER)-report/ \
        && cd release/horace-$(HORVER)-report/ &&\
        make verydeepclean && cd .. &&\
	tar -cjvf horace-$(HORVER)-report.tar.bz2 horace-$(HORVER)-report/ &&\
	rm -rf horace-$(HORVER)-report/

$(LIBSV): $(LIBSVDIR)/*.[Fh] compilers.in
	cd $(LIBSVDIR) && $(MAKE) squared_me.a

$(LIBSVUU): $(LIBZSVDIR)/uubar/*.[Fh] compilers.in
	cd $(LIBZSVDIR)/uubar/ && $(MAKE) squared_me.a

$(LIBSVDD): $(LIBZSVDIR)/ddbar/*.[Fh] compilers.in
	cd $(LIBZSVDIR)/ddbar/ && $(MAKE) squared_me.a

$(LOOPTOOLS): compilers.in
	cd $(LOOPTOOLSDIR) && ./configure --prefix=. && make install

$(CUBALIB): compilers.in
	cd $(CUBADIR) && ./configure && make libcuba.a

matching.o: matching.f $(EXTRADEPS)
	$(F77) -fPIC -c matching.f
main.o: main.f $(EXTRADEPS)
	$(F77) -fPIC -c main.f
routines.o: routines.f $(EXTRADEPS)
	$(F77) -fPIC -c routines.f
subtraction.o: subtraction.f $(EXTRADEPS)
	$(F77) -fPIC -c subtraction.f
PDF.o: PDF.f $(EXTRADEPS)
	$(F77) -fPIC -c PDF.f
pdfsampling.o: pdfsampling.f $(EXTRADEPS)
	$(F77) -fPIC -c pdfsampling.f
QEDPS.o: QEDPS.f $(EXTRADEPS)
	$(F77) -fPIC -c QEDPS.f
event.o: event.f $(EXTRADEPS)
	$(F77) -fPIC -c event.f
distributions.o: distributions.f $(EXTRADEPS) nbin.inc
	$(F77) -fPIC -c distributions.f
cuts.o: cuts.F $(EXTRADEPS) Makefile
	$(F77) -fPIC $(REPORT) -c cuts.F
storage.o: storage.f $(EXTRADEPS)
	$(F77) -fPIC -c storage.f
init.o: init.F $(EXTRADEPS) Makefile
	$(F77) -I$(LOOPTOOLSDIR)/include -fPIC $(REPORT) -c init.F

hadr5n16.o: hadr5n16.f $(EXTRADEPS)
	$(F77) -fPIC -c hadr5n16.f
sampling.o: sampling.f $(EXTRADEPS)
	$(F77) -fPIC -c sampling.f
ranlux.o: ranlux.f $(EXTRADEPS)
	$(F77) -fPIC -c ranlux.f
sq.o: sq.f $(EXTRADEPS) form/formme.f
	$(F77) -fPIC -c sq.f
sqz.o: sqz.f $(EXTRADEPS) form/bornz.f form/zrad.f form/ggborn.f
	$(F77) -fPIC -c sqz.f
phind.o: phind.f $(EXTRADEPS) form/zrad.f form/formme.f
	$(F77) -fPIC -c phind.f
virtual.o: virtual.F $(EXTRADEPS)
	$(F77) -fPIC -ffixed-line-length-132 -I$(LOOPTOOLSDIR)/include -I$(LIBSVDIR) -I. -c virtual.F
softphotonintegral.o: softphotonintegral.f $(EXTRADEPS)
	$(F77) -fPIC -c softphotonintegral.f

# the following are not released
# smear.o: smear.f $(EXTRADEPS)
# 	$(F77) -fPIC -c smear.f
# smearandrecombine.o: smearandrecombine.f $(EXTRADEPS)
# 	$(F77) -fPIC -c smearandrecombine.f
# PDFsets.o: PDFsets.f $(EXTRADEPS)
#	$(F77) -fPIC -c PDFsets.f

libHORACE.a: $(OBJECTS)
	ar cr libHORACE.a $(OBJECTS)

forceLT:
	cd $(LOOPTOOLSDIR)/ && make force && make install

$(EXE): $(LOOPTOOLS) $(CUBALIB) $(LIBSV) $(LIBSVUU) $(LIBSVDD) libHORACE.a main.o
	$(F77) main.o -L. -lHORACE $(eMELALIB) $(LHAPDFLIB) $(LIBSV) $(LIBSVUU) $(LIBSVDD) $(LOOPTOOLS) $(CUBA) -o $(EXE)

ifeq ($(REPORT),-DREPORT)
	@echo " "
	@echo " "
	@echo "Horace compiled with **report setup** (hard wired cuts and some input parameters)."
	@echo "If this is not what you intended, undefine the REPORT variable in Makefile."
	@echo " "
	@echo $(LDMSG1)
	@echo $(LDMSG2)
	@echo " "
	@echo " "
else
	@echo " "
	@echo " "
	@echo "Horace compiled."
	@echo " "
	@echo $(LDMSG1)
	@echo $(LDMSG2)
	@echo " "
	@echo " "
endif
