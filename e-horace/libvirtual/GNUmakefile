# --- variables defined by configure ---

LTDIR = ../LoopTools/
LOOPTOOLS = -L$(LTDIR) -looptools
LIBS = $(LOOPTOOLS) $(CERNLIB)

FFLAGS = -O3 -ffixed-line-length-132
### FFLAGS = -g -ffixed-line-length-132
include ../compilers.in
FC = $(F77) $(FFLAGS) -I. -I$(LTDIR)include 

ifdef DEBUG
FC += -DDEBUG
endif

# --- end defs by configure ---


OBJS = \
  abbr_s.o \
  abbr_angle.o \
  abbr_hel.o \
  born.o \
  sself.o \
  vert.o \
  box.o \
  cts.o

RENCONST = renconst.o

ALLOBJS = $(OBJS) $(RENCONST) squared_me.o

STDDEPS = prefix.h model.h renconst.h

POSSIBLEDEPS = \
  2to2.F 2to2.h gauss.F \
  2to3.F 2to3.h multigauss.F vegas.F \
  sm_ini.F mssm_ini.F


default: run

clean:
	$(RM) $(ALLOBJS) squared_me.a

squared_me.a: squared_me.a($(ALLOBJS))
squared_me.a($(RENCONST)): $(STDDEPS)
squared_me.a($(OBJS) squared_me.o): vars.h $(STDDEPS)

renconst.h:
	touch renconst.h

(%.o): %.F
	$(FC) -c $<
	$(AR) cr $@ $%
##	$(RM) $%

%:: %.F num.F process.h $(STDDEPS) $(POSSIBLEDEPS) squared_me.a
	$(FC) -o $@ $< squared_me.a $(LIBS)
