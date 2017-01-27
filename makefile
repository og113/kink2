# makefile for compiling programs to do kink calculations
#------------------------------------------------------------------------------------------------------------------------
# config file

include config.mk

#------------------------------------------------------------------------------------------------------------------------
# local variables

HEADERS 		= $(wildcard $(HDIR)/*.h)
FNSRC			= $(wildcard $(FNSDIR)/*.cc)
FNOBJS 			= $(patsubst $(FNSDIR)/%.cc,$(FNODIR)/%.o,$(FNSRC))
COMMONSRC		= $(wildcard $(CSDIR)/*.cc) 							# this needs fixing to exclude the fnobjs
COMMONOBJS 		= $(patsubst $(CSDIR)/%.cc,$(CODIR)/%.o,$(COMMONSRC))	# this needs fixing to exclude the fnobjs
SRC				= $(wildcard $(SDIR)/*.cc)
EXE				= $(patsubst $(SDIR)/%.cc,%,$(SRC))
MPISRC			= $(wildcard $(MPISDIR)/*.cc)
MPIEXE			= $(patsubst $(MPISDIR)/%.cc,%,$(MPISRC))

.PHONY : variables
variables :
	@echo HEADERS: $(HEADERS)
	@echo COMMONSRC: $(COMMONSRC)
	@echo COMMONOBJS: $(COMMONOBJS)
	@echo FNSRC: $(FNSRC)
	@echo FNOBJS: $(FNOBJS)
	@echo SRC: $(SRC)
	@echo EXE: $(EXE)
	@echo MPISRC: $(MPISRC)
	@echo MPIEXE: $(MPIEXE)

#------------------------------------------------------------------------------------------------------------------------
# some useful PHONYs

.PHONY: all
all: $(EXE) $(MPIEXE) common fns

.PHONY: common
common: $(COMMONOBJS)

.PHONY: fns
fns: $(FNOBJS)

#------------------------------------------------------------------------------------------------------------------------
# targets, dependencies and rules for executables

analysis: $(ODIR)/analysis.o $(COMMONOBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)
	@echo Simple compiler named analysis has been compiled
	
cleanData: $(ODIR)/cleanData.o $(COMMONOBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)
	@echo Simple compiler named cleanData has been compiled
	
changeInputs: $(ODIR)/changeInputs.o $(COMMONOBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)
	@echo Simple compiler named changeInputs has been compiled

common: $(COMMONOBJS)
	@echo made common objects $(COMMONOBJS)
	
fns: $(FNOBJS)
	@echo made common objects $(FNOBJS)
	
main: $(ODIR)/main.o $(FNODIR)/main_fn.o $(COMMONOBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)
	@echo Simple compiler named main has been compiled
	
<<<<<<< HEAD
main3: $(ODIR)/main3.o $(FNODIR)/main_fn3.o $(COMMONOBJS)
=======
main3: $(ODIR)/main3.o  $(FNODIR)/main_fn3.o $(COMMONOBJS)
>>>>>>> 0e49d8046b835f27d461ba8917ba472e6150864b
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)
	@echo Simple compiler named main3 has been compiled
	
plot: $(ODIR)/plot.o $(COMMONOBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)
	@echo Simple compiler named plot has been compiled
	
pi: $(ODIR)/pi.o $(COMMONOBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)
	@echo Simple compiler named pi has been compiled
	
piEvolve: $(ODIR)/piEvolve.o $(COMMONOBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)
	@echo Simple compiler named piEvolve has been compiled
	
sphaleron: $(ODIR)/sphaleron.o $(FNODIR)/sphaleron_fns.o $(COMMONOBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)
	@echo Simple compiler named sphaleron has been compiled
	
sphaleronPi: $(ODIR)/sphaleronPi.o $(COMMONOBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)
	@echo Simple compiler named sphaleronPi has been compiled
	
staticNewton: $(ODIR)/staticNewton.o $(COMMONOBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)
	@echo Simple compiler named staticNewton has been compiled
	
staticShooting: $(ODIR)/staticShooting.o $(COMMONOBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)
	@echo Simple compiler named staticShooting has been compiled
	
wrapper: $(MPIODIR)/wrapper.o $(FNODIR)/main_fn.o $(COMMONOBJS)
	$(MPICC) -o $@ $^ $(MPICFLAGS) $(INCLUDES) $(MPILIBS)
	@echo Simple compiler named wrapper has been compiled
	
wrapper2: $(MPIODIR)/wrapper2.o $(FNODIR)/main_fn.o $(COMMONOBJS)
	$(MPICC) -o $@ $^ $(MPICFLAGS) $(INCLUDES) $(MPILIBS)
	@echo Simple compiler named wrapper2 has been compiled
	
#------------------------------------------------------------------------------------------------------------------------
	
$(MPIODIR)/%.o: $(MPISDIR)/%.cc
	$(MPICC) -c -o $@ $< $(MPICFLAGS) $(INCLUDES) $(MPILIBS)
	
$(TSDIR)/%: $(ODIR)/%.o $(COMMONOBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)
	
$(SDIR)/%: $(ODIR)/%.o $(COMMONOBJS) 
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)
	
$(CSDIR)/%: $(CSDIR)/%.cc
	$(CC) -c -o $@ $< $(CFLAGS) $(INCLUDES) $(LIBS)
	
$(CODIR)/%.o: $(CSDIR)/%.cc
	$(CC) -c -o $@ $< $(CFLAGS) $(INCLUDES) $(LIBS)
	
$(FNSDIR)/%: $(FNSDIR)/%.cc $(COMMONOBJS)
	$(CC) -c -o $@ $< $(CFLAGS) $(INCLUDES) $(LIBS)
	
$(FNODIR)/%.o: $(FNSDIR)/%.cc $(COMMONOBJS)
	$(CC) -c -o $@ $< $(CFLAGS) $(INCLUDES) $(LIBS)
	
$(ODIR)/%.o: $(SDIR)/%.cc
	$(CC) -c -o $@ $< $(CFLAGS) $(INCLUDES) $(LIBS)
	
$(ODIR)/%.o: $(TSDIR)/%.cc
	$(CC) -c -o $@ $< $(CFLAGS) $(INCLUDES) $(LIBS)

#------------------------------------------------------------------------------------------------------------------------
# clean
	
.PHONY: clean
clean:
	rm -f $(ODIR)/*.o
	rm -f $(CODIR)/*.o
	rm -f $(FNODIR)/*.o
	rm -f $(MPIODIR)/*.o
	rm -f $(TSDIR)/*.o
	rm -f *~ core
	rm -f $(HDIR)/*~
	rm -f data/temp/*
	rm -f temp/*

#------------------------------------------------------------------------------------------------------------------------
# makedepend, NOT WORKING

.PHONY: depend	
depend: $(COMMONSRC) $(SRC) $(FNSRC)
	makedepend -- $(CFLAGS) -- $^ -f- > Makefile.deps

# DO NOT DELETE THIS LINE -- make depend needs it
