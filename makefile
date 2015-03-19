# N.B. the makefile will not work if there are spaces or tabs after variables, in their defintion
SDIR		   	= src
HDIR			= include
ODIR			= objs
CSDIR			= csrc
CODIR			= cobjs
TSDIR			= tests
CC 				= g++
OPTIM 			= 
CFLAGS 			= -Wall -g -std=c++0x -static
#CFLAGS EXPLAINED:
#-std=c++0x 		: added so that auto lambda functions can be used
#-std=c++11 		: for more modern std c++
#-g 				: for using gdb to debug
#-ggdb 				: for using gdb to debug - not sure how this differs from above
#-Wall 				: turns on most compiler warnings
#-O 				: does some optimizations for running (compiling takes longer though)
#-Os 				: optimises as long as code size isn't increased
#-O2 				: does some more optimizations that -O
#-O3 				: does all the optimizations for running
#-static 			: linking to libraries statically
#-ansi -pedantic 	: turns off some extensions of g++ which are incompatible with the ANSI language
#-fopenmp 			: so can use openmp parallelisation
#-pg 				: also known as gprof, the gcc profiling tool 
LFLAGS 			= 
INCLUDES		= -I$(HDIR) -I/home/og/Documents/c++/eigen-eigen-1306d75b4a21/ -I/home/og/Documents/c++/gnuplot/gnuplot-cpp
LIBS 			= -lm -lgsl -lgslcblas

_HEADERS 		= check.h error.h fnptrs.h folder.h gsl_extras.h lattice.h omega.h parameters.h potentials.h print.h simple.h\
					sphaleron_fns.h stepper.h
HEADERS 		= $(patsubst %,$(HDIR)/%,$(_HEADERS))

_COMMONSRC		= check.cc error.cc folder.cc gsl_extras.cc lattice.cc omega.cc parameters.cc potentials.cc print.cc simple.cc\
					sphaleron_fns.cc stepper.cc
_COMMONOBJS		= $(_COMMONSRC:.cc=.o)
COMMONSRC		= $(patsubst %,$(CSDIR)/%,$(_COMMONSRC))
COMMONOBJS 		= $(patsubst %,$(CODIR)/%,$(_COMMONOBJS))

#------------------------------------------------------------------------------------------------------------------------

all: common main pi sphaleron sphaleron_pi	

common: $(COMMONOBJS)
	@echo made common objects $(COMMONOBJS)
	
main: $(ODIR)/main.o $(COMMONOBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)
	@echo Simple compiler named $(MAIN) has been compiled
	
pi: $(ODIR)/pi.o $(COMMONOBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)
	@echo Simple compiler named $(PI) has been compiled
	
sphaleron: $(ODIR)/sphaleron.o $(COMMONOBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)
	@echo Simple compiler named $(PI) has been compiled
	
sphaleron_pi: $(ODIR)/sphaleron_pi.o $(COMMONOBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)
	@echo Simple compiler named $(PI) has been compiled
	
#------------------------------------------------------------------------------------------------------------------------
	
$(TSDIR)/%: $(ODIR)/%.o $(COMMONOBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)
	
$(SDIR)/%: $(ODIR)/%.o $(COMMONOBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)
	
$(CSDIR)/%: $(CSDIR)/%.cc
	$(CC) -c -o $@ $< $(CFLAGS) $(INCLUDES) $(LIBS)
	
$(CODIR)/%.o: $(CSDIR)/%.cc
	$(CC) -c -o $@ $< $(CFLAGS) $(INCLUDES) $(LIBS)
	
$(ODIR)/%.o: $(SDIR)/%.cc
	$(CC) -c -o $@ $< $(CFLAGS) $(INCLUDES) $(LIBS)
	
$(ODIR)/%.o: $(TSDIR)/%.cc
	$(CC) -c -o $@ $< $(CFLAGS) $(INCLUDES) $(LIBS)

#------------------------------------------------------------------------------------------------------------------------
	
.PHONY: clean

clean:
	rm -f $(ODIR)/*.o
	rm -f $(CODIR)/*.o
	rm -f $(TSDIR)/*.o
	rm -f *~ core
	rm -f $(HDIR)/*~
	rm -f $(TSDIR)/testFolder
	rm -f $(TSDIR)/testGsl_extras
	rm -f $(TSDIR)/testLattice
	rm -f $(TSDIR)/testParameters
	rm -f $(TSDIR)/testPotentials
	rm -f $(TSDIR)/testPrint

#------------------------------------------------------------------------------------------------------------------------
	
depend: $(SRC)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it
