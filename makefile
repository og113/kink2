SDIR		   	= src
HDIR			= include
ODIR			= objects
TDIR			= tests
CC 				= g++
OPTIM 			= 
CFLAGS 			= -Wall -g
#CFLAGS EXPLAINED:
#-std=c++0x 		: added so that to_string and similar functions can be used
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
INCLUDES		= -I$(HDIR) -I/home/og/Documents/c++/eigen-eigen-1306d75b4a21 -I/home/og/Documents/c++/gnuplot/gnuplot-cpp
LIBS 			= -lm -lgsl -lgslcblas

_HEADERS 		= error.h fnptrs.h folder.h gsl_extras.h parameters.h potentials.h simple.h thetaT.h 
HEADERS 		= $(patsubst %,$(HDIR)/%,$(_HEADERS))

_COMMONSRC		= error.cc folder.cc gsl_extras.cc parameters.cc potentials.cc simple.cc  thetaT.cc
_COMMONOBJS		= $(_COMMONSRC:.cc=.o)
COMMONSRC		= $(patsubst %,$(SDIR)/%,$(_COMMONSRC))
COMMONOBJS 		= $(patsubst %,$(ODIR)/%,$(_COMMONOBJS))

_MAINSRC		= main.cc
_MAINOBJS		= $(_MAINSRC:.cc=.o)
MAINSRC			= $(patsubst %,$(SDIR)/%,$(_MAINSRC))
MAINOBJS 		= $(patsubst %,$(ODIR)/%,$(_MAINOBJS))
MAIN			= $(_MAINSRC:.cc=)

_PISRC			= pi.cc
_PIOBJ			= $(_PISRC:.cc=.o)
PISRC			= $(patsubst %,$(SDIR)/%,$(_PISRC))
PIOBJS 			= $(patsubst %,$(ODIR)/%,$(_PIOBJS))
PI				= $(_PISRC:.cc=)

_TSRC			= testFolder.cc testPotentials.cc testThetaT.cc testGsl_extras.cc
_TOBJS			= $(_TSRC:.cc=.o)
TSRC			= $(patsubst %,$(TDIR)/%,$(_TSRC))
TOBJS	 		= $(patsubst %,$(TDIR)/%,$(_TOBJS))
T				= $(_TSRC:.cc=)

#------------------------------------------------------------------------------------------------------------------------
	
main: $(MAINOBJS) $(COMMONOBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)
	@echo Simple compiler named $(MAIN) has been compiled
	
pi: $(PIOBJS) $(COMMONOBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)
	@echo Simple compiler named $(PI) has been compiled
	
$(ODIR)/%.o: $(SDIR)/%.cc
	$(CC) -c -o $@ $< $(CFLAGS) $(INCLUDES) $(LIBS)
	
#------------------------------------------------------------------------------------------------------------------------
	
$(TDIR)/%.o: $(TDIR)/%.cc
	$(CC) -c -o $@ $< $(CFLAGS) $(INCLUDES) $(LIBS)
	
#------------------------------------------------------------------------------------------------------------------------

testFolder: $(TDIR)/testFolder
	@echo Simple compiler named $(TDIR)/testFolder has been compiled
	
$(TDIR)/testFolder: $(TDIR)/testFolder.o $(COMMONOBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)

#------------------------------------------------------------------------------------------------------------------------
	
testPotentials: $(TDIR)/testPotentials
	@echo Simple compiler named $(TDIR)/testPotentials has been compiled
	
$(TDIR)/testPotentials: $(TDIR)/testPotentials.o $(COMMONOBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)
	
#------------------------------------------------------------------------------------------------------------------------
	
testThetaT: $(TDIR)/testThetaT
	@echo Simple compiler named $(TDIR)/testThetaT has been compiled

$(TDIR)/testThetaT: $(TDIR)/testThetaT.o $(COMMONOBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)
	
#------------------------------------------------------------------------------------------------------------------------
	
testGsl_extras: $(TDIR)/testGsl_extras
	@echo Simple compiler named $(TDIR)/testGsl_extras has been compiled

$(TDIR)/testGsl_extras: $(TDIR)/testGsl_extras.o $(COMMONOBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)	
	
#------------------------------------------------------------------------------------------------------------------------
testParameters: $(TDIR)/testParameters
	@echo Simple compiler named $(TDIR)/testParameters has been compiled

$(TDIR)/testParameters: $(TDIR)/testParameters.o $(COMMONOBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)	
#------------------------------------------------------------------------------------------------------------------------
	
.PHONY: clean

clean:
	rm -f $(ODIR)/*.o
	rm -f *~ core
	rm -f $(HDIR)/*~
	rm -f $(TDIR)/testFolder
	rm -f $(TDIR)/testPotentials
	rm -f $(TDIR)/testThetaT

#------------------------------------------------------------------------------------------------------------------------
	
depend: $(SRC)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it
