SDIR		   	= src
HDIR			= include
ODIR			= objects
TDIR			= tests
CC 				= g++
OPTIM 			= 
CFLAGS 			= -Wall -g
LFLAGS 			= 
INCLUDES		= -I$(HDIR)
LIBS 			= 
_SRC			= simple.cc error.cc bag.cc test.cc
_HEADERS 		= simple.h error.h bag.h
_OBJS 			= $(_SRC:.cc=.o)
SRC 			= $(patsubst %,$(SDIR)/%,$(_SRC))
HEADERS 		= $(patsubst %,$(HDIR)/%,$(_HEADERS))
OBJS 			= $(patsubst %,$(ODIR)/%,$(_OBJS))
MAIN 			= main

#------------------------------------------------------------------------------------------------------------------------

all: $(MAIN)
	@echo Simple compiler named $(MAIN) has been compiled
	
$(ODIR)/%.o: $(SDIR)/%.cc
	$(CC) -c -o $@ $< $(CFLAGS) $(INCLUDES)

$(MAIN): $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)
	
#------------------------------------------------------------------------------------------------------------------------	
testFolder: $(TDIR)/testFolder
	@echo Simple compiler named $(TDIR)/testfolder has been compiled
	
$(TDIR)/testFolder: 
	$(CC) -o $(TDIR)/testFolder $(TDIR)/testFolder.cc $(SDIR)/simple.cc $(SDIR)/error.cc $(SDIR)/folder.cc $(CFLAGS) $(INCLUDES)
	
#------------------------------------------------------------------------------------------------------------------------	
	
.PHONY: clean

clean:
	rm -f $(ODIR)/*.o
	rm -f *~ core
	rm -f $(HDIR)/*~
	rm -f $(TDIR)/testFolder
