# special variables
CC     = g++ 
OFLAGS  = -O3
LFLAGS  = 
LIBRARIES = -framework veclib -lm
#LIBRARIES = -lfftw3 -lm
#LIBRARIES = -lmpi -lscs_mp

# other variables
MAINFILE = checkcopy
OTHERS   = qobjects lapack_interface
OBJS = $(MAINFILE:%=%.o) $(OTHERS:%=%.o) 

# pattern rules to define implicit rule
%.o: %.cc
	g++ -c $(OFLAGS) $< 
        
%.o: %.cpp
	g++ -c $(OFLAGS) $<
        
%.o: %.f
	f77 -c $<


# default rule(s)
all: $(MAINFILE)
$(MAINFILE): $(OBJS)
	$(LINK.c) $(LFLAGS) -o exec $^ $(LIBRARIES)
        
depend:
	mkdep -f depend $(OFLAGS) $(MAINFILE:%=%.cc) $(OTHERS:%=%.cc)
.PHONY: depend
include depend

