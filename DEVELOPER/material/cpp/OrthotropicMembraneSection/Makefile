include ../../../Makefile.def

ifeq ($(OS_FLAG), -D_LINUX)

LIBRARY 	= OrthotropicMembraneSection.so

OBJS	= OrthotropicMembraneSection.o 

all:   $(OBJS)
	$(CC++) -shared -Wl,-soname,$(LIBRARY) -o $(LIBRARY) $(OBJS) -lc

OrthotropicMembraneSection.o: OrthotropicMembraneSection.cpp
	$(CC++) -fPIC $(INCLUDES) -g -c -Wall OrthotropicMembraneSection.cpp



else

LIBRARY 	= TensionDamage1d.dylib


all: $(OBJS)
	$(CC++) $(OS_FLAG) -dynamiclib $(INCLUDES) -Wl,-undefined,suppress,-flat_namespace OrthotropicMembraneSection.cpp $(OUTSIDE_OBJS)  -current_version 1.0 -compatibility_version 1.0 -fvisibility=hidden -o $(LIBRARY)


endif

# Miscellaneous
tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean:  tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o core *.out *.so *.dylib

spotless: clean
	@$(RM) $(RMFLAGS) $(PROGRAM) fake core

wipe: spotless
	@$(RM) $(RMFLAGS) $(PROGRAM) fake core $(LIBRARY)

# DO NOT DELETE THIS LINE -- make depend depends on it.
