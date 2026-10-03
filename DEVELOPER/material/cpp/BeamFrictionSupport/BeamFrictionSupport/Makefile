include ../../../../Makefile.def

ifeq ($(OS_FLAG), -D_LINUX)

LIBRARY 	= BeamFrictionSupport.so

OBJS	= BeamFrictionSupport.o 

all:   $(OBJS)
	$(CC++) -shared -Wl,-soname,$(LIBRARY) -o $(LIBRARY) $(OBJS) -lc

BeamFrictionSupport.o: BeamFrictionSupport.cpp
	$(CC++) -fPIC $(INCLUDES) -g -c -Wall BeamFrictionSupport.cpp



else

LIBRARY 	= BeamFrictionSupport.dylib


all: $(OBJS)
	$(CC++) $(OS_FLAG) -dynamiclib $(INCLUDES) -Wl,-undefined,suppress,-flat_namespace BeamFrictionSupport.cpp $(OUTSIDE_OBJS)  -current_version 1.0 -compatibility_version 1.0 -fvisibility=hidden -o $(LIBRARY)


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
