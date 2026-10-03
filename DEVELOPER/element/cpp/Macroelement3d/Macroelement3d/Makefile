include ../../../../Makefile.def

ifeq ($(OS_FLAG), -D_LINUX)

LIBRARY 	= Macroelement3d.so

OBJS       = Macroelement3d.o \
	CohesiveSurface.o \
	DamageShearInterface.o \
	GambarottaLagomarsinoModel.o \
	GenericDamagePlasticityShear.o \
	NoTensionSection3d.o \
	WrappedMaterial.o 

all:         $(OBJS)
	$(CC++) -shared -Wl,-soname,$(LIBRARY) -o $(LIBRARY) $(OBJS) -lc

Macroelement3d.o: Macroelement3d.cpp
	$(CC++) -fPIC $(INCLUDES) -g -c -Wall Macroelement3d.cpp

DamageShearInterface.o: DamageShearInterface.cpp
	$(CC++) -fPIC $(INCLUDES) -g -c -Wall DamageShearInterface.cpp	

CohesiveSurface.o: CohesiveSurface.cpp
	$(CC++) -fPIC $(INCLUDES) -g -c -Wall CohesiveSurface.cpp	
	
GambarottaLagomarsinoModel.o: GambarottaLagomarsinoModel.cpp
	$(CC++) -fPIC $(INCLUDES) -g -c -Wall GambarottaLagomarsinoModel.cpp	
	
GenericDamagePlasticityShear.o: GenericDamagePlasticityShear.cpp
	$(CC++) -fPIC $(INCLUDES) -g -c -Wall GenericDamagePlasticityShear.cpp	
	
NoTensionSection3d.o: NoTensionSection3d.cpp
	$(CC++) -fPIC $(INCLUDES) -g -c -Wall NoTensionSection3d.cpp	
	
WrappedMaterial.o: WrappedMaterial.cpp
	$(CC++) -fPIC $(INCLUDES) -g -c -Wall WrappedMaterial.cpp					



else

LIBRARY 	= Macroelement3d.dylib


all: $(OBJS)
	$(CC++) $(OS_FLAG) -dynamiclib $(INCLUDES) -Wl,-undefined,suppress,-flat_namespace Macroelement3d.cpp $(OUTSIDE_OBJS)  -current_version 1.0 -compatibility_version 1.0 -fvisibility=hidden -o $(LIBRARY)


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
