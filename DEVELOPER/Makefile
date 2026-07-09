
include ./Makefile.def
include ../Makefile.def

DEVELOPER_LIBS	= ./element/cpp/Macroelement3d/Macroelement3d/CohesiveSurface.o \
	./element/cpp/Macroelement3d/Macroelement3d/DamageShearInterface.o \
	./element/cpp/Macroelement3d/Macroelement3d/GambarottaLagomarsinoModel.o \
	./element/cpp/Macroelement3d/Macroelement3d/GenericDamagePlasticityShear.o \
	./element/cpp/Macroelement3d/Macroelement3d/Macroelement3d.o \
	./element/cpp/Macroelement3d/Macroelement3d/NoTensionSection3d.o \
	./element/cpp/Macroelement3d/Macroelement3d/WrappedMaterial.o \
	./material/cpp/BeamFrictionSupport/BeamFrictionSupport/BeamFrictionSupport.o \
	./material/cpp/CompressionDamage1d/CompressionDamage1d/CompressionDamage1d.o \
	./material/cpp/TensionDamage1d/TensionDamage1d/TensionDamage1d.o \
	./material/cpp/OrthotropicMembraneSection/OrthotropicMembraneSection.o 

all:
	@cd ./element; $(MAKE);
	@cd ./material; $(MAKE);
#	@cd ./recorder; $(MAKE);
#	@cd ./system; $(MAKE);
#	@cd ./integrator; $(MAKE);
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(DEVELOPER_LIBS) 

tidy:
	@rm -f Makefile.bak *~ #*# core

clean:  tidy
	@rm -f *.o

spotless: clean

wipe: spotless
	@cd ./element; $(MAKE) wipe;
	@cd ./material; $(MAKE) wipe;
#	@cd ./recorder; $(MAKE) wipe;
#	@cd ./system; $(MAKE) wipe;
#	@cd ./integrator; $(MAKE) wipe;

