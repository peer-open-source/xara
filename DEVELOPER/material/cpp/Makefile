include ../../Makefile.def

all:
	echo @pwd;
	@$(CD) CompressionDamage1d/CompressionDamage1d; $(MAKE);
	@$(CD) TensionDamage1d/TensionDamage1d; $(MAKE);
	@$(CD) OrthotropicMembraneSection; $(MAKE);
	@$(CD) BeamFrictionSupport/BeamFrictionSupport; $(MAKE);

# Miscellaneous
tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean:  tidy
	@$(CD) ./CompressionDamage1d/CompressionDamage1d; $(MAKE) wipe;
	@$(CD) ./TensionDamage1d/TensionDamage1d; $(MAKE) wipe;
	@$(CD) ./OrthotropicMembraneSection; $(MAKE) wipe;
	@$(CD) ./BeamFrictionSupport/BeamFrictionSupport; $(MAKE) wipe;
	@$(RM) $(RMFLAGS) $(OBJS) *.o core *.out *.so *.dylib

spotless: clean

wipe: spotless


# DO NOT DELETE THIS LINE -- make depend depends on it.
