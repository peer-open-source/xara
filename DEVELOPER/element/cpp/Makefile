include ../../Makefile.def

all:
	echo @pwd;
	@$(CD) Macroelement3d/Macroelement3d; $(MAKE);

# Miscellaneous
tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean:  tidy
	@$(CD) ./Macroelement3d/Macroelement3d; $(MAKE) wipe;
	@$(RM) $(RMFLAGS) $(OBJS) *.o core *.out *.so *.dylib

spotless: clean

wipe: spotless


# DO NOT DELETE THIS LINE -- make depend depends on it.
