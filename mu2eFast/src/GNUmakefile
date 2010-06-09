# Please see SoftRelTools/HOWTO-GNUmakefile for documentation
# $Id: GNUmakefile 1880 2010-01-25 09:47:55Z stroili $
#################################################################
#++ library products				[build it with 'lib']

LIBREMOVEFILES :=
LIBTMPLFILES :=
LIBDDLORDERED :=

#################################################################
#++ extra binary products	[not in production, build it with extrabin]

EXTRABINS :=

$(addprefix $(bindir),$(EXTRABINS)): $(bindir)% : %.o

#################################################################
#++ binary products				[build it with 'bin']

BINS := 

BINCCFILES := $(BINS:=.cc) $(EXTRABINS:=.cc)

#++ Binary rules		 [in production, build it with 'bin']

$(addprefix $(bindir),$(BINS)): $(bindir)% : %.o

#++ shell script products.. 			[build it with 'bin']
BINSCRIPTS :=

#################################################################
#++ regression test scripts			[build it with 'test']

# $(testdir)mytest.T : mytest.tcl mytesttemp2

#################################################################
#++ include standard makefile from SoftRelTools.
include SoftRelTools/standard.mk
