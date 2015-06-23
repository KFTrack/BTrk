# Please see SoftRelTools/HOWTO-GNUmakefile for documentation
# $Id: GNUmakefile 516 2010-01-15 08:22:00Z stroili $
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
BINCCFILES := 

#++ Binary rules		 [in production, build it with 'bin']

#++ shell script products.. 			[build it with 'bin']
#BINSCRIPTS := testscript

#################################################################
#++ regression test scripts			[build it with 'test']

#################################################################
#++ include standard makefile from SoftRelTools.
include SoftRelTools/standard.mk
