# Please see SoftRelTools/HOWTO-GNUmakefile for documentation
# $Id: GNUmakefile 107 2010-01-21 09:12:24Z stroili $
#################################################################
#++ library products				[build it with 'lib']

LIBREMOVEFILES := 
LIBTMPLFILES := 
LIBDDLORDERED := 

#################################################################
#++ extra binary products	[not in production, build it with extrabin]

EXTRABINS := 

#################################################################
#++ binary products				[build it with 'bin']

BINS := 
BINCCFILES := 

#++ Binary rules		 [in production, build it with 'bin']
$(addprefix $(bindir),$(BINS)): $(bindir)% : %.o

#++ shell script products.. 			[build it with 'bin']
BINSCRIPTS := 

#################################################################
#++ regression test scripts			[build it with 'test']

#################################################################
#++ include standard makefile from SoftRelTools.
include SoftRelTools/standard.mk
