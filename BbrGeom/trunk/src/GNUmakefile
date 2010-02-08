# Please see SoftRelTools/HOWTO-GNUmakefile for documentation
# $Id: GNUmakefile 491 2010-01-13 16:59:16Z stroili $
#################################################################
#++ library products				[build it with 'lib']

LIBREMOVEFILES := 
LIBTMPLFILES := 
LIBDDLORDERED := 
LIBEXTRARDLFILES := 

#################################################################
#++ extra binary products	[not in production, build it with extrabin]

EXTRABINS := 

#################################################################
#++ binary products				[build it with 'bin']

BINS := 
BINCCFILES := $(BINS:=.cc)

#++ Binary rules		 [in production, build it with 'bin']

$(addprefix $(bindir),$(BINS)): $(bindir)% : %.o

#++ shell script products.. 			[build it with 'bin']
BINSCRIPTS := 

#################################################################
#++ regression test scripts			[build it with 'test']

#################################################################
#++ include standard makefile from SoftRelTools.
include SoftRelTools/standard.mk
