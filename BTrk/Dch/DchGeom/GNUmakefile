# Please see SoftRelTools/HOWTO-GNUmakefile for documentation
# $Id: GNUmakefile 118 2010-03-25 14:19:45Z stroili $
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

$(bindir)$(BINS): $(BINCCFILES:.cc=.o)

#++ shell script products.. 			[build it with 'bin']
BINSCRIPTS := 

#################################################################
#++ regression test scripts			[build it with 'test']

SUBDIRS = src
#################################################################
#++ include standard makefile from SoftRelTools.
include SoftRelTools/standard.mk
