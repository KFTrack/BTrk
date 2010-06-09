# Please see SoftRelTools/HOWTO-GNUmakefile for documentation
# $Id: GNUmakefile.template,v 1.2 2005/08/09 23:32:18 chee Exp $
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

SUBDIRS = src bins
#################################################################
#++ include standard makefile from SoftRelTools.
include SoftRelTools/standard.mk
