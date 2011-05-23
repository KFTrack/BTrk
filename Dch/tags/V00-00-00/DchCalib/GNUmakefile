# Please see SoftRelTools/HOWTO-GNUmakefile for documentation
# $Id: GNUmakefile 88 2010-01-14 12:32:57Z stroili $
#################################################################
#++ library products				[build it with 'lib']

LIBREMOVEFILES := 
LIBTMPLFILES := DchLinFun.cc DchCheby.cc DchPoly.cc DchCheby47.cc
LIBDDLORDERED := 

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

SUBDIRS = src
#################################################################
#++ include standard makefile from SoftRelTools.
include SoftRelTools/standard.mk
