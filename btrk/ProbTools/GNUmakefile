# Please see SoftRelTools/HOWTO-GNUmakefile for documentation
# $Id: GNUmakefile 560 2010-03-29 23:44:43Z brownd $
#################################################################
#++ library products				[build it with 'lib']

LIBREMOVEFILES := 
LIBTMPLFILES := 
LIBDDLORDERED := 
LIBEXTRARDLFILES := 

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

SUBDIRS = src
#################################################################
#++ include standard makefile from SoftRelTools.
include SoftRelTools/standard.mk
