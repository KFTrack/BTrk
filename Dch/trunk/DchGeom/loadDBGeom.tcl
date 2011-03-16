#
# $Id: loadDBGeom.tcl,v 1.1 1998/03/20 01:27:37 stroili Exp $
#
#----------------------------------------------------------

path enable AllPath

module talk DchBuildEnv
  usedb set true
  exit

module talk FileInput
  input file bbsim.xdr
  exit

#module talk GenDchGeom
#  verbose set false
#  exit

ev begin -nev 3

exit
