#
# $Id: intersect.tcl,v 1.7 2000/06/28 15:46:36 stroili Exp $
#
#----------------------------------------------------------

#path enable AllPath

global useL1Sim
set useL1Sim true

# path create sample GenBuildEnv BdbSetTime RandomCTL FCClockCTL PdtInit  \
#  DchBuildEnv L1FctMakeDigis InitBunchT0 DchMakeDigi DchMakeHits \
#  DchTrackFinder GenDchGeom

path enable AllPath

module talk FileInput
  input file bbsim.xdr
  exit

sourceFoundFile FcsSim/FcsInputSequence.tcl 

# module talk  L1FctMakeDigis
#   SimuLevelL1A set 2
#   DistWidthL1A set 500
#   exit

# module talk DchInputFilter
#   filterParm set true
#   l1aParm set $useL1Sim
#   windowVal set 2.2e-6
#  exit

# module talk InitBunchT0
#   useSimL1A set $useL1Sim
#   exit

# module talk DchMakeDigi
#   l1aParm set $useL1Sim
#   smearParm set true
#   methodParm set 0
#   # efficiency cuts
#   efficParm set false
#   efficVal set 1.0
#   # crosstalk control
#   xtalkParm set false
#   xtalkVal set 0.01
#   show
#   exit

mod disable BdbCreateCM
mod disable BdbEventUpdate
mod disable BdbInspector
# mod disable EmcBuildEnv
# mod disable L3TConfigProxies
# mod disable L3TBuildEnv
mod disable G3BdbLoad
mod disable DchBdbLoad
mod disable PepBdbLoad
mod disable DrcBdbLoad
mod disable EmcBdbLoad
mod disable IfrBdbLoad
mod disable L1DctBdbLoad
mod disable L1EmtBdbLoad
mod disable L1FctBdbLoad
mod disable L1GltBdbLoad
mod disable L3TBdbLoad
mod disable SvtBdbLoad
mod disable TrkBdbLoad
mod disable BtaBdbLoad
mod disable StdHepBdbLoad
mod disable TagBdbLoad
mod disable BdbEventInput
mod disable DchSegBunchT0
mod disable BunchT0BdbLoad
mod disable EidBdbLoad
mod disable OepBdbLoad
mod disable RdsBdbLoad
mod disable NeutralHadBdbLoad
#
#  use FileInput as input module
#
module input FileInput

path create mytest FcsClockCtl RandomCtl GenBuildEnv HbkTupleEnv \
 PdtInit EvtCounter L1FctMakeDigis InitBunchT0 DchBuildEnv \
 DchTestGeom

path delete AllPath

module talk GenDchGeom
  verbose set false
  exit

module talk DchBuildEnv
  verbose set false
  exit

#action enable ErrLogger
#action off TimeAction
#action off NameAction

mod disable DchCondMonitor
mod disable GenDchGeom
mod disable DchQC
mod list

mod talk DchTestGeom
  IntersectDch set true
  ValidateNextIntersection set false
  NCosThetaScanPoints set 2000
  NPhiScanPoints set 20
  TrajectoryLength set 500
  CosThetaLowRange set -.98
  CosThetaHighRange set .98
  verbose set false
  exit

ev begin -nev 1

exit
