#path create sample PdtInit DchMakeDigi DchDataHist DchDumpData
path enable AllPath

global useL1Sim
set useL1Sim true

sourceFoundFile FcsSim/FcsInputSequence.tcl 

module talk FileInput
  input file bbsim.xdr
  exit

module talk HbkTupleEnv
  ntupleFile set true
  histFileName set digi_bdb.hbook
  exit

# disable unwanted bdb modules
#mod disable BdbSetTime
#mod disable BdbCreateCM
#mod disable BdbEventUpdate
#mod disable BdbInspector
mod disable EmcBuildEnv
mod disable L3TConfigProxies
mod disable L3TBuildEnv
#mod disable G3BdbLoad
#mod disable DchBdbLoad
mod disable PepBdbLoad
mod disable DrcBdbLoad
mod disable EmcBdbLoad
mod disable IfrBdbLoad
mod disable L1DctBdbLoad
mod disable L1EmtBdbLoad`
#mod disable L1FctBdbLoad
mod disable L1GltBdbLoad
mod disable L3TBdbLoad
mod disable SvtBdbLoad
mod disable TrkBdbLoad
mod disable BtaBdbLoad
mod disable StdHepBdbLoad
mod disable TagBdbLoad
mod disable DchSegBunchT0
#mod disable DchTrackFinder
#mod disable DcxHitAdder
#mod disable DcxTrackFinder
#mod disable BdbEventInput

mod disable FcsClockCtl
mod disable RandomCtl

module disable DchDumpData
module input FileInput
module disable DummyInput
module disable FileInputRW
module disable BdbEventInput
mod enable FileInput
# by default, do no output
module output BdbEventOutput
module disable FileOutput
module enable BdbEventOutput
module disable FileOutputRW

path create mybdbtest EvtCounter GenBuildEnv HbkTupleEnv \
 PdtInit BdbSetTime BdbCreateCM DchBdbLoad G3BdbLoad \
 L1FctBdbLoad BdbEventUpdate DchBuildEnv DchDataHist
#  DchMakeHits \
#  InitBunchT0 DchTrackFinder DcxHitAdder DcxTrackFinder DchDataHist


path disable AllPath
path delete AllPath
path enable mybdbtest

# path list

# mod talk DchBuildEnv
#   show
#   verbose set t
#   exit

# mod talk DchTrackFinder
#   verbose set true
#   exit

#
#  set reading of needed data
#
mod talk G3BdbLoad
  writeGTracks set true
  writeGVertex set true
  writeGEvents set true
  exit

mod talk L1FctBdbLoad
  writeDigis set true
  exit

# mod talk L1DctBdbLoad
#   writeTsfDigis set true
#   writeBltDigis set true
#   writePtdDigis set true
#   exit

mod talk DchBdbLoad
  writeDigis set true
  writeGHits set true
  writeMCDigis set true
  exit

events begin -nev 3

exit
