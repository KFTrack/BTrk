#path create sample PdtInit DchMakeDigi DchDataHist DchDumpData
# define sourceFoundFile proc to locate .tcl files


# run the script to configure reconstruction

# Turn off creation of StdHep data as this is in database
echo Disabling StdHepLoad and MakeXReference
catch { module disable StdHepLoad }


global CollectionName
global HFile

set CollectionName events
set HFile dchgeom.hbook

if [ info exists env(DchCollectionName) ] {
  set CollectionName $env(DchCollectionName)
  set HFile dch_$env(DchCollectionName).hbook
}

if [ info exists env(DchHistoFile) ] {
  set HFile $env(DchHistoFile)
}

module talk HbkTupleEnv
  ntupleFile set true
  histFileName set $HFile 
  exit

# get input from the database
sourceFoundFile BdbSequences/BdbSequence.tcl
sourceFoundFile BdbSequences/BdbReadSim.tcl
sourceFoundFile BdbSequences/BdbReadRaw.tcl
sourceFoundFile BdbSequences/BdbReadRec.tcl
sourceFoundFile BdbSequences/BdbReadAod.tcl
sourceFoundFile BdbSequences/BdbReadEsd.tcl
sourceFoundFile BdbSequences/BdbReadTag.tcl

sourceFoundFile  BetaSequences/bdbRecoSetup.tcl 

# path enable AllPath

global useL1Sim
set useL1Sim true

sourceFoundFile FcsSim/FcsInputSequence.tcl 

module talk FileInput
  input file bbsim.xdr
  exit

# disable unwanted bdb modules
#mod disable BdbSetTime
#mod disable BdbCreateCM
#mod disable BdbEventUpdate
#mod disable BdbInspector
#mod disable EmcBuildEnv
#mod disable L3TConfigProxies
#mod disable L3TBuildEnv
#mod disable G3BdbLoad
#mod disable DchBdbLoad
# mod disable PepBdbLoad
# mod disable DrcBdbLoad
mod disable EmcBdbLoad
# mod disable IfrBdbLoad
#mod disable L1DctBdbLoad
#mod disable L1EmtBdbLoad
#mod disable L1FctBdbLoad
# mod disable L1GltBdbLoad
# mod disable L3TBdbLoad
# mod disable SvtBdbLoad
# mod disable TrkBdbLoad
# mod disable BtaBdbLoad
# mod disable StdHepBdbLoad
# mod disable TagBdbLoad
# mod disable DchSegBunchT0
#mod disable DchTrackFinder
#mod disable DcxHitAdder
#mod disable DcxTrackFinder
#mod disable BdbEventInput

# mod disable FcsClockCtl
# mod disable RandomCtl

module disable DchDumpData
module input BdbEventInput
module disable FileInput
module disable DummyInput
# module disable FileInputRW
module enable BdbEventInput
mod disable FileInput
# by default, do no output
module output FileOutput
module disable FileOutput
module disable BdbEventOutput
#module disable FileOutputRW

seq list

path create mybdbtest BdbSequence DchBuildEnv DchTestGeom DchCondMonitor \
  GenDchGeom

seq list

path list
# EvtCounter GenBuildEnv HbkTupleEnv \
#  PdtInit BdbCreateCM DchBdbLoad G3BdbLoad L1DctBdbLoad \
#  SvtBdbLoad TrkBdbLoad TagBdbLoad \
# L1FctBdbLoad BdbEventUpdate DchBuildEnv DchTestGeom DchCondMonitor \
# GenDchGeom
#  DchMakeHits \
#  InitBunchT0 DchTrackFinder DcxHitAdder DcxTrackFinder DchDataHist


# path disable AllPath
# path delete AllPath
# path enable mybdbtest
path disable mybdbtest

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
# mod talk G3BdbLoad
#   readGTracks set true
#   readGVertex set true
#   readGEvents set true
#   exit

# mod talk L1FctBdbLoad
#   readDigis set true
#   exit

# mod talk L1DctBdbLoad
#   readTsfDigis set true
#   readBltDigis set true
#   readPtdDigis set true
#   exit

# mod talk DchBdbLoad
#   readDigis set true
#   readGHits set true
#   readMCDigis set true
#   exit

mod talk BdbEventInput
  collectionName set $CollectionName 
  exit

module talk DchCondMonitor
  verbose set true
  exit

mod talk GenDchGeom
  verbose set true
  trkRecoTrkList set Default
  exit

# action enable ErrLogger
action off TimeAction
action off NameAction


events begin 

exit











