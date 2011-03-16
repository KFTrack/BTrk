#
# $Id: testBdbDchData.tcl,v 1.8 2000/09/07 13:19:27 stroili Exp $
#
#   tcl script to read Dch digi's from the database
#
#
#------------------------------------------------------------------------
#

sourceFoundFile FcsSim/FcsInputSequence.tcl

global CollectionName
global HFile

set CollectionName events
set HFile digi_bdb.hbook

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

# disable unwanted bdb modules
mod disable PepBdbLoad
mod disable DrcBdbLoad
mod disable EmcBdbLoad
mod disable IfrBdbLoad
mod disable L1GltBdbLoad
mod disable L3TBdbLoad
mod disable SvtBdbLoad
mod disable TrkBdbLoad
mod disable BtaBdbLoad
mod disable StdHepBdbLoad
mod disable TagBdbLoad
mod disable DchSegBunchT0
mod disable FcsClockCtl
mod disable RandomCtl
module disable DchDumpData
module disable FileInput
module disable DummyInput
mod disable FileInput

# input module
module input BdbEventInput
# module disable FileInputRW
module enable BdbEventInput

# by default, do no output
module output FileOutput
module disable FileOutput
module disable BdbEventOutput

path create mybdbtest EvtCounter GenBuildEnv MatBuildEnv HbkTupleEnv \
 PdtInit BdbCreateCM DchBdbLoad G3BdbLoad L1DctBdbLoad \
 L1FctBdbLoad BdbEventUpdate DchBuildEnv DchDataHist

path delete AllPath
path enable mybdbtest


set num ""

sourceFoundFile DchDataP/DchBdbReadRaw.tcl
#
# read MC info
#
sourceFoundFile DchDataP/DchBdbReadSim.tcl
sourceFoundFile G3DataP/G3BdbReadSim.tcl

mod talk BdbEventInput
  collectionName set $CollectionName 
  exit


events begin -nev 3

exit
