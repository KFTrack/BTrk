#----------------------------------------------------------------------------
#                                                                        
#    tcl script used to run DchDataApp application on .xdr files
#
# $Id: testDchData.tcl,v 1.14 2006/08/01 16:25:40 kelsey Exp $
#
#----------------------------------------------------------------------------

#path create sample PdtInit DchMakeDigi DchDataHist DchDumpData
path enable AllPath

global useL1Sim
set useL1Sim true

sourceFoundFile FcsSim/FcsInputSequence.tcl

global XDRfile
global HFile

set XDRfile bbsim.xdr
set HFile dch.hbook

if [ info exists env(DchXDRfile) ] {
  set FileName $env(DchXDRfile)
  set HFile dch_$env(DchXDRfile).hbook
}

if [ info exists env(DchHistoFile) ] {
  set HFile $env(DchHistoFile)
}
  
#
#  xdr input file
#
module talk FileInput
  input file $FileName
  exit

module talk HbkTupleEnv
  ntupleFile set true
  histFileName set $HFile
  exit
#
#   disable unwanted bdb modules
#
# mod disable BdbSetTime
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

module disable DchDumpData
#
#  use FileInput as input module
#
module input FileInput

#
#  save digi's to database
#
if [ info exist env(WriteDchDigi) ] {
  echo WriteDchDigi set, saving digi's to database
  module output BdbEventOutput
  module disable FileOutput
  module disable FileOutputRW

  mod talk BdbEventOutput
    collectionName set test
    statistics set true
    exit

  mod talk G3BdbLoad
    writeGTracks set true
    writeGVertex set true
    writeGEvents set true
    exit

  mod talk L1FctBdbLoad
    writeDigis set true
    exit

  mod talk L1DctBdbLoad
    writeTsfDigis set true
    writeBltDigis set true
    writePtdDigis set true
    exit

  mod talk DchBdbLoad
    writeDigis set true
    writeGHits set true
    writeMCDigis set true
    exit

} 

#
#  create execution path
#
path create mytest FcsClockCtl RandomCtl GenBuildEnv HbkTupleEnv \
 PdtInit EvtCounter L1FctDigiMaker InitBunchT0 DchBuildEnv \
 DchMakeDigi DchDataHist
#  DchMakeDigi DchMakeHits InitBunchT0 \
#  DchTrackFinder DcxHitAdder DcxTrackFinder DchDataHist

#
#  disable default path
#
path delete AllPath

mod talk DchBuildEnv
  verbose set t
  exit

mod disable DchMakeDigi

# mod talk DchTrackFinder
#   verbose set true
#   exit

events begin -nev 1000

exit
