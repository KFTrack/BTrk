/*
** ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
**  Package:
**	DCH Dataflow
**
**  Abstract:
**      Defines the status flags used during feature extraction to identify
**      suspicious features of channel data.
**
**  Author:
**      Michael Kelsey, Princeton University, kelsey@slac.stanford.edu
**
**  Creation Date:
**	000 - August 28, 1998
**
**  Revision History:
**	980828  Extracted from DchChanFeature.hh
**      990419  Make dtor public to suppress ODF compiler warnings
**	020116  Add TDC and FADC conversion factors here, from DchDigi
**	021011  Define type-name for flag mask enumerators
** ----------------------------------------------------------------------------
*/
#ifndef DCH_FEATURE_FLAGS_HH
#define DCH_FEATURE_FLAGS_HH

class DchFeatureFlags {
public:
  enum Mask { badWaveform = 0x200,
	      badPedestal = 0x400,
	      clippedWaveform = 0x800,
	      rawDataAppended = 0x1000,
              unstablePedestal = 0x2000,
              chargeCorrected = 0x4000,
              allFlags = 0xFE00 };

  enum { flagShift = 9 };

  ~DchFeatureFlags();		// Public to allow usage of class namespace

  // Conversion factors between electronics counts and real units
  inline static const double TDCtoTime() { return 1000./59.5/256.; }
  inline static const double FADCtoChg() { return 1.; }

private:
  // Make constructor private to prevent instantiation
  DchFeatureFlags();
};

#endif  /* DCH_FEATURE_FLAGS_HH */
