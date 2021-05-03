# transit_tools

A set of Python tools to aid the analysis of exoplanet transits in photometric data. This package is functional but still heavily under construction and is being modified to make it more flexible for a number of exoplanet transit applications. Please contact Ben Hord at benhord@astro.umd.edu about this repo.

Current supported functionality (non-exhaustive):
- fetch TESS light curves via MAST and eleanor
- import of custom light curves
- search for known planets given a target name/RA,Dec and fetch parameters
- fetch stellar parameters
- light curve processing powered via lightkurve
- transit search using Transit Least Squares
- output vetting sheet to show transit diagnostic plots and values
- vetting of signals with TRICERATOPS
- basic vespa functionality

Upcoming:
- TTV search capability
- BLS search capability
- DAVE (implemented but not currently working)
- exoplanet capability
- REBOUND capability
- easy install

All software incoporated here is open source and credit goes to those who developed these useful tools.