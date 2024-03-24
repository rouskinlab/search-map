# Remaining Tasks for SEARCH-MaP / SEISMIC-RNA Paper

## SEISMIC-RNA

- End-aware bias correction
  - Make tables work with new bias correction.
  - Rerun clustering on the 3kb no-ASO ctrl-2 sample to ensure the results are similar.
  - Test clustering algorithm with unit tests.
  - Release version 0.15.0


## Experiments

- rRNA controls 
  - Design experiment.
  - Replicate 1
    - Probe structure with and without ASOs.
    - Sequence.
  - Replicate 2
    - Probe structure with and without ASOs.
    - Sequence.


## Analyses

- Analyze rRNA controls.
- Re-analyze SARS-CoV-2
  - Re-cluster the 3kb no-ASO ctrl2 for FSE and LRI.
  - Re-fold the 1.8kb based on those clusters (and hope for identical results).
- Re-analyze TGEV
  - Re-cluster the FSE and LRI (hope for similar results).
  - Recompute the ensemble average including the bias correction.
  - 