# ocean-obs

This directory contains software for pre-processing (including download) ocean observations 
and datasets. It is organized by different levels of processing:
- Level-1 (in-situ; _spectral_ for satellite: brightness temperature): _raw_ measurement.
- Level-2 (retrievals; not _actual_ measurements).
- Level-4 (gap free gridded datasets, such as optimally interpolated ones).

### Note:
- One should not confuse with other _meanings_ of `levels`, such as this one at the [OceanPredict OS-Eval TT](https://oceanpredict.org/observations-use/#section-argo-profiling-floats)
  that is meant for _how_ the data is being used in specific ocean assimilation systems; not how it is derived in the first place.
- It is advised that every contribution to this repository is accompanied by a [jupyter notebook](https://jupyter.org/). Please `clear` all output so that the notebook is of minimal size.
- Remember to append `utils/pycommon` to the env variable (e.g., `export PYTHONPATH=../../pycommon`) before running the python3 scripts under this directory

