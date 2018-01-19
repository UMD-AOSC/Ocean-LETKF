# Model-specific source files:

The model-specific routines are organized by:
(a) using the model type in the filename (i.e. mom4, mom6, hycom, roms)
(b) storing the model-specific versions of each file in its own directory
(c) using the generic term 'oceanmodel' anywhere in the code that would reference the specific model
(d) ignoring the convention to make the filename and module name the same

Upon compilation, the model can be specified, and assuming appropriate encapsulation,
the rest of the code will be able to compile with any of the supported models.

In progress: nemo, roms
Upcoming model support: wavewatch3

~Steve
