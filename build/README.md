This is the primary directory for building letkf and the observation operators.

Update the make_letkf.<YOUR_MODEL>.sh and make_obsop.<YOUR_MODEL>.sh with a new experiment name.

Run make_letkf.<YOUR_MODEL>.sh to build letkf. The result will be in the build_letkf/ directory.

Run make_obsop.<YOUR_MODEL>.sh to build the observation operators. The result will be in the build_obsop/ directory.


---

The build directory currently uses .sh scripts to perform the build, it will be transitioned to a makefile in the future.

