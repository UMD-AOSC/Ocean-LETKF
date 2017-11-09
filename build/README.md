This is the primary directory for building letkf and the observation operators.

(1) Update the make_letkf.<YOUR_MODEL>.sh and make_obsop.<YOUR_MODEL>.sh with a new experiment name.

(2) Run make_letkf.<YOUR_MODEL>.sh to build letkf. The result will be in the build_letkf/ directory.

(3) Run make_obsop.<YOUR_MODEL>.sh to build the observation operators. The result will be in the build_obsop/ directory.

---

Any changes to the required files used by the model or observation operator will require a change in:
lnkcommon.sh
lnkcommon_obsop.sh 

Any changes to the basic system configuration can be done by updating files in:
../config/machine.sh
../config/<YOUR_MACHINE>.*.sh

---

The build directory currently uses .sh scripts to perform the build, it will be transitioned to a makefile in the future.

