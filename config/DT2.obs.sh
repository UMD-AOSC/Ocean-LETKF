#!/bin/bash

#STEVE: needs to be updated for DT2!!!

# Individual observation directories in original format:
OBSDIR_TPROF=~/lf1u/OBS/historical/TMP_profs
OBSDIR_SPROF=~/lf1u/OBS/historical/SAL_profs_O #observed only
OBSDIR_SPRF2=~/lf1u/OBS/historical/SAL_profs_M #observed and synthetic, for potential temperature computation
OBSDIR_ALTIM=/lustre/f1/unswept/Steve.Penny/OBS/historical/AVISO   #AVISO altimetry data
OBSDIR_TXDIR=/autofs/na1_home1/Steve.Penny/letkf/mom4/obs/NCEP_SFC #datetime table

# Individual observation directories in godas input format:
OBSDIR_TMPA=~/lf1u/OBS/historical/TMP_profs_tmpa_deep
OBSDIR_SALA=~/lf1u/OBS/historical/SAL_profs_O_sala_deep
OBSDIR_SAL2=~/lf1u/OBS/historical/SAL_profs_M_sala_deep
OBSDIR_ALTA=

# godas obs processing scripts
PTMPSAL_SRC=$ROOT/mom4/support/PTMP_SAL/SRC
oaTexe=cmbDLstTh4nc_obsa
oaSexe=cmbDLstPs4nc_obsa

# Combined observation directory in letkf format:
OBSDIR_LETKF_TPROF=~/lf1u/OBS/historical/letkf_fmt/PROFS_gerr_Tonly_deep
OBSDIR_LETKF_SPROF=
OBSDIR_LETKF_TSPRF=~/lf1u/OBS/historical/letkf_fmt/PROFS_gerr_TS_deep
OBSDIR_LETKF_ALTIM=/lustre/f1/unswept/Steve.Penny/OBS/historical/letkf_fmt/SFC_gerr_ALT
