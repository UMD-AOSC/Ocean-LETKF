#!/bin/ksh

DTGSTART=2013010100
DTGEND=2013010200

export PYTHONPATH=/u/prob/interface/rowley/nflux/python

export NFLUXVAR=/u/prob/boundary/jmay/NFLUX_realtime/NFLUX_VAR
export REGION=global_ng_bckg

alias DTG=/home/rowley/linux64/bin/dtg

export CRDATE=$DTGSTART
CRDATE=`DTG -h -6 2>/dev/null`

while [[ $CRDATE -lt $DTGEND ]]; do
  CRDATE=`DTG -h 6 2>/dev/null`
  print Run make_var_plots for $CRDATE ...
  ./make_var_plots.py
  print \ \ \ \ ... done
done
  
