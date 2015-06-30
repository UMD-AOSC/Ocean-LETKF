#!/bin/bash
set -e
date=/bin/date #location of GNU date
CWD=`pwd`

INDIR=/lustre/f1/unswept/Steve.Penny/OBS/synthetic/AVISO
OUTDIR=/lustre/f1/unswept/Steve.Penny/OBS/synthetic/letkf_fmt/SFC_gerr_ALT
#TXDIR=/lustre/f1/unswept/Steve.Penny/Synthetic
TXDIR=/autofs/na1_home1/Steve.Penny/letkf/mom4/obs/NCEP_SFC
#INDIR=.
#OUTDIR=OUT

IY=1992 #2002 #1992
IM=09 #08 #09
ID=24 #10 #24
IH=00

EY=2000 #2002 #1992
EM=01
ED=01
EH=00

inc=1
inc_units='days'

mkdir -p $OUTDIR

# MAIN LOOP
bday=`$date -d 1985-01-01 +%s`
while test $IY$IM$ID$IH -le $EY$EM$ED$EH
do
  cd $CWD
  eday=`$date -d $IY-$IM-$ID +%s`
  yday=`echo "1 + ($eday - $bday) / 86400" | bc`
  echo "yday = $yday"  
  TY=`$date -d "$IY-$IM-$ID $inc $inc_units" +%Y`
  TM=`$date -d "$IY-$IM-$ID $inc $inc_units" +%m`
  TD=`$date -d "$IY-$IM-$ID $inc $inc_units" +%d`
  TH=`$date -d "$IY-$IM-$ID $inc $inc_units" +%H`

# DD=`expr $ID`
# MM=`expr $IM`
# DD=`printf "%02d" $DD`
# MM=`printf "%02d" $MM`

  echo "Computing for $IY$IM$ID"
  rm -f aid.dat
  ./cd.x -y $IY -m $IM -d $ID
  aid=`head -n 1 aid.dat`
  echo "Using aid = $aid"
# afile=`grep -l $aid $INDIR/*.txt`
  line=0
  while test $line -ge 0
  do
    line=`expr $line + 1`
    # Read the index file, one line at a time : txj1j2_cTbl
    linestr=`sed -n -e "${line}{p;q}" $TXDIR/txj1j2_cTbl`

    # Parse the string to get the date range, number of records, and data filename
    date0=${linestr:0:10}
    date1=${linestr:12:10}

    Y0=${date0:0:4}
    M0=${date0:5:2}
    D0=${date0:8:2}
    Y1=${date1:0:4}
    M1=${date1:5:2}
    D1=${date1:8:2}

    if [ $Y0$M0$D0 -le ${IY}${IM}${ID} ] && [ ${IY}${IM}${ID} -le $Y1$M1$D1 ]; then
      records=${linestr:22:7}
      infile=${linestr:29}
      records=`echo $records | sed 's/^ *//'`   #(Remove leading blanks)
      infile=`echo $infile | sed 's/^ *//'`
      echo "line = $line"
      echo "linestr = $linestr"
      echo "parsed line ::"
      echo "date0 = $date0"
      echo "date1 = $date1"
      echo "records = $records"
      echo "infile = $infile"
      afile=$infile

      # Check for an overlap between files
      line=`expr $line + 1`
      linestr=`sed -n -e "${line}{p;q}" $TXDIR/txj1j2_cTbl`
      date02=${linestr:0:10}
      date12=${linestr:12:10}
      Y2=${date02:0:4}
      M2=${date02:5:2}
      D2=${date02:8:2}
      Y3=${date12:0:4}
      M3=${date12:5:2}
      D3=${date12:8:2}
      if [ $Y2$M2$D2 -le ${IY}${IM}${ID} ] && [ ${IY}${IM}${ID} -le $Y3$M3$D3 ]; then
        records2=${linestr:22:7}
        infile2=${linestr:29}
        records2=`echo $records2 | sed 's/^ *//'`   #(Remove leading blanks)
        infile2=`echo $infile2 | sed 's/^ *//'`
        afile="$infile $infile2"
      fi
      line=-1
    fi
  done
  echo "Using altimetry files: afile = $afile"

# exit 1

  idx=0
  for file in $afile; do
    echo "For file: $file"
    if [ $idx -eq 0 ]; then
      ./s2l.x -y $IY -m $IM -d $ID -afile $file -aid $aid -indir $INDIR -outdir $OUTDIR > $OUTDIR/$IY$IM$ID.out
    else
      mv $OUTDIR/$IY$IM$ID.dat $OUTDIR/$IY$IM$ID.dat0
      ./s2l.x -y $IY -m $IM -d $ID -afile $file -aid $aid -indir $INDIR -outdir $OUTDIR > $OUTDIR/$IY$IM$ID.out${idx}
      mv $OUTDIR/$IY$IM$ID.dat $OUTDIR/$IY$IM$ID.dat1
      cat $OUTDIR/$IY$IM$ID.dat0 $OUTDIR/$IY$IM$ID.dat1 > $OUTDIR/$IY$IM$ID.dat
      rm -f $OUTDIR/$IY$IM$ID.dat[0-1]
    fi 
    idx=`expr $idx + 1`
  done

# exit 1

  IY=$TY
  IM=$TM
  ID=$TD
  IH=$TH

done
