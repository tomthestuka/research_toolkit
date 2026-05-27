#!/usr/bin/env bash

# Use by $nohup ./nicer_red.sh > nicer_red.out 2>&1 &

# ---------- Runtime environment ----------- #
# HEASoft
export HEADAS=/opt/heasoft/heasoft-6.34/x86_64-pc-linux-gnu-libc2.34
. "$HEADAS/headas-init.sh"
export HEADASNOQUERY=YES

# Conda env
source /opt/Anaconda/anaconda3/bin/activate
conda activate /opt/anaconda/anaconda3/envs/astropy

# CALDB
export CALDB=/man1/CALDB
export CALDBCONFIG="$CALDB/software/tools/caldb.config"
export CALDBALIAS="$CALDB/software/tools/alias_config.fits"

# Geomagnetic path
export GEOMAG_PATH=/man1/CALDB/data/gen/pcf/geomag

# ------------- Logging setup -------------- #
#LOG_L2="nicer_l2.log"
LOG_L3="nicer_l3.log"
#: > "$LOG_L2"   # truncate existing
: > "$LOG_L3"

# ------------- Main loop ------------------ #
for obsid in [0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]; do
  [[ -d "$obsid" ]] || continue

  # Minimal progress messages to STDOUT
  #echo "processing nicerl2 for $obsid"
  #nicerl2 "$obsid" clobber=YES   >> "$LOG_L2" 2>&1

  echo "processing nicerl3-spect for $obsid"
  nicerl3-spect "$obsid" bkgmodeltype=3c50 bkgformat=file clobber=YES >> "$LOG_L3" 2>&1

  # Spacer lines in logs only
  echo "" >> "$LOG_L3"
  echo ""
done

echo "Done — see $LOG_L2 and $LOG_L3 for details."
