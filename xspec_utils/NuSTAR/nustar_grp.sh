#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# nustar_grp.sh — NuSTAR Level 3 grouping script
# -----------------------------------------------------------------------------
# Groups spectra for all OBSIDs using ftgrouppha.
#
# Usage:
#   nohup ./nustar_grp.sh > nustar_grp.out 2>&1 &
# -----------------------------------------------------------------------------

for obs in [0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]; do
  echo "=== Grouping spectra for OBSID $obs ==="
  
  # Define input and output files for FPMA
  infile_a="/data/fzhou/XSPEC/Cyg_X-1/NuSTAR/$obs/products/nu${obs}A01_sr.pha"
  outfile_a="/data/fzhou/XSPEC/Cyg_X-1/NuSTAR/$obs/products/nu${obs}A01_sr_grp.pha"
  rmf_a="/data/fzhou/XSPEC/Cyg_X-1/NuSTAR/$obs/products/nu${obs}A01_sr.rmf"
  
  # Define input and output files for FPMB
  infile_b="/data/fzhou/XSPEC/Cyg_X-1/NuSTAR/$obs/products/nu${obs}B01_sr.pha"
  outfile_b="/data/fzhou/XSPEC/Cyg_X-1/NuSTAR/$obs/products/nu${obs}B01_sr_grp.pha"
  rmf_b="/data/fzhou/XSPEC/Cyg_X-1/NuSTAR/$obs/products/nu${obs}B01_sr.rmf"
  
  # Run ftgrouppha for FPMA in the background
  echo "Grouping FPMA spectrum for OBSID $obs"
  ftgrouppha \
      infile="$infile_a" \
      outfile="$outfile_a" \
      grouptype=opt \
      respfile="$rmf_a" \
      clobber=yes \
      > "/data/fzhou/XSPEC/Cyg_X-1/NuSTAR/$obs/nu${obs}_grp.log" 2>&1 &
  
  # Run ftgrouppha for FPMB in the background
  echo "Grouping FPMB spectrum for OBSID $obs"
  ftgrouppha \
      infile="$infile_b" \
      outfile="$outfile_b" \
      grouptype=opt \
      respfile="$rmf_b" \
      clobber=yes \
      > "/data/fzhou/XSPEC/Cyg_X-1/NuSTAR/$obs/nu${obs}_grp.log" 2>&1 &
done