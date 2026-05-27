#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# nustar_red.sh — NuSTAR L2/L3 batch driver
# -----------------------------------------------------------------------------
# Processes OBSIDs sequentially; FPMA and FPMB run in parallel within each OBSID.
#
# Usage:
#   nohup ./nustar_red.sh > nustar_red.out 2>&1 &
# -----------------------------------------------------------------------------

# Define function to process one OBSID
process_obsid() {
  obsid=$1
  
  echo "=== Processing OBSID $obsid ==="
  
  # fresh logs with pattern nu<obsid>_lx.log
  L3A_LOG="$obsid/nu${obsid}_l3A.log"
  L3B_LOG="$obsid/nu${obsid}_l3B.log"
  : > "$L3A_LOG"
  : > "$L3B_LOG"

  # --- gunzip ----------------------------------------------------------------
  gunzip -r "$obsid/event_cl/" 2>/dev/null || true
  gunzip -r "$obsid/auxil/" 2>/dev/null || true
  gunzip -r "$obsid/hk/" 2>/dev/null || true

  # --- Generate src.reg and bkg.reg using find_src.py -----------------------
  echo "Generating source and background regions for OBSID $obsid"
  evt_file="$obsid/event_cl/nu${obsid}A01_cl.evt"
  python3 /data/fzhou/XSPEC/Cyg_X-1/lib/find_src.py "$evt_file"

  # Check if regions were created
  if [ -f "$obsid/event_cl/src.reg" ] && [ -f "$obsid/event_cl/bkg.reg" ]; then
    echo "✓ Regions created successfully"
    echo "  - src.reg: $(cat $obsid/event_cl/src.reg | tail -1)"
    echo "  - bkg.reg: $(cat $obsid/event_cl/bkg.reg | tail -1)"
  else
    echo "✗ Failed to create region files"
    exit 1
  fi

  # --- Level 3 (nuproducts + ftgrouppha) -------------------------------------
  echo "Running nuproducts FPMA" | tee -a "$L3A_LOG"
  nuproducts \
      indir="/data/fzhou/XSPEC/Cyg_X-1/NuSTAR/$obsid/event_cl" \
      instrument=FPMA \
      steminputs="nu$obsid" \
      outdir="/data/fzhou/XSPEC/Cyg_X-1/NuSTAR/$obsid/products" \
      srcregionfile="/data/fzhou/XSPEC/Cyg_X-1/NuSTAR/$obsid/event_cl/src.reg" \
      bkgregionfile="/data/fzhou/XSPEC/Cyg_X-1/NuSTAR/$obsid/event_cl/bkg.reg" \
      attfile="/data/fzhou/XSPEC/Cyg_X-1/NuSTAR/$obsid/auxil/nu${obsid}_att.fits" \
      hkfile="/data/fzhou/XSPEC/Cyg_X-1/NuSTAR/$obsid/hk/nu${obsid}A_fpm.hk" \
      clobber=yes \
      >> "$L3A_LOG" 2>&1

  echo "Running nuproducts FPMB" | tee -a "$L3B_LOG"
  nuproducts \
      indir="/data/fzhou/XSPEC/Cyg_X-1/NuSTAR/$obsid/event_cl" \
      instrument=FPMB \
      steminputs="nu$obsid" \
      outdir="/data/fzhou/XSPEC/Cyg_X-1/NuSTAR/$obsid/products" \
      srcregionfile="/data/fzhou/XSPEC/Cyg_X-1/NuSTAR/$obsid/event_cl/src.reg" \
      bkgregionfile="/data/fzhou/XSPEC/Cyg_X-1/NuSTAR/$obsid/event_cl/bkg.reg" \
      attfile="/data/fzhou/XSPEC/Cyg_X-1/NuSTAR/$obsid/auxil/nu${obsid}_att.fits" \
      hkfile="/data/fzhou/XSPEC/Cyg_X-1/NuSTAR/$obsid/hk/nu${obsid}B_fpm.hk" \
      clobber=yes \
      >> "$L3B_LOG" 2>&1

  # Wait for both FPMA and FPMB to complete
  wait
  
  echo "=== Completed OBSID $obsid ==="
}

# Find all 11-digit OBSID directories
obsids=($(ls -d [0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9] 2>/dev/null))

if [ ${#obsids[@]} -eq 0 ]; then
  echo "No OBSID directories found!"
  exit 1
fi

echo "Found ${#obsids[@]} OBSIDs to process (sequentially; FPMA/FPMB parallel within each)"
echo "Start time: $(date)"

# Process OBSIDs sequentially
for obsid in "${obsids[@]}"; do
  process_obsid "$obsid"
done

echo "All OBSIDs processed."
echo "End time: $(date)"
