# This script performs unzipping and instrumental background filtering for IXPE data pulled by the download_IXPE.sh.
# Place in /path/to/IXPE/filter_instrumental_bkg_IXPE.sh
# Automatically recognizes all obsid folders to perform unzipping and filtering.
# Refer to https://github.com/aledimarco/IXPE-background

##############################################################################
#                                                                            #
#   This script skips the obsid 03002599, which has weird event_l1 files.    #
#                                                                            #
##############################################################################


# Use by: $nohup ./filter_inst_bkg_IXPE.sh > filter_inst_bkg_IXPE.out 2>&1 &

#!/usr/bin/env bash

LOGFILE="filter_inst_bkg_IXPE.out"
: > "$LOGFILE"

# Runtime environment setup
source /opt/Anaconda/anaconda3/bin/activate
conda activate /opt/anaconda/anaconda3/envs/astropy

# Loop through each obsid directory in the base directory
for obsid_dir in [0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]; do
    [[ -d "$obsid_dir" ]] || continue
        echo "Processing $obsid_dir..."

    if [[ "$obsid_dir" == "03002599" ]]; then
        echo "Skipping $obsid_dir"
        continue
    fi

    for detid in {1..3}; do

        # Unzip the data files
        gunzip -dkf "${obsid_dir}/event_l1/ixpe${obsid_dir}_det${detid}_evt1_"v[0-9][0-9].fits.gz
        gunzip -dkf "${obsid_dir}/event_l2/ixpe${obsid_dir}_det${detid}_evt2_"v[0-9][0-9].fits.gz
        echo "Unzipped det${detid} files in $obsid_dir."

        # Filter the instrumental background
        python3 filter_background.py \
        "${obsid_dir}/event_l2/ixpe${obsid_dir}_det${detid}_evt2_"v[0-9][0-9].fits \
        "${obsid_dir}/event_l1/ixpe${obsid_dir}_det${detid}_evt1_"v[0-9][0-9].fits --output rej
        echo ""
        echo "---------------------------------------------------------------------------------------------------------"
        echo ""
    done
done
        