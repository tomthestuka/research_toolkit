# HEASoft (PyXspec, FTOOLS) and CALDB for xspec_utils. Absolute paths: the runtime sources a copy.
export HEADAS=/opt/heasoft/heasoft-6.34/x86_64-pc-linux-gnu-libc2.34
. "$HEADAS/headas-init.sh"

export CALDB=/man1/CALDB
export CALDBCONFIG="$CALDB/software/tools/caldb.config"
export CALDBALIAS="$CALDB/software/tools/alias_config.fits"
export GEOMAG_PATH=/man1/CALDB/data/gen/pcf/geomag
