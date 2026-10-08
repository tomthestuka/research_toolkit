# XSPEC session setup for an analyst REPL in this environment (run once at REPL start).
# A failure is reported on stderr; the REPL still starts, and a tool needing XSPEC fails visibly.
import sys

try:
    from xspec import AllModels, Fit, Plot, Xset

    Fit.query = "no"
    Xset.allowPrompting = False
    Xset.chatter = 5
    Plot.xAxis = "keV"
except Exception as exc:
    print(f"boot.py: XSPEC session setup failed: {exc!r}", file=sys.stderr)
else:
    try:
        AllModels.lmod("relxill", "/data/fzhou/XSPEC/relxill/")
    except Exception as exc:
        print(f"boot.py: relxill local models not loaded: {exc!r}", file=sys.stderr)
