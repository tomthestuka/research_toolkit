#!/usr/bin/env python3
import argparse
import hashlib
import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import center_of_mass

def get_event_extension(hdul: list[fits.HDUList]) -> int:
    """Finds the HDU containing event data based on column names."""
    for i, hdu in enumerate(hdul):
        if hdu.data is not None and 'X' in hdu.columns.names and 'Y' in hdu.columns.names:
            return i
    raise ValueError("Could not find an extension with 'X' and 'Y' columns.")

def get_table_wcs(header, colnames):
    """Extract WCS from table column headers (for NuSTAR, etc.)"""
    # Find RA and DEC columns by searching TCTYP keywords
    ra_col, dec_col = None, None
    for key in header.keys():
        if key.startswith('TCTYP'):
            colnum = int(key[5:])
            if 'RA' in header[key]:
                ra_col = colnum
            elif 'DEC' in header[key]:
                dec_col = colnum
    
    if ra_col is None or dec_col is None:
        return None
    
    # Make sure X and Y correspond to these columns
    if colnames[ra_col-1] != 'X' or colnames[dec_col-1] != 'Y':
        return None
    
    # Extract WCS parameters
    ra_refpix = header.get(f'TCRPX{ra_col}')
    ra_refval = header.get(f'TCRVL{ra_col}')
    ra_delta = header.get(f'TCDLT{ra_col}')
    dec_refpix = header.get(f'TCRPX{dec_col}')
    dec_refval = header.get(f'TCRVL{dec_col}')
    dec_delta = header.get(f'TCDLT{dec_col}')
    
    return {
        'ra_refpix': ra_refpix, 'ra_refval': ra_refval, 'ra_delta': ra_delta,
        'dec_refpix': dec_refpix, 'dec_refval': dec_refval, 'dec_delta': dec_delta
    }

def run_extraction(evt_path: str, outdir: str, src_rad: float, bkg_inner_rad: float, bkg_outer_rad: float) -> dict:
    """Extract source and background regions.

    The regions depend on the chosen radii and the PNG shows them, so they are
    science artifacts: ``src.reg``, ``bkg.reg`` and ``verify.png`` go into
    *outdir*, never beside the event file.

    evt_path: [str]
        path to the FITS event file.

    outdir: [str]
        folder the region files (DS9 format) and the verification PNG are written to.

    src_rad: [float]
        source radius in arcminutes.

    bkg_inner_rad: [float]
        inner radius of the background annulus in arcminutes.

    bkg_outer_rad: [float]
        outer radius of the background annulus in arcminutes.

    Returns ``{'ra', 'dec', 'products': [src.reg, bkg.reg], 'figures': [verify.png], 'stdout',
    'exit_code'}``, each product ``{'path', 'caption'}`` and the figure a FigureRef dict.
    """
    os.makedirs(outdir, exist_ok=True)
    src_reg_out = os.path.join(outdir, "src.reg")
    bkg_reg_out = os.path.join(outdir, "bkg.reg")
    png_out = os.path.join(outdir, "verify.png")
    try:
        with fits.open(evt_path) as hdul:
            idx = get_event_extension(hdul)
            header = hdul[idx].header
            data = hdul[idx].data
            colnames = hdul[idx].columns.names
            
            x, y = data['X'], data['Y']
            
            # Mission detection for Plot Title
            mission = header.get('TELESCOP', 'Unknown Mission')
            
            # Try table column WCS first (NuSTAR, etc.)
            table_wcs = get_table_wcs(header, colnames)

        # Centroiding Logic
        bins = 1000
        hist, xedges, yedges = np.histogram2d(x, y, bins=bins)
        # NOTE: np.histogram2d(x,y) returns hist[i,j] where i is x-bin, j is y-bin
        # But for image display, we need hist[row, col] = hist[y, x], so transpose
        hist = hist.T  # Now hist[i,j] has i as y-bin, j as x-bin
        yi, xi = np.unravel_index(np.argmax(hist), hist.shape)
        
        # Local refine (window size 20 is safe for IXPE/NuSTAR PSF)
        win = 20
        sub_hist = hist[yi-win:yi+win, xi-win:xi+win]
        local_y, local_x = center_of_mass(sub_hist)
        
        final_x_px = xedges[xi-win] + local_x * (xedges[1]-xedges[0])
        final_y_px = yedges[yi-win] + local_y * (yedges[1]-yedges[0])
        
        # Sky Conversion
        if table_wcs:
            # Use table column WCS (simple linear transformation)
            ra = table_wcs['ra_refval'] + (final_x_px - table_wcs['ra_refpix']) * table_wcs['ra_delta']
            dec = table_wcs['dec_refval'] + (final_y_px - table_wcs['dec_refpix']) * table_wcs['dec_delta']
        else:
            # Fallback to standard WCS (for IXPE, etc.)
            wcs = WCS(header, keysel=['pixel'], naxis=2)
            ra, dec = wcs.all_pix2world(final_x_px, final_y_px, 0)
            ra, dec = float(ra), float(dec)

        # Output src.reg (convert arcmin to arcsec for DS9 format)
        src_rad_arcsec = src_rad * 60.0
        with open(src_reg_out, 'w') as f:
            f.write(f"fk5\ncircle({ra:.6f}, {dec:.6f}, {src_rad_arcsec:.2f}\")\n")

        # Output bkg.reg (convert arcmin to arcsec for DS9 format)
        bkg_inner_arcsec = bkg_inner_rad * 60.0
        bkg_outer_arcsec = bkg_outer_rad * 60.0
        with open(bkg_reg_out, 'w') as f:
            f.write(f"fk5\nannulus({ra:.6f}, {dec:.6f}, {bkg_inner_arcsec:.2f}\", {bkg_outer_arcsec:.2f}\")\n")

        # Visualization
        if png_out:
            fig, ax = plt.subplots(figsize=(7, 7))
            ax.imshow(np.log10(hist + 1), origin='lower', 
                      extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
            
            # Estimate pixel scale for circle overlay (degrees per pixel)
            if table_wcs:
                pix_scale = abs(table_wcs['ra_delta'])  # degrees per pixel
            else:
                pix_scale = abs(header.get('CDELT1', header.get('CD1_1', 1/3600)))
            
            # Source region
            radius_px = (src_rad / 60.0) / pix_scale  # src_rad in arcmin -> degrees -> pixels
            ax.add_patch(plt.Circle((final_x_px, final_y_px), radius_px, color='red', fill=False, lw=1.5, label='Source'))
            
            # Background annulus
            bkg_inner_px = (bkg_inner_rad / 60.0) / pix_scale
            bkg_outer_px = (bkg_outer_rad / 60.0) / pix_scale
            ax.add_patch(plt.Circle((final_x_px, final_y_px), bkg_inner_px, color='cyan', fill=False, lw=1.5, linestyle='--', label='Bkg'))
            ax.add_patch(plt.Circle((final_x_px, final_y_px), bkg_outer_px, color='cyan', fill=False, lw=1.5, linestyle='--'))
            
            ax.legend(loc='upper right')
            ax.set_title(f"{mission} | RA: {ra:.4f} Dec: {dec:.4f}")
            plt.savefig(png_out)
            plt.close(fig)

        stdout = f"[{mission}] RA: {ra}, Dec: {dec} -> {src_reg_out}, {bkg_reg_out}, {png_out}"
        print(stdout)
        return {
            "ra": float(ra),
            "dec": float(dec),
            "products": [
                {"path": src_reg_out, "caption": f"DS9 source region, {src_rad} arcmin circle"},
                {"path": bkg_reg_out, "caption": f"DS9 background annulus, {bkg_inner_rad}-{bkg_outer_rad} arcmin"},
            ],
            "figures": [{
                "path": png_out,
                "caption": f"{mission} event image with the source and background regions",
                "mime_type": "image/png",
                "produced_by": "xspec_utils.lib.find_src.run_extraction",
                "created_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                "sha256": _sha256(png_out),
            }],
            "stdout": stdout,
            "exit_code": 0,
        }

    except Exception as e:
        message = f"Error processing {evt_path}: {e}"
        print(message)
        return {"stdout": "", "stderr": message, "exit_code": 1}


def _sha256(path: str) -> str:
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("evt", help="FITS event file")
    parser.add_argument("--outdir", required=True,
                        help="folder for src.reg, bkg.reg and verify.png (apart from the observation database)")
    parser.add_argument("--src_rad", type=float, default=1.0, help="Source radius in arcminutes")
    parser.add_argument("--bkg_inner_rad", type=float, default=2.0, help="Background inner radius in arcminutes")
    parser.add_argument("--bkg_outer_rad", type=float, default=3.0, help="Background outer radius in arcminutes")
    args = parser.parse_args()
    sys.exit(run_extraction(args.evt, args.outdir, args.src_rad, args.bkg_inner_rad, args.bkg_outer_rad)["exit_code"])