import glob

def pha(obs: str, nustar_path: str, chatter=False) -> tuple[list[str], list[str]]:
    """
    Get the paths to the pha files for a given observation.

    Parameters
    ----------
    obs : str
        The observation ID.
    nustar_path : str
        The path to the NuSTAR data directory.

    Returns
    -------
    tuple[list[str], list[str]]
        The paths to the source and background pha files. First list is source pha files (FPMA/FPMB), second list is background pha files (FPMA/FPMB).
    """
    src_pha = sorted(glob.glob(f"{nustar_path}/{obs}/products/*sr_grp.pha"))
    bkg_pha = sorted(glob.glob(f"{nustar_path}/{obs}/products/*bk.pha"))
    if len(src_pha) != 2:
        raise ValueError(
            f'Expected 2 source pha files (FPMA/FPMB) for obs: {obs}, but found {len(src_pha)}\n'
            f'the search was: {nustar_path}/{obs}/products/*sr_grp.pha')
    if len(bkg_pha) != 2:
        raise ValueError(
            f'Expected 2 background pha files (FPMA/FPMB) for obs: {obs}, but found {len(bkg_pha)}\n'
            f'the search was: {nustar_path}/{obs}/products/*bk.pha')
    
    if chatter:
        print(f"NuSTAR PHA query results for obs: {obs}")
        print("NuSTAR Source PHA files:")
        for f in src_pha:
            print(f"NuSTAR Source PHA file: {f}")
        print("NuSTAR Background PHA files:")
        for f in bkg_pha:
            print(f"NuSTAR Background PHA file: {f}")

    return src_pha, bkg_pha

def rmf(obs: str, nustar_path: str, chatter=False) -> list[str]:
    """
    Get the paths to the rmf files for a given observation, in order of FPMA, FPMB.

    Parameters
    ----------
    obs : str
        The observation ID.
    nustar_path : str
        The path to the NuSTAR data directory.

    Returns
    -------
    list[str]
        The paths to the rmf files.
    """
    rmf_files = sorted(glob.glob(f"{nustar_path}/{obs}/products/*rmf"))
    if len(rmf_files) != 2:
        raise ValueError(
            f'Expected 2 NuSTAR rmf files (FPMA/FPMB) for obs: {obs}, but found {len(rmf_files)}\n'
            f'the search was: {nustar_path}/{obs}/products/*rmf')
    
    if chatter:
        print(f"NuSTAR RMF query results for obs: {obs}")
        print("NuSTAR RMF files:")
        for f in rmf_files:
            print(f"NuSTAR RMF file: {f}")

    return rmf_files

def arf(obs: str, nustar_path: str, chatter=False) -> list[str]:
    """
    Get the paths to the arf files for a given observation, in order of FPMA, FPMB.

    Parameters
    ----------
    obs : str
        The observation ID.
    nustar_path : str
        The path to the NuSTAR data directory.

    Returns
    -------
    list[str]
        The paths to the arf files.
    """
    arf_files = sorted(glob.glob(f"{nustar_path}/{obs}/products/*arf"))
    if len(arf_files) != 2:
        raise ValueError(
            f'Expected 2 NuSTAR arf files (FPMA/FPMB) for obs: {obs}, but found {len(arf_files)}\n'
            f'the search was: {nustar_path}/{obs}/products/*arf')
    
    if chatter:
        print(f"NuSTAR ARF query results for obs: {obs}")
        print("NuSTAR ARF files:")
        for f in arf_files:
            print(f"NuSTAR ARF file: {f}")

    return arf_files