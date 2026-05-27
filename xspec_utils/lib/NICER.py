import glob

def pha(obs: str, nicer_path: str, chatter=False) -> tuple[str, str]:
    """
    Get the path to the pha file for a given observation.

    Parameters
    ----------
    obs : str
        The observation ID.
    nicer_path : str
        The path to the NICER data directory.

    Returns
    -------
    tuple[str, str]
        The paths to the source and background pha files.
    """
    src_pha = glob.glob(f"{nicer_path}/{obs}/xti/event_cl/*sr.pha")
    bkg_pha = glob.glob(f"{nicer_path}/{obs}/xti/event_cl/*bg.pha")
    if len(src_pha) == 0:
        raise ValueError(f'No source pha file found for obs: {obs}\n'
                         f'the search was: {nicer_path}/{obs}/xti/event_cl/*sr.pha')
    elif len(src_pha) > 1:
        raise ValueError(f'Multiple source pha files found for obs: {obs}\n'
                         f'the search was: {nicer_path}/{obs}/xti/event_cl/*sr.pha')
    if len(bkg_pha) == 0:
        raise ValueError(f'No background pha file found for obs: {obs}\n'
                         f'the search was: {nicer_path}/{obs}/xti/event_cl/*bg.pha')
    elif len(bkg_pha) > 1:
        raise ValueError(f'Multiple background pha files found for obs: {obs}\n'
                         f'the search was: {nicer_path}/{obs}/xti/event_cl/*bg.pha')
    
    if chatter:
        print(f"NICER PHA query results for obs: {obs}")
        print(f"NICER Source PHA files")
        for f in src_pha:
            print(f"NICER Source PHA file: {f}")
        print(f"NICER Background PHA files")
        for f in bkg_pha:
            print(f"NICER Background PHA file: {f}")

    return src_pha[0], bkg_pha[0]

def rmf(obs: str, nicer_path: str, chatter=False) -> str:
    """
    Get the path to the rmf file for a given observation.

    Parameters
    ----------
    obs : str
        The observation ID.
    nicer_path : str
        The path to the NICER data directory.

    Returns
    -------
    str
        The path to the rmf file.
    """
    rmf_file = glob.glob(f"{nicer_path}/{obs}/xti/event_cl/*mpu7.rmf")
    if len(rmf_file) == 0:
        raise ValueError(f'No NICER rmf file found for obs: {obs}\n'
                         f'the search was: {nicer_path}/{obs}/xti/event_cl/*mpu7.rmf')
    elif len(rmf_file) > 1:
        raise ValueError(f'Multiple NICER rmf files found for obs: {obs}\n'
                         f'the search was: {nicer_path}/{obs}/xti/event_cl/*mpu7.rmf')
    
    if chatter:
        print(f"NICER RMF query results for obs: {obs}")
        print(f"NICER RMF files")
        for f in rmf_file:
            print(f"NICER RMF file: {f}")

    return rmf_file[0]

def arf(obs: str, nicer_path: str, chatter=False) -> str:
    """
    Get the path to the arf file for a given observation.

    Parameters
    ----------
    obs : str
        The observation ID.
    nicer_path : str
        The path to the NICER data directory.

    Returns
    -------
    str
        The path to the arf file.
    """
    arf_file = glob.glob(f"{nicer_path}/{obs}/xti/event_cl/*mpu7.arf")
    if len(arf_file) == 0:
        raise ValueError(f'No NICER arf file found for obs: {obs}\n'
                         f'the search was: {nicer_path}/{obs}/xti/event_cl/*mpu7.arf')
    elif len(arf_file) > 1:
        raise ValueError(f'Multiple NICER arf files found for obs: {obs}\n'
                         f'the search was: {nicer_path}/{obs}/xti/event_cl/*mpu7.arf')
    
    if chatter:
        print(f"NICER ARF query results for obs: {obs}")
        print(f"NICER ARF files")
        for f in arf_file:
            print(f"NICER ARF file: {f}")

    return arf_file[0]