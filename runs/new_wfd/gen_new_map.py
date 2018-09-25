import numpy as np
import healpy as hp
from astropy.coordinates import SkyCoord
from astropy import units as u
from lsst.sims.utils import hpid2RaDec
import matplotlib.pylab as plt


def generate_goal(nside, full_dec_min=-90., full_dec_max=+32., wfd_dec_min=-74.,
                  wfd_dec_max=12.5, galactic_avoid=15., full_val=0.2, wfd_val=1.,
                  mask_uy_north = True,
                  wfd_weights={'u': 0.31, 'g': 0.44, 'r': 1., 'i': 1.,
                               'z': 0.9, 'y': 0.9}):
    """
    Generate a goal map that takes out a lot of the galactic plane

    Parameters
    ----------
    nside : int
        A valid healpix nside
    full_dec_min : float (-90.)
        How far to extend the entire survey to the south
    full_dec_max

    Useful notes:
    Observatory is at -30.23 degrees latitude.
    LMC is at dec = -69. SMC at -73. 
    """

    # The final map
    base_map = np.zeros(hp.nside2npix(nside), dtype=float)

    ra, dec = hpid2RaDec(nside, np.arange(hp.nside2npix(nside)))

    coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
    g_long, g_lat = coord.galactic.l.deg, coord.galactic.b.deg

    in_full_area = np.where((dec >= full_dec_min) & (dec <= full_dec_max))
    base_map[in_full_area] = full_val

    in_wfd = np.where((dec >= wfd_dec_min) & (dec <= wfd_dec_max))
    base_map[in_wfd] = wfd_val

    gp = np.where((np.abs(g_lat) < galactic_avoid) & (dec <= full_dec_max))
    base_map[gp] = full_val

    filters = ['u', 'g', 'r', 'i', 'z', 'y']
    wfd_pix = np.where(base_map == wfd_val)
    result = {}
    for filtername in filters:
        result[filtername] = base_map.copy()
        result[filtername][wfd_pix] *= wfd_weights[filtername]

    # Let's take the north out of u and y
    if mask_uy_north:
        far_north = np.where(dec > wfd_dec_max)
        result['u'][far_north] = 0
        result['y'][far_north] = 0

    return result
