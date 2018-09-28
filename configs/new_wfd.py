"""
This is a configuration for a baseline simulation the Feature Based Scheduler. The if statement is used to bypass the
configuration expect when instantiated by the Feature Scheduler Driver. Note that we import and fill in a
SurveyTopology class that informs the (S)OCS about the projects defined on the configuration.

The only things that cannot be changed here are the names of the variables survey_topoly and scheduler. It is possible
as those are expected by the Driver. The way those objects are configured are entirely up to the user though.

09/07/2018 - Ribeiro, T.
"""
import numpy as np
import healpy as hp
import lsst.sims.featureScheduler as fs
from lsst.ts.scheduler.kernel import SurveyTopology
from astropy.coordinates import SkyCoord
from astropy import units as u
from lsst.sims.utils import hpid2RaDec


def generate_target_maps(nside, full_dec_min=-90., full_dec_max=+32., wfd_dec_min=-74.,
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
    ratio_maps = {}
    for filtername in filters:
        ratio_maps[filtername] = base_map.copy()
        ratio_maps[filtername][wfd_pix] *= wfd_weights[filtername]

    # Let's take the north out of u and y
    if mask_uy_north:
        far_north = np.where(dec > wfd_dec_max)
        ratio_maps['u'][far_north] = 0
        ratio_maps['y'][far_north] = 0

    norm_factor = fs.calc_norm_factor(ratio_maps)

    return ratio_maps, norm_factor


def generate_blob_surveys(nside):

    target_maps, norm_factor = generate_target_maps(nside)

    # set up a cloud map
    cloud_map = target_maps['r'][0]*0 + 0.7

    # Set up observations to be taken in blocks
    filter1s = ['u', 'g', 'r', 'i', 'z', 'y']
    filter2s = [None, 'g', 'r', 'i', None, None]

    pair_surveys = []
    for filtername, filtername2 in zip(filter1s, filter2s):
        bfs = []
        bfs.append(fs.M5_diff_basis_function(filtername=filtername, nside=nside))
        if filtername2 is not None:
            bfs.append(fs.M5_diff_basis_function(filtername=filtername2, nside=nside))
        bfs.append(fs.Target_map_basis_function(filtername=filtername,
                                                target_map=target_maps[filtername],
                                                out_of_bounds_val=hp.UNSEEN, nside=nside,
                                                norm_factor=norm_factor))
        if filtername2 is not None:
            bfs.append(fs.Target_map_basis_function(filtername=filtername2,
                                                    target_map=target_maps[filtername2],
                                                    out_of_bounds_val=hp.UNSEEN, nside=nside,
                                                    norm_factor=norm_factor))
        bfs.append(fs.Slewtime_basis_function(filtername=filtername, nside=nside))
        bfs.append(fs.Strict_filter_basis_function(filtername=filtername))
        # Masks, give these 0 weight
        bfs.append(fs.Zenith_shadow_mask_basis_function(nside=nside, shadow_minutes=60., max_alt=76.))
        bfs.append(fs.Moon_avoidance_basis_function(nside=nside, moon_distance=40.))
        bfs.append(fs.Bulk_cloud_basis_function(max_cloud_map=cloud_map, nside=nside))
        weights = np.array([3.0, 3.0, .3, .3, 3., 3., 0., 0., 0.])
        if filtername2 is None:
            # Need to scale weights up so filter balancing still works properly.
            weights = np.array([6.0, 0.6, 3., 3., 0., 0., 0.])
        if filtername2 is None:
            survey_name = 'blob, %s' % filtername
        else:
            survey_name = 'blob, %s%s' % (filtername, filtername2)
        pair_surveys.append(fs.Blob_survey(bfs, weights, filtername=filtername, filter2=filtername2,
                                           survey_note=survey_name, ignore_obs='DD'))

    return pair_surveys


def generate_greedy(nside):
    target_maps, norm_factor = generate_target_maps(nside)

    cloud_map = target_maps['r'][0]*0 + 0.7

    # filters = ['u', 'g', 'r', 'i', 'z', 'y']
    filters = ['g', 'r', 'i', 'z', 'y']
    surveys = []

    for filtername in filters:
        bfs = list()
        bfs.append(fs.M5_diff_basis_function(filtername=filtername, nside=nside))
        bfs.append(fs.Target_map_basis_function(filtername=filtername,
                                                target_map=target_maps[filtername],
                                                out_of_bounds_val=hp.UNSEEN, nside=nside,
                                                norm_factor=norm_factor))
        bfs.append(fs.Aggressive_Slewtime_basis_function(filtername=filtername, nside=nside, order=6., hard_max=120.))
        bfs.append(fs.Strict_filter_basis_function(filtername=filtername))
        bfs.append(fs.Zenith_shadow_mask_basis_function(nside=nside, shadow_minutes=0., max_alt=76.))
        bfs.append(fs.Bulk_cloud_basis_function(max_cloud_map=cloud_map, nside=nside))
        bfs.append(fs.Moon_avoidance_basis_function(nside=nside, moon_distance=40.))
        weights = np.array([3., 1., 3., 3., 0., 0, 0])
        surveys.append(fs.Greedy_survey_fields(bfs, weights, block_size=1,
                                               filtername=filtername, dither=True,
                                               nside=nside,
                                               ignore_obs='DD'))
    return surveys


def generate_dd(nside):

    # XXX-TODO:  Not sure what this should actually be.
    filter_prop = None

    # Set up the DD
    # ELAIS S1
    dd_surveys = list()
    dd_surveys.append(fs.Deep_drilling_survey(9.45, -44., sequence='rgizy',
                                              nvis=[20, 10, 20, 26, 20],
                                              survey_name='DD:ELAISS1', reward_value=100, moon_up=None,
                                              fraction_limit=0.148, ha_limits=([0., 0.5], [23.5, 24.]),
                                              nside=nside,
                                              avoid_same_day=True,
                                              filter_goals=filter_prop))
    dd_surveys.append(fs.Deep_drilling_survey(9.45, -44., sequence='u',
                                              nvis=[7],
                                              survey_name='DD:u,ELAISS1', reward_value=100, moon_up=False,
                                              fraction_limit=0.0012, ha_limits=([0., 0.5], [23.5, 24.]),
                                              nside=nside))

    # XMM-LSS
    dd_surveys.append(fs.Deep_drilling_survey(35.708333, -4 - 45 / 60., sequence='rgizy',
                                              nvis=[20, 10, 20, 26, 20],
                                              survey_name='DD:XMM-LSS', reward_value=100, moon_up=None,
                                              fraction_limit=0.148, ha_limits=([0., 0.5], [23.5, 24.]),
                                              nside=nside,
                                              avoid_same_day=True,
                                              filter_goals=filter_prop))
    dd_surveys.append(fs.Deep_drilling_survey(35.708333, -4 - 45 / 60., sequence='u',
                                              nvis=[7],
                                              survey_name='DD:u,XMM-LSS', reward_value=100, moon_up=False,
                                              fraction_limit=0.0012, ha_limits=([0., 0.5], [23.5, 24.]),
                                              nside=nside))

    # Extended Chandra Deep Field South
    # XXX--Note, this one can pass near zenith. Should go back and add better planning on this.
    dd_surveys.append(fs.Deep_drilling_survey(53.125, -28. - 6 / 60., sequence='rgizy',
                                              nvis=[20, 10, 20, 26, 20],
                                              survey_name='DD:ECDFS', reward_value=100, moon_up=None,
                                              fraction_limit=0.148, ha_limits=[[0.5, 1.0], [23., 22.5]],
                                              nside=nside,
                                              avoid_same_day=True,
                                              filter_goals=filter_prop))
    dd_surveys.append(fs.Deep_drilling_survey(53.125, -28. - 6 / 60., sequence='u',
                                              nvis=[7],
                                              survey_name='DD:u,ECDFS', reward_value=100, moon_up=False,
                                              fraction_limit=0.0012, ha_limits=[[0.5, 1.0], [23., 22.5]],
                                              nside=nside))
    # COSMOS
    dd_surveys.append(fs.Deep_drilling_survey(150.1, 2. + 10. / 60. + 55 / 3600., sequence='rgizy',
                                              nvis=[20, 10, 20, 26, 20],
                                              survey_name='DD:COSMOS', reward_value=100, moon_up=None,
                                              fraction_limit=0.148, ha_limits=([0., 0.5], [23.5, 24.]),
                                              nside=nside,
                                              avoid_same_day=True,
                                              filter_goals=filter_prop))
    dd_surveys.append(fs.Deep_drilling_survey(150.1, 2. + 10. / 60. + 55 / 3600., sequence='u',
                                              nvis=[7], ha_limits=([0., .5], [23.5, 24.]),
                                              survey_name='DD:u,COSMOS', reward_value=100, moon_up=False,
                                              fraction_limit=0.0012,
                                              nside=nside))

    # Extra DD Field, just to get to 5. Still not closed on this one
    dd_surveys.append(fs.Deep_drilling_survey(349.386443, -63.321004, sequence='rgizy',
                                              nvis=[20, 10, 20, 26, 20],
                                              survey_name='DD:290', reward_value=100, moon_up=None,
                                              fraction_limit=0.148, ha_limits=([0., 0.5], [23.5, 24.]),
                                              nside=nside,
                                              avoid_same_day=True,
                                              filter_goals=filter_prop))
    dd_surveys.append(fs.Deep_drilling_survey(349.386443, -63.321004, sequence='u',
                                              nvis=[7],
                                              survey_name='DD:u,290', reward_value=100, moon_up=False,
                                              fraction_limit=0.0012, ha_limits=([0., 0.5], [23.5, 24.]),
                                              nside=nside,
                                              filter_goals=filter_prop))
    return dd_surveys

if __name__ == 'config':
    survey_topology = SurveyTopology()
    survey_topology.num_general_props = 4
    survey_topology.general_propos = ["NorthEclipticSpur", "SouthCelestialPole",
                                      "WideFastDeep", "GalacticPlane"]
    survey_topology.num_seq_props = 1
    survey_topology.sequence_propos = ["DeepDrillingCosmology1"]

    nside = fs.set_default_nside(nside=32)  # Required

    greedy_surveys = generate_greedy(nside)
    blob_surveys = generate_blob_surveys(nside)
    dd_surveys = generate_dd(nside)

    scheduler = fs.Core_scheduler([dd_surveys, blob_surveys, greedy_surveys], nside=nside)  # Required
