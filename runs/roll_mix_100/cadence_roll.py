import numpy as np
import lsst.sims.featureScheduler as fs
from lsst.sims.speedObservatory import Speed_observatory
import matplotlib.pylab as plt
import healpy as hp
import time
from rolling_slair import Target_map_modulo_basis_function
from drive_cadence import Cadence_enhance_basis_function


# Try out the rolling cadence

survey_length = 365.25*10  # days
nside = fs.set_default_nside(nside=32)
years = np.round(survey_length/365.25)
t0 = time.time()

target_map = fs.standard_goals(nside=nside)
# Set up a map of where to drive the cadence
cadence_area = target_map['r'] * 0
cadence_area[np.where(target_map['r'] == 1)] = 1.
hp.mollview(cadence_area, title='Where to drive cadence')


even_year_target = fs.standard_goals(nside=nside)
odd_year_target = fs.standard_goals(nside=nside)

up = 2.0
down = 0.


# Let's find the healpix that divides the WFD area in half
wfd = even_year_target['r'] * 0
wfd[np.where(even_year_target['r'] == 1)] = 1
wfd_accum = np.cumsum(wfd)
hp.mollview(wfd_accum)
split_indx = np.max(np.where(wfd_accum < wfd_accum.max()/2.))

indx = np.arange(even_year_target['r'].size)
top_half_wfd = np.where((even_year_target['r'] == 1) & (indx <= split_indx))
bottom_half_wfd = np.where((even_year_target['r'] == 1) & (indx > split_indx))
top_half_wfd[0].size, bottom_half_wfd[0].size


for filtername in even_year_target:
    even_year_target[filtername][top_half_wfd] *= up
    even_year_target[filtername][bottom_half_wfd] *= down

    odd_year_target[filtername][top_half_wfd] *= down
    odd_year_target[filtername][bottom_half_wfd] *= up

even_norm = fs.calc_norm_factor(even_year_target)
odd_norm = fs.calc_norm_factor(odd_year_target)

surveys = []
mod_year = 2
offset = 1
# Set up observations to be taken in blocks
filter1s = ['u', 'g', 'r', 'i', 'z', 'y']
filter2s = [None, 'r', 'i', 'z', None, None]
for filtername, filtername2 in zip(filter1s, filter2s):
    bfs = []
    bfs.append(fs.M5_diff_basis_function(filtername=filtername, nside=nside))
    if filtername2 is not None:
        bfs.append(fs.M5_diff_basis_function(filtername=filtername2, nside=nside))
    bfs.append(Target_map_modulo_basis_function(filtername=filtername,
                                                target_map=even_year_target[filtername],
                                                mod_year=mod_year, offset=0,
                                                out_of_bounds_val=hp.UNSEEN, nside=nside,
                                                norm_factor=even_norm))
    if filtername2 is not None:
        bfs.append(Target_map_modulo_basis_function(filtername=filtername2,
                                                    target_map=even_year_target[filtername2],
                                                    mod_year=mod_year, offset=0,
                                                    out_of_bounds_val=hp.UNSEEN, nside=nside,
                                                    norm_factor=even_norm))

    bfs.append(Target_map_modulo_basis_function(filtername=filtername,
                                                target_map=odd_year_target[filtername],
                                                mod_year=mod_year, offset=offset,
                                                out_of_bounds_val=hp.UNSEEN, nside=nside,
                                                norm_factor=odd_norm))
    if filtername2 is not None:
        bfs.append(Target_map_modulo_basis_function(filtername=filtername2,
                                                    target_map=odd_year_target[filtername2],
                                                    mod_year=mod_year, offset=offset,
                                                    out_of_bounds_val=hp.UNSEEN, nside=nside,
                                                    norm_factor=odd_norm))
    bfs.append(fs.Slewtime_basis_function(filtername=filtername, nside=nside))
    bfs.append(fs.Strict_filter_basis_function(filtername=filtername))
    bfs.append(fs.Zenith_shadow_mask_basis_function(nside=nside, shadow_minutes=60., max_alt=76.))
    bfs.append(Cadence_enhance_basis_function(nside=nside, apply_area=cadence_area, filtername='gri'))
    weights = np.array([3.0, 3.0, .3, .3, 0.3, 0.3, 3., 3., 0., 3.])
    if filtername2 is None:
        # Need to scale weights up so filter balancing still works properly.
        weights = np.array([6.0, 0.6, 0.6, 3., 3., 0., 3.])
    # XXX-
    # This is where we could add a look-ahead basis function to include m5_diff in the future.
    # Actually, having a near-future m5 would also help prevent switching to u or g right at twilight?
    # Maybe just need a "filter future" basis function?
    if filtername2 is None:
        survey_name = 'blob, %s' % filtername
    else:
        survey_name = 'blob, %s%s' % (filtername, filtername2)
    surveys.append(fs.Blob_survey(bfs, weights, filtername=filtername, filter2=filtername2,
                                  survey_note=survey_name))

# Set up the greedy surveys for filling time when can't take pairs.
filters = ['u', 'g', 'r', 'i', 'z', 'y']
greedy_surveys = []
for filtername in filters:
    bfs = []
    bfs.append(fs.M5_diff_basis_function(filtername=filtername, nside=nside))
    bfs.append(Target_map_modulo_basis_function(filtername=filtername,
                                                mod_year=mod_year, offset=0,
                                                target_map=even_year_target[filtername],
                                                norm_factor=even_norm))

    bfs.append(Target_map_modulo_basis_function(filtername=filtername,
                                                mod_year=mod_year, offset=offset,
                                                target_map=even_year_target[filtername],
                                                out_of_bounds_val=hp.UNSEEN, nside=nside,
                                                norm_factor=odd_norm))

    bfs.append(fs.North_south_patch_basis_function(zenith_min_alt=50., nside=nside))
    bfs.append(fs.Slewtime_basis_function(filtername=filtername, nside=nside))
    bfs.append(fs.Strict_filter_basis_function(filtername=filtername))
    bfs.append(fs.Zenith_shadow_mask_basis_function(nside=nside, shadow_minutes=60., max_alt=76.))
    bfs.append(Cadence_enhance_basis_function(nside=nside, apply_area=cadence_area, filtername='gri'))
    weights = np.array([3.0, 0.3, 0.3, 1., 3., 3., 0., 3.])
    # Might want to try ignoring DD observations here, so the DD area gets covered normally--DONE
    sv = fs.Greedy_survey_fields(bfs, weights, block_size=1, filtername=filtername,
                                 dither=True, nside=nside, ignore_obs='DD')
    greedy_surveys.append(sv)

# Set up the DD surveys
dd_surveys = fs.generate_dd_surveys()

survey_list_o_lists = [dd_surveys, surveys, greedy_surveys]

scheduler = fs.Core_scheduler(survey_list_o_lists, nside=nside)
n_visit_limit = None
observatory = Speed_observatory(nside=nside, quickTest=True)
observatory, scheduler, observations = fs.sim_runner(observatory, scheduler,
                                                     survey_length=survey_length,
                                                     filename='rolling_mix_%iyrs.db' % years,
                                                     delete_past=True, n_visit_limit=n_visit_limit)
t1 = time.time()
delta_t = t1-t0
print('ran in %.1f min = %.1f hours' % (delta_t/60., delta_t/3600.))
