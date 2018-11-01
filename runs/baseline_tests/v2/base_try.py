import numpy as np
import lsst.sims.featureScheduler as fs
from lsst.sims.featureScheduler.utils import (sim_runner, set_default_nside, standard_goals,
                                              calc_norm_factor, generate_goal_map)
import lsst.sims.featureScheduler.surveys as survey
import lsst.sims.featureScheduler.basis_functions as basis_functions
from lsst.sims.speedObservatory import Speed_observatory
import matplotlib.pylab as plt
import healpy as hp
import time

# Let's add a zenith mask to see if that brings the slewtime down

t0 = time.time()

survey_length = 365.25*10  # days
nside = set_default_nside(nside=32)
# Define what we want the final visit ratio map to look like
years = np.round(survey_length/365.25)
# get rid of silly northern strip.
target_map = standard_goals(nside=nside)
norm_factor = calc_norm_factor(target_map)

# set up a cloud map
cloud_map = target_map['r']*0 + 0.7

# List to hold all the surveys (for easy plotting later)
surveys = []

# Set up observations to be taken in blocks
filter1s = ['u', 'g', 'r', 'i', 'z', 'y']
filter2s = [None, 'g', 'r', 'i', None, None]

pair_surveys = []
times_needed = [22., 44.]
for filtername, filtername2 in zip(filter1s, filter2s):
    bfs = []
    bfs.append(basis_functions.M5_diff_basis_function(filtername=filtername, nside=nside))
    if filtername2 is not None:
        bfs.append(basis_functions.M5_diff_basis_function(filtername=filtername2, nside=nside))
    bfs.append(basis_functions.Target_map_basis_function(filtername=filtername,
                                                         target_map=target_map[filtername],
                                                         out_of_bounds_val=hp.UNSEEN, nside=nside,
                                                         norm_factor=norm_factor))
    if filtername2 is not None:
        bfs.append(basis_functions.Target_map_basis_function(filtername=filtername2,
                                                             target_map=target_map[filtername2],
                                                             out_of_bounds_val=hp.UNSEEN, nside=nside,
                                                             norm_factor=norm_factor))
    bfs.append(basis_functions.Slewtime_basis_function(filtername=filtername, nside=nside))
    bfs.append(basis_functions.Strict_filter_basis_function(filtername=filtername))
    # Masks, give these 0 weight
    bfs.append(basis_functions.Zenith_shadow_mask_basis_function(nside=nside, shadow_minutes=60., max_alt=76.))
    bfs.append(basis_functions.Moon_avoidance_basis_function(nside=nside, moon_distance=30.))
    bfs.append(basis_functions.Bulk_cloud_basis_function(max_cloud_map=cloud_map, nside=nside))
    filternames = [fn for fn in [filtername, filtername2] if fn is not None]
    bfs.append(basis_functions.Filter_loaded_basis_function(filternames=filternames))
    if filtername2 is None:
        time_needed = times_needed[0]
    else:
        time_needed = times_needed[1]
    bfs.append(basis_functions.Time_to_twilight_basis_function(time_needed=time_needed))
    bfs.append(basis_functions.Not_twilight_basis_function(sun_alt_limit=-18.5))  #XXX--possible bug in pyephem?
    weights = np.array([3.0, 3.0, .3, .3, 3., 3., 0., 0., 0., 0., 0., 0.])
    if filtername2 is None:
        # Need to scale weights up so filter balancing still works properly.
        weights = np.array([6.0, 0.6, 3., 3., 0., 0., 0., 0., 0., 0.])
    # XXX-
    # This is where we could add a look-ahead basis function to include m5_diff in the future.
    # Actually, having a near-future m5 would also help prevent switching to u or g right at twilight?
    # Maybe just need a "filter future" basis function?
    if filtername2 is None:
        survey_name = 'blob, %s' % filtername
    else:
        survey_name = 'blob, %s%s' % (filtername, filtername2)
    surveys.append(survey.Blob_survey(bfs, weights, filtername1=filtername, filtername2=filtername2,
                                      survey_note=survey_name, ignore_obs='DD', dither=True))
    pair_surveys.append(surveys[-1])


# Let's set up some standard surveys as well to fill in the gaps. This is my old silly masked version.
# It would be good to put in Tiago's verion and lift nearly all the masking. That way this can also
# chase sucker holes.
#filters = ['u', 'g', 'r', 'i', 'z', 'y']
filters = ['i', 'z', 'y']

greedy_target_map = standard_goals(nside=nside)
# Let's take out the NES area on the target maps. This way we won't
# take images in the NES that aren't paired.
temp_map = generate_goal_map(nside=nside, NES_fraction=1.,
                             WFD_fraction=0, SCP_fraction=0,
                             GP_fraction=0, WFD_upper_edge_fraction=0.)
nes_pix = np.where(temp_map == 1)
for filtername in greedy_target_map:
    greedy_target_map[filtername][nes_pix] = 0

greedy_surveys = []
for filtername in filters:
    bfs = []
    bfs.append(basis_functions.M5_diff_basis_function(filtername=filtername, nside=nside))
    bfs.append(basis_functions.Target_map_basis_function(filtername=filtername,
                                                         target_map=greedy_target_map[filtername],
                                                         out_of_bounds_val=hp.UNSEEN, nside=nside,
                                                         norm_factor=norm_factor))

    #bfs.append(fs.North_south_patch_basis_function(zenith_min_alt=50., nside=nside))
    bfs.append(basis_functions.Slewtime_basis_function(filtername=filtername, nside=nside))
    bfs.append(basis_functions.Strict_filter_basis_function(filtername=filtername))
    bfs.append(basis_functions.Zenith_shadow_mask_basis_function(nside=nside, shadow_minutes=60., max_alt=76.))
    bfs.append(basis_functions.Moon_avoidance_basis_function(nside=nside, moon_distance=30.))
    bfs.append(basis_functions.Bulk_cloud_basis_function(max_cloud_map=cloud_map, nside=nside))

    bfs.append(basis_functions.Filter_loaded_basis_function(filternames=filtername))
    weights = np.array([3.0, 0.3, 3., 3., 0., 0., 0., 0.])
    # Might want to try ignoring DD observations here, so the DD area gets covered normally--DONE
    surveys.append(survey.Greedy_survey(bfs, weights, block_size=1, filtername=filtername,
                                        dither=True, nside=nside, ignore_obs='DD'))
    greedy_surveys.append(surveys[-1])

# Set up the DD surveys
dd_surveys = survey.generate_dd_surveys()
survey_list_o_lists = [dd_surveys, pair_surveys, greedy_surveys]

# Debug to stop at a spot if needed
n_visit_limit = None

# put in as list-of-lists so pairs get evaluated first.
scheduler = fs.Core_scheduler(survey_list_o_lists, nside=nside)
observatory = Speed_observatory(nside=nside, quickTest=True)
observatory, scheduler, observations = sim_runner(observatory, scheduler,
                                                  survey_length=survey_length,
                                                  filename='baseline_test%iyrs.db' % years,
                                                  delete_past=True, n_visit_limit=n_visit_limit)
t1 = time.time()
delta_t = t1-t0
print('ran in %.1f min = %.1f hours' % (delta_t/60., delta_t/3600.))
