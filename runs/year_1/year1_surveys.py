import numpy as np
import lsst.sims.featureScheduler as fs
import lsst.sims.featureScheduler.utils as utils
from lsst.sims.speedObservatory import Speed_observatory
import matplotlib.pylab as plt
import healpy as hp
from lsst.sims.skybrightness_pre import M5percentiles
from lsst.sims.utils import hpid2RaDec
import lsst.sims.featureScheduler.features as features

class Time_limit_basis_function(fs.Base_basis_function):
    """Limit how long a survey can run with a basis function
    """
    def __init__(self, mjd_start=None, day_max=365.25, day_min=0,
                 survey_features=None, condition_features=None, **kwargs):
        """
        Parameters
        ----------
        day_max : float (365.25)

        """

        self.mjd_start = mjd_start
        self.day_max = day_max
        self.day_min = day_min

        if condition_features is None:
            self.condition_features = {}
            self.condition_features['Current_mjd'] = fs.features.Current_mjd()

    def update_conditions(self, conditions):
        for feature in self.condition_features:
            self.condition_features[feature].update_conditions(conditions)
        # if we haven't set a start date, use the current conditions.
        # XXX-This might not be cold start compatible
        if self.mjd_start is None:
            self.mjd_start = self.condition_features['Current_mjd'].feature

    def check_feasibility(self):

        day = self.condition_features['Current_mjd'].feature - self.mjd_start
        if (day > self.day_max) | (day < self.day_min):
            result = False
        else:
            result = True
        return result

    def __call__(self):
        return 0


class Seeing_limit_basis_function(fs.Base_basis_function):
    def __init__(self, max_seeing=1.2, filtername='r', nside=None,
                 survey_features=None, condition_features=None, **kwargs):
        """
        Parameters
        ----------
        max_seeing : float (365.25)

        """
        if nside is None:
            nside = fs.utils.set_default_nside()
        self.max_seeing = max_seeing
        self.filtername = filtername

        if condition_features is None:
            self.condition_features = {}
            self.condition_features['Current_seeing'] = fs.features.Current_seeing(filtername=filtername,
                                                                                   nside=nside)
        self.result_map = np.zeros(hp.nside2npix(self.nside))

    def check_feasibility(self):
        if np.max(self()) == hp.UNSEEN:
            return False
        else:
            return True

    def __call__(self):
        result = self.result_map.copy()
        poor_seeing = np.where(self.condition_features['Current_seeing'] > self.max_seeing)
        result[poor_seeing] = hp.UNSEEN
        return result


class Nvis_limit_basis_function(Seeing_limit_basis_function):
    """Shut off observations after a given number of visits
    """
    def __init__(self, ):
        pass

    def __call__(self):
        pass


class Limit_m5_percentile_basis_function(Seeing_limit_basis_function):
    """
    """
    def __init__(self, percentile_limit=60., nside=None, filtername='r',
                 survey_features=None, condition_features=None, **kwargs):
        """
        Parameters
        ----------
        percentile_limit : float (60.)
            The lowest acceptable percentile 5-sigma depth to consider for observing.
            Anything below the limit will be masked.
        """
        if nside is None:
            nside = utils.set_default_nside()
        self.filtername = filtername
        self.percentile_limit = percentile_limit
        self.m5p = M5percentiles()
        if self.condition_features is None:
            self.condition_features = {}
            self.condition_features['m5depth'] = fs.M5Depth(filtername=filtername)
        self.result_map = np.zeros(hp.nside2npix(nside))

    def __call__(self):
        result = self.result_map.copy()
        per_map = self.m5p.m5map2percentile(self.condition_features['m5depth'])
        below_limit = np.where(per_map < self.percentile_limit)
        result[below_limit] = hp.UNSEEN
        return result


class N_obs_good_conditions_feature(fs.BaseSurveyFeature):
    """
    Track the number of observations that have been made accross the sky.
    """
    def __init__(self, filtername='r', nside=None, mask_indx=None):
        """
        Parameters
        ----------
        filtername : str ('r')
            String or list that has all the filters that can count.
        nside : int (32)
            The nside of the healpixel map to use
        mask_indx : list of ints (None)
            List of healpixel indices to mask and interpolate over
        """

        if nside is None:
            nside = utils.set_default_nside()

        self.feature = np.zeros(hp.nside2npix(nside), dtype=float)
        self.filtername = filtername
        self.mask_indx = mask_indx

        self.extra_features = {}
        self.extra_features['last_observed'] = 

    def _conditions_good(self, observation):
        """Check if the observation counts as "good"
        """

        result = True
        # Are the conditions good enough?

        # has it been long enough since it was last observed?

        return result

    def add_observation(self, observation, indx=None):
        """
        Parameters
        ----------
        indx : ints
            The indices of the healpixel map that have been observed by observation
        """

        if self.filtername is None or observation['filter'][0] in self.filtername:
            if self._conditions_good(observation):
                self.feature[indx] += 1


def large_target_map(nside, dec_max=34.3):
    """Create a large target map
    """
    result = np.ones(hp.nside2npix(nside), dtype=float)
    hpids = np.arange(result.size)
    ra, dec = hpid2RaDec(nside, hpids)
    out_of_bounds = np.where(dec > dec_max)
    result[out_of_bounds] = hp.UNSEEN
    return result


def year_1_surveys(nside=32, mjd0=None):
    """
    Generate a list of surveys for executing in year 1
    """

    nside = 32
    filters = ['u', 'g', 'r', 'i', 'z', 'y']

    target_map = large_target_map(nside, dec_max=34.3)
    norm_factor = fs.calc_norm_factor({'r': target_map})

    # set up a cloud map
    cloud_map = target_map*0 + 0.7

    surveys = []

    for filtername in filters:
        bfs = []
        bfs.append(fs.M5_diff_basis_function(filtername=filtername, nside=nside))
        bfs.append(fs.Target_map_basis_function(filtername=filtername,
                                                target_map=target_map,
                                                out_of_bounds_val=hp.UNSEEN, nside=nside,
                                                norm_factor=norm_factor))

        bfs.append(fs.Slewtime_basis_function(filtername=filtername, nside=nside))
        bfs.append(fs.Strict_filter_basis_function(filtername=filtername))
        bfs.append(fs.Zenith_shadow_mask_basis_function(nside=nside, shadow_minutes=0., max_alt=76.))
        bfs.append(fs.Moon_avoidance_basis_function(nside=nside, moon_distance=40.))
        bfs.append(fs.Bulk_cloud_basis_function(max_cloud_map=cloud_map, nside=nside))
        # add in some constriants to make sure we only observe in good conditions and shut off after 3 good ones

        weights = np.array([3.0, 0.3, 3., 3., 0, 0., 0.])
        # Might want to try ignoring DD observations here, so the DD area gets covered normally--DONE
        surveys.append(fs.Greedy_survey_fields(bfs, weights, block_size=1, filtername=filtername,
                                               dither=True, nside=nside, ignore_obs='DD'))


    # Do we want to cover all the potential area LSST could observe? In case a GW goes off
    # in the north not in the NES.
    return surveys


    # Don't let observations be taken in the same night--maybe just set to 13 hours or something.
    #fs.Avoid_Fast_Revists

    # Maybe this would be a good time to invoke the percentile limit again! That way we won't take
    # images in really poor depth conditions.
