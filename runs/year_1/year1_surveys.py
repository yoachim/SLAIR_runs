import numpy as np
import lsst.sims.featureScheduler as fs
import lsst.sims.featureScheduler.utils as utils
from lsst.sims.speedObservatory import Speed_observatory
import matplotlib.pylab as plt
import healpy as hp
from lsst.sims.skybrightness_pre import M5percentiles
from lsst.sims.utils import hpid2RaDec
import lsst.sims.featureScheduler.features as features
import lsst.sims.skybrightness_pre as sb


class Current_seeing(fs.BaseConditionsFeature):
    def __init__(self, filtername='r', nside=None):
        if nside is None:
            nside = fs.utils.set_default_nside()
        self.nside = nside
        self.feature = None
        self.filtername = filtername

    def update_conditions(self, conditions, **kwargs):
        self.feature = conditions['FWHMeff_%s' % self.filtername]
        self.feature = hp.ud_grade(self.feature, nside_out=self.nside)


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

        self.condition_features = condition_features
        self.survey_features = survey_features

        if self.condition_features is None:
            self.condition_features = {}
            self.condition_features['Current_mjd'] = fs.features.Current_mjd()
        if self.survey_features is None:
            self.survey_features = {}

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

    def __call__(self, **kwards):
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
        self.nside = nside
        self.filtername = filtername

        if condition_features is None:
            self.condition_features = {}
            self.condition_features['Current_seeing'] = Current_seeing(filtername=filtername,
                                                                       nside=nside)
        if survey_features is None:
            self.survey_features = {}

        self.result_map = np.zeros(hp.nside2npix(self.nside))

    def check_feasibility(self):
        if np.max(self.__call__()) == hp.UNSEEN:
            return False
        else:
            return True

    def __call__(self, indx=None):
        result = self.result_map.copy()
        poor_seeing = np.where(self.condition_features['Current_seeing'].feature > self.max_seeing)
        result[poor_seeing] = hp.UNSEEN
        return result


class Nvis_limit_basis_function(Seeing_limit_basis_function):
    """Shut off observations after a given number of visits that meet criteria.
    """
    def __init__(self, nside=None, filtername='r', n_limit=3,
                 seeing_limit=1.2, time_lag=0.45,
                 m5_limit_map=None,
                 survey_features=None, condition_features=None,
                 **kwargs):
        """

        """
        if nside is None:
            nside = utils.set_default_nside()
        self.filtername = filtername
        self.n_limit = n_limit

        self.condition_features = condition_features
        self.survey_features = survey_features
        if self.survey_features is None:
            self.survey_features = {}
            self.survey_features['N_good'] = N_obs_good_conditions_feature(filtername=filtername,
                                                                           nside=nside,
                                                                           seeing_limit=seeing_limit,
                                                                           time_lag=time_lag,
                                                                           m5_limit_map=m5_limit_map,
                                                                           **kwargs)
        if self.condition_features is None:
            self.condition_features = {}
        self.result = np.zeros(hp.nside2npix(nside), dtype=float)

    def __call__(self, indx=None):
        result = self.result.copy()
        over_count = np.where(self.survey_features['N_good'].feature >= self.n_limit)
        result[over_count] = hp.UNSEEN
        return result


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
        self.percentile_limit = percentile_limit/100. #  Saved percentiles are between 0-1
        self.m5p = M5percentiles()
        self.condition_features = condition_features
        self.survey_features = survey_features
        if self.condition_features is None:
            self.condition_features = {}
            self.condition_features['m5depth'] = fs.M5Depth(filtername=filtername)
        if self.survey_features is None:
            self.survey_features = {}

        self.result_map = np.zeros(hp.nside2npix(nside))

    def __call__(self, indx=None):
        result = self.result_map.copy()
        per_map = self.m5p.m5map2percentile(self.condition_features['m5depth'].feature)
        below_limit = np.where(per_map < self.percentile_limit)
        result[below_limit] = hp.UNSEEN
        return result


class Limit_m5_map_basis_function(Seeing_limit_basis_function):
    """
    """
    def __init__(self, m5_limit, nside=None, filtername='r',
                 survey_features=None, condition_features=None, **kwargs):
        """
        Parameters
        ----------
        m5_limit : healpy array
            The faintest acceptable percentile 5-sigma depth to consider for observing.
            Anything below the limit will be masked.
        """
        if nside is None:
            nside = utils.set_default_nside()
        if hp.npix2nside(np.size(m5_limit)) != nside:
            raise ValueError('m5_limit map nside does not match basis function nside')
        self.filtername = filtername
        self.m5_limit = m5_limit
        self.condition_features = condition_features
        self.survey_features = survey_features
        if self.condition_features is None:
            self.condition_features = {}
            self.condition_features['m5depth'] = fs.M5Depth(filtername=filtername, nside=nside)
        if self.survey_features is None:
            self.survey_features = {}

        self.result_map = np.zeros(hp.nside2npix(nside))

    def __call__(self, indx=None):
        result = self.result_map.copy()
        diff = self.condition_features['m5depth'].feature - self.m5_limit
        below_limit = np.where(diff < 0)
        result[below_limit] = hp.UNSEEN
        return result


class N_obs_good_conditions_feature(fs.BaseSurveyFeature):
    """
    Track the number of observations that have been made accross the sky.
    """
    def __init__(self, filtername='r', seeing_limit=1.2, time_lag=0.45,
                 nside=None, m5_limit_map=None, mask_indx=None):
        """
        Parameters
        ----------
        filtername : str ('r')
            String or list that has all the filters that can count.
        seeing_limit : float (1.2)
            Only count an observation if the seeing is less than seeing_limit (arcsec).
            Uses the FWHMeff to compare.
        time_lag : float (0.45)
            Only count an observation if at least time_lag has elapsed (days).
        nside : int (32)
            The nside of the healpixel map to use
        m5_limit_map : healpix array (None)
            The 5-sigma limiting depth the observation must have to count as "good".
        mask_indx : list of ints (None)
            List of healpixel indices to mask and interpolate over
        """

        if nside is None:
            nside = utils.set_default_nside()

        self.feature = np.zeros(hp.nside2npix(nside), dtype=float)
        self.filtername = filtername
        self.mask_indx = mask_indx
        self.time_lag = time_lag
        self.seeing_limit = seeing_limit
        if m5_limit_map is None:
            self.m5_limit_map = np.zeros(hp.nside2npix(nside), dtype=float)
        else:
            self.m5_limit_map = m5_limit_map

        self.extra_features = {}
        # Note, this will only count last_observed in good conditions.
        self.extra_features['last_observed'] = fs.Last_observed(filtername=filtername, nside=nside)

    def _conditions_good(self, observation, indx):
        """Check if the observation counts as "good"
        """

        # How long has it been?
        dt = observation['mjd'] - self.extra_features['last_observed'].feature[indx]
        good_time = np.where(dt > self.time_lag)[0]
        # Is the observation deep enough?
        d_mag = observation['fivesigmadepth'] - self.m5_limit_map[indx]
        good_depth = np.where(d_mag > 0)[0]

        if (good_time.size > 0) & (good_depth.size > 0) & ((observation['FWHMeff'] - self.seeing_limit) < 0):
            result = True
        else:
            result = False

        return result

    def add_observation(self, observation, indx=None):
        """
        Parameters
        ----------
        indx : ints
            The indices of the healpixel map that have been observed by observation
        """

        if self.filtername is None or observation['filter'][0] in self.filtername:
            if self._conditions_good(observation, indx):
                self.feature[indx] += 1

                for feature in self.extra_features:
                        self.extra_features[feature].add_observation(observation, indx=indx)


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

    nside = nside
    filters = ['u', 'g', 'r', 'i', 'z', 'y']

    target_map = large_target_map(nside, dec_max=34.3)
    norm_factor = fs.calc_norm_factor({'r': target_map})

    # set up a cloud map
    cloud_map = target_map*0 + 0.7

    # Set up map m5-depth limits:
    m5_limits = {}
    percentile_cut = 0.7
    m52per = sb.M5percentiles()
    for filtername in filters:
        m5_limits[filtername] = m52per.percentile2m5map(percentile_cut, filtername=filtername, nside=nside)

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
        weights = [3.0, 0.3, 3., 3., 0, 0., 0.]
        # add in some constriants to make sure we only observe in good conditions and shut off after 3 good ones
        bfs.append(Limit_m5_map_basis_function(m5_limits[filtername], nside=nside, filtername=filtername))
        bfs.append(Seeing_limit_basis_function(nside=nside, filtername=filtername))
        bfs.append(Time_limit_basis_function(day_max=365.25))
        # XXX--Do I need a m5-depth limit on here too?
        bfs.append(Nvis_limit_basis_function(nside=nside, filtername=filtername, n_limit=3,
                                             seeing_limit=1.2, time_lag=0.45,
                                             m5_limit_map=m5_limits[filtername]))
        weights.extend([0, 0, 0, 0])
        #weights.extend([0., 0., 0.])

        weights = np.array(weights)
        # Might want to try ignoring DD observations here, so the DD area gets covered normally--DONE
        surveys.append(fs.Greedy_survey_fields(bfs, weights, block_size=1, filtername=filtername,
                                               dither=True, nside=nside, ignore_obs='DD',
                                               survey_name='templates'))


    # Do we want to cover all the potential area LSST could observe? In case a GW goes off
    # in the north not in the NES.
    return surveys


    # Don't let observations be taken in the same night--maybe just set to 13 hours or something.
    #fs.Avoid_Fast_Revists

    # Maybe this would be a good time to invoke the percentile limit again! That way we won't take
    # images in really poor depth conditions.
