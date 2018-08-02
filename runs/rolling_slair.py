import numpy as np
from lsst.sims.featureScheduler import BaseSurveyFeature, Base_basis_function, utils, features
import lsst.sims.featureScheduler as fs
import healpy as hp


default_nside = None


class N_observations_mod(BaseSurveyFeature):
    """
    Track the number of observations that have been made accross the sky, modulo year
    """
    def __init__(self, filtername=None, nside=default_nside,
                 mask_indx=None, mjd0=59580.035, mod_year=2, offset=0):
        """
        Parameters
        ----------
        filtername : str ('r')
            String or list that has all the filters that can count.
        nside : int (32)
            The nside of the healpixel map to use
        mjd0 : float
            The start of the survey
        mod_year : int
            Only record observations on these years.
        offset : int
            The offset to apply when calculating mod
        """
        if nside is None:
            nside = utils.set_default_nside()

        self.feature = np.zeros(hp.nside2npix(nside), dtype=float)
        self.filtername = filtername
        self.mjd0 = mjd0
        self.offset = offset
        self.mod_year = mod_year

    def add_observation(self, observation, indx=None):
        """
        Parameters
        ----------
        indx : ints
            The indices of the healpixel map that have been observed by observation
        """
        year = np.floor((observation['mjd'] - self.mjd0)/365.25)

        if (year + self.offset) % self.mod_year == 0:
            if observation['filter'][0] in self.filtername:
                self.feature[indx] += 1


class N_obs_count_mod(BaseSurveyFeature):
    """Count the number of observations.
    """
    def __init__(self, mod_year=2, offset=0, mjd0=59580.035,
                 filtername=None):
        self.feature = 0
        self.filtername = filtername
        self.mod_year = mod_year
        self.offset = offset
        self.mjd0 = mjd0

    def add_observation(self, observation, indx=None):
        # Track all observations if the year is correct
        year = np.floor((observation['mjd'] - self.mjd0)/365.25)

        if (year + self.offset) % self.mod_year == 0:
            if self.filtername is None:
                self.feature += 1
            else:
                if observation['filter'][0] in self.filtername:
                    self.feature += 1


class Target_map_modulo_basis_function(Base_basis_function):
    """Like target map, but modulo a year
    """
    def __init__(self, filtername='r', nside=default_nside, target_map=None,
                 survey_features=None, condition_features=None, norm_factor=0.00010519,
                 out_of_bounds_val=-10., mod_year=2, offset=0, mjd0=59580.035):
        """
        Parameters
        ----------
        filtername: (string 'r')
            The name of the filter for this target map.
        nside: int (default_nside)
            The healpix resolution.
        target_map : numpy array (None)
            A healpix map showing the ratio of observations desired for all points on the sky
        norm_factor : float (0.00010519)
            for converting target map to number of observations. Should be the area of the camera
            divided by the area of a healpixel divided by the sum of all your goal maps. Default
            value assumes LSST foV has 1.75 degree radius and the standard goal maps.
        out_of_bounds_val : float (-10.)
            Point value to give regions where there are no observations requested
        """
        if nside is None:
            nside = utils.set_default_nside()

        self.norm_factor = norm_factor
        if survey_features is None:
            survey_features = {}
            # Map of the number of observations in filter
            survey_features['N_obs'] = N_observations_mod(filtername=filtername,
                                                                        mod_year=mod_year,
                                                                        offset=offset,
                                                                        mjd0=mjd0, nside=nside)
            # Count of all the observations
            survey_features['N_obs_count_all'] = N_obs_count_mod(filtername=None,
                                                                               mod_year=mod_year,
                                                                               offset=offset,
                                                                               mjd0=mjd0)
        if condition_features is None:
            condition_features = {}
            condition_features['Current_mjd'] = features.Current_mjd()
        super(Target_map_modulo_basis_function, self).__init__(survey_features=survey_features,
                                                               condition_features=condition_features)
        self.nside = nside
        if target_map is None:
            self.target_map = utils.generate_goal_map(filtername=filtername)
        else:
            self.target_map = target_map
        self.out_of_bounds_area = np.where(self.target_map == 0)[0]
        self.out_of_bounds_val = out_of_bounds_val
        self.mjd0 = mjd0
        self.mod_year = mod_year
        self.offset = offset

    def __call__(self, indx=None):
        """
        Parameters
        ----------
        indx : list (None)
            Index values to compute, if None, full map is computed
        Returns
        -------
        Healpix reward map
        """

        # Check if the current year is one we should be calculating for
        year = np.floor((self.condition_features['Current_mjd'].feature - self.mjd0)/365.25)
        if (year + self.offset) % self.mod_year == 0:
            result = np.zeros(hp.nside2npix(self.nside), dtype=float)
            if indx is None:
                indx = np.arange(result.size)

            # Find out how many observations we want now at those points
            goal_N = self.target_map[indx] * self.survey_features['N_obs_count_all'].feature * self.norm_factor

            result[indx] = goal_N - self.survey_features['N_obs'].feature[indx]
            result[self.out_of_bounds_area] = self.out_of_bounds_val
        else:
            result = 0

        return result
