import numpy as np
from lsst.sims.featureScheduler import BaseSurveyFeature, Base_basis_function, utils, features, Last_observed
import lsst.sims.featureScheduler as fs
import healpy as hp


default_nside = None

# Let's make it possible to drive a particular cadence that we like.


class Cadence_enhance_basis_function(Base_basis_function):
    """Drive a certain cadence"""
    def __init__(self, filtername='gri', nside=default_nside,
                 supress_window=[0, 1.8], supress_val=-0.5,
                 enhance_window=[2.1, 3.2], enhance_val=1.,
                 apply_area=None,
                 survey_features=None, condition_features=None):
        """
        Parameters
        ----------
        filtername : str ('gri')
            The filter(s) that should be grouped together
        supress_window : list of float
            The start and stop window for when observations should be repressed (days)
        apply_area : healpix map
            The area over which to try and drive the cadence. Good values as 1, no candece drive 0.
            Probably works as a bool array too.
        """
        self.nside = nside

        self.supress_window = np.sort(supress_window)
        self.supress_val = supress_val
        self.enhance_window = np.sort(enhance_window)
        self.enhance_val = enhance_val

        if survey_features is None:
            survey_features = {}
            survey_features['last_observed'] = Last_observed(filtername=filtername)
        if condition_features is None:
            condition_features = {}
            condition_features['Current_mjd'] = features.Current_mjd()

        super(Cadence_enhance_basis_function, self).__init__(survey_features=survey_features,
                                                             condition_features=condition_features)
        self.empty = np.zeros(hp.nside2npix(self.nside), dtype=float)
        # No map, try to drive the whole area
        if apply_area is None:
            self.apply_indx = np.arange(self.empty.size)
        else:
            self.apply_indx = np.where(apply_area != 0)[0]

    def __call__(self, indx=None):
        # copy an empty array
        result = self.empty + 0
        if indx is not None:
            ind = np.intersect1d(indx, self.apply_indx)
        else:
            ind = self.apply_indx
        if np.size(ind) == 0:
            result = 0
        else:
            mjd_diff = self.condition_features['Current_mjd'].feature - self.survey_features['last_observed'].feature[ind]
            to_supress = np.where((mjd_diff > self.supress_window[0]) & (mjd_diff < self.supress_window[1]))
            result[ind[to_supress]] = self.supress_val
            to_enhance = np.where((mjd_diff > self.enhance_window[0]) & (mjd_diff < self.enhance_window[1]))
            result[ind[to_enhance]] = self.enhance_val
        return result




