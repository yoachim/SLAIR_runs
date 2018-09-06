import numpy as np
import lsst.sims.featureScheduler as fs
from lsst.sims.speedObservatory import Speed_observatory
import matplotlib.pylab as plt
import healpy as hp
from lsst.sims.skybrightness_pre import M5percentiles


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


class Nvis_limit_basis_function(fs.Seeing_limit_basis_function):
    """Shut off observations after a given number of visits
    """
    def __init__(self, ):
        pass


class Limit_m5_percentile_basis_function(fs.Seeing_limit_basis_function):
    """
    """
    def __init__(self, percentile_limit=60., nside=None,
                 survey_features=None, condition_features=None, **kwargs):
        
        self.m5p =  M5percentiles()

    def __call__(self):
        pass




def year_1_surveys(nside=32, mjd0=None):
    """
    Generate a list of surveys for executing in year 1
    """

    filters = ['u', 'g', 'r', 'i', 'z', 'y']


    # Don't let observations be taken in the same night--maybe just set to 13 hours or something.
    fs.Avoid_Fast_Revists

    # Maybe this would be a good time to invoke the percentile limit again! That way we won't take
    # images in really poor depth conditions.
