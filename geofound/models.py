from sfsimodels import models
from sfsimodels.models import SoilProfile, SoilStructureSystem

from geofound.exceptions import DesignError


class Soil(models.Soil):
    gwl = 100000.0  # ground water level


class FoundationRaft(models.RaftFoundation):
    q_ult = 0.0
    nc_factor = 0.0
    ng_factor = 0.0
    nq_factor = 0.0

    def __repr__(self):
        return self.length, self.width, self.height

class RaftFoundation(models.RaftFoundation):
    q_ult = 0.0
    nc_factor = 0.0
    ng_factor = 0.0
    nq_factor = 0.0

    def __repr__(self):
        return self.length, self.width, self.height


class Foundation(FoundationRaft):

    def __repr__(self):
        return self.length, self.width, self.height


class PadFoundation(models.PadFoundation):
    q_ult = 0.0
    nc_factor = 0.0
    ng_factor = 0.0
    nq_factor = 0.0

    def __repr__(self):
        return self.length, self.width, self.height


def create_foundation(length, width, depth=0.0, height=0.0):
    """
    Can define a Foundation Object from dimensions.
    :param length: Foundation length
    :param width: Foundation width
    :param depth: Foundation depth
    :param height: Foundation height
    :return: A Foundation object
    """
    a_foundation = FoundationRaft()
    a_foundation.length = length
    a_foundation.width = width
    a_foundation.depth = depth
    a_foundation.height = height
    return a_foundation


def create_soil(phi=0.0, cohesion=0.0, unit_dry_weight=0.0, pw=9800):
    """
    Can define a Soil object.
    :param phi: Internal friction angle
    :param cohesion: Cohesion of soil
    :param unit_dry_weight: The dry unit weight of the soil.
    :param pw: specific weight of water
    :return: A Soil object.
    """
    soil = Soil(pw=pw)
    soil.phi = phi
    soil.cohesion = cohesion
    soil.unit_dry_weight = unit_dry_weight
    return soil


def check_required(obj, required_parameters):
    """
    Check if a parameter is available on an object

    :param obj: Object
    :param required_parameters: list of parameters
    :return:
    """
    for parameter in required_parameters:
        if not hasattr(obj, parameter) or getattr(obj, parameter) is None:
            raise DesignError("parameter '%s' must be set for '%s' object." % (parameter, obj.base_type))
