from sfsimodels import models


class Soil(models.Soil):
    gwl = 100000.0  # ground water level


class RaftFoundation(models.RaftFoundation):
    q_ult = 0.0
    nc_factor = 0.0
    ng_factor = 0.0
    nq_factor = 0.0

    def __repr__(self):
        return self.length, self.width, self.height


class Foundation(RaftFoundation):

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
    a_foundation = RaftFoundation()
    a_foundation.length = length
    a_foundation.width = width
    a_foundation.depth = depth
    a_foundation.height = height
    return a_foundation


def create_soil(phi=0.0, cohesion=0.0, unit_dry_weight=0.0):
    """
    Can define a Soil object.
    :param phi: Internal friction angle
    :param cohesion: Cohesion of soil
    :param unit_dry_weight: The dry unit weight of the soil.
    :return: A Soil object.
    """
    soil = Soil()
    soil.phi = phi
    soil.cohesion = cohesion
    soil.unit_dry_weight = unit_dry_weight
    return soil
