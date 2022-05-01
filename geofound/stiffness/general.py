from geofound.stiffness.gazetas_1991 import calc_rot_via_gazetas_1991
from geofound.stiffness.pais_1988 import calc_rot_via_pais_1988


def rotational_stiffness(sl, fd, ip_axis="width", axis=None, a0=0.0, method='gazetas_1991', **kwargs):
    """
    Rotation stiffness of foundation.

    Parameters
    ----------
    fd: Foundation object
    sl: Soil Object.
    ip_axis: str
        The axis that is in the plane of deformation
    axis: str
        The axis which it should be computed around (if not None, then ip_axis is ignored)

    k_f: float
        Rotational stiffness of the foundation
    """
    if method == 'gazetas_1991':
        return calc_rot_via_gazetas_1991(sl, fd, ip_axis=ip_axis, axis=axis, a0=a0, **kwargs)
    else:
        return calc_rot_via_pais_1988(sl, fd, ip_axis=ip_axis, axis=axis, a0=a0, **kwargs)