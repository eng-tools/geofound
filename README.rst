***********
foundations
***********

A Python package to assess the settlement and bearing capacity of foundations.

How to Use
==========

The geofound package relies on two main Objects, the Foundation object and the Soil object.

These objects are inherited from the sfsimodels python module.

Examples
--------

.. code-block:: python

    import geofound
    import geofound.settlement
    import geofound.capacity
    import geofound.models as gm

    length = 21.7
    width = 1.5  # metres
    depth = 1.5  # metres
    phi = 38
    cohesion = 0
    unit_weight = 8.5  # submerged
    youngs_modulus_soil = 30000


    fd = gm.create_foundation(length, width, depth)
    sl = gm.create_soil(phi, cohesion, unit_weight)

    sl.unit_sat_weight = 18.5

    q_lim = geofound.capacity.vesics_capacity(sl, fd)
    Q_lim = q_lim*length*width

    print(' ')
    print('Ultimate bearing stress is q_lim = ' + str(round(q_lim,0)) + ' kPa')
    print('Ultimate load is Q_lim = ' + str(round(Q_lim,0)) + ' kN')


and for settlement.

.. code-block:: python
    load = 21484
    s = geofound.settlement.schmertmann_settlement(sl, fd, load, youngs_modulus_soil)
    print(' ')
    print('Settlement is si = ' + str(round(s,2)) + ' m')

Useful material
===============

http://geo.cv.nctu.edu.tw/foundation/download/BearingCapacityOfFoundations.pdf

http://geo.cv.nctu.edu.tw/foundation/download/FoundationsSettlements.pdf

http://geo.cv.nctu.edu.tw/foundation/download/SchmertmannMethod.pdf

http://www.civil.utah.edu/~bartlett/CVEEN5305/Handout%206%20-%20Schmertmann%20CPT%20method%20Example.pdf