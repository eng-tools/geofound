.. image:: https://travis-ci.org/eng-tools/geofound.svg?branch=master
   :target: https://travis-ci.org/eng-tools/geofound
   :alt: Testing Status

.. image:: https://img.shields.io/pypi/v/geofound.svg
   :target: https://pypi.python.org/pypi/geofound
   :alt: PyPi version

.. image:: https://coveralls.io/repos/github/eng-tools/geofound/badge.svg
   :target: https://coveralls.io/github/eng-tools/geofound

.. image:: https://img.shields.io/badge/license-MIT-blue.svg
    :target: https://github.com/eng-tools/geofound/blob/master/LICENSE
    :alt: License

.. image:: https://circleci.com/gh/eng-tools/geofound.svg?style=svg
    :target: https://circleci.com/gh/eng-tools/geofound
    :alt: Build Status

.. image:: https://zenodo.org/badge/112601559.svg
    :target: https://zenodo.org/badge/latestdoi/112601559
    :alt: DOI

********
geofound
********

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


    length = 21.7
    width = 1.5  # metres
    depth = 1.5  # metres
    phi = 38
    cohesion = 0
    unit_weight = 18.5  # submerged
    youngs_modulus_soil = 30000


    fd = geofound.create_foundation(length, width, depth)
    sl = geofound.create_soil(phi, cohesion, unit_weight)

    sl.unit_sat_weight = 18.5

    q_lim = geofound.capacity.capacity_vesics_1975(sl, fd)
    p_max = q_lim * length * width

    print(' ')
    print('Ultimate bearing stress is q_lim = ' + str(round(q_lim,0)) + ' kPa')
    print('Ultimate load is Q_lim = ' + str(round(p_max, 0)) + ' kN')


and for settlement.

.. code-block:: python

    load = 21484
    s = geofound.settlement.settlement_schmertmann(sl, fd, load, youngs_modulus_soil)
    print(' ')
    print('Settlement is si = ' + str(round(s,2)) + ' m')

Useful material
===============

http://geo.cv.nctu.edu.tw/foundation/download/BearingCapacityOfFoundations.pdf

http://geo.cv.nctu.edu.tw/foundation/download/FoundationsSettlements.pdf

http://geo.cv.nctu.edu.tw/foundation/download/SchmertmannMethod.pdf

http://www.civil.utah.edu/~bartlett/CVEEN5305/Handout%206%20-%20Schmertmann%20CPT%20method%20Example.pdf