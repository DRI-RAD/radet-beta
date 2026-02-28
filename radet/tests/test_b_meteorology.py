from datetime import datetime

import ee
import pytest

import radet.meteorology as meteorology
import radet.utils as utils
# TODO: import utils from openet.core
# import openet.core.utils as utils


@pytest.mark.parametrize(
    'source, xy, expected',
    [
        ['IDAHO_EPSCOR/GRIDMET', [-120.113, 36.336], 90.0],
        ['GRIDMET', [-120.113, 36.336], 90.0],
        # ['IDAHO_EPSCOR/GRIDMET', [-106.03249, 37.17777], 2402.9],
        # ['GRIDMET', [-106.03249, 37.17777], 2402.9],
    ]
)
def test_elevation(source, xy, expected, scale=4000, tol=0.1):
    output = meteorology.elevation(source=source).rename('elevation')
    assert abs(utils.point_image_value(output, xy, scale)['elevation'] - expected) <= tol


@pytest.mark.parametrize(
    'date, xy, variable, expected',
    [
        ['2017-07-16', [-120.113, 36.336], 'tmin', 292.7],
        ['2017-07-16', [-120.113, 36.336], 'tmax', 312.2],
        ['2017-07-16', [-120.113, 36.336], 'qa', 0.00723],
        ['2017-07-16', [-120.113, 36.336], 'u10', 4.3],
        ['2017-07-16', [-120.113, 36.336], 'srad', 351.2],
    ]
)
def test_gridmet(variable, date, xy, expected, scale=4000, tol=0.001):
    output = meteorology.gridmet(
        variable=variable, time_start=utils.millis(datetime.strptime(date, '%Y-%m-%d'))
    )
    assert abs(utils.point_image_value(output, xy, scale)[variable] - expected) <= tol
