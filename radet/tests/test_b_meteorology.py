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
    'variable, date, xy, expected',
    [
        ['tmin', '2017-07-16', [-120.113, 36.336], 292.7],
        ['tmax', '2017-07-16', [-120.113, 36.336], 312.2],
        ['qa', '2017-07-16', [-120.113, 36.336], 0.00723],
        ['u10', '2017-07-16', [-120.113, 36.336], 4.3],
        ['srad', '2017-07-16', [-120.113, 36.336], 351.2],
    ]
)
def test_gridmet(variable, date, xy, expected, scale=4000, tol=0.001):
    output = meteorology.gridmet(
        variable=variable, time_start=utils.millis(datetime.strptime(date, '%Y-%m-%d'))
    )
    assert abs(utils.point_image_value(output, xy, scale)[variable] - expected) <= tol


@pytest.mark.parametrize(
    'source, variable, date, xy, expected',
    [
        ['IDAHO_EPSCOR/GRIDMET', 'tmin', '2017-07-16', [-120.113, 36.336], 292.7],
    ]
)
def test_get_source_variable(source, variable, date, xy, expected, scale=4000, tol=0.1):
    output = meteorology.get_source_variable(
        source=source, variable=variable, time_start=utils.millis(datetime.strptime(date, '%Y-%m-%d'))
    )
    assert abs(utils.point_image_value(output, xy, scale)[variable] - expected) <= tol
