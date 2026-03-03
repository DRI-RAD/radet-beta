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
        # LC08_044033_20170716 centroid
        ['IDAHO_EPSCOR/GRIDMET', [-121.668, 38.905], 7.08],
        ['GRIDMET', [-121.668, 38.905], 7.08],
        # LC08_042035_20150713 centroid?
        ['IDAHO_EPSCOR/GRIDMET', [-120.113, 36.336], 90.0],
        ['GRIDMET', [-120.113, 36.336], 90.0],
        #
        ['100', [-120.113, 36.336], 100.0],
        [100, [-120.113, 36.336], 100.0],
    ]
)
def test_elevation(source, xy, expected, scale=4000, tol=0.1):
    output = meteorology.elevation(source=source).rename('elevation')
    assert abs(utils.point_image_value(output, xy, scale)['elevation'] - expected) <= tol


@pytest.mark.parametrize(
    'variable, date, xy, expected',
    [
        ['tmin', '2017-07-16', [-121.668, 38.905], 289.6],
        ['tmax', '2017-07-16', [-121.668, 38.905], 312.4],
        ['qa', '2017-07-16', [-121.668, 38.905], 0.00937],
        ['u10', '2017-07-16', [-121.668, 38.905], 1.8],
        ['srad', '2017-07-16', [-121.668, 38.905], 353.7],
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
        ['IDAHO_EPSCOR/GRIDMET', 'tmin', '2017-07-16', [-121.668, 38.905], 289.6],
    ]
)
def test_get_source_variable(source, variable, date, xy, expected, scale=4000, tol=0.1):
    output = meteorology.get_source_variable(
        source=source, variable=variable, time_start=utils.millis(datetime.strptime(date, '%Y-%m-%d'))
    )
    assert abs(utils.point_image_value(output, xy, scale)[variable] - expected) <= tol
