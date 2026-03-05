from datetime import datetime
import math
import pprint

import ee
import pytest

import radet.model as model
import radet.utils as utils
# TODO: import utils from openet.core
# import openet.core.utils as utils


# TODO: Try moving to conftest and/or make a fixture
SCENE_ID = 'LC08_044033_20170716'
SCENE_DT = datetime.strptime(SCENE_ID[-8:], '%Y%m%d')
SCENE_DATE = SCENE_DT.strftime('%Y-%m-%d')
SCENE_TIME = utils.millis(SCENE_DT)
SCENE_TIME = 1500230731090
TEST_POINT = [-121.668, 38.905]


@pytest.mark.parametrize(
    'landcover, lai, wet, water',
    [
        [11, 0, 1, 1],
        [81, 0, 1, 0],
        [82, 0, 1, 0],
        [95, 0, 1, 0],
        [90, 0, 1, 0],
        [90, 0.9, 1, 0],
        [90, 1.0, 0, 0],
        [71, 0, 0, 0],
        [71, 4, 0, 0],
    ]
)
def test_wet_mask(landcover, lai, wet, water):
    output = utils.constant_image_value(model.wet_mask(
        landcover=ee.Image.constant(landcover), lai=ee.Number(lai)
    ))
    assert output['remapped'] == wet
    assert output['remapped_1'] == water


@pytest.mark.parametrize(
    'lai, tauL, tauS, fc',
    [
        # Check a low and high value
        [0.1, 0.909373, 0.945539, 0.039211],
        [3.5, 0.035973, 0.140858, 0.753403],
        # Check for low value clamping at LAI values outside clamped range
        [0.0, 1.0, 1.0, 0.01],
        [8.0, 0.01, 0.0113334, 0.959238],
        [10, 0.01, 0.01, 0.981684],
    ]
)
def test_transmissivities(lai, tauL, tauS, fc, tol=0.000001):
    output = utils.constant_image_value(model.transmissivities(ee.Image.constant(lai)))
    assert abs(output['tauL'] - tauL) < tol
    assert abs(output['tauS'] - tauS) < tol
    assert abs(output['fc'] - fc) < tol


@pytest.mark.parametrize(
    'lai, lst, albedo, emissivity, landcover, elevation, tmin, tmax, qa, u10, srad, '
    'meteo_elevation, time_start, latitude, longitude, expected',
    [
        # High ET (-121.668, 38.905)
        [
            3.5, 304, 0.15, 0.97, 82, 10, 290, 312, 0.007, 2, 350,
            10, SCENE_TIME, 38.910, -121.66, 7.6365
        ],
        # Low ET (-121.650, 38.913)
        [
            0.1, 326, 0.175, 0.97, 82, 10, 290, 312, 0.007, 2, 350,
            10, 1500230731090, 38.910, -121.66, 1.2952
        ],
        # Check ET if tmin=tmax to mimic test_c_image.py defaults
        [
            3.5, 304, 0.15, 0.97, 82, 10, 301, 301, 0.007, 2, 350,
            10, 1500230731090, 38.910, -121.66, 7.6685
        ],
    ]
)
def test_et(
        lai, lst, albedo, emissivity, landcover, elevation, tmin, tmax, qa, u10, srad,
        meteo_elevation, time_start, latitude, longitude, expected, tol=0.0001
):
    output = utils.constant_image_value(model.et(
        lai=ee.Image.constant(lai),
        lst=ee.Image.constant(lst),
        albedo=ee.Number(albedo),
        emissivity=ee.Number(emissivity),
        landcover=ee.Image.constant(landcover),
        elevation=ee.Number(elevation),
        tmin=ee.Image.constant(tmin),
        tmax=ee.Image.constant(tmax),
        qa=ee.Number(qa),
        u10=ee.Number(u10),
        srad=ee.Image.constant(srad),
        meteo_elevation=ee.Number(meteo_elevation),
        time_start=ee.Number(time_start),
        latitude=ee.Number(latitude),
        longitude=ee.Number(longitude),
    ))
    assert abs(output['ET'] - expected) < tol


@pytest.mark.parametrize(
    'lai, lst, albedo, emissivity, landcover, elevation, tmin, tmax, qa, u10, srad, '
    'meteo_elevation, time_start, latitude, longitude, expected',
    [
        # High ET site
        [
            3.5, 304, 0.15, 0.97, 82, 10, 290, 312, 0.007, 2, 350,
            10, 1500230731090, 38.910, -121.66, 7.6365
        ],
        # Low ET site
        [
            0.1, 326, 0.175, 0.97, 82, 10, 290, 312, 0.007, 2, 350,
            10, 1500230731090, 38.910, -121.66, 1.2952
        ],
    ]
)
def test_et_positional_arguments(
        lai, lst, albedo, emissivity, tmin, tmax, qa, u10, srad, landcover, elevation,
        meteo_elevation, time_start, latitude, longitude, expected, tol=0.0001
):
    output = utils.constant_image_value(model.et(
        ee.Image.constant(lai),
        ee.Image.constant(lst),
        ee.Number(albedo),
        ee.Number(emissivity),
        ee.Image.constant(landcover),
        ee.Number(elevation),
        ee.Image.constant(tmin),
        ee.Image.constant(tmax),
        ee.Number(qa),
        ee.Number(u10),
        ee.Image.constant(srad),
        ee.Number(meteo_elevation),
        ee.Number(time_start),
        ee.Number(latitude),
        ee.Number(longitude),
    ))
    assert abs(output['ET'] - expected) < tol
