from datetime import datetime

import ee
import pytest

import radet.meteorology as meteorology
import radet.utils as utils
# TODO: import utils from openet.core
# import openet.core.utils as utils


# # TODO: Try moving to conftest and/or make a fixture
# SCENE_ID = 'LC08_044033_20170716'
# # SCENE_ID = 'LC08_042035_20150713'
# SCENE_DT = datetime.strptime(SCENE_ID[-8:], '%Y%m%d')
# SCENE_DATE = SCENE_DT.strftime('%Y-%m-%d')
# SCENE_TIME = utils.millis(SCENE_DT)


@pytest.mark.parametrize(
    'date, xy, scale, expected',
    [
        # CGM - Arbitrary test values just to make sure function is returning something reasonable
        [
            '2017-07-16', [-106.03249, 37.17777], 4000,
            {'srad': 196.8, 'tmmn': 281.8, 'tmmx': 297.3, 'sph': 0.01041, 'vs': 2.4}
        ],
    ]
)
def test_meteorology_gridmet(date, xy, scale, expected, tol=0.001):
    time_start = utils.millis(datetime.strptime(date, '%Y-%m-%d'))
    srad, tmin, tmax, qa, u10 = meteorology.gridmet(time_start=time_start)
    assert abs(utils.point_image_value(srad, xy, scale)['srad'] - expected['srad']) <= tol
    assert abs(utils.point_image_value(tmin, xy, scale)['tmmn'] - expected['tmmn']) <= tol
    assert abs(utils.point_image_value(tmax, xy, scale)['tmmx'] - expected['tmmx']) <= tol
    assert abs(utils.point_image_value(qa, xy, scale)['sph'] - expected['sph']) <= tol
    assert abs(utils.point_image_value(u10, xy, scale)['vs'] - expected['vs']) <= tol


# def test_meteorology_era5land():
#     # time_start, meteorology_source
#     output = meteorology.era5land()
#     # [tmin, tmax, tair_c, wind_med, rh, rso_inst, swdown24h, tfac]
#     assert False
