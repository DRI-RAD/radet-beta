import datetime

import ee
import pytest

import radet.landsat as landsat
import radet.utils as utils
# TODO: import utils from openet.core
# import openet.core.utils as utils


# TODO: Try moving to conftest and/or make a fixture
SCENE_ID = 'LC08_044033_20170716'
# SCENE_ID = 'LC08_042035_20150713'
SCENE_DT = datetime.datetime.strptime(SCENE_ID[-8:], '%Y%m%d')
SCENE_DATE = SCENE_DT.strftime('%Y-%m-%d')
SCENE_TIME = utils.millis(SCENE_DT)


def sr_image(ub=0.0, blue=0.2, green=0.2, red=0.2, nir=0.7, swir1=0.2, swir2=0.2, bt=300, qa_pixel=0):
    """Construct a fake Landsat 8 image with renamed bands"""
    return (
        ee.Image.constant([ub, blue, green, red, nir, swir1, swir2, bt, qa_pixel])
        .rename(['ultra_blue', 'blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'tir', 'QA_PIXEL'])
        .set({
            'system:time_start': ee.Date(SCENE_DATE).millis(),
            'k1_constant': ee.Number(607.76),
            'k2_constant': ee.Number(1260.56),
        })
    )


# Test the static methods of the class first
# Do these need to be inside the TestClass?
def test_Image_evi2_band_name():
    output = landsat.evi2(sr_image()).getInfo()['bands'][0]['id']
    assert output == 'evi2'


@pytest.mark.parametrize(
    'red, nir, expected',
    [
        [0.0, 0.0, 0.0],
        [0.2, 0.7, 0.57],
        [0.2, 0.2, 0.0],
        [0.1, 0.9, 0.93],
        [0.001, 0.999, 1.25],
    ]
)
def test_Image_evi2_calculation(red, nir, expected, tol=0.01):
    output = utils.constant_image_value(landsat.evi2(sr_image(red=red, nir=nir)))
    assert abs(output['evi2'] - expected) <= tol


def test_Image_ndmi_scaled_band_name():
    output = landsat.ndmi_scaled(sr_image()).getInfo()['bands'][0]['id']
    assert output == 'ndmi_scaled'


@pytest.mark.parametrize(
    'nir, swir1, expected',
    [
        [0.0, 0.0, 1.0],
        [0.1, 0.9, 0.0],
        [0.7, 0.1, 1.0],
        [0.9, 0.1, 1.0],
        # [0.0, 0.1, ],
    ]
)
def test_Image_ndmi_scaled_calculation(nir, swir1, expected, tol=0.01):
    output = utils.constant_image_value(landsat.ndmi_scaled(sr_image(nir=nir, swir1=swir1)))
    assert abs(output['ndmi_scaled'] - expected) <= tol


def test_Image_lai_band_name():
    output = landsat.lai(sr_image()).getInfo()['bands'][0]['id']
    assert output == 'lai'


@pytest.mark.parametrize(
    'red, nir, swir1, expected',
    [
        [0.2, 0.7, 0.2, 3.17],
        [0.2, 0.2, 0.2, 0.18],
        # Check high values and clamping
        [0.1, 0.9, 0.1, 5.73],
        [0.01, 0.99, 0.01, 7.79],
        [0.001, 0.999, 0.001, 8.0],
        [0.1, 0.1, 0.9, 0.0],
    ]
)
def test_Image_lai_calculation(red, nir, swir1, expected, tol=0.01):
    output = utils.constant_image_value(landsat.lai(sr_image(red=red, nir=nir, swir1=swir1)))
    assert abs(output['lai'] - expected) <= tol


def test_Image_ndvi_band_name():
    output = landsat.ndvi(sr_image()).getInfo()['bands'][0]['id']
    assert output == 'ndvi'


@pytest.mark.parametrize(
    'red, nir, expected',
    [
        [0.2, 9.0 / 55, -0.1],
        [0.2, 0.2, 0.0],
        [0.1, 11.0 / 90,  0.1],
        [0.2, 0.3, 0.2],
        [0.1, 13.0 / 70, 0.3],
        [0.3, 0.7, 0.4],
        [0.2, 0.6, 0.5],
        [0.2, 0.8, 0.6],
        [0.1, 17.0 / 30, 0.7],
        [0.2, 0.7, 0.55555555],
        # TODO: Modify NDVI to vetter handle negative and very low values
        # # Check that negative values are not masked
        # [-0.01, 0.1, 1.0],
        # [0.1, -0.01, -1.0],
        # # Check that low values are set to 0
        # [-0.1, -0.1, 0.0],
        # [0.0, 0.0, 0.0],
        # [0.009, 0.009, 0.0],
        # [0.009, -0.01, 0.0],
        # [-0.01, 0.009, 0.0],
        # # Don't adjust NDVI if only one reflectance value is low
        # [0.005, 0.1, 0.9047619104385376],
    ]
)
def test_Image_ndvi_calculation(red, nir, expected, tol=0.000001):
    output = utils.constant_image_value(landsat.ndvi(sr_image(red=red, nir=nir)))
    assert abs(output['ndvi'] - expected) <= tol


def test_Image_ndwi_band_name():
    output = landsat.ndwi(sr_image()).getInfo()['bands'][0]['id']
    assert output == 'ndwi'


@pytest.mark.parametrize(
    'green, nir, expected',
    [
        [0.2, 0.2, 0.0],
        [0.2, 0.7, -0.56],
        [0.7, 0.2, 0.56],
    ]
)
def test_Image_ndwi_calculation(green, nir, expected, tol=0.01):
    output = utils.constant_image_value(landsat.ndwi(sr_image(green=green, nir=nir)))
    assert abs(output['ndwi'] - expected) <= tol


# TODO: Add test
# def test_landsat_c2_qa_water_mask():
#     output = landsat_c2_qa_water_mask()
#     return False


def test_Image_albedo_disalexi_band_name():
    output = landsat.albedo_disalexi(sr_image()).getInfo()['bands'][0]['id']
    assert output == 'albedo'


@pytest.mark.parametrize(
    'blue, green, red, nir, swir1, swir2, expected',
    [
        # [0.2, 0.2, 0.2, 0.7, 0.2, 0.2, 0.3555],
        [0.2, 0.2, 0.2, 0.7, 0.2, 0.2, 0.3879],
        [0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2],
    ]
)
def test_Image_albedo_disalexi_calculation(blue, green, red, nir, swir1, swir2, expected, tol=0.01):
    output = utils.constant_image_value(landsat.albedo_disalexi(
        sr_image(blue=blue, green=green, red=red, nir=nir, swir1=swir1, swir2=swir2)
    ))
    assert abs(output['albedo'] - expected) <= tol


def test_Image_albedo_l89_band_name():
    output = landsat.albedo_l89(sr_image()).getInfo()['bands'][0]['id']
    assert output == 'albedo'


@pytest.mark.parametrize(
    'ub, blue, green, red, nir, swir1, swir2, expected',
    [
        [0.2, 0.2, 0.2, 0.2, 0.7, 0.2, 0.2, 0.3403],
    ]
)
def test_Image_albedo_l89_calculation(ub, blue, green, red, nir, swir1, swir2, expected, tol=0.000001):
    output = utils.constant_image_value(landsat.albedo_l89(
        sr_image(ub=ub, blue=blue, green=green, red=red, nir=nir, swir1=swir1, swir2=swir2)
    ))
    assert abs(output['albedo'] - expected) <= tol


# def test_Image_albedo_l457_band_name():
#     output = landsat.albedo_l457(sr_image()).getInfo()['bands'][0]['id']
#     assert output == 'albedo'


# TODO: Add test
# def test_cloud_mask_sr_l457():
#     output = cloud_mask_sr_l457()
#     return False


# TODO: Add test
# def test_cloud_mask_sr_l8():
#     output = cloud_mask_sr_l8()
#     return False


# TODO: Add test
# def test_water_mask():
#     output = water_mask()
#     return False
