from datetime import datetime

import ee
import pytest

import radet.landsat as landsat
import radet.utils as utils
# TODO: import utils from openet.core
# import openet.core.utils as utils


# TODO: Try moving to conftest and/or make a fixture
SCENE_ID = 'LC08_044033_20170716'
SCENE_DT = datetime.strptime(SCENE_ID[-8:], '%Y%m%d')
SCENE_DATE = SCENE_DT.strftime('%Y-%m-%d')
SCENE_TIME = utils.millis(SCENE_DT)
# SCENE_TIME = 1500230731090


def sr_image(blue=0.2, green=0.2, red=0.2, nir=0.7, swir1=0.2, swir2=0.2, bt=300, qa_pixel=0):
    # Construct a fake Landsat 8 image with renamed bands
    # Excluding the ultra_blue band since it is not currently being used in the model
    return (
        ee.Image.constant([blue, green, red, nir, swir1, swir2, bt, qa_pixel])
        .rename(['blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'tir', 'QA_PIXEL'])
        .set({'system:time_start': ee.Date(SCENE_DATE).millis()})
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
        # Check negative reflectance values
        [-0.01, 0.1, 0.26],
        [0.1, -0.01, 0],
        [-0.01, 0.0, 0.03],
        [0.0, -0.01, 0],
        [-0.1, -0.1, 0],
        # Check saturated reflectance values
        # TODO: Check calculation for saturated reflectance values
        [1.1, 1.1, 0.0],
        [0.9, 1.1, 0.12],
        [1.1, 0.9, 0],
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
        [0.2, 0.4, 0.0],
        [0.2, 0.3, 0.33],
        [0.2, 0.25, 0.63],
        [0.7, 0.2, 1.0],
        [0.0, 0.0, 1.0],
        # Check negative reflectance values
        [-0.01, 0.1, 1.0],
        [0.1, -0.01, 1.0],
        [-0.01, 0.0, 1.0],
        [0.0, -0.01, 1.0],
        # Check saturated reflectance values
        [1.1, 1.1, 1.0],
        [0.9, 1.1, 0.67],
        [1.1, 0.9, 1.0],
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
        # Check negative reflectance values
        [-0.01, 0.1, 0.1, 1.09],
        [0.1, -0.01, 0.1, 0.18],
        [0.1, 0.1, -0.01, 0.18],
        # Check saturated reflectance values
        [1.1, 0.9, 0.9, 0.18],
        [0.9, 1.1, 0.9, 0.33],
        [0.9, 0.9, 1.1, 0.12],
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
        # NDVI is not being used in the model so extreme values are okay for now
        # Check negative reflectance values
        [-0.01, 0.1, 0],
        [0.1, -0.01, 0],
        # Check saturated reflectance values
        [1.1, 1.1, 0],
        [0.9, 1.1, 0.1],
        [1.1, 0.9, -0.1],
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
        # NDWI is not being used in the model so extreme values are okay for now
        # Check negative reflectance values
        [-0.01, 0.1, 0],
        [0.1, -0.01, 0],
        # Check saturated reflectance values
        [1.1, 1.1, 0],
        [0.9, 1.1, -0.1],
        [1.1, 0.9, 0.1],
    ]
)
def test_Image_ndwi_calculation(green, nir, expected, tol=0.01):
    output = utils.constant_image_value(landsat.ndwi(sr_image(green=green, nir=nir)))
    assert abs(output['ndwi'] - expected) <= tol


@pytest.mark.parametrize(
    'qa_pixel,  expected',
    [
        ['0000000010000000', 1],
        ['0000000000000000', 0],
    ]
)
def test_landsat_c2_qa_water_mask(qa_pixel, expected):
    output = utils.constant_image_value(
        landsat.landsat_c2_qa_water_mask(sr_image(qa_pixel=int(qa_pixel, 2)))
    )
    assert output['qa_water'] == expected


def test_Image_albedo_disalexi_band_name():
    assert landsat.albedo_disalexi(sr_image()).getInfo()['bands'][0]['id'] == 'albedo'


@pytest.mark.parametrize(
    'blue, green, red, nir, swir1, swir2, expected',
    [
        [0.2, 0.2, 0.2, 0.7, 0.2, 0.2, 0.3879],
        [0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2014],
    ]
)
def test_Image_albedo_disalexi_calculation(blue, green, red, nir, swir1, swir2, expected, tol=0.0001):
    output = utils.constant_image_value(landsat.albedo_disalexi(
        sr_image(blue=blue, green=green, red=red, nir=nir, swir1=swir1, swir2=swir2)
    ))
    assert abs(output['albedo'] - expected) <= tol


def test_Image_albedo_l89_band_name():
    # The test sr_image object is built without the ultra_blue band
    #   that is needed for this albedo calculation
    output = (
        landsat.albedo_l89(sr_image().addBands(ee.Image.constant([0]).rename('ultra_blue')))
        .getInfo()['bands'][0]['id']
    )
    assert output == 'albedo'


@pytest.mark.parametrize(
    'ub, blue, green, red, nir, swir1, swir2, expected',
    [
        [0.2, 0.2, 0.2, 0.2, 0.7, 0.2, 0.2, 0.3403],
    ]
)
def test_Image_albedo_l89_calculation(ub, blue, green, red, nir, swir1, swir2, expected, tol=0.0001):
    # The test sr_image object is built without the ultra_blue band
    #   that is needed for this albedo calculation
    output = utils.constant_image_value(landsat.albedo_l89(
        sr_image(blue=blue, green=green, red=red, nir=nir, swir1=swir1, swir2=swir2)
        .addBands(ee.Image.constant([ub]).rename('ultra_blue'))
    ))
    assert abs(output['albedo'] - expected) <= tol


def test_Image_albedo_l457_band_name():
    assert landsat.albedo_l457(sr_image()).getInfo()['bands'][0]['id'] == 'albedo'


@pytest.mark.parametrize(
    'ub, blue, green, red, nir, swir1, swir2, expected',
    [
        [0.2, 0.2, 0.2, 0.2, 0.7, 0.2, 0.2, 0.3555],
    ]
)
def test_Image_albedo_l457_calculation(ub, blue, green, red, nir, swir1, swir2, expected, tol=0.0001):
    output = utils.constant_image_value(landsat.albedo_l457(
        sr_image(blue=blue, green=green, red=red, nir=nir, swir1=swir1, swir2=swir2)
    ))
    assert abs(output['albedo'] - expected) <= tol


@pytest.mark.parametrize(
    'qa_pixel, expected',
    [
        # This cloud mask function returns 1 for non-masked pixels
        #   so that the output can be passed directly to an updateMask call
        # These are the two values that are not masked
        # All other combinations of bits will be masked
        ['0001010101000000', 1],  # Clear sky (5440)
        ['0001010110000000', 1],  # Water (5504)
        # Check the individual bits
        ['0000000000000000', 0],
        ['0000000000000001', 0],
        ['0000000000000010', 0],  # Dilate bit
        ['0000000000000100', 0],  # Cirrus bit
        ['0000000000001000', 0],  # Cloud bit
        ['0000000000010000', 0],  # Shadow bit
        ['0000000001000000', 0],  # Snow bit`
        ['0000000010000000', 0],  # Water bit`
    ]
)
def test_cloud_mask_C2_l457(qa_pixel, expected):
    output = utils.constant_image_value(
        landsat.cloud_mask_C2_l457(sr_image(qa_pixel=int(qa_pixel, 2)))
    )
    assert output['QA_PIXEL'] == expected


@pytest.mark.parametrize(
    'qa_pixel, expected',
    [
        # This cloud mask function returns 1 for non-masked pixels
        #   so that the output can be passed directly to an updateMask call
        # These are the two values that are not masked
        # All other combinations of bits will be masked
        ['0101010101000000', 1],  # Clear sky (21824)
        ['0101010111000000', 1],  # Water (21952)
        # Check the individual bits
        ['0000000000000000', 0],
        ['0000000000000001', 0],
        ['0000000000000010', 0],  # Dilate bit
        ['0000000000000100', 0],  # Cirrus bit
        ['0000000000001000', 0],  # Cloud bit
        ['0000000000010000', 0],  # Shadow bit
        ['0000000001000000', 0],  # Snow bit`
        ['0000000010000000', 0],  # Water bit`
    ]
)
def test_cloud_mask_C2_l89(qa_pixel, expected):
    output = utils.constant_image_value(
        landsat.cloud_mask_C2_l89(sr_image(qa_pixel=int(qa_pixel, 2)))
    )
    assert output['QA_PIXEL'] == expected


@pytest.mark.parametrize(
    'product, xy, expected',
    [
        # These water masks primarily identify ocean the great lakes
        ['GLO', [-120.0, 39.0], 0],
        ['GLO', [-121.0, 39.0], 0],
        ['GLO', [-125.0, 39.0], 1],
        ['OSM', [-120.0, 39.0], 0],
        ['OSM', [-121.0, 39.0], 0],
        ['OSM', [-125.0, 39.0], 1],
        # Check that lower case values are supported
        ['glo', [-125.0, 39.0], 1],
        ['osm', [-125.0, 39.0], 1],
    ]
)
def test_water_mask(product, xy, expected):
    output = utils.point_image_value(landsat.water_mask(product), xy, scale=10)
    assert output['water_mask'] == expected


def test_Image_meteorology_source_sources_exception():
    with pytest.raises(ValueError):
        utils.getinfo(landsat.water_mask(product=''))
