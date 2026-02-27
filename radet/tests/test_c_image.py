import datetime
# import logging
# import pprint

import ee
import pytest

import radet as radet
import radet.utils as utils
# TODO: import utils from openet.core
# import openet.core.utils as utils


# TODO: Try moving to conftest and/or make a fixture
COLL_ID = 'LANDSAT/LC08/C02/T1_L2/'
SCENE_ID = 'LC08_044033_20170716'
SCENE_DT = datetime.datetime.strptime(SCENE_ID[-8:], '%Y%m%d')
SCENE_DATE = SCENE_DT.strftime('%Y-%m-%d')
SCENE_DOY = int(SCENE_DT.strftime('%j'))
SCENE_0UTC_DT = datetime.datetime.strptime(SCENE_DATE, '%Y-%m-%d')
SCENE_TIME = 1500230731090
# TEST_POINT = (-121.5265, 38.7399)
TEST_POINT = [-120.113, 36.336]
SUN_ELEVATION = 64.27935371


# Should these be test fixtures instead?
# I'm not sure how to make them fixtures and allow input parameters
# CGM - This function is not currently used
# def sr_image(blue=0.2, green=0.2, red=0.2, nir=0.7, swir1=0.2, swir2=0.2, bt=300):
#     """Construct a fake Landsat 8 image with renamed bands"""
#     return (
#         ee.Image.constant([blue, green, red, nir, swir1, swir2, bt])
#         .rename(['blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'lst'])
#         .set({
#             'system:time_start': SCENE_TIME,
#             'k1_constant': ee.Number(607.76),
#             'k2_constant': ee.Number(1260.56)})
#     )


def default_image(albedo=0.2, emissivity=0.99, lai=3, lst=300, ndvi=0.8, ndwi=-0.5):
    # First construct a fake 'prepped' input image
    return (
        ee.Image.constant([albedo, emissivity, lai, lst, ndvi, ndwi])
        .rename(['albedo', 'emissivity', 'lai', 'lst', 'ndvi', 'ndwi'])
        .set({
            'system:index': SCENE_ID,
            'system:time_start': SCENE_TIME,
            'system:id': COLL_ID + SCENE_ID,
            'SUN_ELEVATION': SUN_ELEVATION,
        })
    )


# Setting etr_source and etr_band on the default image to simplify testing
#   but these do not have defaults in the Image class init
def default_image_args(
        albedo=0.2,
        emissivity=0.99,
        lai=3,
        lst=300,
        ndvi=0.8,
        ndwi=-0.5,
        meteorology_source='IDAHO_EPSCOR/GRIDMET',
        landcover=81,
        elevation=100,
        latitude=36,
        # longitude=-120,
        et_reference_source=10,
        et_reference_band='eto',
        et_reference_factor=1.0,
        et_reference_resample='nearest',
):
    return {
        'image': default_image(
            albedo=albedo, emissivity=emissivity, lai=lai, lst=lst, ndvi=ndvi, ndwi=ndwi
        ),
        'meteorology_source': meteorology_source,
        'landcover': landcover,
        'elevation': elevation,
        'latitude': latitude,
        # 'longitude': longitude,
        'et_reference_source': et_reference_source,
        'et_reference_band': et_reference_band,
        'et_reference_factor': et_reference_factor,
        'et_reference_resample': et_reference_resample,
    }


def default_image_obj(
        albedo=0.2,
        emissivity=0.99,
        lai=3,
        lst=300,
        ndvi=0.8,
        ndwi=-0.5,
        meteorology_source='IDAHO_EPSCOR/GRIDMET',
        landcover=81,
        elevation=100,
        latitude=36,
        # longitude=-120,
        et_reference_source=15,
        et_reference_band='etr',
        et_reference_factor=0.85,
        et_reference_resample='nearest',
):
    return radet.Image(**default_image_args(
        albedo=albedo,
        emissivity=emissivity,
        lai=lai,
        lst=lst,
        ndvi=ndvi,
        ndwi=ndwi,
        meteorology_source=meteorology_source,
        landcover=landcover,
        elevation=elevation,
        latitude=latitude,
        # longitude=longitude,
        et_reference_source=et_reference_source,
        et_reference_band=et_reference_band,
        et_reference_factor=et_reference_factor,
        et_reference_resample=et_reference_resample,
    ))


def test_Image_init_default_parameters():
    m = radet.Image(default_image())
    assert m.meteorology_source == 'IDAHO_EPSCOR/GRIDMET'
    # assert m.et_reference_source is None
    # assert m.et_reference_band is None
    # assert m.et_reference_factor is None
    # assert m.et_reference_factor is None
    # assert m.latitude is None
    # assert m.longitude is None


# Todo: Break these up into separate functions?
def test_Image_init_calculated_properties():
    m = radet.Image(default_image())
    assert utils.getinfo(m.time_start) == SCENE_TIME
    # assert m.scene_id.getInfo() == SCENE_ID
    # assert m.wrs2_tile.getInfo() == 'p{}r{}'.format(
    #     SCENE_ID.split('_')[1][:3], SCENE_ID.split('_')[1][3:])


def test_Image_init_date_properties():
    m = radet.Image(default_image())
    assert utils.getinfo(m.date)['value'] == SCENE_TIME
    assert utils.getinfo(m.year) == int(SCENE_DATE.split('-')[0])
    assert utils.getinfo(m.month) == int(SCENE_DATE.split('-')[1])
    assert utils.getinfo(m.start_date)['value'] == utils.millis(SCENE_0UTC_DT)
    assert utils.getinfo(m.end_date)['value'] == utils.millis(
        SCENE_0UTC_DT + datetime.timedelta(days=1))
    assert utils.getinfo(m.doy) == SCENE_DOY


# CGM - scene_id is not currently being used in the model
# def test_Image_init_scene_id_property():
#     """Test that the system:index from a merged collection is parsed"""
#     input_img = default_image()
#     m = radet.Image(input_img.set('system:index', '1_2_' + SCENE_ID))
#     assert utils.getinfo(m.scene_id) == SCENE_ID


@pytest.mark.parametrize('variable', ['albedo', 'emissivity', 'lai', 'lst', 'ndvi', 'ndwi'])
def test_Image_init_variable_properties(variable):
    """Test the band name and if properties are set on the variable images"""
    output = utils.getinfo(getattr(default_image_obj(), variable))
    # This assumes the band name is the lowercase of the variable name
    assert output['bands'][0]['id'] == variable.lower()
    assert output['properties']['system:index'] == SCENE_ID
    assert output['properties']['system:time_start'] == SCENE_TIME
    assert output['properties']['image_id'] == COLL_ID + SCENE_ID


@pytest.mark.parametrize(
    'meteorology_source, xy, expected',
    [
        ['IDAHO_EPSCOR/GRIDMET', TEST_POINT, 7.77],
        # ['projects/openet/assets/meteorology/era5land/na/daily', TEST_POINT, 1080.6605],
    ]
)
def test_Image_meteorology_source_sources(meteorology_source, xy, expected, tol=0.01):
    """Test daily meteorology source ET values for a single date at a real point"""
    m = default_image_obj(meteorology_source=meteorology_source)
    output = utils.point_image_value(ee.Image(m.et), xy, scale=30)
    assert abs(output['et'] - expected) <= tol


def test_Image_meteorology_source_sources_exception():
    with pytest.raises(ValueError):
        utils.getinfo(default_image_obj(meteorology_source='').et)


@pytest.mark.parametrize(
    'source, xy, expected',
    [
        ['USGS/SRTMGL1_003', TEST_POINT, 85],
        # CGM - This one causes ee.Initialization errors
        # [ee.Image('USGS/SRTMGL1_003'), TEST_POINT, 3],
        ['2364.351', TEST_POINT, 2364.351],
        [2364.351, TEST_POINT, 2364.351],
    ]
)
def test_Image_elevation_source(source, xy, expected, tol=0.001):
    output = utils.point_image_value(
        radet.Image(default_image(), elevation_source=source).elevation, xy
    )
    assert abs(output['elevation'] - expected) <= tol


# def test_Image_elevation_source_exception():
#     with pytest.raises(ValueError):
#         utils.getinfo(disalexi.Image(default_image(), elevation_source='').elevation)


def test_Image_elevation_band_name():
    output = utils.getinfo(radet.Image(default_image()).elevation)['bands'][0]['id']
    assert output == 'elevation'


@pytest.mark.parametrize(
    'source, xy, expected',
    [
        [
            'projects/sat-io/open-datasets/USGS/ANNUAL_NLCD/LANDCOVER/Annual_NLCD_LndCov_2020_CU_C1V1',
            TEST_POINT, 82
        ],
        ['projects/sat-io/open-datasets/USGS/ANNUAL_NLCD/LANDCOVER', TEST_POINT, 82],
        # CGM - Not currently support using the old NLCD collections
        # ['USGS/NLCD_RELEASES/2021_REL/NLCD/2021', TEST_POINT, 82],
        # ['USGS/NLCD_RELEASES/2021_REL/NLCD', TEST_POINT, 82],
        # ['USGS/NLCD_RELEASES/2019_REL/NLCD', TEST_POINT, 82],
        # ['USGS/NLCD_RELEASES/2019_REL/NLCD/2016', TEST_POINT, 82],
        # CGM - This one causes ee.Initialization errors
        # [ee.Image('USGS/SRTMGL1_003').multiply(0).add(82), TEST_POINT, 82],
        ['82', TEST_POINT, 82],
        [82, TEST_POINT, 82],
    ]
)
def test_Image_landcover_source(source, xy, expected, tol=0.001):
    output = utils.point_image_value(
        radet.Image(default_image(), landcover_source=source).landcover, xy
    )
    assert abs(output['landcover'] - expected) <= tol


@pytest.mark.parametrize(
    'source',
    [
        # No source
        '',
        # Fail on collection ID with trailing slash
        'USGS/NLCD_RELEASES/2019_REL/NLCD/',
    ]
)
def test_Image_landcover_source_invalid(source):
    with pytest.raises(ValueError):
        utils.getinfo(radet.Image(default_image(), landcover_source=source).landcover)


def test_Image_landcover_band_name():
    output = utils.getinfo(radet.Image(default_image()).landcover)['bands'][0]['id']
    assert output == 'landcover'


def test_Image_albedo_properties():
    """"""
    output = utils.getinfo(default_image_obj().albedo)
    assert output['bands'][0]['id'] == 'albedo'


def test_Image_albedo_value():
    """Test that a non-zero value is returned for the default inputs"""
    output = utils.constant_image_value(ee.Image(default_image_obj().albedo))
    assert output['albedo'] > 0


def test_Image_lai_properties():
    output = utils.getinfo(ee.Image(default_image_obj().lai))
    assert output['bands'][0]['id'] == 'lai'


def test_Image_lai_values():
    assert utils.constant_image_value(default_image_obj().lai)['lai'] == 3


def test_Image_lst_properties():
    output = utils.getinfo(ee.Image(default_image_obj().lst))
    assert output['bands'][0]['id'] == 'lst'


def test_Image_lst_values():
    assert utils.constant_image_value(default_image_obj().lst)['lst'] == 300


def test_Image_ndvi_properties():
    output = utils.getinfo(ee.Image(default_image_obj().ndvi))
    assert output['bands'][0]['id'] == 'ndvi'


def test_Image_ndvi_values():
    assert utils.constant_image_value(default_image_obj().ndvi)['ndvi'] == 0.8


def test_Image_ndwi_properties():
    output = utils.getinfo(ee.Image(default_image_obj().ndwi))
    assert output['bands'][0]['id'] == 'ndwi'


def test_Image_ndwi_values():
    assert utils.constant_image_value(default_image_obj().ndwi)['ndwi'] == -0.5


def test_Image_mask_properties():
    """Test if properties are set on the time image"""
    output = utils.getinfo(default_image_obj().mask)
    assert output['bands'][0]['id'] == 'mask'
    assert output['properties']['system:index'] == SCENE_ID
    assert output['properties']['system:time_start'] == SCENE_TIME
    assert output['properties']['image_id'] == COLL_ID + SCENE_ID


def test_Image_mask_values():
    assert utils.constant_image_value(default_image_obj().mask)['mask'] == 1


def test_Image_time_properties():
    """Test if properties are set on the time image"""
    output = utils.getinfo(default_image_obj().time)
    assert output['bands'][0]['id'] == 'time'
    assert output['properties']['system:index'] == SCENE_ID
    assert output['properties']['system:time_start'] == SCENE_TIME
    assert output['properties']['image_id'] == COLL_ID + SCENE_ID


def test_Image_time_values():
    """The time band should have the 0 UTC time in it for interpolation"""
    assert utils.constant_image_value(
        default_image_obj().time)['time'] == utils.millis(SCENE_0UTC_DT)


def test_Image_et_properties():
    """Test band name and if properties are set on the image"""
    output = utils.getinfo(default_image_obj().et)
    assert output['bands'][0]['id'] == 'et'
    assert output['properties']['system:index'] == SCENE_ID
    assert output['properties']['system:time_start'] == SCENE_TIME
    assert output['properties']['image_id'] == COLL_ID + SCENE_ID


@pytest.mark.xfail(reason="constant_image_value calls won't work with gridded meteorology")
def test_Image_et_defaults(expected=6.128, tol=0.001):
    output = utils.constant_image_value(ee.Image(default_image_obj().et))
    assert abs(output['et'] - expected) <= tol


# CGM - Test if there are any conditions that should return nodata
# @pytest.mark.parametrize(
#     'albedo, emissivity, lai, lst, ndvi, ndwi, expected',
#     [
#         [0.2, 0.99, 300, 0.80, None],
#     ]
# )
# def test_Image_et_nodata(albedo, emissivity, lst, ndvi, expected):
#     output_img = default_image_obj(albedo=albedo, emissivity=emissivity, lst=lst, ndvi=ndvi)
#     output = utils.constant_image_value(ee.Image(output_img.et))
#     assert output['et'] is None


def test_Image_et_reference_properties():
    """Test if properties are set on the reference ET image"""
    output = utils.getinfo(default_image_obj().et_reference)
    assert output['bands'][0]['id'] == 'et_reference'
    assert output['properties']['system:index'] == SCENE_ID
    assert output['properties']['system:time_start'] == SCENE_TIME
    assert output['properties']['image_id'] == COLL_ID + SCENE_ID


@pytest.mark.parametrize(
    'source, band, factor, xy, expected',
    [
        ['IDAHO_EPSCOR/GRIDMET', 'etr', 1, TEST_POINT, 12.9],
        ['IDAHO_EPSCOR/GRIDMET', 'etr', 0.85, TEST_POINT, 12.9 * 0.85],
        [
            'projects/openet/assets/reference_et/california/cimis/daily/v1',
            'etr', 1, TEST_POINT, 11.7893
        ],
        [10, 'FOO', 1, TEST_POINT, 10.0],
        [10, 'FOO', 0.85, TEST_POINT, 8.5],
    ]
)
def test_Image_et_reference_sources(source, band, factor, xy, expected, tol=0.001):
    """Test getting reference ET values for a single date at a real point"""
    output = utils.point_image_value(default_image_obj(
        et_reference_source=source, et_reference_band=band,
        et_reference_factor=factor).et_reference, xy)
    assert abs(output['et_reference'] - expected) <= tol


def test_Image_et_fraction_properties():
    """Test if properties are set on the ET fraction image"""
    output = utils.getinfo(default_image_obj().et_fraction)
    assert output['bands'][0]['id'] == 'et_fraction'
    assert output['properties']['system:index'] == SCENE_ID
    assert output['properties']['system:time_start'] == SCENE_TIME
    assert output['properties']['image_id'] == COLL_ID + SCENE_ID


def test_Image_calculate_properties():
    """Test if properties are set on the output image"""
    output = utils.getinfo(default_image_obj().calculate(['ndvi']))
    assert output['properties']['system:index'] == SCENE_ID
    assert output['properties']['system:time_start'] == SCENE_TIME
    assert output['properties']['image_id'] == COLL_ID + SCENE_ID


def test_Image_calculate_variables_default():
    output = utils.getinfo(default_image_obj().calculate())
    assert set([x['id'] for x in output['bands']]) == {'et'}


def test_Image_calculate_variables_custom():
    variables = {'ndvi'}
    output = utils.getinfo(default_image_obj().calculate(variables))
    assert set([x['id'] for x in output['bands']]) == variables


def test_Image_calculate_variables_all():
    variables = {'et', 'et_fraction', 'et_reference', 'ndvi'}
    output = utils.getinfo(default_image_obj().calculate(variables=list(variables)))
    assert set([x['id'] for x in output['bands']]) == variables


@pytest.mark.xfail(reason="constant_image_value calls won't work with gridded meteorology")
def test_Image_calculate_values():
    """Test if the calculate method returns values"""
    output_img = default_image_obj().calculate(['et'])
    assert utils.constant_image_value(output_img.select(['et']))['et'] > 0


def test_Image_calculate_variables_valueerror():
    """Test if calculate method raises a valueerror for invalid variables"""
    with pytest.raises(ValueError):
        utils.getinfo(default_image_obj().calculate(['FOO']))


def test_Image_from_landsat_c2_sr_default_image():
    """Test that the classmethod is returning a class object"""
    output = radet.Image.from_landsat_c2_sr(
        ee.Image('LANDSAT/LC08/C02/T1_L2/LC08_038031_20130828')
        # 'LANDSAT/LC08/C02/T1_L2/LC08_038031_20130828'
    )
    assert type(output) == type(default_image_obj())


@pytest.mark.parametrize(
    'image_id',
    [
        'LANDSAT/LT04/C02/T1_L2/LT04_044033_19830812',
        'LANDSAT/LT05/C02/T1_L2/LT05_044033_20110716',
        'LANDSAT/LE07/C02/T1_L2/LE07_044033_20170708',
        'LANDSAT/LC08/C02/T1_L2/LC08_044033_20170716',
        'LANDSAT/LC09/C02/T1_L2/LC09_044033_20220127',
    ]
)
def test_Image_from_landsat_c2_sr_landsat_image(image_id):
    """Test instantiating the class from a real Landsat images"""
    output = utils.getinfo(radet.Image.from_landsat_c2_sr(ee.Image(image_id)).ndvi)
    assert output['properties']['system:index'] == image_id.split('/')[-1]


def test_Image_from_landsat_c2_sr_exception():
    """Test instantiating the class for an invalid image ID"""
    with pytest.raises(Exception):
        # Intentionally using .getInfo() since utils.getinfo() might catch the exception
        radet.Image.from_landsat_c2_sr(ee.Image('FOO')).ndvi.getInfo()


@pytest.mark.xfail(reason="constant_image_value calls won't work with gridded meteorology")
def test_Image_from_landsat_c2_sr_scaling():
    """Test if Landsat SR images images are being scaled"""
    sr_img = ee.Image('LANDSAT/LC08/C02/T1_L2/LC08_044033_20170716')
    input_img = (
        ee.Image.constant([10909, 10909, 10909, 10909, 10909, 10909, 10909, 44177.6, 21824, 0, 9800])
        .rename(['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7',
                 'ST_B10', 'QA_PIXEL', 'QA_RADSAT', 'ST_EMIS'])
        .set({'SPACECRAFT_ID': ee.String(sr_img.get('SPACECRAFT_ID')),
              'SUN_ELEVATION': ee.String(sr_img.get('SUN_ELEVATION')),
              'system:id': ee.String(sr_img.get('system:id')),
              'system:index': ee.String(sr_img.get('system:index')),
              'system:time_start': ee.Number(sr_img.get('system:time_start'))})
    )

    # LST correction and cloud score masking do not work with a constant image
    #   and must be explicitly set to False
    output = utils.constant_image_value(radet.Image.from_landsat_c2_sr(
        input_img, c2_lst_correct=False,
        # cloudmask_args={'cloud_score_flag': False, 'filter_flag': False}
    ).albedo)
    assert abs(output['albedo'] - 0.1) <= 0.01

    output = utils.constant_image_value(radet.Image.from_landsat_c2_sr(
        input_img, c2_lst_correct=False,
        # cloudmask_args={'cloud_score_flag': False, 'filter_flag': False}
    ).lst)
    assert abs(output['lst'] - 300) <= 0.1


# def test_Image_from_landsat_c2_sr_cloud_mask_args():
#     """Test if the cloud_mask_args parameter can be set (not if it works)"""
#     output = radet.Image.from_landsat_c2_sr(
#         'LANDSAT/LC08/C02/T1_L2/LC08_044033_20170716',
#         cloudmask_args={'snow_flag': True, 'cirrus_flag': True}
#     )
#     assert type(output) == type(default_image_obj())
#
#
# def test_Image_from_landsat_c2_sr_cloud_score_mask_arg():
#     """Test if the cloud_score_flag parameter can be set in cloudmask_args"""
#     output = radet.Image.from_landsat_c2_sr(
#         'LANDSAT/LC08/C02/T1_L2/LC08_044033_20170716',
#         cloudmask_args={'cloud_score_flag': True, 'cloud_score_pct': 50})
#     assert type(output) == type(default_image_obj())


def test_Image_from_landsat_c2_sr_c2_lst_correct_arg():
    """Test if the c2_lst_correct parameter can be set (not if it works)"""
    output = radet.Image.from_landsat_c2_sr(
        'LANDSAT/LC08/C02/T1_L2/LC08_031034_20160702', c2_lst_correct=True)
    assert type(output) == type(default_image_obj())


def test_Image_from_landsat_c2_sr_c2_lst_correct_fill():
    """Test if the c2_lst_correct fills the LST holes in Nebraska"""
    image_id = 'LANDSAT/LC08/C02/T1_L2/LC08_031034_20160702'
    xy = [-102.08284, 37.81728]

    # CGM - Is the uncorrected test needed?
    uncorrected = utils.point_image_value(
        radet.Image.from_landsat_c2_sr(image_id, c2_lst_correct=False).lst, xy)
    assert uncorrected['lst'] is None
    corrected = utils.point_image_value(
        radet.Image.from_landsat_c2_sr(image_id, c2_lst_correct=True).lst, xy)
    assert corrected['lst'] > 0
    # # Exact test values copied from openet-core
    # assert abs(corrected['lst'] - 306.83) <= 0.25


@pytest.mark.parametrize(
    'image_id',
    [
        'LANDSAT/LC08/C02/T1_L2/LC08_044033_20170716',
        'LANDSAT/LC09/C02/T1_L2/LC09_044033_20220127',
    ]
)
def test_Image_from_image_id(image_id):
    """Test instantiating the class using the from_image_id method"""
    output = utils.getinfo(radet.Image.from_image_id(image_id).ndvi)
    assert output['properties']['system:index'] == image_id.split('/')[-1]
    assert output['properties']['image_id'] == image_id


# def test_Image_from_method_kwargs():
#     """Test that the init parameters can be passed through the helper methods"""
#     assert radet.Image.from_landsat_c2_sr(
#         'LANDSAT/LC08/C02/T1_L2/LC08_042035_20150713',
#         ea_source='FOO').ea_source == 'FOO'
