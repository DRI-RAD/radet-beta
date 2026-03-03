from datetime import datetime, timedelta
import pprint

import ee
import pytest

import radet as model
import radet.utils as utils
# TODO: import utils from openet.core
# import openet.core.utils as utils


# TODO: Try moving to conftest and/or make a fixture
COLL_ID = 'LANDSAT/LC08/C02/T1_L2/'
SCENE_ID = 'LC08_044033_20170716'
SCENE_DT = datetime.strptime(SCENE_ID[-8:], '%Y%m%d')
SCENE_DATE = SCENE_DT.strftime('%Y-%m-%d')
SCENE_DOY = int(SCENE_DT.strftime('%j'))
SCENE_0UTC_DT = datetime.strptime(SCENE_DATE, '%Y-%m-%d')
SCENE_TIME = 1500230731090
TEST_POINT = [-121.668, 38.905]


def default_image(albedo=0.2, emissivity=0.99, lai=3, lst=300, ndvi=0.8, ndwi=-0.5):
    # First construct a fake 'prepped' input image
    return (
        ee.Image.constant([albedo, emissivity, lai, lst, ndvi, ndwi])
        .rename(['albedo', 'emissivity', 'lai', 'lst', 'ndvi', 'ndwi'])
        .set({
            'system:index': SCENE_ID,
            'system:time_start': SCENE_TIME,
            'system:id': COLL_ID + SCENE_ID,
            # 'SUN_ELEVATION': SUN_ELEVATION,
        })
    )


# Setting etr_source and etr_band on the default image to simplify testing
#   but these do not have defaults in the Image class init
def default_image_args(
        albedo=0.15,
        emissivity=0.97,
        lai=3.5,
        lst=304,
        ndvi=0.8,
        ndwi=-0.5,
        landcover_source=82,
        elevation_source=10,
        temperature_source=301,
        humidity_source=0.007,
        windspeed_source=2,
        solar_radiation_source=350,
        meteo_elevation_source=10,
        latitude=38.91,
        longitude=-121.66,
        et_reference_source=10,
        et_reference_band='eto',
        et_reference_factor=1.0,
        et_reference_resample='nearest',
        # Currently only a single temperature source is allowed,
        #   so keep track of tmin and tmax separately
        tmin_source=290,
        tmax_source=312,
):
    return {
        'image': default_image(
            albedo=albedo, emissivity=emissivity, lai=lai, lst=lst, ndvi=ndvi, ndwi=ndwi
        ),
        'landcover_source': landcover_source,
        'elevation_source': elevation_source,
        'temperature_source': temperature_source,
        'humidity_source': humidity_source,
        'windspeed_source': windspeed_source,
        'solar_radiation_source': solar_radiation_source,
        'meteo_elevation_source': meteo_elevation_source,
        'latitude': latitude,
        'longitude': longitude,
        'et_reference_source': et_reference_source,
        'et_reference_band': et_reference_band,
        'et_reference_factor': et_reference_factor,
        'et_reference_resample': et_reference_resample,
        'tmin_source': tmin_source,
        'tmax_source': tmax_source,
    }


def default_image_obj(
        albedo=0.15,
        emissivity=0.97,
        lai=3.5,
        lst=304,
        ndvi=0.8,
        ndwi=-0.5,
        landcover_source=82,
        elevation_source=10,
        temperature_source=301,
        humidity_source=0.007,
        windspeed_source=2,
        solar_radiation_source=350,
        meteo_elevation_source=10,
        latitude=38.91,
        longitude=-121.66,
        et_reference_source=10,
        et_reference_band='eto',
        et_reference_factor=1.0,
        et_reference_resample='nearest',
        # Currently only a single temperature source is allowed,
        #   so keep track of tmin and tmax separately
        tmin_source=290,
        tmax_source=312,
):
    return model.Image(**default_image_args(
        albedo=albedo,
        emissivity=emissivity,
        lai=lai,
        lst=lst,
        ndvi=ndvi,
        ndwi=ndwi,
        landcover_source=landcover_source,
        elevation_source=elevation_source,
        temperature_source=temperature_source,
        humidity_source=humidity_source,
        windspeed_source=windspeed_source,
        solar_radiation_source=solar_radiation_source,
        # meteo_elevation_source=meteo_elevation_source,
        latitude=latitude,
        longitude=longitude,
        et_reference_source=et_reference_source,
        et_reference_band=et_reference_band,
        et_reference_factor=et_reference_factor,
        et_reference_resample=et_reference_resample,
    ))


def test_Image_init_default_parameters():
    m = model.Image(default_image())
    assert m.temperature_source == 'IDAHO_EPSCOR/GRIDMET'
    assert m.humidity_source == 'IDAHO_EPSCOR/GRIDMET'
    assert m.windspeed_source == 'IDAHO_EPSCOR/GRIDMET'
    assert m.solar_radiation_source == 'IDAHO_EPSCOR/GRIDMET'
    # assert m.et_reference_source is None
    # assert m.et_reference_band is None
    # assert m.et_reference_factor is None
    # assert m.et_reference_factor is None
    # assert m.latitude is None
    # assert m.longitude is None


# Todo: Break these up into separate functions?
def test_Image_init_calculated_properties():
    m = model.Image(default_image())
    assert utils.getinfo(m.time_start) == SCENE_TIME
    # assert m.scene_id.getInfo() == SCENE_ID
    # assert m.wrs2_tile.getInfo() == 'p{}r{}'.format(
    #     SCENE_ID.split('_')[1][:3], SCENE_ID.split('_')[1][3:])


def test_Image_init_date_properties():
    m = model.Image(default_image())
    assert utils.getinfo(m.date)['value'] == SCENE_TIME
    assert utils.getinfo(m.year) == int(SCENE_DATE.split('-')[0])
    assert utils.getinfo(m.month) == int(SCENE_DATE.split('-')[1])
    assert utils.getinfo(m.start_date)['value'] == utils.millis(SCENE_0UTC_DT)
    assert utils.getinfo(m.end_date)['value'] == utils.millis(SCENE_0UTC_DT + timedelta(days=1))
    assert utils.getinfo(m.doy) == SCENE_DOY


def test_Image_init_scene_id_property():
    """Test that the system:index from a merged collection is parsed"""
    input_img = default_image()
    m = model.Image(input_img.set('system:index', '1_2_' + SCENE_ID))
    assert utils.getinfo(m.scene_id) == SCENE_ID


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
    'source, xy, scale, expected',
    [
        ['IDAHO_EPSCOR/GRIDMET', TEST_POINT, 4000, {'tmin': 293.3, 'tmax': 314.6}],
        ['300', TEST_POINT, 4000, {'tmin': 300, 'tmax': 300}],
        [300, TEST_POINT, 4000, {'tmin': 300, 'tmax': 300}],
    ]
)
def test_Image_temperature_source(source, xy, scale, expected, tol=0.001):
    output = utils.point_image_value(
        model.Image(default_image(), temperature_source=source).tmin, xy, scale
    )
    assert abs(output['tmin'] - expected['tmin']) <= tol
    output = utils.point_image_value(
        model.Image(default_image(), temperature_source=source).tmax, xy, scale
    )
    assert abs(output['tmax'] - expected['tmax']) <= tol


def test_Image_temperature_source_exception():
    with pytest.raises(ValueError):
        utils.getinfo(model.Image(default_image(), temperature_source='').tmin)


@pytest.mark.parametrize(
    'source, xy, scale, expected',
    [
        ['IDAHO_EPSCOR/GRIDMET', TEST_POINT, 4000, 0.00878],
        ['0.001', TEST_POINT, 4000, 0.001],
        [0.001, TEST_POINT, 4000, 0.001],
    ]
)
def test_Image_humidity_source_values(source, xy, scale, expected, tol=0.0001):
    output = utils.point_image_value(
        model.Image(default_image(), humidity_source=source).qa, xy, scale
    )
    assert abs(output['qa'] - expected) <= tol


def test_Image_humidity_source_exception():
    with pytest.raises(ValueError):
        utils.getinfo(model.Image(default_image(), humidity_source='FOO').qa)


@pytest.mark.parametrize(
    'source, xy, scale, expected',
    [
        ['IDAHO_EPSCOR/GRIDMET', TEST_POINT, 4000, 2.3],
        [2, TEST_POINT, 4000, 2],
    ]
)
def test_Image_windspeed_source_values(source, xy, scale, expected, tol=0.001):
    output = utils.point_image_value(
        model.Image(default_image(), windspeed_source=source).u10, xy, scale
    )
    assert abs(output['u10'] - expected) <= tol


def test_Image_windspeed_source_exception():
    with pytest.raises(ValueError):
        utils.getinfo(model.Image(default_image(), windspeed_source='FOO').u10)


@pytest.mark.parametrize(
    'source, xy, scale, expected',
    [
        ['IDAHO_EPSCOR/GRIDMET', TEST_POINT, 4000, 351.6],
        [300, TEST_POINT, 4000, 300],
    ]
)
def test_Image_solar_radiation_source_values(source, xy, scale, expected, tol=0.01):
    output = utils.point_image_value(
        model.Image(default_image(), solar_radiation_source=source).srad, xy, scale
    )
    assert abs(output['srad'] - expected) <= tol


def test_Image_solar_radiation_source_exception():
    with pytest.raises(ValueError):
        utils.getinfo(model.Image(default_image(), solar_radiation_source='FOO').srad)


@pytest.mark.parametrize(
    'source, xy, expected',
    [
        ['IDAHO_EPSCOR/GRIDMET', TEST_POINT, 7.08],
        [10, TEST_POINT, 10],
    ]
)
def test_Image_meteo_elevation_source_not_set(source, xy, expected, tol=0.01):
    # If the meteorology elevation source is not set, default to the temperature source
    output = utils.point_image_value(
        model.Image(default_image(), temperature_source=source).meteo_elevation, xy
    )
    assert abs(output['elevation'] - expected) <= tol


@pytest.mark.parametrize(
    'source, xy, expected',
    [
        ['IDAHO_EPSCOR/GRIDMET', TEST_POINT, 7.08],
        ['2364.351', TEST_POINT, 2364.35],
        [2364.351, TEST_POINT, 2364.35],
    ]
)
def test_Image_meteo_elevation_source_values(source, xy, expected, tol=0.01):
    output = utils.point_image_value(
        model.Image(default_image(), meteo_elevation_source=source).meteo_elevation, xy
    )
    assert abs(output['elevation'] - expected) <= tol


@pytest.mark.parametrize(
    'source, xy, expected',
    [
        ['USGS/SRTMGL1_003', TEST_POINT, 4],
        # CGM - This one causes ee.Initialization errors
        # [ee.Image('USGS/SRTMGL1_003'), TEST_POINT, 3],
        ['2364.351', TEST_POINT, 2364.35],
        [2364.351, TEST_POINT, 2364.35],
    ]
)
def test_Image_elevation_source_values(source, xy, expected, tol=0.01):
    output = utils.point_image_value(
        model.Image(default_image(), elevation_source=source).elevation, xy
    )
    assert abs(output['elevation'] - expected) <= tol


def test_Image_elevation_source_exception():
    with pytest.raises(ValueError):
        utils.getinfo(model.Image(default_image(), elevation_source=None).elevation)


def test_Image_elevation_band_name():
    output = utils.getinfo(model.Image(default_image()).elevation)['bands'][0]['id']
    assert output == 'elevation'


@pytest.mark.parametrize(
    'source, xy, expected',
    [
        [
            'projects/sat-io/open-datasets/USGS/ANNUAL_NLCD/LANDCOVER/Annual_NLCD_LndCov_2020_CU_C1V1',
            TEST_POINT, 82
        ],
        ['projects/sat-io/open-datasets/USGS/ANNUAL_NLCD/LANDCOVER', TEST_POINT, 82],
        # CGM - This one causes ee.Initialization errors
        # [ee.Image('USGS/SRTMGL1_003').multiply(0).add(82), TEST_POINT, 82],
        ['82', TEST_POINT, 82],
        [82, TEST_POINT, 82],
    ]
)
def test_Image_landcover_source_values(source, xy, expected, tol=0.001):
    output = utils.point_image_value(
        model.Image(default_image(), landcover_source=source).landcover, xy
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
def test_Image_landcover_source_exception(source):
    with pytest.raises(ValueError):
        utils.getinfo(model.Image(default_image(), landcover_source=source).landcover)


def test_Image_landcover_band_name():
    output = utils.getinfo(model.Image(default_image()).landcover)['bands'][0]['id']
    assert output == 'landcover'


def test_Image_albedo_default_value_set():
    # Test that the lazy property returns the default test value and the band name is set
    assert utils.constant_image_value(ee.Image(default_image_obj().albedo))['albedo'] > 0


def test_Image_emissivity_default_value_set():
    # Test that the lazy property returns the default test value and the band name is set
    assert utils.constant_image_value(ee.Image(default_image_obj().emissivity))['emissivity'] > 0


def test_Image_lai_default_value_set():
    # Test that the lazy property returns the default test value and the band name is set
    assert utils.constant_image_value(ee.Image(default_image_obj().lai))['lai'] > 0


def test_Image_lst_default_value_set():
    # Test that the lazy property returns the default test value and the band name is set
    assert utils.constant_image_value(ee.Image(default_image_obj().lst))['lst'] > 0


def test_Image_ndvi_default_value_set():
    # Test that the lazy property returns the default test value and the band name is set
    assert utils.constant_image_value(ee.Image(default_image_obj().ndvi))['ndvi'] > 0


def test_Image_ndwi_default_value_set():
    # Test that the lazy property returns the default test value and the band name is set
    assert utils.constant_image_value(ee.Image(default_image_obj().ndwi))['ndwi'] >= -1


def test_Image_mask_properties():
    """Test if properties are set on the time image"""
    output = utils.getinfo(default_image_obj().mask)
    assert output['bands'][0]['id'] == 'mask'
    assert output['properties']['system:index'] == SCENE_ID
    assert output['properties']['system:time_start'] == SCENE_TIME
    assert output['properties']['image_id'] == COLL_ID + SCENE_ID


def test_Image_mask_default_value():
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


def test_Image_et_defaults(expected=7.6685, tol=0.0001):
    # Note that tmin and tmax are the same t_avg value in this test
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
        ['IDAHO_EPSCOR/GRIDMET', 'eto', 1, TEST_POINT, 8.2],
        ['IDAHO_EPSCOR/GRIDMET', 'etr', 1, TEST_POINT, 10.8],
        ['IDAHO_EPSCOR/GRIDMET', 'etr', 0.85, TEST_POINT, 10.8 * 0.85],
        [
            'projects/openet/assets/reference_et/california/cimis/daily/v1',
            'etr', 1, TEST_POINT, 9.8313
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
    output = model.Image.from_landsat_c2_sr(
        ee.Image('LANDSAT/LC08/C02/T1_L2/LC08_038031_20130828')
        # 'LANDSAT/LC08/C02/T1_L2/LC08_038031_20130828'
    )
    assert type(output) == type(default_image_obj())


@pytest.mark.parametrize(
    'image_id',
    [
        # 'LANDSAT/LT04/C02/T1_L2/LT04_044033_19830812',
        'LANDSAT/LT05/C02/T1_L2/LT05_044033_20110716',
        'LANDSAT/LE07/C02/T1_L2/LE07_044033_20170708',
        'LANDSAT/LC08/C02/T1_L2/LC08_044033_20170716',
        'LANDSAT/LC09/C02/T1_L2/LC09_044033_20220127',
    ]
)
def test_Image_from_landsat_c2_sr_landsat_image(image_id):
    """Test instantiating the class from a real Landsat images"""
    output = utils.getinfo(model.Image.from_landsat_c2_sr(ee.Image(image_id)).lai)
    assert output['properties']['system:index'] == image_id.split('/')[-1]


def test_Image_from_landsat_c2_sr_exception():
    """Test instantiating the class for an invalid image ID"""
    with pytest.raises(Exception):
        # Intentionally using .getInfo() since utils.getinfo() might catch the exception
        model.Image.from_landsat_c2_sr(ee.Image('FOO')).lai.getInfo()


def test_Image_from_landsat_c2_sr_scaling():
    """Test if Landsat SR images images are being scaled"""
    sr_img = ee.Image('LANDSAT/LC08/C02/T1_L2/LC08_044033_20170716')
    input_img = (
        ee.Image.constant([10909, 10909, 10909, 10909, 10909, 10909, 10909, 44177.6, 21824, 0, 9800])
        .rename(['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7',
                 'ST_B10', 'QA_PIXEL', 'QA_RADSAT', 'ST_EMIS'])
        .set({'SPACECRAFT_ID': ee.String(sr_img.get('SPACECRAFT_ID')),
              'system:id': ee.String(sr_img.get('system:id')),
              'system:index': ee.String(sr_img.get('system:index')),
              'system:time_start': ee.Number(sr_img.get('system:time_start'))})
    )

    # LST correction, ocean masking, and cloud score masking
    #   do not work with a constant image and must be explicitly set to False
    # Other source inputs must be set to constant values
    output = utils.constant_image_value(model.Image.from_landsat_c2_sr(
        input_img, c2_lst_correct=False, mask_ocean_flag=False,
        latitude=36, longitude=-120, landcover_source=81, elevation_source=100,
        temperature_source=300, humidity_source=0.007, windspeed_source=4,
        solar_radiation_source=350,
        # cloudmask_args={'cloud_score_flag': False, 'filter_flag': False}
    ).albedo)
    assert abs(output['albedo'] - 0.1) <= 0.01

    output = utils.constant_image_value(model.Image.from_landsat_c2_sr(
        input_img, c2_lst_correct=False, mask_ocean_flag=False,
        latitude=36, longitude=-120, landcover_source=81, elevation_source=100,
        temperature_source=300, humidity_source=0.007, windspeed_source=4,
        solar_radiation_source=350,
        # cloudmask_args={'cloud_score_flag': False, 'filter_flag': False}
    ).lst)
    assert abs(output['lst'] - 300) <= 0.1


# def test_Image_from_landsat_c2_sr_cloud_mask_args():
#     """Test if the cloud_mask_args parameter can be set (not if it works)"""
#     output = model.Image.from_landsat_c2_sr(
#         'LANDSAT/LC08/C02/T1_L2/LC08_044033_20170716',
#         cloudmask_args={'snow_flag': True, 'cirrus_flag': True}
#     )
#     assert type(output) == type(default_image_obj())
#
#
# def test_Image_from_landsat_c2_sr_cloud_score_mask_arg():
#     """Test if the cloud_score_flag parameter can be set in cloudmask_args"""
#     output = model.Image.from_landsat_c2_sr(
#         'LANDSAT/LC08/C02/T1_L2/LC08_044033_20170716',
#         cloudmask_args={'cloud_score_flag': True, 'cloud_score_pct': 50})
#     assert type(output) == type(default_image_obj())


def test_Image_from_landsat_c2_sr_c2_lst_correct_arg():
    """Test if the c2_lst_correct parameter can be set (not if it works)"""
    output = model.Image.from_landsat_c2_sr(
        'LANDSAT/LC08/C02/T1_L2/LC08_031034_20160702', c2_lst_correct=True)
    assert type(output) == type(default_image_obj())


def test_Image_from_landsat_c2_sr_c2_lst_correct_fill():
    """Test if the c2_lst_correct fills the LST holes in Nebraska"""
    image_id = 'LANDSAT/LC08/C02/T1_L2/LC08_031034_20160702'
    xy = [-102.08284, 37.81728]

    # CGM - Is the uncorrected test needed?
    uncorrected = utils.point_image_value(
        model.Image.from_landsat_c2_sr(image_id, c2_lst_correct=False).lst, xy
    )
    assert uncorrected['lst'] is None
    corrected = utils.point_image_value(
        model.Image.from_landsat_c2_sr(image_id, c2_lst_correct=True).lst, xy
    )
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
    output = utils.getinfo(model.Image.from_image_id(image_id).lai)
    assert output['properties']['system:index'] == image_id.split('/')[-1]
    assert output['properties']['image_id'] == image_id


@pytest.mark.parametrize(
    'image_id, xy, expected',
    [
        [
            'LANDSAT/LC08/C02/T1_L2/LC08_044033_20170716', [-121.668, 38.905],
            {'lai': 3.5250, 'lst': 304.2305, 'albedo': 0.1514, 'emissivity': 0.9713}
        ],
        [
            'LANDSAT/LC08/C02/T1_L2/LC08_044033_20170716', [-121.650, 38.913],
            {'lai': 0.1052, 'lst': 326.5032, 'albedo': 0.1773, 'emissivity': 0.9752}
        ],
        # [ 'LANDSAT/LC09/C02/T1_L2/LC09_044033_20220127', ]
    ]
)
def test_Image_from_landsat_c2_sr_image_values(image_id, xy, expected, tol=0.0001):
    """Test instantiating the class using the from_image_id method"""
    lai = utils.point_image_value(model.Image.from_landsat_c2_sr(image_id).lai, xy, 30)
    assert abs(lai['lai'] - expected['lai']) < tol

    lst = utils.point_image_value(model.Image.from_landsat_c2_sr(image_id).lst, xy, 30)
    assert abs(lst['lst'] - expected['lst']) < tol

    albedo = utils.point_image_value(model.Image.from_landsat_c2_sr(image_id).albedo, xy, 30)
    assert abs(albedo['albedo'] - expected['albedo']) < tol

    emissivity = utils.point_image_value(model.Image.from_landsat_c2_sr(image_id).emissivity, xy, 30)
    assert abs(emissivity['emissivity'] - expected['emissivity']) < tol

    # print(lai, lst, albedo, emissivity)

    # landcover = utils.point_image_value(model.Image.from_landsat_c2_sr(image_id).landcover, xy, 30)
    # assert abs(landcover['landcover'] - 82) < tol
    # elevation = utils.point_image_value(model.Image.from_landsat_c2_sr(image_id).elevation, xy, 30)
    # assert abs(elevation['elevation'] - 10) < tol

