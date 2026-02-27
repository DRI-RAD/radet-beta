from datetime import datetime

import ee
import pytest

import radet.model as model
import radet.utils as utils
# TODO: import utils from openet.core
# import openet.core.utils as utils


# TODO: Try moving to conftest and/or make a fixture
SCENE_ID = 'LC08_044033_20170716'
# SCENE_ID = 'LC08_042035_20150713'
SCENE_DT = datetime.strptime(SCENE_ID[-8:], '%Y%m%d')
SCENE_DATE = SCENE_DT.strftime('%Y-%m-%d')
SCENE_TIME = utils.millis(SCENE_DT)



# def test_transmissivities():
#     lai =
#     output = model.transmissivities(lai)
#     # tauL, tauS, fc
#     assert False
#
#
# def test_wet_mask():
#     time_start =
#     lai =
#     output = model.wet_mask(time_start, lai)
#     # wet.Or(wet_ww2), water
#     assert False
#
#
# def test_clear_sky_terms():
#     # time_start, lat, elev, srad
#     output = model.clear_sky_terms(time_start, lat, elev, srad)
#     # Rs_MJ, Ra_MJ, fcd, sunrise_ts
#     assert False
#
#
# def test_daily_avg_lst():
#     # tminK, lst, sunrise_ts, t_avg
#     output = model.daily_avg_lst(tminK, lst, sunrise_ts, t_avg)
#     # LST_max.add(LST_min).multiply(0.5).max(t_avg)
#     assert False
#
#
# def test_canopy_and_soil_LST():
#     # LST_avg, t_avg, Rs_MJ, Rld_MJ, fc, tauS, tauL, mu_c, mu_s, RHs, DELTA, gamma, emissivity, albedo
#     output = model.canopy_and_soil_LST(
#         LST_avg, t_avg, Rs_MJ, Rld_MJ, fc, tauS, tauL, mu_c, mu_s, RHs, DELTA, gamma, emissivity, albedo
#     )
#     # LST_canopy, LST_soil
#     assert False
#
#
# def test_RHs_model():
#     # ea, esat, DELTA, LST_soil, t_avg, mu_s, water
#     output = model.RHs_model(ea, esat, DELTA, LST_soil, t_avg, mu_s, water)
#     # RHs.where(water, 1).rename("RHs")
#     assert False
#
#
# def test_net_radiation():
#     # emissivity, LST_canopy, LST_soil, Rs_MJ, Rld_MJ, albedo, tauS, tauL
#     output = model.net_radiation(emissivity, LST_canopy, LST_soil, Rs_MJ, Rld_MJ, albedo, tauS, tauL)
#     # Rnc.add(Rns), Rnc, Rns, G, Rns.subtract(G)
#     assert False
#
#
# def test_isothermal_net_radiation():
#     # Rnc, AEs, tauL, emissivity, t_avg, gg, LST_canopy, LST_soil
#     output = model.isothermal_net_radiation(Rnc, AEs, tauL, emissivity, t_avg, gg, LST_canopy, LST_soil)
#     # Rnci, AEsi
#     assert False
#
#
# def test_mu_terms():
#     # Rnc, Rnci, DELTA, gamma, AEs, AEsi, RHs
#     output = model.mu_terms(Rnc, Rnci, DELTA, gamma, AEs, AEsi, RHs)
#     #  mu_c, mu_s
#     assert False


@pytest.mark.parametrize(
    'gamma, DELTA, Rnc, AEs, RHs, mu_c, mu_s, expected',
    [
        # Test values pulled from an arbitrary image
        [0.06122, 0.12630, 4.28263, 5.31659, 0.21313, 1.11053, 4.64080, 1.32411],
    ]
)
def test_DIF_model(gamma, DELTA, Rnc, AEs, RHs, mu_c, mu_s, expected, tol=0.0001):
    # gamma, DELTA, Rnc, AEs, RHs, mu_c, mu_s
    output = utils.constant_image_value(model.DIF_model(
        ee.Image.constant(gamma),
        ee.Image.constant(DELTA),
        ee.Image.constant(Rnc),
        ee.Image.constant(AEs),
        ee.Image.constant(RHs),
        ee.Image.constant(mu_c),
        ee.Image.constant(mu_s)
    ).rename('constant'))
    assert abs(output['constant'] - expected) <= tol


# def test_aerodynamic_term():
#     # del_LC, fc, LST_soil, RHs, gamma, DELTA, u2, esat, ea
#     output = model.aerodynamic_term(del_LC, fc, LST_soil, RHs, gamma, DELTA, u2, esat, ea)
#     assert False
#
#
# def test_Rld_atm_ASCE():
#     # emissivity, fcd, ea, t_avg
#     output = model.Rld_atm_ASCE(emissivity, fcd, ea, t_avg)
#     assert False
#
#
# def test_add_lapse_correction():
#     # ta, elev, elev_gridmet
#     output = model.add_lapse_correction(ta, elev, elev_gridmet)
#     # t_corrected
#     assert False
#
#
# def test_terrain_shade_correct_srad():
#     # Rs_MJ, Ra_MJ, elev, time_start, albedo
#     output = model.terrain_shade_correct_srad()
#     assert False
#
#
# def test_et():
#     output = model.et(
#         # image,
#         lai,
#         lst,
#         albedo,
#         emissivity,
#         meteorology_source,
#         landcover,
#         elevation,
#         time_start,
#         proj,
#         latitude,
#         # geometry_image,
#         # coords
#     )
#     assert False


# TODO: Check that a ValueError is raised for invalid daily meteorology sources
# def test_Image_meteorology_source_sources_exception():
#     with pytest.raises(ValueError):
#         utils.getinfo(model.et(meteorology_source='').et)
