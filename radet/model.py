import logging
import math

import ee

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

pi = math.pi
SIG = 4.901e-9


def et(
    image,
    lai,
    lst,
    albedo,
    emissivity,
    meteorology_source_daily,
    # elev_product,
    time_start,
    geometry_image,
    proj,
    coords,
):
    """
    Daily Evapotranspiration [mm day-1].

    Parameters
    ----------
    image : ee.Image
        Landsat image.
    lai : ee.Image
        Leaf Area index.
    lst : ee.Image
        Land Surface Temperature [K].
    albedo : ee.Image
        Surface albedo.
    emissivity : ee.Image
        Broad-band surface emissivity.
    meteorology_source_inst : ee.ImageCollection
        Meteorological dataset [inst]
    meteorology_source_daily : ee.ImageCollection
        Meteorological dataset [daily]
    elev_product : ee.Image
        Elevation image [m].
    time_start : str
        Image property: time start of the image.
    geometry_image : ee.Geometry
        Image geometry.
    proj : ee.Image
        Landsat image projection.
    coords : ee.Image
        Landsat image Latitude and longitude.

    Returns
    -------
    ee.Image

    Reference
    ----------
    TODO: Update reference
    [Kim2026] Yeonuk Kim, Justin L Huntington, Bruno Comini de Andrade,
        John M Volk, Sayantan Majumdar, Charles Morton, ...

    """

    # TODO: check if time_start local time affects calculations. (BC: Don't think so)
    # TODO: check with CM if ancillary elev and lat images are really needed.

    # Meteorological data
    # Only accepts GridMET for now
    # TODO: add ERA5-Land compatibility
    if meteorology_source_daily == "IDAHO_EPSCOR/GRIDMET":
        srad, tminK, tmaxK, qa, u10 = meteorology_gridmet(
            time_start,
            meteorology_source_daily,
        )
        lat = ee.Image("projects/openet/assets/meteorology/gridmet/ancillary/latitude")
        elev = elev = ee.Image("projects/openet/assets/meteorology/gridmet/ancillary/elevation")

    elif (meteorology_source_daily == "projects/openet/assets/meteorology/era5land/na/daily") or (
        meteorology_source_daily == "projects/openet/assets/meteorology/era5land/sa/daily"
    ):

        tmin, tmax, tair, ux, rh, rso24h, tfac = meteorology_era5land(
            time_start,
            meteorology_source_daily,
        )

    else:
        raise Exception("Error: wrong daily or instant met data source assigned.")

    #################################
    # Meteorological variables
    #################################
    # Elevation image
    elev = ee.Image("USGS/SRTMGL1_003")

    # Air pressure (kPa)
    PA = ee.Image().expression("101.3 * ((293 - 0.0065*elev)/293)**5.26", {"elev": elev})

    # Psychrometric constant
    gamma = ee.Image().expression("PA * 0.000665", {"PA": PA})

    # 2m wind speed
    u2 = ee.Image().expression("u10 * 4.87 / log(67.8*10 - 5.42)", {"u10": u10})

    # Vapor pressure (kPa)
    ea = ee.Image().expression("qa * PA/((1 - 0.622) * qa + 0.622)", {"PA": PA, "qa": qa})

    ###################################################
    # temperature variables after lapse rate correction
    ##################################################
    # GridMET Elevation image
    elev_gridmet = ee.Image("projects/openet/assets/meteorology/gridmet/ancillary/elevation")
    
    tmaxK_cor = add_lapse_correction(tmaxK, elev, elev_gridmet)
    tmninK_cor = add_lapse_correction(tminK, elev, elev_gridmet)

    # average air temperature (K)
    t_avg = tmninK_cor.add(tmaxK_cor).multiply(0.5)

    # Saturation vapor pressure (kPa)
    esat = ee.Image().expression("0.6108 * exp(17.27*(t-273.15)/((t-273.15)+237.3))", {"t": t_avg})

    # Slope of vapor pressure curve
    DELTA = ee.Image().expression(
        "2503*exp(17.27*(t-273.15)/((t-273.15)+237.3))/((t-273.15)+237.3)**2", {"t": t_avg})


    ############################################
    # Clear Sky terms and slope shade correction
    ############################################
    # Clear-Sky terms
    Rs_MJ, Ra_MJ, fcd, sunrise_ts = clear_sky_terms(time_start, lat, elev, srad)

    # surface shortwave correction based on shade effect (Allen et al., 2006)
    Rs_MJ_cor = terrain_shade_correct_srad(Rs_MJ, Ra_MJ, elev, time_start, albedo)
    
    ############################################
    # Other variables indepedent to mu terms
    ############################################
    # land cover mask
    del_LC, water = wet_mask(time_start,lai)

    # transmissivity 
    tauL, tauS, fc = transmissivities(lai)
    
    # donward long wave from atmosphere
    Rld_MJ = Rld_atm_ASCE(emissivity, fcd, ea, t_avg)

    # Daily mean LST 
    LST_avg = daily_avg_lst(tmninK_cor, lst, sunrise_ts, t_avg)

    # soil conductive exchange coefficient
    gg = ee.Image().expression('1000*((pi/86400)**0.5)*86400/(10**6)',{"pi":pi})

    ##########################################
    ### initial guess: mu_s = 1, mu_c = 1
    ##########################################    
    mu_c = ee.Number(1)
    mu_s = ee.Number(1)

    # RHs (mu_s=1 leads to RHs = RHa)
    RHs = ea.divide(esat)
        
    # soil and canopy LST
    LST_canopy, LST_soil= canopy_and_soil_LST(LST_avg, t_avg, Rs_MJ_cor, Rld_MJ, fc, tauS, tauL, mu_c, mu_s, RHs, DELTA, gamma, emissivity, albedo)

    # Daily mean net radiation (total, soil and canopy, soil heat flux) (MJ m-2 d-1)
    Rn, Rnc, Rns, G, AEs = net_radiation(emissivity, LST_canopy, LST_soil, Rs_MJ_cor, Rld_MJ, albedo, tauS, tauL)

    # Isothermal canopy net radiation and isothermal soil available energy
    Rnci, AEsi = isothermal_net_radiation(Rnc, AEs, tauL, emissivity, t_avg, gg, LST_canopy, LST_soil)

    # mu terms
    mu_c, mu_s = mu_terms(Rnc, Rnci, DELTA, gamma, AEs, AEsi, RHs)

    ##########################################
    ### Second run
    ##########################################    
    # RHs 
    RHs = RHs_model(ea, esat, DELTA, LST_soil, t_avg, mu_s, water)
    
    # soil and canopy LST
    LST_canopy, LST_soil= canopy_and_soil_LST(LST_avg, t_avg, Rs_MJ_cor, Rld_MJ, fc, tauS, tauL, mu_c, mu_s, RHs, DELTA, gamma, emissivity, albedo)

    # Daily mean net radiation (total, soil and canopy, soil heat flux) (MJ m-2 d-1)
    Rn, Rnc, Rns, G, AEs = net_radiation(emissivity, LST_canopy, LST_soil, Rs_MJ_cor, Rld_MJ, albedo, tauS, tauL)

    # Isothermal canopy net radiation and isothermal soil available energy
    Rnci, AEsi = isothermal_net_radiation(Rnc, AEs, tauL, emissivity, t_avg, gg, LST_canopy, LST_soil)

    # mu terms
    mu_c, mu_s = mu_terms(Rnc, Rnci, DELTA, gamma, AEs, AEsi, RHs)

    ##########################################
    # ET model (mm/day)
    ##########################################
    # ET_DIF (mm/day)
    ET_DIF = DIF_model(gamma, DELTA, Rnc, AEs, RHs, mu_c, mu_s)

    # Penman aerodynamic term
    Ea = aerodynamic_term(del_LC, fc, LST_soil, RHs, gamma, DELTA, u2, esat, ea)

    return ee.Image(ET_DIF.add(Ea).rename("ET")).setDefaultProjection(proj)


def transmissivities(lai):

    tauL = lai.expression(
              'exp(-0.95*lai)',
              {'lai': lai}
            ).rename('tauL').clamp(0.01,1)
    
    tauS = lai.expression(
              'exp(-0.56*lai)',
              {'lai': lai}
            ).rename('tauS').clamp(0.01,1)
    
    fc = lai.expression(
              '1-exp(-0.4*lai)',
              {'lai': lai}
            ).rename('fc').clamp(0.01,1)
    
    return tauL, tauS, fc

    
def wet_mask(time_start, lai):
    # USGS NLCD land cover 
    # Load land cover image

    YEAR = ee.Number(ee.Date(time_start).get('year'))
    landcover = ee.ImageCollection("projects/sat-io/open-datasets/USGS/ANNUAL_NLCD/LANDCOVER")
    landcover_img = landcover.filter(ee.Filter.eq("year", YEAR)).first()
    if landcover_img.getInfo() is None:
        landcover_img = landcover.sort("year", False).first()
        land_cover_year = landcover_img.get("year").getInfo()
        logger.info(f"No NLCD data for year {YEAR.getInfo()}. Using closest available year: {land_cover_year}.")
    water = landcover_img.remap([11], [1], 0).eq(1)
    # Define wet surfaces
    wet = landcover_img.remap([11, 81, 82, 95], [1, 1, 1, 1], 0).eq(1)
    wet_ww = landcover_img.remap([90], [1], 0).eq(1)
    wet_ww2 = wet_ww.And(lai.lt(1))
    return wet.Or(wet_ww2), water


def clear_sky_terms(time_start, lat, elev, srad):
    J = ee.Number(ee.Date(time_start).getRelative("day", "year")).add(1)
    phi = lat.multiply(pi / 180)

    dr = ee.Image().expression("1 + 0.033 * cos(2*PI/365 * J)", {"PI": pi, "J": J})
    delta = ee.Image().expression("0.409 * sin(2*PI/365 * J - 1.39)", {"PI": pi, "J": J})
    cosW = ee.Image().expression("-tan(phi) * tan(delta)", {"phi": phi, "delta": delta})
    omega = cosW.clamp(-1, 1).acos()

    Ra_MJ = ee.Image().expression(
        "(24/PI) * 4.92 * dr * (omega*sin(phi)*sin(delta) + cos(phi)*cos(delta)*sin(omega))",
        {"PI": pi, "dr": dr, "omega": omega, "phi": phi, "delta": delta},
    )

    Rso_MJ = ee.Image().expression("(0.75 + 2e-5*elev) * Ra", {"elev": elev, "Ra": Ra_MJ})

    Rs_MJ = srad.multiply(86400).divide(1e6)
    ratio = Rs_MJ.divide(Rso_MJ.where(Rso_MJ.eq(0), 1)).clamp(0.3, 1)
    fcd = ee.Image(1.35).multiply(ratio).subtract(0.35)
    day_length = ee.Image().expression("24/PI*omega", {"PI": pi, "omega": omega})
    sunrise_ts = ee.Image().expression("12 - day_length/2", {"day_length": day_length})

    return Rs_MJ, Ra_MJ, fcd, sunrise_ts


def daily_avg_lst(tminK, lst, sunrise_ts, t_avg):
    LST_min = tminK.subtract(1)
    LST_max = lst.expression(
        "LST_min + (LST - LST_min)/(cos((10-12.5)/(12.5 - sunrise_ts)/2*PI))",
        {"PI": pi, "LST_min": LST_min, "LST": lst, "sunrise_ts": sunrise_ts},
    )
    return LST_max.add(LST_min).multiply(0.5).max(t_avg)


    
def  canopy_and_soil_LST(LST_avg, t_avg, Rs_MJ,Rld_MJ, fc, tauS, tauL, mu_c, mu_s, RHs, DELTA, gamma, emissivity, albedo):
    
    ## canopy LST
    beta = ee.Image().expression('fc/(fc + mu_s/mu_c*(1-fc)*(DELTA+gamma*mu_c)/(DELTA*RHs+gamma*mu_s))',
          {"fc":fc, "DELTA":DELTA, "gamma":gamma, "RHs": RHs, "mu_c": mu_c, "mu_s": mu_s})
    
    LST_canopy = ee.Image().expression('t_avg + beta*(LST_avg-t_avg)',
              {"t_avg":t_avg, "LST_avg":LST_avg, "beta":beta});
    
    ## maximum soil LST
    LST_soil_max = ee.Image().expression('((tauS*(1-albedo)*Rs_MJ + tauL*Rld + (1-tauL)*emis*SIG *LST_canopy**4)/emis/SIG)**(1/4)',
             {"tauS":tauS, "tauL":tauL, "Rs_MJ":Rs_MJ, "albedo":albedo, "Rld":Rld_MJ, "LST_canopy":LST_canopy,"emis":emissivity, "SIG":SIG})
    
    ## soil LST
    LST_soil = ee.Image().expression('(((LST_avg)**4 - (1-tauL)*(LST_canopy)**4)/tauL)**(1/4)',
              {"LST_canopy":LST_canopy, "LST_avg":LST_avg, "tauL":tauL});
    LST_soil = LST_soil.min(LST_soil_max)    
    
    return LST_canopy, LST_soil

def RHs_model(ea, esat, DELTA, LST_soil, t_avg, mu_s, water):

    RHs = ee.Image().expression('ea/(esat+DELTA*(LST_soil - ta)*(mu_s-1)/mu_s)',
          {"ea":ea, "esat":esat, "mu_s":mu_s, "DELTA":DELTA, "LST_soil":LST_soil, "ta": t_avg})
    
    return RHs.where(water,1).rename("RHs")
    
def net_radiation(emissivity, LST_canopy, LST_soil, Rs_MJ, Rld_MJ, albedo, tauS, tauL):

    ## net radiation at soil 
    Rns  = ee.Image().expression(
        'tauS*(1-albedo)*Rs_MJ + tauL*Rld + (1-tauL)*emis*SIG *LST_canopy**4 - emis*SIG *LST_soil**4', 
        {"tauS":tauS, "tauL":tauL, "Rs_MJ":Rs_MJ, "albedo":albedo, "Rld":Rld_MJ, 
         "emis":emissivity, "SIG":SIG, "LST_canopy":LST_canopy, "LST_soil":LST_soil}).rename('Rns').max(0)
    ## net radiation at canopy 
    Rnc  = ee.Image().expression(
        '(1-tauS)*(1-albedo)*Rs_MJ + (1-tauL)*(Rld + emis*SIG *LST_soil**4 - 2*emis*SIG *LST_canopy**4)', 
        {"tauS":tauS, "tauL":tauL, "Rs_MJ":Rs_MJ, "albedo":albedo, "Rld":Rld_MJ, 
         "emis":emissivity, "SIG":SIG, "LST_canopy":LST_canopy, "LST_soil":LST_soil}).rename('Rnc').max(0)
    ## soil heat flux  
    G = ee.Image().expression('0.35*Rns - 1.5', {"Rns":Rns}).rename('G')
    
    return Rnc.add(Rns), Rnc, Rns, G, Rns.subtract(G)


def isothermal_net_radiation(Rnc, AEs, tauL, emissivity, t_avg, gg, LST_canopy, LST_soil):
    # Isothermal canopy net radiation and isothermal soil available energy
    Rnci  = ee.Image().expression(
        'Rnc + 4*2*(1-tauL)*emis*SIG*ta**3*(LST-ta)', 
        {"Rnc":Rnc, "tauL":tauL, "emis":emissivity, "SIG":SIG, "ta":t_avg, "LST": LST_canopy}).rename('Rnci')

    AEsi  = ee.Image().expression(
        'AEs + (4*emis*SIG*ta**3 + gg)*(LST-ta)', 
        {"AEs":AEs, "emis":emissivity, "SIG":SIG, "ta":t_avg, "LST": LST_soil, "gg":gg }).rename('AEsi')

    return Rnci, AEsi


def mu_terms(Rnc, Rnci, DELTA, gamma, AEs, AEsi, RHs):
    
    mu_c = ee.Image().expression(
        '(Rnci + sqrt(Rnci**2 + 4*DELTA/gamma*Rnc*(Rnci-Rnc)))/(2*Rnc)', 
        {"Rnc": Rnc, "Rnci": Rnci, "DELTA": DELTA, "gamma": gamma}).max(1).rename('mu_c')

    mu_s = ee.Image().expression(
        '(AEsi + sqrt(AEsi**2 + 4*RHs*DELTA/gamma*AEs*(AEsi-AEs)))/(2*AEs)', 
        {"AEs": AEs, "AEsi": AEsi, "DELTA": DELTA, "gamma": gamma, "RHs":RHs}).max(1).rename('mu_s')

    return mu_c, mu_s


def DIF_model(gamma, DELTA, Rnc, AEs, RHs, mu_c, mu_s):
    
    return (
        ee.Image().expression(
            "(DELTA/(DELTA+gamma*mu_c)*Rnc + RHs*DELTA/(RHs*DELTA+gamma*mu_s)*AEs)/2.45",
            {"gamma": gamma, "DELTA": DELTA, "Rnc": Rnc, "AEs": AEs, "RHs": RHs, "mu_c": mu_c, "mu_s": mu_s},
        )
        .rename("ET_DIF")
        .max(0)
    )


def aerodynamic_term(del_LC, fc, LST_soil, RHs, gamma, DELTA, u2, esat, ea):

    esat_LST_soil = ee.Image().expression("0.6108 * exp(17.27*(t-273.15)/((t-273.15)+237.3))", {"t": LST_soil})

    del_water = ee.Image().expression(
        'fc + (1-fc)*(RHs**(esat_soil - esat_soil*RHs))/(1+exp(10-LST_soil))',
        {"fc":fc, "RHs":RHs, "esat_soil":esat_LST_soil, "LST_soil":LST_soil})
    
    Ea = ee.Image().expression(
        'del_LC*del_water*gamma/(DELTA + gamma)*(esat-ea)*2.6*(1+0.54*u2)', 
        {"gamma":gamma, "DELTA":DELTA, "u2": u2, "esat": esat, "ea":  ea, "del_water":del_water, "del_LC":del_LC}).rename('Ea')

    return Ea


def Rld_atm_ASCE(emissivity, fcd, ea, t_avg):

    Rld = ee.Image().expression('emissivity*SIG * (1 - fcd*(0.34 - 0.14 * ea ** 0.5)) * t_avg **4',
              {"emissivity":emissivity, "SIG":SIG, "fcd":fcd, "ea":ea, "t_avg":t_avg}).rename("Rld") 
    return Rld

def add_lapse_correction(ta, elev, elev_gridmet):
    zGrid = elev_gridmet
    zPix  = elev

    lapse_corr = zPix.subtract(zGrid).multiply(-0.0065)
    t_corrected = ta.add(lapse_corr).rename('t_corrected')

    return t_corrected

def meteorology_gridmet(time_start, meteorology_source_daily):
    """
    Parameters
    ----------
    time_start : str
        Image property: time start of the image.
    meteorology_source_inst: ee.ImageCollection, str
        Instantaneous meteorological data.
    meteorology_source_daily :  ee.ImageCollection, str
        Daily meteorological data.

    Returns
    -------
    ee.Image

    Notes
    -----
    Accepted collections:
    Daily : IDAHO_EPSCOR/GRIDMET

    References
    ----------

    """
    # Get date information
    time_start = ee.Number(time_start)

    # Filtering Daily data
    meteorology_daily = (
        ee.ImageCollection(meteorology_source_daily)
        .filterDate(ee.Date(time_start).advance(-1, "day"), ee.Date(time_start))
        .first()
    )

    srad = meteorology_daily.select("srad")
    tmaxK = meteorology_daily.select("tmmx")
    tminK = meteorology_daily.select("tmmn")
    qa = meteorology_daily.select("sph")
    u10 = meteorology_daily.select("vs")
    
    return srad, tminK, tmaxK, qa, u10

def meteorology_era5land(time_start, meteorology_source_inst, meteorology_source_daily):
    """
    Parameters
    ----------
    time_start : str
        Image property: time start of the image.
    meteorology_source_inst: ee.ImageCollection, str
        Instantaneous meteorological data.
    meteorology_source_daily :  ee.ImageCollection, str
        Daily meteorological data.

    Returns
    -------
    ee.Image

    Notes
    -----
    Accepted collections:
    Inst : ECMWF/ERA5_LAND/HOURLY
    Daily : projects/openet/assets/meteorology/era5land/na/daily
            projects/openet/assets/meteorology/era5land/sa/daily

    References
    ----------

    """

    # Get date information
    time_start = ee.Number(time_start)

    # Filtering Daily data
    meteorology_daily = (
        ee.ImageCollection(meteorology_source_daily)
        .filterDate(ee.Date(time_start).advance(-1, "day"), ee.Date(time_start).advance(1, "day"))
        .first()
    )

    # Instantaneous data
    meteorology_inst_collection = ee.ImageCollection(meteorology_source_inst)

    # Linear interpolation
    previous_time = time_start.subtract(1 * 60 * 60 * 1000)
    next_time = time_start.add(1 * 60 * 60 * 1000)

    previous_image = (
        meteorology_inst_collection.filterDate(previous_time, time_start)
        .limit(1, "system:time_start", False)
        .first()
    )

    next_image = (
        meteorology_inst_collection.filterDate(time_start, next_time)
        .limit(1, "system:time_start", True)
        .first()
    )

    image_previous_time = ee.Number(previous_image.get("system:time_start"))
    image_next_time = ee.Number(next_image.get("system:time_start"))

    delta_time = time_start.subtract(image_previous_time).divide(
        image_next_time.subtract(image_previous_time)
    )

    # Incoming shorwave down [W m-2]
    swdown24h = meteorology_daily.select("surface_solar_radiation_downwards").divide(1 * 60 * 60 * 24)

    # Minimum air temperature [K]
    tmin = meteorology_daily.select("temperature_2m_min").rename("tmin")

    # Maximum air temperature [K]
    tmax = meteorology_daily.select("temperature_2m_max").rename("tmax")

    # Instantaneous incoming shortwave radiation [W m-2]
    rso_inst = (
        ee.ImageCollection(meteorology_source_inst)
        .filterDate(ee.Date(time_start), ee.Date(time_start).advance(1, "hour"))
        .select("surface_solar_radiation_downwards_hourly")
        .mean()
        .divide(1 * 60 * 60)
        .rename("rso_inst")
    )

    # Air temperature [C]
    # TODO: LL- Change all temperatures to K ?
    tair_c = (
        next_image.select("temperature_2m")
        .subtract(previous_image.select("temperature_2m"))
        .multiply(delta_time)
        .add(previous_image.select("temperature_2m"))
        .subtract(273.15)
        .rename("tair")
    )

    # Wind speed [ m/s]
    wind_u = (
        next_image.select("u_component_of_wind_10m")
        .subtract(previous_image.select("u_component_of_wind_10m"))
        .multiply(delta_time)
        .add(previous_image.select("u_component_of_wind_10m"))
    )

    wind_v = (
        next_image.select("v_component_of_wind_10m")
        .subtract(previous_image.select("v_component_of_wind_10m"))
        .multiply(delta_time)
        .add(previous_image.select("v_component_of_wind_10m"))
    )

    wind_med = wind_u.expression(
        "sqrt(ux_u ** 2 + ux_v ** 2)",
        {"ux_u": wind_u, "ux_v": wind_v},
    ).rename("ux")

    wind_med = wind_med.expression("ux * (4.87) / log(67.8 * z - 5.42)", {"ux": wind_med, "z": 10.0}).rename(
        "ux"
    )

    # Dew point temperature [°K]
    tdp = (
        next_image.select("dewpoint_temperature_2m")
        .subtract(previous_image.select("dewpoint_temperature_2m"))
        .multiply(delta_time)
        .add(previous_image.select("dewpoint_temperature_2m"))
        .rename("tdp")
    )

    # Actual vapour pressure [kPa]
    ea = tdp.expression("0.6108 * (exp((17.27 * T_air) / (T_air + 237.3)))", {"T_air": tdp.subtract(273.15)})

    # SATURATED VAPOR PRESSURE [kPa]
    esat = tair_c.expression("0.6108 * (exp((17.27 * T_air) / (T_air + 237.3)))", {"T_air": tair_c})

    # RELATIVE HUMIDITY (%)
    rh = ea.divide(esat).multiply(100).rename("RH")

    # Surface temperature correction based on precipitation and reference ET

    # Accumulation time period
    accum_period = -60

    # Accum meteo data
    gridmet_accum = ee.ImageCollection(meteorology_source_daily).filterDate(
        ee.Date(time_start).advance(accum_period, "days"), ee.Date(time_start)
    )

    # Reference ET
    etr_accum = gridmet_accum.select("etr_asce").sum()

    # Precipitation
    precipt_accum = gridmet_accum.select("total_precipitation").sum()

    # Ratio between precipt/etr
    ratio = precipt_accum.divide(etr_accum)

    # Temperature adjustment offset (Allen2013 Eqn 8)

    tfac = etr_accum.expression("2.6 - 13 * ratio", {"ratio": ratio})

    tfac = ee.Image(tfac.where(ratio.gt(0.2), 0)).rename("tfac")

    # Resample
    tmin = tmin.subtract(273.15).resample("bilinear")
    tmax = tmax.subtract(273.15).resample("bilinear")
    rso_inst = rso_inst.resample("bilinear")
    tair_c = tair_c.resample("bilinear")
    wind_med = wind_med.resample("bilinear")
    rh = rh.resample("bilinear")
    swdown24h = swdown24h.resample("bilinear")

    return [tmin, tmax, tair_c, wind_med, rh, rso_inst, swdown24h, tfac]


    #############################################################################
    ############ TERRAIN SHADE CORRECTION of SRAD (based on Allen et al., 2006)
    #############################################################################
    
def terrain_shade_correct_srad(Rs_MJ, Ra_MJ, elev, time_start, albedo):
    """
    Terrain shade correction of daily solar radiation (Allen et al., 2006-based).
    
    Args:
        Rs_MJ (ee.Image): incoming solar radiation (MJ m-2 d-1)
        Ra_MJ (ee.Image): extraterrestrial radiation (MJ m-2 d-1)
        elev  (ee.Image): elevation DEM (m)
        time_start 
        albedo (float): surface albedo used in Eq.38 term

    Returns: 'Rs_MJ_corr'   (ee.Image)  # band name 'Rs_MJ_shade'
          
    """

    # 0) Julian date
    J = ee.Number(ee.Date(time_start).getRelative("day", "year")).add(1)

    # 1) tau = Rs / Ra (Eq.39)
    tau = Rs_MJ.divide(Ra_MJ.where(Ra_MJ.eq(0), 1)).clamp(0, 1)

    # 2) KB, KD (Eq.41)
    KB1 = ee.Image().expression("1.56*tau - 0.55", {"tau": tau})
    KB2 = ee.Image().expression("0.022 - 0.28*tau + 0.828*tau**2 + 0.765*tau**3", {"tau": tau})
    KB3 = ee.Image().expression("0.016*tau", {"tau": tau})

    KB = KB1.where(tau.lt(0.42), KB2).where(tau.lt(0.175), KB3)
    KD = tau.subtract(KB)

    # 3) Integral of cosine theta (Eq.5 + Appendix A)
    def cosThetaIntegral(slope, aspect, doy, lat, lon):
        slope  = ee.Image(slope)
        aspect = ee.Image(aspect)
        lat    = ee.Image(lat).multiply(pi/180.0)  # rad
        doy    = ee.Image(doy)

        delta = ee.Image().expression(
            "0.409 * sin(2 * pi / 365 * J - 1.39)",
            {"pi": pi, "J": doy}
        )

        sin_d   = delta.sin()
        cos_d   = delta.cos()
        sin_phi = lat.sin()
        cos_phi = lat.cos()
        sin_s   = slope.sin()
        cos_s   = slope.cos()
        cos_g   = aspect.cos()
        sin_g   = aspect.sin()

        polarDay_NH   = delta.add(lat).gt(pi/2)
        polarDay_SH   = delta.add(lat).lt(-pi/2)
        condPolarDay  = polarDay_NH.Or(polarDay_SH)

        polarNight_NH  = delta.subtract(lat).gt(pi/2)
        polarNight_SH  = delta.subtract(lat).lt(-pi/2)
        condPolarNight = polarNight_NH.Or(polarNight_SH)

        v124_polarDay   = ee.Image(-pi)
        v224_polarDay   = ee.Image( pi)
        v124_polarNight = ee.Image(0)
        v224_polarNight = ee.Image(0)

        a = ee.Image().expression(
            "sd*cp*ss*cg - sd*sp*cs",
            {"sd": sin_d, "sp": sin_phi, "cp": cos_phi, "ss": sin_s, "cs": cos_s, "cg": cos_g}
        )
        b = ee.Image().expression(
            "cd*cp*cs + cd*sp*ss*cg",
            {"cd": cos_d, "sp": sin_phi, "cp": cos_phi, "ss": sin_s, "cs": cos_s, "cg": cos_g}
        )
        c = ee.Image().expression(
            "cd * sg * ss",
            {"cd": cos_d, "sg": sin_g, "ss": sin_s}
        )

        vS = ee.Image().expression(
            "acos(-tan(phi) * tan(delta))",
            {"phi": lat, "delta": delta}
        )
        vS_neg = vS.multiply(-1)

        def cosu(v):
            return ee.Image().expression(
                "-a + b*cos(v) + c*sin(v)",
                {"a": a, "b": b, "c": c, "v": v}
            )

        cosu_sunrise = cosu(vS_neg)
        cosu_sunset  = cosu(vS)

        disc = b.pow(2).add(c.pow(2)).subtract(a.pow(2)).max(0.0001)
        sqrt_disc = disc.sqrt()

        sin_v1_cand = ee.Image().expression(
            "(a*c - b*sqrtD) / (b*b + c*c)",
            {"a": a, "b": b, "c": c, "sqrtD": sqrt_disc}
        ).clamp(-1, 1)

        sin_v2_cand = ee.Image().expression(
            "(a*c + b*sqrtD) / (b*b + c*c)",
            {"a": a, "b": b, "c": c, "sqrtD": sqrt_disc}
        ).clamp(-1, 1)

        v1_raw = sin_v1_cand.asin()
        v2_raw = sin_v2_cand.asin()

        v1_mirror = v1_raw.multiply(-1).add(-pi)
        v2_mirror = v2_raw.multiply(-1).add( pi)

        cosu_v1  = cosu(v1_raw)
        cosu_v1m = cosu(v1_mirror)
        cosu_v2  = cosu(v2_raw)
        cosu_v2m = cosu(v2_mirror)

        # Step B: sunrise on slope
        v124 = vS_neg
        condB1 = cosu_sunrise.lte(cosu_v1).And(cosu_v1.lt(0.001))
        v124 = v124.where(condB1, v1_raw)

        condB3 = condB1.Not().And(cosu_v1m.lte(0.001)).And(v1_mirror.gt(vS_neg))
        v124 = v124.where(condB3, v1_mirror)
        v124 = v124.where(v124.lt(vS_neg), vS_neg)

        # Step C: sunset on slope
        v224 = vS
        condC1 = cosu_sunset.lte(cosu_v2).And(cosu_v2.lt(0.001))
        v224 = v224.where(condC1, v2_raw)

        condC3 = condC1.Not().And(cosu_v2m.lte(0.001)).And(v2_mirror.lt(vS))
        v224 = v224.where(condC3, v2_mirror)

        v224 = v224.where(v224.gt(vS), vS)

        noBeam = v224.lt(v124)

        # Step D: two beam period logic
        doStepD = sin_s.gt(sin_phi.multiply(cos_d).add(cos_phi.multiply(sin_d)))

        A = v2_raw
        B = v1_raw
        v224b = A.min(B)
        v124b = A.max(B)

        cosu_v224b = cosu(v224b)
        cosu_v124b = cosu(v124b)

        cond224_adj = cosu_v224b.lt(-0.001).Or(cosu_v224b.gt(0.001))
        v224b = v224b.where(cond224_adj, ee.Image().expression("-pi - v", {"pi": pi, "v": v224b}))

        cond124_adj = cosu_v124b.lt(-0.001).Or(cosu_v124b.gt(0.001))
        v124b = v124b.where(cond124_adj, ee.Image().expression("pi - v", {"pi": pi, "v": v124b}))

        constraintsOk = v224b.gte(v124).And(v124b.lte(v224))
        doStepD = doStepD.And(constraintsOk).And(noBeam.Not())

        def I(v1, v2):
            return ee.Image().expression(
                " sd*sp*cs*(v2-v1)"
                "- sd*cp*ss*cg*(v2-v1)"
                "+ cd*cp*cs*( sin(v2)-sin(v1) )"
                "+ cd*sp*ss*cg*( sin(v2)-sin(v1) )"
                "- cd*ss*sg*( cos(v2)-cos(v1) )",
                {
                    "sd": sin_d, "cd": cos_d,
                    "sp": sin_phi, "cp": cos_phi,
                    "ss": sin_s, "cs": cos_s,
                    "cg": cos_g, "sg": sin_g,
                    "v1": v1, "v2": v2
                }
            )

        I_single = I(v124, v224).where(noBeam, 0)

        X = I(v224b, v124b)
        middayShade = X.lt(0)
        doStepD = doStepD.And(middayShade)

        I_total = I_single.where(doStepD, I(v124, v224b).add(I(v124b, v224)))

        return (ee.Image(0)
                .where(condPolarDay, I(v124_polarDay, v224_polarDay))
                .where(condPolarNight, 0)
                .where(condPolarDay.Not().And(condPolarNight.Not()), I_total)
                .max(0.0001))

    # 4) daily_ratio (Eq.38)
    def daily_ratio(J, KB, KD, slope, aspect, lat, lon, albedo):
        doy = J

        cos_theta_flat  = cosThetaIntegral(ee.Image(0), ee.Image(0), doy, lat, lon)
        cos_theta_slope = cosThetaIntegral(slope, aspect, doy, lat, lon)
        cos_theta_slope = cos_theta_slope.where(cos_theta_flat.lte(0.0001), 0.0001)

        fB = cos_theta_slope.divide(cos_theta_flat)

        fi = (ee.Image(0.75)
              .add(slope.cos().multiply(0.25))
              .subtract(slope.multiply(0.5 / pi)))

        sin_s_2 = slope.multiply(0.5).sin()

        fia = ee.Image().expression(
            "(1 - KB) * (1 + ((KB/(KB+KD))**0.5) * (sin_s_2**3))*fi + fB*KB",
            {"KB": KB, "KD": KD, "fi": fi, "fB": fB, "sin_s_2": sin_s_2}
        )

        ratio = ee.Image().expression(
            "fB*KB/(KB+KD) + fia*KD/(KB+KD) + albedo * (1 - fi)",
            {"fB": fB, "fia": fia, "fi": fi, "KD": KD, "KB": KB, "albedo": albedo}
        )
        return ratio

    slope  = ee.Terrain.slope(elev).multiply(pi/180.0)
    aspect = ee.Terrain.aspect(elev).subtract(180).multiply(pi/180.0)
    ll = ee.Image.pixelLonLat()
    lat = ll.select("latitude")
    lon = ll.select("longitude")

    shade_coeff2 = daily_ratio(J, KB, KD, slope, aspect, lat, lon, albedo)

    # corrected SRAD
    return Rs_MJ.multiply(shade_coeff2).rename("Rs_MJ_shade")