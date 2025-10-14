import math

import ee

pi = math.pi
SIG = 4.901e-9


def et(
    image,
    ndvi,
    ndwi,
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
    ndvi : ee.Image
        Normalized difference vegetation index.
    ndwi : ee.Image
        Normalized difference water index.
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
        srad, tminK, tmaxK, qa, t_avg, u2, gamma, ea, esat, DELTA = meteorology_gridmet(
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

    # Clear-Sky terms
    Rs_MJ, fcd, sunrise_ts = clear_sky_terms(time_start, lat, elev, srad)

    # LAI <- 10 if open water
    lai_corrected = lai_correct_water(lai)

    # Penman aerodynamic term
    Ea = aerodynamic_term(lai_corrected, gamma, DELTA, u2, esat, ea)

    # Daily mean LST and surface humidity
    LST_avg = daily_avg_lst(tminK, lst, sunrise_ts, t_avg)

    esat_LST = ee.Image().expression("0.6108 * exp(17.27*(t-273.15)/((t-273.15)+237.3))", {"t": LST_avg})
    RHs = ee.Image().expression("ea/esat_LST", {"ea": ea, "esat_LST": esat_LST})

    # Daily mean net radiation (total, soil and canopy)
    Rn, Rns, Rnc = net_radiation(emissivity, LST_avg, ea, esat_LST, fcd, t_avg, Rs_MJ, albedo, lai_corrected)

    # Soil heat flux
    G = ee.Image().expression("0.35*Rns - 2", {"Rns": Rns}).rename("G")

    # Isothermal net radiation and available energy
    Rni, AEs, AEsi = Rni_AE(Rn, emissivity, t_avg, LST_avg, Rns, G)

    # mu terms
    mu_c, mu_s = mu_terms(Rn, Rni, DELTA, gamma, AEs, AEsi, RHs)

    # ET model (mm/day)
    ET_DIF = DIF_model(gamma, DELTA, Rnc, AEs, RHs, mu_c, mu_s)

    return ee.Image(ET_DIF.add(Ea).rename("ET")).setDefaultProjection(proj)

def lai_correct_water(lai):
    # USGS NLCD land cover (2020)
    # Load land cover image
    # TODO: make the search for land cover images dynamic

    landcover = ee.ImageCollection("projects/sat-io/open-datasets/USGS/ANNUAL_NLCD/LANDCOVER")
    landcover_2020 = landcover.filter(ee.Filter.eq("year", 2020)).first()

    # Define open water mask and set LAI = 10
    water = landcover_2020.remap([11], [1], 0).eq(1)
    lai_corrected = lai.where(water, 10)
    
    return lai_corrected


def wet_mask(lai):
    # USGS NLCD land cover (2020)
    # Load land cover image
    # TODO: make the search for land cover images dynamic

    landcover = ee.ImageCollection("projects/sat-io/open-datasets/USGS/ANNUAL_NLCD/LANDCOVER")
    landcover_2020 = landcover.filter(ee.Filter.eq("year", 2020)).first()

    # Define wet surfaces
    wet = landcover_2020.remap([11, 81, 82, 95], [1, 1, 1, 1], 0).eq(1)
    wet_ww = landcover_2020.remap([90], [1], 0).eq(1)
    wet_ww2 = wet_ww.And(lai.lt(1))
    return ee.Image(wet.Or(wet_ww2))


def wind_10m_2m(u):
    # Wind speed at 2 m
    return ee.Image().expression("u * 4.87 / log(67.8*10 - 5.42)", {"u": u})


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
    Inst : NASA/NLDAS/FORA0125_H002
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
    t_avg = tmaxK.add(tminK).multiply(0.5).rename("t_avg")
    u2 = wind_10m_2m(u10)

    # Elevation image
    elev = ee.Image("projects/openet/assets/meteorology/gridmet/ancillary/elevation")

    # Air pressure (kPa)
    PA = ee.Image().expression("101.3 * ((293 - 0.0065*elev)/293)**5.26", {"elev": elev})

    # Psychrometric constant
    gamma = ee.Image().expression("PA * 0.000665", {"PA": PA})

    # Vapor pressure (kPa)
    ea = ee.Image().expression("qa * PA/((1 - 0.622) * qa + 0.622)", {"PA": PA, "qa": qa})

    # Saturation vapor pressure (kPa)
    esat = ee.Image().expression("0.6108 * exp(17.27*(t-273.15)/((t-273.15)+237.3))", {"t": t_avg})

    # Slope of vapor pressure curve
    DELTA = ee.Image().expression(
        "2503*exp(17.27*(t-273.15)/((t-273.15)+237.3))/((t-273.15)+237.3)**2", {"t": t_avg}
    )

    return srad, tminK, tmaxK, qa, t_avg, u2, gamma, ea, esat, DELTA


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

    return Rs_MJ, fcd, sunrise_ts


def daily_avg_lst(tminK, lst, sunrise_ts, t_avg):
    LST_min = tminK
    LST_max = lst.expression(
        "LST_min + (LST - LST_min)/(cos((10-12.5)/(12.5 - sunrise_ts)/2*PI))",
        {"PI": pi, "LST_min": LST_min, "LST": lst, "sunrise_ts": sunrise_ts},
    )
    return LST_max.add(LST_min).multiply(0.5).max(t_avg)


def net_radiation(emissivity, LST_avg, ea, esat_LST, fcd, t_avg, Rs_MJ, albedo, lai_corrected):
    Rlu = emissivity.multiply(LST_avg.pow(4).multiply(SIG)).rename("Rlup")

    # Incoming longwave radiation
    Rld = ee.Image().expression(
        "emis*SIG * (1 - fcd*(0.34 - 0.14 * ea ** 0.5)) * t_avg **4",
        {"emis": emissivity, "SIG": SIG, "fcd": fcd, "ea": ea, "t_avg": t_avg},
    )

    # Net radiation
    Rn = (
        ee.Image()
        .expression(
            "(1-albedo)*Rs_MJ - Rlu + Rld",
            {"Rs_MJ": Rs_MJ, "albedo": albedo, "Rld": Rld, "Rlu": Rlu},
        )
        .rename("Rn")
    )

    Rns = ee.Image().expression("Rn*exp(-0.6*lai)", {"Rn": Rn, "lai": lai_corrected}).rename("Rns")
    Rnc = ee.Image().expression("Rn - Rns", {"Rn": Rn, "Rns": Rns}).rename("Rnc")

    return Rn, Rns, Rnc


def Rni_AE(Rn, emissivity, t_avg, LST_avg, Rns, G):
    # Isothermal net radiation and available energy
    Rni = (
        ee.Image()
        .expression(
            "Rn + 4*emis*SIG*ta**3*(LST-ta)",
            {"Rn": Rn, "emis": emissivity, "SIG": SIG, "ta": t_avg, "LST": LST_avg},
        )
        .rename("Rni")
    )

    AEs = ee.Image().expression("Rns - G", {"Rns": Rns, "G": G}).rename("AEs")

    AEsi = (
        ee.Image()
        .expression(
            "AEs + (4*emis*SIG*ta**3 + 1.3)*(LST-ta)",
            {"AEs": AEs, "emis": emissivity, "SIG": SIG, "ta": t_avg, "LST": LST_avg},
        )
        .rename("AEsi")
    )

    return Rni, AEs, AEsi


def mu_terms(Rn, Rni, DELTA, gamma, AEs, AEsi, RHs):
    mu_c = (
        ee.Image()
        .expression(
            "(Rni + sqrt(Rni**2 + 4*DELTA/gamma*Rn*(Rni-Rn)))/(2*Rn)",
            {"Rn": Rn, "Rni": Rni, "DELTA": DELTA, "gamma": gamma},
        )
        .max(1)
        .rename("mu_c")
    )

    mu_s = (
        ee.Image()
        .expression(
            "(AEsi + sqrt(AEsi**2 + 4*RHs*DELTA/gamma*AEs*(AEsi-AEs)))/(2*AEs)",
            {"AEs": AEs, "AEsi": AEsi, "DELTA": DELTA, "gamma": gamma, "RHs": RHs},
        )
        .max(1)
        .rename("mu_s")
    )

    return mu_c, mu_s


def aerodynamic_term(lai_corrected, gamma, DELTA, u2, esat, ea):
    wet = wet_mask(lai_corrected)

    return (
        ee.Image()
        .expression(
            "wet_mask*(1-exp(-0.7*lai))*gamma/(DELTA + gamma)*(esat-ea)*2.6*(1+0.54*u2)",
            {"gamma": gamma, "DELTA": DELTA, "u2": u2, "esat": esat, "ea": ea, "lai": lai_corrected, "wet_mask": wet},
        )
        .rename("Ea")
    )


def DIF_model(gamma, DELTA, Rnc, AEs, RHs, mu_c, mu_s):
    return (
        ee.Image()
        .expression(
            "(DELTA/(DELTA+gamma*mu_c)*Rnc + RHs*DELTA/(RHs*DELTA+gamma*mu_s)*AEs)/2.45",
            {"gamma": gamma, "DELTA": DELTA, "Rnc": Rnc, "AEs": AEs, "RHs": RHs, "mu_c": mu_c, "mu_s": mu_s},
        )
        .rename("ET_DIF")
        .max(0)
    )
