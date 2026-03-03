import ee

from radet import utils


# TODO: Come up with a better name for this function
def get_source_variable(source, variable, time_start):
    """Helper function for selecting the meteorology variable from the target source"""
    if utils.is_number(source):
        meteo_img = ee.Image.constant(float(source))
    elif isinstance(source, ee.computedobject.ComputedObject):
        meteo_img = ee.Image(source)
    elif source in ["IDAHO_EPSCOR/GRIDMET", "GRIDMET"]:
        meteo_img = gridmet(variable=variable, time_start=time_start)
    else:
        raise ValueError(f"Unsupported source: {source}\n")

    return meteo_img.rename([variable])


def elevation(source):
    """Get the meteorology elevation based on the temperature source"""
    if utils.is_number(source):
        return ee.Image.constant(float(source)).rename('elevation')
    elif source in ["IDAHO_EPSCOR/GRIDMET", "GRIDMET"]:
        return ee.Image("projects/openet/assets/meteorology/gridmet/ancillary/elevation")
    else:
        raise ValueError("Unsupported temperature source for selecting meteorology elevation: {variable}")


# TODO: Should the collection ID be an input to the function?
#   It would make it easier for the user to reuse one meteorology function
#   for a different similar one (like URMA and RTMA) but is probably not needed
def gridmet(variable, time_start):
    """GRIDMET daily meteorology

    Parameters
    ----------
    variable : {"tmin", "tmax", "qa", "u10", "srad"}
        Standard meteorology variable names used in RADET model.
    time_start : int, ee.Number
        Image property: time start of the image.

    Returns
    -------
    ee.Image

    """
    # Get date information
    time_start = ee.Number(time_start)

    # Filtering daily data
    meteorology_daily = (
        ee.ImageCollection("IDAHO_EPSCOR/GRIDMET")
        .filterDate(ee.Date(time_start).advance(-1, "day"), ee.Date(time_start))
        .first()
    )

    # Map the standardized meteorology variable names to the GRIDMET bands
    if variable == "tmin":
        return meteorology_daily.select(["tmmn"], [variable])
    elif variable == "tmax":
        return meteorology_daily.select(["tmmx"], [variable])
    elif variable == "qa":
        return meteorology_daily.select(["sph"], [variable])
    elif variable == "u10":
        return meteorology_daily.select(["vs"], [variable])
    elif variable == "srad":
        return meteorology_daily.select(["srad"], [variable])
    else:
        raise ValueError("Unsupported GRIDMET meteorology variable: {variable}")


# def era5land(time_start, meteorology_source_inst, meteorology_source_daily):
#     """
#
#     Parameters
#     ----------
#     time_start : str
#         Image property: time start of the image.
#     meteorology_source_inst: ee.ImageCollection, str
#         Instantaneous meteorological data.
#     meteorology_source_daily :  ee.ImageCollection, str
#         Daily meteorological data.
#
#     Returns
#     -------
#     ee.Image
#
#     Notes
#     -----
#     Accepted collections:
#     Inst : ECMWF/ERA5_LAND/HOURLY
#     Daily : projects/openet/assets/meteorology/era5land/na/daily
#             projects/openet/assets/meteorology/era5land/sa/daily
#
#     """
#
#     # Get date information
#     time_start = ee.Number(time_start)
#
#     # Filtering Daily data
#     meteorology_daily = (
#         ee.ImageCollection(meteorology_source_daily)
#         .filterDate(ee.Date(time_start).advance(-1, "day"), ee.Date(time_start).advance(1, "day"))
#         .first()
#     )
#
#     # Instantaneous data
#     meteorology_inst_collection = ee.ImageCollection(meteorology_source_inst)
#
#     # Linear interpolation
#     previous_time = time_start.subtract(1 * 60 * 60 * 1000)
#     next_time = time_start.add(1 * 60 * 60 * 1000)
#
#     previous_image = (
#         meteorology_inst_collection.filterDate(previous_time, time_start)
#         .limit(1, "system:time_start", False)
#         .first()
#     )
#
#     next_image = (
#         meteorology_inst_collection.filterDate(time_start, next_time)
#         .limit(1, "system:time_start", True)
#         .first()
#     )
#
#     image_previous_time = ee.Number(previous_image.get("system:time_start"))
#     image_next_time = ee.Number(next_image.get("system:time_start"))
#
#     delta_time = time_start.subtract(image_previous_time).divide(
#         image_next_time.subtract(image_previous_time)
#     )
#
#     # Incoming shorwave down [W m-2]
#     swdown24h = meteorology_daily.select("surface_solar_radiation_downwards").divide(1 * 60 * 60 * 24)
#
#     # Minimum air temperature [K]
#     tmin = meteorology_daily.select("temperature_2m_min").rename("tmin")
#
#     # Maximum air temperature [K]
#     tmax = meteorology_daily.select("temperature_2m_max").rename("tmax")
#
#     # Instantaneous incoming shortwave radiation [W m-2]
#     rso_inst = (
#         ee.ImageCollection(meteorology_source_inst)
#         .filterDate(ee.Date(time_start), ee.Date(time_start).advance(1, "hour"))
#         .select("surface_solar_radiation_downwards_hourly")
#         .mean()
#         .divide(1 * 60 * 60)
#         .rename("rso_inst")
#     )
#
#     # Air temperature [C]
#     # TODO: LL- Change all temperatures to K ?
#     tair_c = (
#         next_image.select("temperature_2m")
#         .subtract(previous_image.select("temperature_2m"))
#         .multiply(delta_time)
#         .add(previous_image.select("temperature_2m"))
#         .subtract(273.15)
#         .rename("tair")
#     )
#
#     # Wind speed [ m/s]
#     wind_u = (
#         next_image.select("u_component_of_wind_10m")
#         .subtract(previous_image.select("u_component_of_wind_10m"))
#         .multiply(delta_time)
#         .add(previous_image.select("u_component_of_wind_10m"))
#     )
#
#     wind_v = (
#         next_image.select("v_component_of_wind_10m")
#         .subtract(previous_image.select("v_component_of_wind_10m"))
#         .multiply(delta_time)
#         .add(previous_image.select("v_component_of_wind_10m"))
#     )
#
#     wind_med = wind_u.expression("sqrt(ux_u ** 2 + ux_v ** 2)", {"ux_u": wind_u, "ux_v": wind_v})
#
#     wind_med = wind_med.expression("ux * (4.87) / log(67.8 * z - 5.42)", {"ux": wind_med, "z": 10.0}).rename("ux")
#
#     # Dew point temperature [°K]
#     tdp = (
#         next_image.select("dewpoint_temperature_2m")
#         .subtract(previous_image.select("dewpoint_temperature_2m"))
#         .multiply(delta_time)
#         .add(previous_image.select("dewpoint_temperature_2m"))
#         .rename("tdp")
#     )
#
#     # Actual vapour pressure [kPa]
#     ea = tdp.expression("0.6108 * (exp((17.27 * T_air) / (T_air + 237.3)))", {"T_air": tdp.subtract(273.15)})
#
#     # SATURATED VAPOR PRESSURE [kPa]
#     esat = tair_c.expression("0.6108 * (exp((17.27 * T_air) / (T_air + 237.3)))", {"T_air": tair_c})
#
#     # RELATIVE HUMIDITY (%)
#     rh = ea.divide(esat).multiply(100).rename("RH")
#
#     # Surface temperature correction based on precipitation and reference ET
#
#     # Accumulation time period
#     accum_period = -60
#
#     # Accum meteo data
#     gridmet_accum = ee.ImageCollection(meteorology_source_daily).filterDate(
#         ee.Date(time_start).advance(accum_period, "days"), ee.Date(time_start)
#     )
#
#     # Reference ET
#     etr_accum = gridmet_accum.select("etr_asce").sum()
#
#     # Precipitation
#     precipt_accum = gridmet_accum.select("total_precipitation").sum()
#
#     # Ratio between precipt/etr
#     ratio = precipt_accum.divide(etr_accum)
#
#     # Temperature adjustment offset (Allen2013 Eqn 8)
#
#     tfac = etr_accum.expression("2.6 - 13 * ratio", {"ratio": ratio})
#
#     tfac = ee.Image(tfac.where(ratio.gt(0.2), 0)).rename("tfac")
#
#     # Resample
#     tmin = tmin.subtract(273.15).resample("bilinear")
#     tmax = tmax.subtract(273.15).resample("bilinear")
#     rso_inst = rso_inst.resample("bilinear")
#     tair_c = tair_c.resample("bilinear")
#     wind_med = wind_med.resample("bilinear")
#     rh = rh.resample("bilinear")
#     swdown24h = swdown24h.resample("bilinear")
#
#     return [tmin, tmax, tair_c, wind_med, rh, rso_inst, swdown24h, tfac]
