import ee


def evi2(landsat_image):
    """two band enhanced vegetation index (Jiang et al., 2008)

    Parameters
    ----------
    landsat_image : ee.Image
        Landsat image with standardized band names.

    Returns
    -------
    ee.Image

    References
    ----------

    """
    return landsat_image.expression(
        "2.5 * (nir - red) / (nir + red * 2.4 + 1)",
        {"red": landsat_image.select("red"), "nir": landsat_image.select("nir")}
    ).rename("evi2").clamp(0, 1.25)


def ndmi_scaled(landsat_image):
    """Normalized difference moisture index

    Parameters
    ----------
    landsat_image : ee.Image
        Landsat image with standardized band names.

    Returns
    -------
    ee.Image

    References
    ----------
    .. [Gao1996] Gao 1996.

    """
    return (
        ee.Image(landsat_image).normalizedDifference(["nir","swir1"])
        .unmask(0).add(0.3).divide(0.3).clamp(0, 1)
        .rename("ndmi_scaled")
    )


def lai(landsat_image):
    """Leaf area index (Kang et al., 2016 + NDMI correction)

    Parameters
    ----------
    landsat_image : ee.Image
        Landsat image with standardized band names.

    Returns
    -------
    ee.Image

    References
    ----------
    .. [Kang2016] Kang et. al. (2016). How Universal Is the Relationship
        between Remotely Sensed Vegetation Indices and Crop Leaf Area Index?
        A Global Assessment. Remote Sensing, https://doi.org/10.3390/rs8070597

    """
    return ee.Image(landsat_image).expression(
        "NDMI_scaled * (2.92 * sqrt(EVI2) - 0.43) ** 2",
        {"EVI2": evi2(ee.Image(landsat_image)), "NDMI_scaled": ndmi_scaled(ee.Image(landsat_image))}
    ).clamp(0, 8).rename("lai")


def ndvi(landsat_image):
    """Normalized difference vegetation index

    Parameters
    ----------
    landsat_image : ee.Image
        Landsat image with standardized band names.

    Returns
    -------
    ee.Image

    References
    ----------

    """
    return ee.Image(landsat_image).normalizedDifference(["nir", "red"]).rename(["ndvi"]).unmask(0)
    

def ndwi(landsat_image):
    """Normalized difference water index

    Parameters
    ----------
    landsat_image : ee.Image
        Landsat image with standardized band names.

    Returns
    -------
    ee.Image

    References
    ----------

    """
    return ee.Image(landsat_image).normalizedDifference(["green", "nir"]).rename("ndwi").unmask(0)


def landsat_c2_qa_water_mask(landsat_image):
    """Extract water mask from the Landsat Collection 2 SR QA_PIXEL band.

    Parameters
    ----------
    landsat_image : ee.Image
        Landsat image with standardized band names.

    Returns
    -------
    ee.Image

    """
    return (
        ee.Image(landsat_image).select(["QA_PIXEL"])
        .rightShift(7).bitwiseAnd(1).neq(0)
        .rename(["qa_water"])
    )


def albedo_disalexi(landsat_image):
    """Total shortwave broadband albedo following [Liang2001]

    Parameters
    ----------
    landsat_image : ee.Image
        "Prepped" Landsat image with standardized band names.

    Returns
    -------
    ee.Image

    Notes
    -----
    [Liang2001] indicate that the green band is intentionally not used.
    Coefficients were derived for Landsat 7 ETM+, but were found to be
        "suitable" to Landsat 4/5 TM also.

    References
    ----------
    .. [Liang2001] Shunlin Liang (2001).  Narrowband to broadband conversions
        of land surface albedo - I Algorithms,
        Remote Sensing of Environment, Volume 76, Issue2, Pages 213-238,
        http://doi.org/10.1016/S0034-4257(00)00205-4

    """
    albedo_img = (
        landsat_image
        .select(["blue", "red", "nir", "swir1", "swir2"])
        .multiply([0.356, 0.130, 0.373, 0.085, 0.072])
    )
    return (
        albedo_img.select([0])
        .add(albedo_img.select([1])).add(albedo_img.select([2]))
        .add(albedo_img.select([3])).add(albedo_img.select([4]))
        .subtract(0.0018)
        .rename(["albedo"])
    )


def albedo_l457(landsat_image):
    """Albedo (Landsat 4/5/7)

    Parameters
    ----------
    landsat_image : ee.Image
        Landsat image with standardized band names.

    Returns
    -------
    ee.Image

    References
    ----------
    .. [Tasumi2008] M. Tasumi, R. Allen, R; Trezza,
       At-Surface Reflectance and Albedo from Satellite for Operational
       Calculation of Land Surface Energy Balance, Journal of Hydrology,
       https://doi.org/10.1061/(ASCE)1084-0699(2008)13:2(51)

    """
    return landsat_image.expression(
        "(0.254 * B1) + (0.149 * B2) + (0.147 * B3) + (0.311 * B4) + (0.103 * B5) + (0.036 * B7)",
        {
            "B1": landsat_image.select(["blue"]),
            "B2": landsat_image.select(["green"]),
            "B3": landsat_image.select(["red"]),
            "B4": landsat_image.select(["nir"]),
            "B5": landsat_image.select(["swir1"]),
            "B7": landsat_image.select(["swir2"]),
        }
    ).rename("albedo")


def albedo_l89(landsat_image):
    """Albedo (Landsat 8/9)

    Parameters
    ----------
    landsat_image : ee.Image
        Landsat image with standardized band names.

    Returns
    -------
    ee.Image

    References
    ----------
    .. [Ke2016] Y. Ke, J. Im, S. Park, H. Gong,
       Downscaling of MODIS One Kilometer Evapotranspiration
       Using Landsat-8 Data and Machine Learning Approaches,
       Remote Sensing, https://doi.org/10.3390/rs8030215

    """

    # CGM - These coefficients don"t sum to 1.0?
    return landsat_image.expression(
        "(0.130 * B1) + (0.115 * B2) + (0.143 * B3) + (0.180 * B4) + "
        "(0.281 * B5) + (0.108 * B6) + (0.042 * B7)", {
            "B1": landsat_image.select(["ultra_blue"]),
            "B2": landsat_image.select(["blue"]),
            "B3": landsat_image.select(["green"]),
            "B4": landsat_image.select(["red"]),
            "B5": landsat_image.select(["nir"]),
            "B6": landsat_image.select(["swir1"]),
            "B7": landsat_image.select(["swir2"]),
        }
    ).rename("albedo")


def cloud_mask_sr_l457(landsat_image):
    """Cloud mask (Landsat 4/5/7)

    Parameters
    ----------
    landsat_image : ee.Image
        Landsat image with standardized band names.

    Returns
    -------
    ee.Image

    Notes
    -----
    clear sky value = 66  (01000010)
    water value = 68   (01000100)

    References
    ----------

    """
    quality = landsat_image.select("pixel_qa")
    c01 = quality.eq(66)  # Clear (01000010)
    c02 = quality.eq(68)  # Water (01000100)
    mask = c01.Or(c02)

    return mask


def cloud_mask_sr_l8(landsat_image):
    """Cloud mask (Landsat 8)

    Parameters
    ----------
    landsat_image : ee.Image
        Landsat image with standardized band names.

    Returns
    -------
    ee.Image

    Notes
    -----
    clear sky value = 322  (00101000010)
    water value = 324   (00101000100)

    References
    ----------

    """
    quality = landsat_image.select("pixel_qa")
    c01 = quality.eq(322)
    c02 = quality.eq(324)
    c03 = quality.eq(1346)  # (10101000010)
    mask = c01.Or(c02).Or(c03)

    return mask


def cloud_mask_C2_l457(landsat_image):
    """Cloud mask (Landsat 4/5/7)

    Parameters
    ----------
    landsat_image : ee.Image
        Landsat image with standardized band names.

    Returns
    -------
    ee.Image

    Notes
    -----
    clear sky value = 5440  (0001010101000000)
    water value = 5504   (0001010110000000)

    References
    ----------

    """
    quality = landsat_image.select("QA_PIXEL")
    c01 = quality.eq(5440)
    c02 = quality.eq(5504)
    mask = c01.Or(c02)

    return mask


def cloud_mask_C2_l89(landsat_image):
    """Cloud mask (Landsat 8/9)

    Parameters
    ----------
    landsat_image : ee.Image
        Landsat image with standardized band names.

    Returns
    -------
    ee.Image

    Notes
    -----
    clear sky value = 21824  (0001010101000000)
    water value = 21952   (0001010110000000)

    References
    ----------

    """
    quality = landsat_image.select("QA_PIXEL")
    c01 = quality.eq(21824)
    c02 = quality.eq(21952)
    mask = c01.Or(c02)

    return mask


def water_mask(product="GLO"):
    """Water maskfor precomputing

    Parameters
    ----------
    product : string
        Product name. 
            GLO: Water mask from the Copernicus GLO DEM.
            OSM: Water mask from OpenStreetMap.

    Returns
    -------
    ee.Image

    Notes
    -----
   

    References
    ----------
    GLO 30 DGED: https://doi.org/10.5270/ESA-c5d3d65 
    OSM: https://gee-community-catalog.org/projects/osm_water/

    """
    if product =="GLO":

        water_mask = ee.Image("projects/openet/assets/features/water_mask_glo30_0p001")

    elif product =="OSM":

        water_mask = ee.Image("projects/openet/assets/features/water_mask_osm_0p001")

    return water_mask
