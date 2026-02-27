import re

import ee
import openet.core.common

from radet import landsat
from radet import model
from radet import utils


def lazy_property(fn):
    """Decorator that makes a property lazy-evaluated
    https://stevenloria.com/lazy-properties/
    """
    attr_name = "_lazy_" + fn.__name__

    @property
    def _lazy_property(self):
        if not hasattr(self, attr_name):
            setattr(self, attr_name, fn(self))
        return getattr(self, attr_name)

    return _lazy_property


class Image:
    """Google Earth Engine - RADET for Landsat image"""

    _C2_LST_CORRECT = True  # Enable (True) C2 LST correction to recalculate LST

    def __init__(
        self,
        image,
        meteorology_source="IDAHO_EPSCOR/GRIDMET",
        landcover_source="projects/sat-io/open-datasets/USGS/ANNUAL_NLCD/LANDCOVER",
        elevation_source="USGS/SRTMGL1_003",
        latitude=None,
        **kwargs,
    ):
        """Construct a generic RADET Image

        Parameters
        ----------
        image : ee.Image
            A "prepped" RADET input image.
            Image must have bands:
                albedo, emissivity, lai, lst, emissivity, ndvi, ndwi
            Image must have properties:
                SUN_ELEVATION, system:id, system:index, system:time_start
        meteorology_source : {"IDAHO_EPSCOR/GRIDMET"}
            Daily meteorology source collection ID.
            Meteorology collection must have bands for:
                min temp, max temp, solar rad, humidity, wind speed
        landcover_source : {"projects/sat-io/open-datasets/USGS/ANNUAL_NLCD/LANDCOVER",
                            "projects/sat-io/open-datasets/USGS/ANNUAL_NLCD/LANDCOVER/Annual_NLCD_LndCov_2023_CU_C1V1"}
            Land cover source collection or image ID.
            The default is "projects/sat-io/open-datasets/USGS/ANNUAL_NLCD/LANDCOVER".
        elevation_source : str, ee.Image
            Elevation source keyword or asset.
            The default is "USGS/SRTMGL1_003".
            Units must be in meters.
        kwargs : dict, optional
            et_reference_source : str, float
                Reference ET source (the default is None).
                Parameter is required if computing "et_fraction" or "et_reference".
            et_reference_band : str
                Reference ET band name (the default is None).
                Parameter is required if computing "et_fraction" or "et_reference".
            et_reference_factor : float, None
                Reference ET scaling factor.  The default is None which is
                equivalent to 1.0 (or no scaling).
            et_reference_resample : {"nearest", "bilinear", "bicubic", None}
                Reference ET resampling.  The default is None which is
                equivalent to nearest neighbor resampling.
            latitude : ee.Image, ee.Number, float, optional
                Latitude [deg].  If not set will default to ee.Image.pixelLonLat().

        Notes
        -----
        Standard percentiles are from Allen et al. (2013)

        """
        # Image
        self.image = image

        # Copy system properties
        self.id = self.image.get("system:id")
        self.index = self.image.get("system:index")
        self.time_start = self.image.get("system:time_start")
        self.properties = {
            "system:index": self.index,
            "system:time_start": self.time_start,
            "image_id": self.id,
        }
        # Build SCENE_ID from the (possibly merged) system:index
        scene_id = ee.List(ee.String(self.index).split("_")).slice(-3)
        self.scene_id = (
            ee.String(scene_id.get(0))
            .cat("_")
            .cat(ee.String(scene_id.get(1)))
            .cat("_")
            .cat(ee.String(scene_id.get(2)))
        )

        # Build WRS2_TILE from the scene_id
        self.wrs2_tile = (
            ee.String("p").cat(self.scene_id.slice(5, 8)).cat("r").cat(self.scene_id.slice(8, 11))
        )

        # Set server side date/time properties using the "system:time_start"
        self.date = ee.Date(self.time_start)
        self.year = ee.Number(self.date.get("year"))
        self.month = ee.Number(self.date.get("month"))
        self.start_date = ee.Date(utils.date_to_time_0utc(self.date))
        self.end_date = self.start_date.advance(1, "day")
        self.doy = ee.Number(self.date.getRelative("day", "year")).add(1).int()

        # Model input parameters
        self.meteorology_source = meteorology_source
        self.landcover_source = landcover_source
        self.elevation_source = elevation_source

        # Reference ET parameters
        try:
            self.et_reference_source = kwargs["et_reference_source"]
        except:
            self.et_reference_source = None
        try:
            self.et_reference_band = kwargs["et_reference_band"]
        except:
            self.et_reference_band = None
        try:
            self.et_reference_factor = kwargs["et_reference_factor"]
        except:
            self.et_reference_factor = None
        try:
            self.et_reference_resample = kwargs["et_reference_resample"]
        except:
            self.et_reference_resample = None

        # Check reference ET parameters
        if self.et_reference_factor and not utils.is_number(self.et_reference_factor):
            raise ValueError("et_reference_factor must be a number")
        if self.et_reference_factor and (self.et_reference_factor < 0):
            raise ValueError("et_reference_factor must be greater than zero")
        resample_methods = ["nearest", "bilinear", "bicubic"]
        if self.et_reference_resample and (self.et_reference_resample.lower() not in resample_methods):
            raise ValueError("unsupported et_reference_resample method")

        # Image projection and geotransform
        self.proj = self.image.select([0]).projection()
        self.crs = self.image.select([0]).projection().crs()
        self.transform = ee.List(ee.Dictionary(ee.Algorithms.Describe(image.projection())).get("transform"))
        self.geometry = self.image.select([0]).geometry()
        # self.latlon = ee.Image.pixelLonLat().reproject(self.proj)
        # self.coords = self.latlon.select(["longitude", "latitude"])

        # CGM - Needed for running the tests
        if ("latitude" not in kwargs.keys()) or (not kwargs["latitude"]):
            self.latitude = self.lai.multiply(0).add(ee.Image.pixelLonLat().select(["latitude"]))
        elif utils.is_number(kwargs["latitude"]):
            self.latitude = ee.Image.constant(latitude)
        elif isinstance(kwargs["latitude"], ee.computedobject.ComputedObject):
            self.latitude = latitude
        else:
            raise ValueError("invalid latitude parameter")

    @classmethod
    def from_image_id(cls, image_id, **kwargs):
        """Constructs a RADET Image instance from an image ID

        Parameters
        ----------
        image_id : str
            An earth engine image ID.
            (i.e. "LANDSAT/LC08/C02/T1_L2/LC08_044033_20170716")
        kwargs
            Keyword arguments to pass through to model init.

        Returns
        -------
        new instance of Image class

        """

        collection_methods = {
            "LANDSAT/LT04/C02/T1_L2": "from_landsat_c2_sr",
            "LANDSAT/LT05/C02/T1_L2": "from_landsat_c2_sr",
            "LANDSAT/LE07/C02/T1_L2": "from_landsat_c2_sr",
            "LANDSAT/LC08/C02/T1_L2": "from_landsat_c2_sr",
            "LANDSAT/LC09/C02/T1_L2": "from_landsat_c2_sr",
        }

        try:
            method_name = collection_methods[image_id.rsplit("/", 1)[0]]
        except KeyError:
            raise ValueError(f"unsupported collection ID: {image_id}")
        except Exception as e:
            raise Exception(f"unhandled exception: {e}")

        method = getattr(Image, method_name)

        return method(ee.Image(image_id), **kwargs)

    @classmethod
    def from_landsat_c2_sr(cls, sr_image, cloudmask_args={}, **kwargs):
        """Returns a RADET Image instance from a Landsat Collection 2 SR image

        Parameters
        ----------
        sr_image : ee.Image, str
            A raw Landsat Collection 2 SR image or image ID.
        cloudmask_args : dict
            keyword arguments to pass through to cloud mask function
        kwargs : dict
            Keyword arguments to pass through to Image init function

        Returns
        -------
        Image

        """

        sr_image = ee.Image(sr_image)

        # Use the SPACECRAFT_ID property identify each Landsat type
        spacecraft_id = ee.String(sr_image.get("SPACECRAFT_ID"))

        # Rename bands to generic names
        # CGM - Intentionally letting these lines be long to improve readability
        #   If the plan long term is to use the DisALEXI albedo calculation
        #   then we could remove the ultra_blue band and greatly simplify this section
        input_bands = ee.Dictionary(
            {
                "LANDSAT_4": ["SR_B1", "SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B7", "ST_B6", "QA_PIXEL", "ST_EMIS"],
                "LANDSAT_5": ["SR_B1", "SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B7", "ST_B6", "QA_PIXEL", "ST_EMIS"],
                "LANDSAT_7": ["SR_B1", "SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B7", "ST_B6", "QA_PIXEL", "ST_EMIS"],
                "LANDSAT_8": [
                    "SR_B1", "SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B6", "SR_B7", "ST_B10", "QA_PIXEL", "ST_EMIS"
                ],
                "LANDSAT_9": [
                    "SR_B1", "SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B6", "SR_B7", "ST_B10", "QA_PIXEL", "ST_EMIS"
                ],
            }
        )
        output_bands = ee.Dictionary(
            {
                "LANDSAT_4": ["blue", "green", "red", "nir", "swir1", "swir2", "lst", "QA_PIXEL", "ASTER_GED_emissivity"],
                "LANDSAT_5": ["blue", "green", "red", "nir", "swir1", "swir2", "lst", "QA_PIXEL", "ASTER_GED_emissivity"],
                "LANDSAT_7": ["blue", "green", "red", "nir", "swir1", "swir2", "lst", "QA_PIXEL", "ASTER_GED_emissivity"],
                "LANDSAT_8": [
                    "ultra_blue", "blue", "green", "red", "nir", "swir1", "swir2", "lst", "QA_PIXEL", "ASTER_GED_emissivity"
                ],
                "LANDSAT_9": [
                    "ultra_blue", "blue", "green", "red", "nir", "swir1", "swir2", "lst", "QA_PIXEL", "ASTER_GED_emissivity"
                ],
            }
        )
        scalars = ee.Dictionary(
            {
                "LANDSAT_4": [0.0000275, 0.0000275, 0.0000275, 0.0000275, 0.0000275, 0.0000275, 0.00341802, 1, 0.0001],
                "LANDSAT_5": [0.0000275, 0.0000275, 0.0000275, 0.0000275, 0.0000275, 0.0000275, 0.00341802, 1, 0.0001],
                "LANDSAT_7": [0.0000275, 0.0000275, 0.0000275, 0.0000275, 0.0000275, 0.0000275, 0.00341802, 1, 0.0001],
                "LANDSAT_8": [
                    0.0000275, 0.0000275, 0.0000275, 0.0000275, 0.0000275, 0.0000275, 0.0000275, 0.00341802, 1, 0.0001
                ],
                "LANDSAT_9": [
                    0.0000275, 0.0000275, 0.0000275, 0.0000275, 0.0000275, 0.0000275, 0.0000275, 0.00341802, 1, 0.0001
                ],
            }
        )
        offsets = ee.Dictionary(
            {
                "LANDSAT_4": [-0.2, -0.2, -0.2, -0.2, -0.2, -0.2, 149.0, 0, 0],
                "LANDSAT_5": [-0.2, -0.2, -0.2, -0.2, -0.2, -0.2, 149.0, 0, 0],
                "LANDSAT_7": [-0.2, -0.2, -0.2, -0.2, -0.2, -0.2, 149.0, 0, 0],
                "LANDSAT_8": [-0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, 149.0, 0, 0],
                "LANDSAT_9": [-0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, 149.0, 0, 0],
            }
        )

        prep_image = (
            sr_image.select(input_bands.get(spacecraft_id), output_bands.get(spacecraft_id))
            .multiply(ee.Number(ee.List(scalars.get(spacecraft_id))))
            .add(ee.Number(ee.List(offsets.get(spacecraft_id))))
        )

        # YK - RADET used Disalexi albedo 
        albedo = landsat.albedo_disalexi(prep_image)
        
        cloud_mask = ee.Algorithms.If(
            ee.List(["LANDSAT_8", "LANDSAT_9"]).contains(spacecraft_id),
            landsat.cloud_mask_C2_l89(sr_image),
            landsat.cloud_mask_C2_l457(sr_image),
        )

        # Water mask
        water_mask = landsat.water_mask(product="GLO")

        # # Default the cloudmask flags to True if they were not
        # # Eventually these will probably all default to True in openet.core
        # if "cirrus_flag" not in cloudmask_args.keys():
        #     cloudmask_args["cirrus_flag"] = True
        # if "dilate_flag" not in cloudmask_args.keys():
        #     cloudmask_args["dilate_flag"] = True
        # if "shadow_flag" not in cloudmask_args.keys():
        #     cloudmask_args["shadow_flag"] = True
        # if "snow_flag" not in cloudmask_args.keys():
        #     cloudmask_args["snow_flag"] = True
        # cloud_mask = openet.core.common.landsat_c2_sr_cloud_mask(
        #     sr_image, **cloudmask_args)

        # Check if passing c2_lst_correct arguments
        if "c2_lst_correct" in kwargs.keys():
            assert isinstance(kwargs["c2_lst_correct"], bool), "selection type must be a boolean"
            # Remove from kwargs since it is not a valid argument for Image init
            c2_lst_correct = kwargs.pop("c2_lst_correct")
        else:
            c2_lst_correct = cls._C2_LST_CORRECT

        if c2_lst_correct:
            lst = openet.core.common.landsat_c2_sr_lst_correct(sr_image, landsat.ndvi(prep_image))
        else:
            lst = prep_image.select(["lst"])

        # Build the input image
        # Don't compute LST since it is being provided
        # YK: Don't compute emissivity 
        input_image = ee.Image(
            [
                albedo,
                prep_image.select(["ASTER_GED_emissivity"], ["emissivity"]),
                landsat.lai(prep_image),
                lst.rename(["lst"]),
                landsat.ndvi(prep_image),
                landsat.ndwi(prep_image),
            ]
        )

        # Apply the cloud mask and add properties
        input_image = (
            input_image.updateMask(cloud_mask)
            .updateMask(water_mask.Not())
            .set(
                {
                    "system:index": sr_image.get("system:index"),
                    "system:time_start": sr_image.get("system:time_start"),
                    "system:id": sr_image.get("system:id"),
                    "SUN_ELEVATION": sr_image.get("SUN_ELEVATION"),
                }
            )
        )

        # Instantiate the class
        return cls(input_image, **kwargs)

    def calculate(self, variables=["et"]):
        """Return a multiband image of calculated variables

        Parameters
        ----------
        variables : list

        Returns
        -------
        ee.Image

        """

        output_images = []
        for v in variables:
            if v.lower() == "et":
                output_images.append(self.et.float())
            elif v.lower() == "et_fraction":
                output_images.append(self.et_fraction.float())
            elif v.lower() == "et_reference":
                output_images.append(self.et_reference.float())
            elif v.lower() == "lst":
                output_images.append(self.lst.float())
            elif v.lower() == "ndvi":
                output_images.append(self.ndvi.float())
            # CGM - The time and mask bands are needed for the interpolation
            elif v.lower() == "mask":
                output_images.append(self.mask)
            elif v.lower() == 'time':
                output_images.append(self.time)
            else:
                raise ValueError(f"unsupported variable: {v}")

        return ee.Image(output_images).set(self.properties)

    @lazy_property
    def et(self):
        """Compute RADET actual ET [mm day-1]"""
        et = model.et(
            albedo=self.albedo,
            emissivity=self.emissivity,
            lai=self.lai,
            lst=self.lst,
            meteorology_source=self.meteorology_source,
            landcover=self.landcover,
            elevation=self.elevation,
            time_start=self.time_start,
            # TODO: Remove after testing
            # proj=self.proj,
            # geometry_image=self.geometry,
            # coords=self.coords,
        )

        return et.rename("et").set(self.properties)

    @lazy_property
    def et_reference(self):
        """Reference ET for the image date"""
        if utils.is_number(self.et_reference_source):
            # Interpret numbers as constant images
            # CGM - Should we use the ee_types here instead?
            #   i.e. ee.ee_types.isNumber(self.et_reference_source)
            et_reference_img = ee.Image.constant(self.et_reference_source)
        elif type(self.et_reference_source) is str:
            # Assume a string source is an image collection ID (not an image ID)
            et_reference_coll = (
                ee.ImageCollection(self.et_reference_source)
                .filterDate(self.start_date, self.end_date)
                .select([self.et_reference_band])
            )
            et_reference_img = ee.Image(et_reference_coll.first())
            if self.et_reference_resample in ["bilinear", "bicubic"]:
                et_reference_img = et_reference_img.reproject(self.image.projection()).resample(
                    self.et_reference_resample
                )
        else:
            raise ValueError(f"unsupported et_reference_source: {self.et_reference_source}")

        if self.et_reference_factor:
            et_reference_img = et_reference_img.multiply(self.et_reference_factor)

        return et_reference_img.rename(["et_reference"]).set(self.properties)

        # Map ETr values directly to the input (i.e. Landsat) image pixels
        # The benefit of this is the ETr image is now in the same crs as the
        #   input image.  Not all models may want this though.
        # Note, doing this will cause the reference ET to be cloud masked.
        # return (
        #     self.ndvi.multiply(0).add(et_reference_img)
        #     .rename(["et_reference"]).set(self.properties)
        # )

    @lazy_property
    def et_fraction(self):
        """Fraction of reference ET (equivalent to the Kc)"""
        return self.et.divide(self.et_reference).rename(["et_fraction"]).set(self.properties)

    @lazy_property
    def albedo(self):
        """Albedo"""
        return self.image.select(["albedo"]).set(self.properties)

    @lazy_property
    def emissivity(self):
        """Emissivity"""
        return self.image.select(["emissivity"]).set(self.properties)

    @lazy_property
    def lai(self):
        """Leaf area index (LAI)"""
        return self.image.select(["lai"]).set(self.properties)

    @lazy_property
    def lst(self):
        """Land surface temperature (LST)"""
        return self.image.select(["lst"]).set(self.properties)

    @lazy_property
    def ndvi(self):
        """Normalized difference vegetation index (NDVI)"""
        return self.image.select(["ndvi"]).set(self.properties)

    @lazy_property
    def ndwi(self):
        """Normalized difference water index (NDWI)"""
        return self.image.select(["ndwi"]).set(self.properties)

    # TODO: If the model does not do any additional masking we might be able to
    #   build the mask from the NDVI or QA band instead of the ET
    # CGM - Had to switch to building the mask from ndvi to get the tests to pass
    @lazy_property
    def mask(self):
        """Mask of all active pixels (based on the final et)"""
        return self.ndvi.multiply(0).add(1).updateMask(1).uint8().rename(["mask"]).set(self.properties)
        # return self.et.multiply(0).add(1).updateMask(1).uint8().rename(["mask"]).set(self.properties)

    # CGM - This band shouldn't be needed but removing it is causing problems
    @lazy_property
    def time(self):
        """Return an image of the 0 UTC time (in milliseconds)"""
        return (
            self.mask
            .double().multiply(0).add(utils.date_to_time_0utc(self.date))
            .rename(['time']).set(self.properties)
        )

    # TODO: Consider adding support for elevation image collections that are tiled (e.g. new 3DEP)
    @lazy_property
    def elevation(self):
        """Elevation"""
        if utils.is_number(self.elevation_source):
            elev_img = ee.Image.constant(float(self.elevation_source))
        elif isinstance(self.elevation_source, ee.computedobject.ComputedObject):
            elev_img = self.elevation_source
        elif type(self.elevation_source) is str:
            elev_img = ee.Image(self.elevation_source)
        else:
            raise ValueError(f"Unsupported elevation_source: {self.elevation_source}\n")

        return elev_img.select([0], ["elevation"])

    @lazy_property
    def landcover(self):
        """Landcover"""
        if utils.is_number(self.landcover_source):
            lc_img = ee.Image.constant(int(self.landcover_source)).rename(["landcover"])
            self.lc_type = "NLCD"
        elif isinstance(self.landcover_source, ee.computedobject.ComputedObject):
            # If the source is an ee.Image assume it is an NLCD image
            lc_img = self.landcover_source.rename(["landcover"])
            self.landcover_type = "NLCD"
        elif re.match("projects/sat-io/open-datasets/USGS/ANNUAL_NLCD/LANDCOVER/Annual_NLCD_LndCov_\\d{4}_CU_\\w+",
                      self.landcover_source, re.I):
            # Assume an annual NLCD image ID was passed in and use it directly
            lc_img = ee.Image(self.landcover_source).rename(["landcover"])
            self.landcover_type = "NLCD"
        elif self.landcover_source == "projects/sat-io/open-datasets/USGS/ANNUAL_NLCD/LANDCOVER":
            # Select the closest year in time from the Annual NLCD image collection
            # Hardcoding the year ranges for now but we might want to change this to
            #   a more dynamic approach to allow for additional years to be added.
            lc_coll = ee.ImageCollection(self.landcover_source)
            lc_year = (
                ee.Number(self.year)
                .max(ee.Date(lc_coll.aggregate_min("system:time_start")).get("year"))
                .min(ee.Date(lc_coll.aggregate_max("system:time_start")).get("year"))
            )
            lc_img = (
                lc_coll.filter(ee.Filter.calendarRange(lc_year, lc_year, "year")).first()
                .rename(["landcover"])
            )
            self.landcover_type = "NLCD"
        else:
            raise ValueError(f"Unsupported landcover_source: {self.landcover_source}\n")

        return lc_img
