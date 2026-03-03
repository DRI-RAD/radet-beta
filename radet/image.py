import re

import ee
import openet.core.common

from radet import landsat
from radet import meteorology
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
        temperature_source="IDAHO_EPSCOR/GRIDMET",
        humidity_source="IDAHO_EPSCOR/GRIDMET",
        windspeed_source="IDAHO_EPSCOR/GRIDMET",
        solar_radiation_source="IDAHO_EPSCOR/GRIDMET",
        landcover_source="projects/sat-io/open-datasets/USGS/ANNUAL_NLCD/LANDCOVER",
        elevation_source="USGS/SRTMGL1_003",
        **kwargs,
    ):
        """Construct a generic RADET Image

        Parameters
        ----------
        image : ee.Image
            A "prepped" RADET input image.
            Image must have bands:
                albedo, emissivity, lai, lst, ndvi, ndwi
            Image must have properties:
                system:id, system:index, system:time_start
        temperature_source : {"IDAHO_EPSCOR/GRIDMET"}
            Temperature source collection ID or keyword.
        humidity_source : {"IDAHO_EPSCOR/GRIDMET"}
            Humidity source collection ID or keyword.
        windspeed_source : {"IDAHO_EPSCOR/GRIDMET"}
            Wind speed source collection ID or keyword.
        solar_radiation_source : {"IDAHO_EPSCOR/GRIDMET"}
            Solar radiation source collection ID or keyword.
        landcover_source : {"projects/sat-io/open-datasets/USGS/ANNUAL_NLCD/LANDCOVER",
                            "projects/sat-io/open-datasets/USGS/ANNUAL_NLCD/LANDCOVER/Annual_NLCD_LndCov_2023_CU_C1V1"}
            Land cover source collection or image ID.
            The default is "projects/sat-io/open-datasets/USGS/ANNUAL_NLCD/LANDCOVER".
        elevation_source : str, ee.Image
            Elevation source keyword or asset. The default is "USGS/SRTMGL1_003".
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
            longitude : ee.Image, ee.Number, float, optional
                Longitude [deg].  If not set will default to ee.Image.pixelLonLat().
            meteo_elevation_source : str, optional
                Meteorology elevation source keyword or asset.
                If not set the temperature_source keyword or collection ID will be used.

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
        self.temperature_source = temperature_source
        self.humidity_source = humidity_source
        self.windspeed_source = windspeed_source
        self.solar_radiation_source = solar_radiation_source
        self.landcover_source = landcover_source
        self.elevation_source = elevation_source

        self.init_kwargs = kwargs

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
    def from_landsat_c2_sr(cls, sr_image, mask_ocean_flag=True, cloudmask_args={}, **kwargs):
        """Returns a RADET Image instance from a Landsat Collection 2 SR image

        Parameters
        ----------
        sr_image : ee.Image, str
            A raw Landsat Collection 2 SR image or image ID.
        mask_ocean_flag : bool
            Mask ocean pixels.
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
        #   then we could remove the ultra_blue band and simplify this section
        #   to a single set of band names and coefficients
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
            .multiply(ee.Image.constant(ee.List(scalars.get(spacecraft_id))))
            .add(ee.Image.constant(ee.List(offsets.get(spacecraft_id))))
        )

        # YK - RADET used Disalexi albedo 
        albedo = landsat.albedo_disalexi(prep_image)
        
        cloud_mask = ee.Algorithms.If(
            ee.List(["LANDSAT_8", "LANDSAT_9"]).contains(spacecraft_id),
            landsat.cloud_mask_C2_l89(sr_image),
            landsat.cloud_mask_C2_l457(sr_image),
        )

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
            lst = openet.core.common.landsat_c2_sr_lst_correct(sr_image)
        else:
            lst = prep_image.select(["lst"])

        # Build the input image from the components
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

        input_image = input_image.updateMask(cloud_mask)

        # Apply the ocean mask
        if mask_ocean_flag:
            # TODO: Consider renaming this variable (and function) to "ocean_mask"
            #   since the two supported products are not general water masks
            input_image = input_image.updateMask(landsat.water_mask(product="GLO").Not())

        input_image = input_image.set(
            {
                "system:index": sr_image.get("system:index"),
                "system:time_start": sr_image.get("system:time_start"),
                "system:id": sr_image.get("system:id"),
            }
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
            landcover=self.landcover,
            elevation=self.elevation,
            tmin=self.tmin,
            tmax=self.tmax,
            qa=self.qa,
            u10=self.u10,
            srad=self.srad,
            meteo_elevation=self.meteo_elevation,
            time_start=self.time_start,
            latitude=self.latitude,
            longitude=self.longitude,
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

    @lazy_property
    def tmin(self):
        """Daily minimum air temperature [K]"""
        return meteorology.get_source_variable(
            self.temperature_source, variable="tmin", time_start=self.time_start
        )

    @lazy_property
    def tmax(self):
        """Daily maximum air temperature [K]"""
        return meteorology.get_source_variable(
            self.temperature_source, variable="tmax", time_start=self.time_start
        )

    # TODO: Consider renaming to something else like sph or specific_humidity
    @lazy_property
    def qa(self):
        """Specific humidity"""
        return meteorology.get_source_variable(
            self.humidity_source, variable="qa", time_start=self.time_start
        )

    @lazy_property
    def u10(self):
        """Daily average wind speed [m s-1]"""
        return meteorology.get_source_variable(
            self.windspeed_source, variable="u10", time_start=self.time_start
        )

    @lazy_property
    def srad(self):
        """Daily incoming solar radiation [W m-2]"""
        return meteorology.get_source_variable(
            self.solar_radiation_source, variable="srad", time_start=self.time_start
        )

    @lazy_property
    def meteo_elevation(self):
        """Meteorology (temperature) elevation"""
        if "meteo_elevation_source" in self.init_kwargs.keys():
            return meteorology.elevation(self.init_kwargs["meteo_elevation_source"])
        else:
            return meteorology.elevation(self.temperature_source)

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
            lc_img = ee.Image.constant(int(self.landcover_source))
        elif isinstance(self.landcover_source, ee.computedobject.ComputedObject):
            # If the source is an ee.Image object assume it is an NLCD image
            lc_img = ee.Image(self.landcover_source)
        elif re.match("projects/sat-io/open-datasets/USGS/ANNUAL_NLCD/LANDCOVER/Annual_NLCD_LndCov_\\d{4}_CU_\\w+",
                      self.landcover_source, re.I):
            # Check if the source is similar to an awesome GEE catalog annual NLCD image ID
            lc_img = ee.Image(self.landcover_source)
        elif self.landcover_source == "projects/sat-io/open-datasets/USGS/ANNUAL_NLCD/LANDCOVER":
            # Check if the source is similar to an awesome GEE catalog annual NLCD collection ID
            # Select the closest image in time to the target Landsat image
            lc_coll = ee.ImageCollection(self.landcover_source)
            lc_year = (
                ee.Number(self.year)
                .max(ee.Date(lc_coll.aggregate_min("system:time_start")).get("year"))
                .min(ee.Date(lc_coll.aggregate_max("system:time_start")).get("year"))
            )
            lc_img = (
                lc_coll.filter(ee.Filter.calendarRange(lc_year, lc_year, "year")).first()
                .set({'landcover_year': lc_year})
            )
        else:
            raise ValueError(f"Unsupported landcover_source: {self.landcover_source}\n")

        return lc_img.rename(["landcover"])

    @lazy_property
    def latitude(self):
        """Longitude [deg]"""
        if ("latitude" not in self.init_kwargs.keys()) or (not self.init_kwargs["latitude"]):
            return self.lai.multiply(0).add(ee.Image.pixelLonLat().select(["latitude"]))
            # return ee.Image.pixelLonLat().select(["latitude"])
        elif utils.is_number(self.init_kwargs["latitude"]):
            return ee.Image.constant(self.init_kwargs["latitude"])
        elif isinstance(self.init_kwargs["latitude"], ee.computedobject.ComputedObject):
            return self.init_kwargs["latitude"]
        else:
            raise ValueError("invalid latitude parameter")

    @lazy_property
    def longitude(self):
        """Latitude [deg]"""
        if ("longitude" not in self.init_kwargs.keys()) or (not self.init_kwargs["longitude"]):
            return self.lai.multiply(0).add(ee.Image.pixelLonLat().select(["longitude"]))
            # return ee.Image.pixelLonLat().select(["longitude"])
        elif utils.is_number(self.init_kwargs["longitude"]):
            return ee.Image.constant(self.init_kwargs["longitude"])
        elif isinstance(self.init_kwargs["longitude"], ee.computedobject.ComputedObject):
            return self.init_kwargs["longitude"]
        else:
            raise ValueError("invalid longitude parameter")

    @lazy_property
    def mask(self):
        """Mask of all active pixels (based on the final et)

        Notes
        -----
        This function is needed for the collection interpolation

        """
        # CGM - Had to switch to building the mask from ndvi to get the tests to pass
        #   and it may be okay to do this if no additional masking is applied within
        #   the model.et() function
        return self.lai.multiply(0).add(1).updateMask(1).uint8().rename(["mask"]).set(self.properties)
        # return self.et.multiply(0).add(1).updateMask(1).uint8().rename(["mask"]).set(self.properties)

    @lazy_property
    def time(self):
        """Image of the 0 UTC time (in milliseconds) for all active pixels

        Notes
        -----
        This function is needed for the collection interpolation

        """
        return (
            self.mask
            .double().multiply(0).add(utils.date_to_time_0utc(self.date))
            .rename(['time']).set(self.properties)
        )