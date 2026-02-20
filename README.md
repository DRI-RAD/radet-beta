# RADET - beta


**WARNING:** This code is in development and may change without notice.

This repository provides a Google Earth Engine (Python API) implementation of the **RADET model (Radiation Advection Diffusivity-independent Evapotranspiration)** for estimating actual evapotranspiration (ET). RADET estimates ET based on the Diffusivity-Independent Flux hypothesis and conditionally incorporates Penman’s aerodynamic term when and where advection is expected to be significant (Kim et al., 2026). The RADET-beta implementation here is designed to be consistent with the OpenET Python pipeline to facilitate interoperability and integration within existing workflows.

## Model Design

The primary component of the RADET model is the Image() class. The Image class can be used to compute a single ET image from a single input image.  The Image class should generally be instantiated from an Earth Engine Landsat image using the collection specific methods listed below.  ET image collections can be built by computing ET in a function that is mapped over a collection of input images. 

## Input Collections

RADET-beta can currently be computed for Landsat Collection 2 Level 2 (SR/ST) images  images from the following Earth Engine image collections:

 * LANDSAT/LT05/C02/T1_L2
 * LANDSAT/LE07/C02/T1_L2
 * LANDSAT/LC08/C02/T1_L2
 * LANDSAT/LC09/C02/T1_L2

### Landsat Collection 2 SR/ST Input Image

To instantiate the class for a Landsat Collection 2 SR/ST image, use the Image.from_landsat_c2_sr method.

The input Landsat image must have the following bands and properties:

| SPACECRAFT_ID | Band Names |
|---------------|------------|
| LANDSAT_5 | SR_B1, SR_B2, SR_B3, SR_B4, SR_B5, SR_B7, ST_B6, ST_EMIS, QA_PIXEL |
| LANDSAT_7 | SR_B1, SR_B2, SR_B3, SR_B4, SR_B5, SR_B7, ST_B6, ST_EMIS, QA_PIXEL |
| LANDSAT_8 | SR_B1, SR_B2, SR_B3, SR_B4, SR_B5, SR_B6, SR_B7, ST_B10, ST_EMIS, QA_PIXEL |
| LANDSAT_9 | SR_B1, SR_B2, SR_B3, SR_B4, SR_B5, SR_B6, SR_B7, ST_B10, ST_EMIS, QA_PIXEL |

| Property          | Description |
|-------------------|-------------|
| `system:index`    | - Landsat Scene ID<br>- Must be in Earth Engine format (e.g. `LC08_044033_20170716`)<br>- Used to lookup the scene-specific c-factor |
| `system:time_start` | Image datetime in milliseconds since 1970 |
| `SPACECRAFT_ID`   | - Used to determine Landsat sensor type<br>- Must be one of: `LANDSAT_5`, `LANDSAT_7`, `LANDSAT_8`, `LANDSAT_9` |

### Model Output

The primary output of the RADET-beta is the actual ET (ETa) in millimeters.

### Example

Examples of how to use RADET model are detailed in examples folder.

## Dependencies

- [earthengine-api](https://github.com/google/earthengine-api)
- [openet-core](https://github.com/Open-ET/openet-core)


## References

Kim, Y., Huntington, J. L., Comini de Andrade, B., Johnson, M. S., Volk, J. M., Majumdar, S., Morton, C., & ReVelle, P. (2026). Thermodynamically constrained closed-form surface energy balance using medium-resolution remote sensing for efficient evapotranspiration mapping. *EarthArXiv (preprint)*. [https://doi.org/10.31223/X51B4P](https://doi.org/10.31223/X51B4P)
