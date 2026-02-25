# RADET - beta

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/License-Apache%202.0-green)](LICENSE)
[![GEE](https://img.shields.io/badge/Google%20Earth%20Engine-4285F4?logo=google-earth&logoColor=white)](https://earthengine.google.com/)
![Status](https://img.shields.io/badge/Status-Beta-yellow)
[![EarthArXiv Preprint](https://img.shields.io/badge/EarthArXiv-10.31223%2FX51B4P-blue)](https://doi.org/10.31223/X51B4P)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18225226.svg)](https://doi.org/10.5281/zenodo.18225226)

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

### Examples

The `examples/` folder contains the following:

- [radet_single_image.ipynb](examples/radet_single_image.ipynb) — Compute RADET for a single Landsat image
- [radet_collection_interpolate.ipynb](examples/radet_collection_interpolate.ipynb) — Build a RADET image collection and interpolate
- [runtime_comparison.ipynb](examples/runtime_comparison.ipynb) — Compare runtimes of OpenET models
- [eecu_analysis.py](examples/eecu_analysis.py) — Analyze Earth Engine Compute Unit (EECU) usage across OpenET models

## Project Structure

```
radet-beta/
├── radet/
│   ├── __init__.py
│   ├── collection.py      # ET image collection builder
│   ├── image.py           # Core Image class for single ET computation
│   ├── interpolate.py     # Temporal interpolation utilities
│   ├── landsat.py         # Landsat-specific preprocessing
│   ├── model.py           # RADET model implementation
│   └── utils.py           # Helper functions
├── examples/
│   ├── README.md
│   ├── radet_single_image.ipynb
│   ├── radet_collection_interpolate.ipynb
│   ├── runtime_comparison.ipynb
│   ├── eecu_analysis.py
│   ├── eecu_data/         # Raw EECU input data
│   └── eecu_output/       # Generated analysis results and plots
├── .gitignore
├── LICENSE
└── README.md
```

## Dependencies

- [earthengine-api](https://github.com/google/earthengine-api) # main RADET model dependency
- [openet-core](https://pypi.org/project/openet-core/) # main RADET model dependency
- [openet-sims](https://pypi.org/project/openet-sims/) # Only for runtime comparisons ([runtime_comparison.ipynb](examples/runtime_comparison.ipynb))
- [openet-ssebop](https://pypi.org/project/openet-ssebop/) # Only for runtime comparisons ([runtime_comparison.ipynb](examples/runtime_comparison.ipynb))
- [openet-ptjpl](https://pypi.org/project/openet-ptjpl/) # Only for runtime comparisons ([runtime_comparison.ipynb](examples/runtime_comparison.ipynb))
- [openet-geesebal](https://pypi.org/project/openet-geesebal/) # Only for runtime comparisons ([runtime_comparison.ipynb](examples/runtime_comparison.ipynb))
- [openet-disalexi](https://pypi.org/project/openet-disalexi/) # Only for runtime comparisons ([runtime_comparison.ipynb](examples/runtime_comparison.ipynb))
- [pandas](https://pypi.org/project/pandas/) # For analysis scripts
- [seaborn](https://seaborn.pydata.org/) # For analysis scripts

## Installation

```
pip install earthengine-api openet-core pandas seaborn openet-sims openet-ptjpl openet-ssebop openet-disalexi openet-geesebal
```

### Google Earth Engine Authentication

This project uses the Google Earth Engine (GEE) Python API for geospatial data extraction.

1. Install [Google Cloud CLI](https://cloud.google.com/sdk/docs/install-sdk)
2. Create a GCloud project (e.g., `gee-radet`) with GEE API enabled at https://console.cloud.google.com/
3. Configure the project:
   ```bash
   gcloud config set project gee-radet
   gcloud auth application-default set-quota-project gee-radet  # if prompted
   earthengine authenticate
   ```

See the [Earth Engine Python installation guide](https://developers.google.com/earth-engine/guides/python_install) for details.

## References

Kim, Y., Huntington, J. L., Comini de Andrade, B., Johnson, M. S., Volk, J. M., Majumdar, S., Morton, C., & ReVelle, P. (2026). Thermodynamically constrained closed-form surface energy balance using medium-resolution remote sensing for efficient evapotranspiration mapping. *EarthArXiv (preprint)*. [https://doi.org/10.31223/X51B4P](https://doi.org/10.31223/X51B4P)
