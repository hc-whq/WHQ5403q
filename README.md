#  WRF-Hydro® <img src=".github/images/wrf_hydro_symbol_logo_2017_09.png" width=100 align="left" />

[![Build Status](https://github.com/NCAR/wrf_hydro_nwm_public/actions/workflows/test-pr.yml/badge.svg?branch=main)](https://github.com/NCAR/wrf_hydro_nwm_public/actions/workflows/test-pr.yml)
[![Doc Build Status](https://readthedocs.org/projects/wrf-hydro/badge/?version=latest&style=flat)](https://wrf-hydro.readthedocs.io)
[![Release](https://img.shields.io/github/release/NCAR/wrf_hydro_nwm_public.svg)](https://github.com/NCAR/wrf_hydro_nwm_public/releases/latest)
[![DOI](https://zenodo.org/badge/121802384.svg)](https://zenodo.org/badge/latestdoi/121802384)

## Description
This is the code repository for [WRF-Hydro®](https://ral.ucar.edu/projects/wrf_hydro).

WRF-Hydro is a community modeling system and framework for hydrologic modeling and model coupling.  In 2016 a configuration of WRF-Hydro was implemented as the [National Water Model](http://water.noaa.gov/about/nwm) (NWM) for the continental United States.

## Documentation
Documentation is available at [wrf-hydro.readthedocs.io](https://wrf-hydro.readthedocs.io) with further guides available at our [project website](https://ral.ucar.edu/projects/wrf_hydro/technical-description-user-guide). For getting started see the [dependencies](https://wrf-hydro.readthedocs.io/en/latest/model-code-config.html#requirements) and [build instructions](https://wrf-hydro.readthedocs.io/en/latest/model-code-config.html#cmake-build) in the documentation.

## WRF-HydroQual (SedCNP) notes
The water quality model (SedCNP: sediment, C, N, P) runs after the hydro model
(`SedCNP_option = 3` in `namelist.hrldas`) and reads its outputs from `LDASOUT_dir`.

**Lateral loads** (`SedCNP_lateral_option` in `whq.namelist`, section `SedCNPmodel_OFFLINE`)
- `1` (default): surface runoff and interflow loads follow the cell water paths of the hydro
  model; groundwater loads are divided among the stream pixels of each gw basin as in the hydro
  model. Requires `RTOUT_DOMAIN = 1`, `rt_option = 1` and `AGGFACTRT = 1` in `hydro.namelist`,
  and an RTOUT written by this code version (cell fluxes `sfc_in`, `sfc_out`, `sfc_dir`,
  `sfc_chan`, `sfc_lake`, `sfc_bdry`, `sfc_rem`, `sub_dir`, `sub_exfil`).
- `0`: previous method, summing the loads of all cells of the gw basin of each link.
- A per-step load budget is written to `debug/lateral_budget.csv`.

**Lakes/reservoirs** (level-pool lakes, `lake_option = 1`)
- Each lake is simulated as a completely mixed reactor at its outlet link. Requires the variable
  `BottomE` (lake bottom elevation, m) in `LAKEPARM.nc` and `outlake = 1` (LAKEOUT files).
- Outputs `LAKESEDOUT`, `LAKECOUT`, `LAKENOUT`, `LAKEPOUT` in `WHQOUT_dir`: lake inflow,
  storage, outflow and the CHSEDOUT/CHCOUT/CHNOUT/CHPOUT variables at the lake outlet links.

**Restart (initial condition) files** (`whq.namelist`, NetCDF)
- Output: `<WHQOUT_dir>/RESTART_<kind>.YYYYMMDDHH_DOMAIN1.nc` for kind = `SOC`, `SON`, `SOP`,
  `GWC`, `GWN`, `GWP`, `CHC`, `CHN`, `CHP`, every `WHQ_RESTART_DT` hours of model time
  (<= 0 or omitted: only at the end of the simulation) and at the end of the simulation.
  YYYYMMDDHH is the valid time of the state (soil kg ha-1, groundwater and channels kg).
- Input: `RESTART_FILENAME_SOC = './RESTART_WHQ/RESTART_SOC.2020083116_DOMAIN1.nc'` etc.; start the
  run with `SedCNP_START_*` = the valid time of the files. A kind without a file name falls back to
  the text files (`CiniSo_file`, ...) or a cold start.

## Resources and Support
For news and updates regarding the WRF-Hydro project please subscribe to our [email list](https://ral.ucar.edu/projects/wrf_hydro/subscribe).

For user support and general inquiries please use our [contact form](https://ral.ucar.edu/projects/wrf_hydro/contact).

If you have found a bug or would like to propose changes to the model code please refer to our [contributing guidelines](.github/CONTRIBUTING.md).

## Contributions
For more information on how to contribute to this project please refer to our [contributing guidelines](.github/CONTRIBUTING.md).

## License
The license and terms of use for this software can be found [here](LICENSE.txt).
The Crocus snowpack module and related files are from the SURFEX surface model developed by Météo-France, the French national meteorological service.
These files are under the [CeCILL-C](http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html) license.

## Acknowledgements
Funding support for the development and application of the WRF-Hydro system has been provided by:
- The National Science Foundation and the National Center for Atmospheric Research
- The National Oceanic and Atmospheric Administration (NOAA) Office of Water Prediction (OWP)
- The U.S. National Weather Service (NWS)
- The U.S. Department of Energy (DOE)
- The Colorado Water Conservation Board
- Baron Advanced Meteorological Services
- National Aeronautics and Space Administration (NASA)
