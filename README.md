# Snowline Altitude Analysis in the Himalayas

This repository contains the exact code used for the analyses presented in the paper:

"Contrasting patterns of change in snowline altitude across five Himalayan catchments"
Sasaki, O., Miles, E., et al. (2025) – _The Cryosphere_  

---

## Overview

The scripts in this repository were used to extract and analyze snowline altitude (SLA) time series across five Himalayan catchments, using satellite imagery (Landsat 5/7/8 and Sentinel-2) via Google Earth Engine (GEE).

---

## Contents

- run_calcSLA.sh: Shell script to run the SLA extraction Python scripts across the five glacier basins.
- SLA_LS5.py: Python script for Extracting Snowline altitude (SLA) using Landsat 5 imagery and Google Earth Engine (GEE).
- SLA_LS7.py: Same as SLA_LS5.py, but for Landsat 7 data.
- SLA_LS8.py: Same as SLA_LS5.py, but for Landsat 8 data.
- SLA_S2.py: Same as SLA_LS5.py, but for Sentinel-2 data.
- data_5glaciers.csv: Latitude and longitude information for the five study glacier basins.

---

## Important Notes

- The code in this repository reflects the exact version used in the publication.
- **Some satellite datasets (e.g., older Landsat TOA collections) are no longer available on GEE**, and the scripts may not run without modification.
- We provide the code as a record of our processing steps and methodology for transparency and reproducibility.
- If you wish to re-run the analysis, we recommend updating the scripts to match currently available GEE datasets.

---

## Citation

If you use this code or build upon it, please cite the original paper:

> Sasaki, O., Miles, E., et al. (2025). Contrasting patterns of change in snowline altitude across five Himalayan catchments. *The Cryosphere*.
