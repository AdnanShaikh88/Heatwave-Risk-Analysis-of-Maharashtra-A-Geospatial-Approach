# Heatwave Risk Analysis of Maharashtra: A Geospatial Approach

**Term Project for ES 216**
**Group 2:** Shaikh Adnan Matiyoddin (23b4207), Anika Saanvi Sahoo (23b4227), Devanshi Khandelwal (23b4233)
**Supervisor:** Prof. Srinidhi Balasubramanian

---

### Project Overview

This project presents a high-resolution (1km) geospatial heatwave risk assessment for the state of Maharashtra, India. The analysis was conducted entirely within the **Google Earth Engine (GEE)** platform for the 2023 pre-monsoon season (March-May).

The objective is to identify and rank all districts in Maharashtra based on their heatwave risk, providing actionable data for policymakers. The final report, `Heatwave_Risk_Analysis_Maharashtra_A...pdf`, is included in this repository.

### Methodology

The analysis follows the standard climate risk framework:
**`Risk = Hazard × Exposure × Vulnerability`**

#### 1. Hazard (H)
* **Data:** ERA5-Land (11km native resolution)
* **Metric:** A health-based **Heat Index (HI)** was calculated using temperature and dewpoint.
* **Process:**
    1.  A 30-year (1991-2020) "Decadal Baseline" was created by calculating the 90th percentile HI for each month, averaged across three decades.
    2.  The number of days in 2023 where the daily HI exceeded this baseline was counted, creating a coarse 11km `hazardImage_coarse`.
    3.  This 11km image was **interpolated to 1km** using `bicubic` resampling to match the other layers.

#### 2. Exposure (E)
* **Data:** GPWv4 (Gridded Population of the World, v4)
* **Metric:** Population density at 1km resolution. This serves as a direct measure of the population exposed to the heat hazard.

#### 3. Vulnerability (V)
* **Data:** MODIS (1km resolution)
* **Metric:** A composite index was created by combining two proxies:
    1.  **Environmental Vulnerability:** Inverted NDVI (`MOD13A3`), where low vegetation = high vulnerability.
    2.  **Urban Heat Island (UHI) Vulnerability:** Nighttime Land Surface Temperature (`MYD11A1`), where high night temps = high vulnerability.
* **Weighting:** `Vulnerability = (v_env * 0.4) + (v_uhi * 0.6)`

#### 4. Final Risk Index
1.  All three 1km layers (Hazard, Exposure, Vulnerability) were normalized to a 0-1 scale.
2.  The final 1km risk map was calculated using the multiplicative model:
    `Risk = normalizedHazard × normalizedExposure × normalizedVulnerability`
3.  Zonal statistics were run to find the mean risk score for every district in Maharashtra.

### Repository Contents

* `Heatwave_Risk_Analysis_Maharashtra_A...pdf`: The final project report in PDF format.
* `Heatwave_Risk_Maharashtra_v2.js`: The **final 1km "smooth" GEE script** that was used for the analysis in the report (this is the one that uses interpolation).
* `Heatwave_code_v1.js`: The initial **11km "blocky" GEE script** that was tested. This version is statistically robust but was not used for the final maps due to computational errors in the GEE console.
* `README.md`: This file.

### How to Use the GEE Scripts
1.  Go to [code.earthengine.google.com](https://code.earthengine.google.com/).
2.  Copy the contents of `Heatwave_Risk_Maharashtra_v2.js`.
3.  Paste the code into a new, blank script in the GEE Code Editor.
4.  Click **"Run"**.
5.  The script will display all map layers (Hazard, Exposure, Vulnerability, and Final Risk) and add legends.
6.  The **Console** tab will print the ranked list of all 35 districts.
7.  The **Tasks** tab will be populated with export tasks to save the final map (GeoTIFF) and the district list (CSV) to your Google Drive. You must click "Run" on these tasks manually.
