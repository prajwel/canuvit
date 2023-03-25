# **CanUVIT**
> To check whether a field can be safely observed with UVIT.

<p align="center">
<img src="https://i.imgur.com/b0hoB04.png" width="400"/>
</p>


You can install the CanUVIT Python package using the following command.

```bash
pip install canuvit --upgrade
``` 
	
> **IMPORTANT:** Even if you have CanUVIT already installed, make sure you use the latest version by running the above command. Current version of CanUVIT is shown on the badge below: <br> <a href="https://pypi.org/project/canuvit/"><img src="https://img.shields.io/pypi/v/canuvit?style=for-the-badge"/></a> <br>

> **Note:** If you don't want to install CanUVIT or are facing problems during installation, You can run it online using Binder: 
<br> [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/prajwel/UVIT_notebooks/main?labpath=notebook2_UVIT_VIS_UV_safety_check.ipynb)



After installation, you can run CanUVIT on a Python command prompt or as a script. For example, if your primary instrument is UVIT and the RA, DEC coordinates of the field are (12:12:12, 12:12:12),
you may run the CanUVIT package as follows.

```python
>>> import canuvit
>>> canuvit.observe('uvit', '12:12:12', '12:12:12')
```

> **Note:** In general, `canuvit.observe(instrument, RA, DEC)` where the `instrument` can be either 'uvit', 'sxt', 'czti', or 'laxpc' and `RA` and `DEC` field coordinates should be in sexagesimal format.

For the above example, you should get an output as shown below. Please also check the working directory for the output GALEX images with sources marked in the primary instrument field of view.

```
=========================================================
Payload: uvit, Coordinates: 12:12:12, 12:12:12
=========================================================

### VIS

sl_no   ra_hms     dec_dms   mag  B-V SpecType  VIS3   VIS2  VIS1 ND1   BK7  
----- ---------- ----------- ---- --- -------- ------ ----- ----- ---- ------
    1 12:11:52.8 +12:07:47.5 11.1 0.9       K1 1333.0 124.4  88.6 29.7 1624.9
    2 12:12:22.9 +12:17:23.9 11.1 0.8       K0 1296.7 121.0  86.2 28.9 1580.6
    3 12:11:35.0 +12:12:04.6 11.4 0.5       F5 1457.4 234.9 180.9 32.8 1915.2
    4 12:11:01.7 +12:08:35.9 11.9 0.7       G5  754.3  89.1  69.2 16.8  950.9
    5 12:11:11.5 +12:03:14.0 12.2 0.3       F0  803.1 143.5 101.4 18.2 1061.7
    6 12:12:05.5 +12:19:09.8 12.3 0.8       K0  452.9  42.3  30.1 10.1  552.1

Safe filters: ['VIS3', 'VIS2', 'VIS1', 'ND1', 'BK7']

FUV observations seem to be absent! Using M_fuv = M_nuv - 1.65.

### NUV

sl_no   ra_hms     dec_dms   Mag  Mag_corrected silica  b4 b13 b15  n2
----- ---------- ----------- ---- ------------- ------ --- --- --- ---
    1 12:12:32.4 +12:07:27.4 19.3          19.3    1.9 0.4 0.5 0.1 0.1
    2 12:11:11.7 +12:03:14.8 16.2          16.2   31.8 7.0 8.6 2.3 1.7
    3 12:12:41.1 +12:14:58.3 16.2          16.2   34.7 7.6 9.4 2.6 1.9
    4 12:12:15.3 +12:29:18.1 19.5          19.5    1.6 0.3 0.4 0.1 0.1
    5 12:11:35.0 +12:12:04.7 16.5          16.5   25.9 5.7 7.0 1.9 1.4

Safe filters in NUV: ['Silica', 'NUV-grating', 'NUV-B4', 'NUV-B13', 'NUV-B15', 'NUV-N2']

### FUV

sl_no   ra_hms     dec_dms   Mag  Mag_corrected caf2 baf2 sapphire silica
----- ---------- ----------- ---- ------------- ---- ---- -------- ------
    1 12:12:32.4 +12:07:27.4 19.3          17.7  1.7  1.4      1.0    0.4
    2 12:11:11.7 +12:03:14.8 16.2          14.6 28.2 23.9     17.7    6.2
    3 12:12:41.1 +12:14:58.3 16.2          14.5 30.7 26.1     19.4    6.8
    4 12:12:15.3 +12:29:18.1 19.5          17.8  1.4  1.2      0.9    0.3
    5 12:11:35.0 +12:12:04.7 16.5          14.8 23.0 19.6     14.5    5.1

Safe filters in FUV: ['CaF2', 'FUV-grating', 'BaF2', 'Sapphire', 'Silica']
```

Please choose the VIS filters such that none of the stars in the field gives >4800 counts/second. Further, for good tracking of the aspect, there should be at least two stars within a 12 arcminute radius of the target with >30 counts/second (for good S/N) and <1000 counts/second (to avoid saturation) in the chosen filter. Please also avoid configuring multiple VIS filters. 

Two additional functions are also available, which takes the same input arguments as `canuvit.observe()`.

* `canuvit.observe_VIS()`: to find safe VIS filters.
* `canuvit.observe_UV()`: to find safe UV filters.

### Command Line Interface

After installation with pip, you can also access CanUVIT from the command line. Here's an example:

```bash
canuvit -i uvit -r "12:12:12" -d "12:12:12"
```
The help page of the command-line tool can be accessed as follows:

```bash
canuvit -h
```

```
Usage: canuvit [OPTIONS]

  Program to check if a given coordinate can be safely observed using UVIT.

  Example usage:
  canuvit -r "13:12:14" -d "-14:15:13" 

Options:
  --all                           Check safety for all filters.  [default:
                                  all]
  --vis                           Check saftey for only visible filters.
  --uv                            Check safety for only UV filters.
  -r, --ra RA                     Right ascension of the coordinate. Format:
                                  hh:mm:ss[.ss] e.g. "00:54:53.45"  [required]
  -d, --dec DEC                   Declination of the coordinate. Format:
                                  [-]dd:mm:ss[.ss] e.g. "-37:41:03.23".
                                  [required]
  -i, --instrument [uvit|sxt|czti|laxpc]
                                  Instrument to check for.  [default: uvit]
  -v, --verbose                   Increase output verbosity.
  --version                       Show the version and exit.
  -h, --help                      Show this message and exit.
```

### Acknowledgements

CanUVIT depends on the following web tools and API for its functioning. 

* Bright Source Warning Tool (https://uvit.iiap.res.in/Software/bswt)
* Exposure Time Calculator (https://uvit.iiap.res.in/Software/etc)
* MAST API (https://mast.stsci.edu/api/v0/) 

### Requirements

CanUVIT works with Python 3.6 or later. CanUVIT depends on the following packages:

* astropy
* astroquery
* beautifulsoup4
* click
* matplotlib
* numpy
* requests

