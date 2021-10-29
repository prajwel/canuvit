# **CanUVIT**
> To check whether a field can be safely observed with UVIT.

<p align="center">
<img src="https://i.imgur.com/b0hoB04.png" width="400"/>
</p>
You can install the "CanUVIT" Python package using the following command.

```bash
pip install canuvit
``` 

After installation, you can run `canuvit` on a Python command prompt or as a script. For example, 
```python
>>> import canuvit
>>> canuvit.observe('uvit', '12:12:12', '12:12:12')
```

> **Note:** In general, `canuvit.observe(instrument, RA, DEC)`. `instrument` can be either 'uvit', 'sxt', 'czti', or 'laxpc'. `RA` and `DEC` are expected to be in sexagesimal format. 

You may get an output as shown below. Also, please check the working directory for the output GALEX images with sources marked. 

```
### VIS

ra_hms	dec_dms	mag	B-V	SpecType	VIS3	VIS2	VIS1	ND1	BK7

12:11:52.7568	+12:07:47.532	11.096	0.864	K1	1333.0	124.4	88.6	29.7	1624.9
12:12:22.944	+12:17:23.856	11.126	0.814	K0	1296.7	121.0	86.2	28.9	1580.6
12:11:35.016	+12:12:04.644	11.426	0.451	F5	1457.4	234.9	180.9	32.8	1915.2
12:11:01.656	+12:08:35.916	11.874	0.694	G5	754.3	89.1	69.2	16.8	950.9
12:11:11.5272	+12:03:14.04	12.177	0.322	F0	803.1	143.5	101.4	18.2	1061.7
12:12:05.5368	+12:19:09.768	12.268	0.787	K0	452.9	42.3	30.1	10.1	552.1


Safe filters: ['VIS3', 'VIS2', 'VIS1', 'ND1', 'BK7']


FUV observations seem to be absent! Using M_fuv = M_nuv - 1.65.


### NUV

sl_no     ra_hms       dec_dms      Mag  Mag_corrected silica  b4  b13  b15   n2 
----- ------------- -------------- ----- ------------- ------ ---- ---- ---- ----
    1 12:12:32.3946 +12:07:27.4144 19.32         19.32   1.86 0.41 0.50 0.14 0.10
    2 12:11:11.6503 +12:03:14.7794 16.25         16.25  31.75 6.99 8.57 2.35 1.75
    3 12:12:41.0882 +12:14:58.2679 16.15         16.15  34.66 7.63 9.36 2.56 1.91
    4 12:12:15.3493 +12:29:18.1277 19.50         19.50   1.59 0.35 0.43 0.12 0.09
    5 12:11:35.0116 +12:12:04.7063 16.47         16.47  25.93 5.70 7.00 1.92 1.43



Safe filters in NUV: ['Silica', 'NUV-grating', 'NUV-B4', 'NUV-B13', 'NUV-B15', 'NUV-N2']


### FUV 

sl_no     ra_hms       dec_dms      Mag  Mag_corrected  caf2  baf2 sapphire silica
----- ------------- -------------- ----- ------------- ----- ----- -------- ------
    1 12:12:32.3946 +12:07:27.4144 19.32         17.67  1.65  1.40     1.04   0.36
    2 12:11:11.6503 +12:03:14.7794 16.25         14.60 28.17 23.94    17.75   6.20
    3 12:12:41.0882 +12:14:58.2679 16.15         14.50 30.75 26.14    19.37   6.76
    4 12:12:15.3493 +12:29:18.1277 19.50         17.85  1.41  1.20     0.89   0.31
    5 12:11:35.0116 +12:12:04.7063 16.47         14.82 23.00 19.55    14.49   5.06




Safe filters in FUV: ['CaF2', 'FUV-grating', 'BaF2', 'Sapphire', 'Silica']

```

Two additional functions are also available which takes the same input arguments as `canuvit.observe()`.

* `canuvit.observe_VIS()`: to find safe VIS filters.
* `canuvit.observe_UV()`: to find safe UV filters.