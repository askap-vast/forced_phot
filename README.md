# forced_phot
Forced radio photometry (for ASKAP)

Contents:
* `forced_phot.py` main code containing classes and methods.  Usage: 
```python
from astropy.table import Table
from astropy import units as u, constants as c
import numpy as np
from astropy.coordinates import SkyCoord
import pandas as pd
import time

import forced_phot
    
# read in a selavy catalog with pandas 
df=pd.read_fwf('selavy-image.i.SB9668.cont.VAST_0341-50A.linmos.taylor.0.restored.islands.txt',skiprows=[1,])

# and convert to astropy Table for easier handling
data_islands=Table.from_pandas(df)
# construct a SkyCoord object from the sources
P_islands=SkyCoord(data_islands['ra_deg_cont']*u.deg,data_islands['dec_deg_cont']*u.deg)

# image, background, and noise maps from ASKAPSoft
image='image.i.SB9668.cont.VAST_0341-50A.linmos.taylor.0.restored.fits'
background='meanMap.image.i.SB9668.cont.VAST_0341-50A.linmos.taylor.0.restored.fits'
noise='noiseMap.image.i.SB9668.cont.VAST_0341-50A.linmos.taylor.0.restored.fits'

# make the Forced Photometry object
FP=forced_phot.ForcedPhot(image, background, noise)

# run the forced photometry
flux_islands,flux_err_islands,chisq_islands,DOF_islands=FP.measure(P_islands,
                                                                   data_islands['maj_axis']*u.arcsec, data_islands['min_axis']*u.arcsec, data_islands['pos_ang']*u.deg,
                                                                   cluster_threshold=3)
```
* `test_forced_phot_inject.py`: test photometry via injection
* `test_clusterphot.py`: test photometry of close source pairs via injection and turning cluster fitting on/off (this also has a couple of plots associated with it [`injection_1.0arcsec.pdf`](injection_1.0arcsec.pdf), [`injection_statistics.pdf`](injection_statistics.pdf)

The cluster fitting was tested with a variety of thresholds, and we found that if the sources were closer than ~1.0 BMAJ there was a significant bias in not using the cluster fitting.  So the default threshold has been set to 1.5.
