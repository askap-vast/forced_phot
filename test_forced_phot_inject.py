from astropy.table import Table
from astropy import units as u, constants as c
import numpy as np
from astropy.coordinates import SkyCoord
import pandas as pd
import time
import astropy.wcs
from astropy.io import fits
from astropy.io.fits.verify import VerifyWarning
from astropy.utils.exceptions import AstropyWarning
import warnings
# suppress FITS verification warnings
warnings.simplefilter('ignore', category=AstropyWarning)

from matplotlib import pyplot as plt

import forced_phot


    
image='image.i.SB9668.cont.VAST_0341-50A.linmos.taylor.0.restored.fits'
background='meanMap.image.i.SB9668.cont.VAST_0341-50A.linmos.taylor.0.restored.fits'
noise='noiseMap.image.i.SB9668.cont.VAST_0341-50A.linmos.taylor.0.restored.fits'


FP=forced_phot.forced_phot(image, background, noise)

n=500
t=time.time()

x=(np.random.random_sample((n,))-0.5)*8000+7046.5
y=(np.random.random_sample((n,))-0.5)*8000+7046.5
P_inj=astropy.wcs.utils.pixel_to_skycoord(x, y, FP.w)
flux_inj=np.random.random_sample((n,))*100e-3+0.5e-3

# inject with a wider kernel than recovery
FP.inject(flux_inj, P_inj, nbeam=15)
flux_recover,flux_err_recover=FP.measure(P_inj,cluster_threshold=None)

print(time.time()-t)
plt.clf()
plt.errorbar(flux_inj,(flux_recover-flux_inj),yerr=flux_err_recover,fmt='o')
plt.plot([0,flux_inj.max()],[0,0],'k--')
plt.xlabel('Injected flux density (Jy)')
plt.ylabel('(Recovered - injected (Jy)')
plt.title('$\chi^2=%.1f$ (%d DOF)' % ((((flux_recover-flux_inj)/flux_err_recover)**2).sum(),
                                     n))

