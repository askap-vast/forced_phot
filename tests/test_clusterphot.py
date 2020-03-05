import time
import warnings
from astropy import units as u, constants as c
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.wcs
from astropy.io import fits
from astropy.io.fits.verify import VerifyWarning
from astropy.utils.exceptions import AstropyWarning
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

import forced_phot

# suppress FITS verification warnings
warnings.simplefilter("ignore", category=AstropyWarning)

image = "image.i.SB9668.cont.VAST_0341-50A.linmos.taylor.0.restored.fits"
background = "meanMap.image.i.SB9668.cont.VAST_0341-50A.linmos.taylor.0.restored.fits"
noise = "noiseMap.image.i.SB9668.cont.VAST_0341-50A.linmos.taylor.0.restored.fits"


FP = forced_phot.ForcedPhot(image, background, noise)

separations = np.logspace(0, 2, 15) * u.arcsec
n = 20
chi2_inj_recover = np.zeros(len(separations))
chi2_inj_cluster = np.zeros(len(separations))


for j in range(len(separations)):
    s = separations[j]
    t = time.time()

    x = (np.random.random_sample((n,)) - 0.5) * 8000 + 7046.5
    y = (np.random.random_sample((n,)) - 0.5) * 8000 + 7046.5
    flux_inj = np.random.random_sample((n,)) * 100e-3 + 0.5e-3

    # add another source nearby to each main source
    theta = np.random.random_sample((n,)) * 2 * np.pi * u.rad
    x2 = x + ((s / FP.pixelscale).decompose() * np.cos(theta)).value
    y2 = y + ((s / FP.pixelscale).decompose() * np.sin(theta)).value

    x = np.r_[x, x2]
    y = np.r_[y, y2]
    P_inj = astropy.wcs.utils.pixel_to_skycoord(x, y, FP.w)
    flux_inj = np.r_[flux_inj, flux_inj + np.random.random_sample((n,)) * 0.5e-3]

    # inject with a wider kernel than recovery
    FP.inject(flux_inj, P_inj, nbeam=15)
    # recover assuming no clusters
    flux_recover, flux_err_recover, chi2_recover, DOF_recover = FP.measure(
        P_inj, cluster_threshold=None
    )
    # recover assuming  clusters
    flux_cluster, flux_err_cluster, chi2_cluster, DOF_cluster = FP.measure(
        P_inj, cluster_threshold=2 * (s / FP.pixelscale).decompose().value
    )

    chi2_inj_recover[j] = (((flux_recover - flux_inj) / flux_err_recover) ** 2).sum()
    chi2_inj_cluster[j] = (((flux_cluster - flux_inj) / flux_err_cluster) ** 2).sum()

    print(time.time() - t)

    if j == 0:
        plt.clf()
        plt.errorbar(
            flux_inj, (flux_recover - flux_inj), yerr=flux_err_recover, fmt="o"
        )
        plt.errorbar(
            flux_inj, (flux_cluster - flux_inj), yerr=flux_err_cluster, fmt="s"
        )
        plt.plot([0, flux_inj.max()], [0, 0], "k--")
        plt.xlabel("Injected flux density (Jy)")
        plt.ylabel("(Recovered - injected (Jy)")
        plt.title(r"$\chi^2=%.1f$ (%d DOF)" % (chi2_recover[j], 2 * n))
        plt.savefig("injection_%.1farcsec.pdf" % s.value)

plt.clf()
plt.plot(separations, chi2_inj_recover / 2 / n, "o", label="single source")
plt.plot(separations, chi2_inj_cluster / 2 / n, "s", label="cluster")
plt.plot(plt.gca().get_xlim(), np.array([1, 1]), "k--", label=r"$\chi^2/$DOF=1")
plt.plot(
    FP.BMAJ.to(u.arcsec) * np.array([1, 1]), plt.gca().get_ylim(), ":", label="BMAJ"
)
plt.xlabel("Source separation (arcsec)")
plt.ylabel(r"$\chi^2$ for recovered vs. injected")
plt.legend()

plt.savefig("injection_statistics.pdf")
