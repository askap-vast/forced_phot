"""
Usage:

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
flux_islands, flux_err_islands, chisq_islands, DOF_islands = FP.measure(
    P_islands,
    data_islands['maj_axis']*u.arcsec,
    data_islands['min_axis']*u.arcsec,
    data_islands['pos_ang']*u.deg,
    cluster_threshold=3
)
"""

from itertools import chain
from typing import Any, Dict, List, Optional, Tuple, Union

import logging
import astropy
import astropy.nddata
import astropy.wcs
import numpy as np
import scipy.spatial
from astropy import units as u
from astropy.io import fits
from astropy.modeling import fitting, models
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales

import numba as nb
from numba import float64, njit
from numba.experimental import jitclass

from timeit import default_timer as timer

logger = logging.getLogger(__name__)

numba_logger = logging.getLogger('numba')
numba_logger.setLevel(logging.WARNING)

pa_offset = 90 * u.deg
pa_offset_deg = pa_offset.value


class ArgumentError(Exception):
    pass

class G2D:
    """2D Gaussian for use as a kernel.

    Example usage:
        create the kernel:
        g = G2D(x0, y0, fwhm_x, fwhm_y, PA)
        and return the kernel:
        g(x, y)

    Args:
        x0 (float): the mean x coordinate (pixels)
        y0 (float): the mean y coordinate (pixels)
        fwhm_x (float): the FWHM in the x coordinate (pixels)
        fwhm_y (float): the FWHM in the y coordinate (pixels)
        pa (float): the position angle of the Gaussian (E of N) as a Quantity or in
            radians.
    """

    def __init__(self, x0: float, y0: float, fwhm_x: float, fwhm_y: float, pa: float):
        self.x0 = x0
        self.y0 = y0
        self.fwhm_x = fwhm_x
        self.fwhm_y = fwhm_y
        # adjust the PA to agree with the selavy convention
        # E of N
        self.pa = pa - pa_offset
        self.sigma_x = self.fwhm_x / 2 / np.sqrt(2 * np.log(2))
        self.sigma_y = self.fwhm_y / 2 / np.sqrt(2 * np.log(2))

        self.a = (
            np.cos(self.pa) ** 2 / 2 / self.sigma_x ** 2
            + np.sin(self.pa) ** 2 / 2 / self.sigma_y ** 2
        )
        self.b = (
            np.sin(2 * self.pa) / 2 / self.sigma_x ** 2
            - np.sin(2 * self.pa) / 2 / self.sigma_y ** 2
        )
        self.c = (
            np.sin(self.pa) ** 2 / 2 / self.sigma_x ** 2
            + np.cos(self.pa) ** 2 / 2 / self.sigma_y ** 2
        )

    def __call__(self, x: float, y: float) -> np.ndarray:
        """Return the kernel evaluated at given pixel coordinates.

        Args:
            x (float): x coordinate for evaluation
            y (float): y coordinate for evaluation

        Returns:
            np.ndarray: the kernel evaluated at the given coordinates
        """
        return np.exp(
            -self.a * (x - self.x0) ** 2
            - self.b * (x - self.x0) * (y - self.y0)
            - self.c * (y - self.y0) ** 2
        )

specs = (
    nb.int64[:,:],
    nb.int64[:,:],
    float64,
    float64,
    float64,
    float64,
    float64
)
@njit(specs)
def get_kernel(
    x: np.ndarray, 
    y: np.ndarray,
    x0: float,
    y0: float,
    fwhm_x: float,
    fwhm_y: float,
    pa: float
    ):
    """Calculate 2D Gaussian kernel, optimised with numba

    Args:
        x0 (float): the mean x coordinate (pixels)
        y0 (float): the mean y coordinate (pixels)
        fwhm_x (float): the FWHM in the x coordinate (pixels)
        fwhm_y (float): the FWHM in the y coordinate (pixels)
        pa (float): the position angle of the Gaussian (E of N) in deg
    """
    
    pa = np.deg2rad(pa - pa_offset_deg)
    
    sigma_fac = 0.5/np.sqrt(2 * np.log(2))
    
    
    sigma_x = fwhm_x * sigma_fac
    sigma_y = fwhm_y * sigma_fac
    
    sigma_x_fac = 0.5 / sigma_x ** 2
    sigma_y_fac = 0.5 /sigma_y ** 2

    a = (
        np.cos(pa) ** 2 * sigma_x_fac
        + np.sin(pa) ** 2 * sigma_y_fac 
    )
    b = (
        np.sin(2 * pa) * sigma_x_fac
        - np.sin(2 * pa) * sigma_y_fac 
    )
    c = (
        np.sin(pa) ** 2 * sigma_x_fac
        + np.cos(pa) ** 2 * sigma_y_fac 
    )
    
    xdiff = x-x0
    ydiff = y-y0

    return np.exp(
        -a * (xdiff) ** 2
        - b * (xdiff) * (ydiff)
        - c * (ydiff) ** 2
    )

@njit
def _convolution(d, n, kernel):
    """
    Calculate the convolution of the data and noise with the kernel using numba.
    
    Args:
        d: numpy array containing the background-subtracted data
        n: numpy array containing the noise map
        k: numpy array containing the kernel
    
    Returns:
        Flux, error and chi-squared
    """
    
    kernel_n2 = kernel / n ** 2

    flux = ((d) * kernel_n2).sum() / (kernel ** 2 / n ** 2).sum()
    flux_err = ((n) * kernel_n2).sum() / (kernel_n2).sum()
    chisq = (((d - kernel * flux) / n) ** 2).sum()

    return flux, flux_err, chisq

@njit
def _meshgrid(xmin, xmax, ymin, ymax):
    """
    Generate a 2x2 meshgrid as required.
    
    With numba, this specific code is faster than the generalised np.meshgrid
    """
    
    x = np.arange(xmin, xmax)
    y = np.arange(ymin, ymax)
    
    xx = np.empty(shape=(y.size, x.size), dtype=x.dtype)
    yy = np.empty_like(xx)
    for i in range(x.size):
        for j in range(y.size):
            xx[j,i] = x[i]
            yy[j,i] = y[j]
    return xx, yy

class ForcedPhot:
    """Create a ForcedPhotometry object for processing an ASKAPSoft image.

    Example usage:
        forced_phot_obj = ForcedPhot(image, background, noise)
        flux_islands, flux_err_islands, chisq_islands, dof_islands = forced_phot_obj.measure(islands)

        where `islands` is an array `astropy.coordinates.SkyCoord` objects.

    Args:
        image (Union[str, fits.HDUList]): name of the primary image or FITS handle.
        background (Union[str, fits.HDUList]): name of the background image or FITS handle.
        noise (Union[str, fits.HDUList]): name of the noise map image or FITS handle.
        verbose (bool, optional): whether to be verbose in output. Defaults to False.
        use_numba (bool, optional): whether to use numba, which may be slower for a single fit. Default to False.

    Raises:
        ArgumentError: an input type is not a supported.
        FileNotFoundError: an input could not be opened.
        KeyError: could not get required header info from the image
    """

    def __init__(
        self,
        image: Union[str, fits.HDUList],
        background: Union[str, fits.HDUList],
        noise: Union[str, fits.HDUList],
        verbose: Optional[bool] = False,
        use_numba: Optional[bool] = False,
    ):
        self.verbose = verbose
        self.use_numba = use_numba
        
        if isinstance(image, str):
            try:
                fi = fits.open(image)
            except FileNotFoundError:
                logger.exception("Unable to open image %s", image)
        elif isinstance(image, fits.HDUList):
            fi = image
        else:
            raise ArgumentError("Do not understand input image")
        if isinstance(background, str):
            try:
                fb = fits.open(background)
            except FileNotFoundError:
                logger.exception("Unable to open background image %s", background)
        elif isinstance(background, fits.HDUList):
            fb = background
        else:
            raise ArgumentError("Do not understand input background image")
        if isinstance(noise, str):
            try:
                fn = fits.open(noise)
            except FileNotFoundError:
                logger.exception("Unable to open noise image %s", noise)
        elif isinstance(noise, fits.HDUList):
            fn = noise
        else:
            raise ArgumentError("Do not understand input noise image")
        
        img_header = fi[0].header
        if not (
            ("BMAJ" in img_header.keys())
            and ("BMIN" in img_header.keys())
            and ("BPA" in img_header.keys())
        ):

            raise KeyError("Image header does not have BMAJ, BMIN, BPA keywords")

        self.BMAJ = img_header["BMAJ"] * u.deg
        self.BMIN = img_header["BMIN"] * u.deg
        self.BPA = img_header["BPA"] * u.deg

        self.NAXIS1 = img_header["NAXIS1"]
        self.NAXIS2 = img_header["NAXIS2"]
        
        self.data = (fi[0].data-fb[0].data).squeeze()
        self.noisedata = fn[0].data.squeeze()
        
        # do this so that 0-regions in the background
        # or noise maps are set to NaN, and then
        # will be handled through other means
        self.noisedata[self.noisedata == 0] = np.nan
        
        if self.use_numba:
            self.data = self._byte_reorder(self.data)
            self.noisedata = self._byte_reorder(self.noisedata)

        self.w = WCS(img_header, naxis=2).celestial
        self.pixelscale = (proj_plane_pixel_scales(self.w)[1] * u.deg).to(u.arcsec)

    def _byte_reorder(self, arr):
        """Convert an array to use native byte ordering as required by numba
        
        Args:
            arr (np.ndarray): the array
            
        Returns:
            The array, with reordering
        """
        if arr.dtype.byteorder != '=':
            logger.warning("Numba requires native byte ordering. Converting.")
            arr = arr.astype(arr.dtype.newbyteorder('='))
        return arr

    def cluster(self, X0: np.ndarray, Y0: np.ndarray, threshold: Optional[float] = 1.5):
        """Identifies clusters among the given X, Y points that are within `threshold` * BMAJ
            of each other using a KDTree algorithm. Results are stored in `self.clusters`
            and `self.in_cluster`:
                - `self.clusters` is a dict mapping cluster indices to a set of their members.
                - `self.in_cluster` is a list of all of the sources in a cluster

        Args:
            X0 (np.ndarray): array of X coordinates of sources.
            Y0 (np.ndarray): array of Y coordinates of sources.
            threshold (float, optional): multiple of BMAJ for finding clusters.
                Set to 0 or None to disable. Defaults to 1.5.
        """
        self.clusters: Dict[int, set]
        self.in_cluster: List[int]
        if threshold is None or threshold == 0:
            self.clusters = {}
            self.in_cluster = []
            return

        threshold_pixels = threshold * (self.BMAJ / self.pixelscale).decompose().value
        t = scipy.spatial.KDTree(np.c_[X0, Y0])

        # this is somewhat convoluted
        # ideally the KDTree query would work on its own
        # but we want to make sure that if sources 3,4,5 should be grouped,
        # they will be grouped no matter whether we query first for 3, 4, or 5
        # but that they are only a single cluster
        self.clusters = {}
        for i in range(len(X0)):
            dists, indices = t.query(
                np.c_[X0[i], Y0[i]], k=10, distance_upper_bound=threshold_pixels
            )
            indices = indices[~np.isinf(dists)]
            if len(indices) > 1:
                # too close to another source: do a simultaneous fit
                n = np.isin(indices, list(self.clusters.keys()))
                if np.any(n):
                    j = indices[n][0]
                    for k in indices:
                        self.clusters[j].add(k)
                else:
                    self.clusters[i] = set(indices)
        self.in_cluster = sorted(
            list(chain.from_iterable([*self.clusters.values()]))
        )

    def measure(
        self,
        positions: "astropy.coordinates.SkyCoord",
        major_axes: Optional["astropy.coordinates.Angle"] = None,
        minor_axes: Optional["astropy.coordinates.Angle"] = None,
        position_angles: Optional["astropy.coordinates.Angle"] = None,
        nbeam: int = 3,
        cluster_threshold: Optional[float] = 1.5,
        allow_nan: bool = True,
        stamps: bool = False,
        edge_buffer: float = 1.0,
        use_clusters: bool = True,
    ) -> Union[Tuple[Any, Any, Any, Any, Any], Tuple[Any, Any, Any, Any, Any, Any, Any]]:
        """Perform the forced photometry returning the flux density and uncertainty.
        Example usage:
            flux, flux_err, chisq, dof, cluster_id = forced_phot_obj.measure(positions, nbeam=3)

            or

            flux, flux_err, chisq, dof, cluster_id, data, model = forced_phot_obj.measure(
                positions, nbeam=3, allow_nan=True, stamps=True)

        Args:
            positions: Coordinates of sources to measure.
            major_axes: FWHMs along major axes of sources to measure. If None, will use
                header BMAJ. Defaults to None.
            minor_axes: FWHMs along minor axes of sources to measure. If None, will use
                header BMIN. Defaults to None.
            position_angles: Position angles of sources to measure. If None, will use
                header BPA. Defaults to None.
            nbeam: Diameter of the square cutout for fitting in units of
                the major axis. Defaults to 3.
            cluster_threshold: Multiple of `major_axes` to use for identifying clusters.
                Set to 0 or None to disable. Defaults to 1.5.
            allow_nan: whether or not to try to measure sources where some RMS values may
                be NaN.  Defaults to True.
            stamps: whether or not to also return a postage stamp. Can only be used on a
                single source. Defaults to False.
            use_clusters: Whether to use clustering, or fit each source individually.
                Defaults to True.

        Raises:
            ArgumentError: stamps were requested for more than one object.
            ArgumentError: an input argument was not a supported type.

        Returns:
            A tuple containing the flux, flux error, chi-squared value, degrees of
            freedom, cluster ID. If `stamps` is True, the data and model are also returned.
        """
        X0, Y0 = map(
            np.atleast_1d, astropy.wcs.utils.skycoord_to_pixel(positions, self.w)
        )
        X0, Y0 = self._filter_out_of_range(X0, Y0, nbeam, edge_buffer)
        if use_clusters:
            self.cluster(X0, Y0, threshold=cluster_threshold)
        
        if stamps:
            if len(X0) > 1 and not (
                len(self.in_cluster) == len(X0) and len(self.clusters.keys()) == 1
            ):
                raise ArgumentError("Cannot output postage stamps for >1 object")

        if major_axes is None:
            a = np.ones(len(X0)) * (self.BMAJ).to(u.arcsec)
        else:
            if not isinstance(major_axes, astropy.units.Quantity):
                raise ArgumentError("Major axes must be a quantity")

            if major_axes.isscalar:
                a = (major_axes * np.ones(len(X0))).to(u.arcsec)
            else:
                a = major_axes.to(u.arcsec)
                a[np.isnan(a)] = (self.BMAJ).to(u.arcsec)
        
        if minor_axes is None:
            b = np.ones(len(X0)) * (self.BMIN).to(u.arcsec)
        else:
            if not isinstance(minor_axes, astropy.units.Quantity):
                raise ArgumentError("Minor axes must be a quantity")

            if minor_axes.isscalar:
                b = (minor_axes * np.ones(len(X0))).to(u.arcsec)
            else:
                b = minor_axes.to(u.arcsec)
                b[np.isnan(b)] = (self.BMIN).to(u.arcsec)
        
        if position_angles is None:
            pa = np.ones(len(X0)) * (self.BPA)
        else:
            if not isinstance(position_angles, astropy.units.Quantity):
                raise ArgumentError("Position angles must be a quantity")

            if position_angles.isscalar:
                pa = position_angles * np.ones(len(X0))
            else:
                pa = position_angles
                pa[np.isnan(pa)] = self.BPA

        # set up the postage stamps for the sources
        # goes from [xmin,xmax) and [ymin,ymax)
        # so add 1 to the maxes to be inclusive
        # and then check against boundaries
        npix = ((nbeam / 2.0) * a / self.pixelscale).value
        xmin = np.int16(np.round(X0 - npix))
        xmax = np.int16(np.round(X0 + npix)) + 1
        ymin = np.int16(np.round(Y0 - npix))
        ymax = np.int16(np.round(Y0 + npix)) + 1
        xmin[xmin < 0] = 0
        ymin[ymin < 0] = 0
        xmax[xmax > self.data.shape[-1]] = self.data.shape[-1]
        ymax[ymax > self.data.shape[-2]] = self.data.shape[-2]

        flux = np.zeros(len(X0))
        flux_err = np.zeros(len(X0))
        chisq = np.zeros(len(X0))
        dof = np.zeros(len(X0), dtype=np.int16)
        iscluster = np.zeros(len(X0), dtype=np.int16)
        
        
        for i in range(len(X0)):
            if use_clusters:
                if i in self.in_cluster:
                    continue
            if self.use_numba:
                out = self._numba_measure(
                        X0[i],
                        Y0[i],
                        xmin[i],
                        xmax[i],
                        ymin[i],
                        ymax[i],
                        a[i],
                        b[i],
                        pa[i],
                        allow_nan,
                        stamps
                    )
            else:
                out = self._measure(
                    X0[i],
                    Y0[i],
                    xmin[i],
                    xmax[i],
                    ymin[i],
                    ymax[i],
                    a[i],
                    b[i],
                    pa[i],
                    allow_nan=allow_nan,
                    stamps=stamps,
                )
            flux[i], flux_err[i], chisq[i], dof[i], *_ = out

        if use_clusters:
            clusters = list(self.clusters.values())
            for j in range(len(clusters)):
                ii = np.array(list(clusters[j]))
                if self.verbose:
                    logger.info("Fitting a cluster of sources %s" % ii)
                xmin = max(int(round((X0[ii] - npix[ii]).min())), 0)
                xmax = min(int(round((X0[ii] + npix[ii]).max())), self.data.shape[-1]) + 1
                ymin = max(int(round((Y0[ii] - npix[ii]).min())), 0)
                ymax = min(int(round((Y0[ii] + npix[ii]).max())), self.data.shape[-2]) + 1

                out = self._measure_cluster(
                    X0[ii],
                    Y0[ii],
                    xmin,
                    xmax,
                    ymin,
                    ymax,
                    a[ii],
                    b[ii],
                    pa[ii],
                    allow_nan=allow_nan,
                    stamps=stamps,
                )
                f, f_err, csq, _dof = out[:4]
                for k in range(len(ii)):
                    flux[ii[k]] = f[k]
                    flux_err[ii[k]] = f_err[k]
                    chisq[ii[k]] = csq[k]
                    dof[ii[k]] = _dof[k]
                    iscluster[ii[k]] = j + 1

        if positions.isscalar:
            if stamps:
                return (
                    flux[0],
                    flux_err[0],
                    chisq[0],
                    dof[0],
                    iscluster[0],
                    out[-2],
                    out[-1],
                )
            else:
                return flux[0], flux_err[0], chisq[0], dof[0], iscluster[0]

        else:
            flux, flux_err, chisq, dof = self.reshape_output(
                [flux, flux_err, chisq, dof],
                self.idx_mask
            )
            if stamps:
                return flux, flux_err, chisq, dof, iscluster, out[-2], out[-1]
            else:
                return flux, flux_err, chisq, dof, iscluster

    def inject(
        self,
        flux: Union[float, np.ndarray],
        positions: Union[float, np.ndarray],
        major_axes: Optional["astropy.coordinates.Angle"] = None,
        minor_axes: Optional["astropy.coordinates.Angle"] = None,
        position_angles: Optional["astropy.coordinates.Angle"] = None,
        nbeam: int = 3,
    ):
        """Inject one or more fake point sources (defined by the header) into `self.data`
        with the flux(es) and position(s) specified.

        Args:
            flux: Flux(es) of source(s) to inject in same units as the image.
            positions: Position(s) of source(s) to inject.
            major_axes: FWHMs along major axes of sources to measure. If None, will use
                header BMAJ. Defaults to None.
            minor_axes: FWHMs along minor axes of sources to measure. If None, will use
                header BMIN. Defaults to None.
            position_angles: Position angles of sources to measure. If None, will use
                header BPA. Defaults to None.
            nbeam: Diameter of the square cutout for injection in units of the major axis.
                Defaults to 3.
        """
        X0, Y0 = map(
            np.atleast_1d, astropy.wcs.utils.skycoord_to_pixel(positions, self.w)
        )
        flux = np.atleast_1d(flux)

        if major_axes is None:
            a = self.BMAJ.to(u.arcsec) * np.ones(len(X0))
        else:
            if not isinstance(major_axes, astropy.units.Quantity):
                raise ArgumentError("Major axes must be a quantity")

            if major_axes.isscalar:
                a = (major_axes * np.ones(len(X0))).to(u.arcsec)
            else:
                a = major_axes.to(u.arcsec)
                a[np.isnan(a)] = (self.BMAJ).to(u.arcsec)
        if minor_axes is None:
            b = self.BMIN.to(u.arcsec) * np.ones(len(X0))
        else:
            if not isinstance(minor_axes, astropy.units.Quantity):
                raise ArgumentError("Minor axes must be a quantity")

            if minor_axes.isscalar:
                b = (minor_axes * np.ones(len(X0))).to(u.arcsec)
            else:
                b = minor_axes.to(u.arcsec)
                b[np.isnan(b)] = (self.BMIN).to(u.arcsec)
        if position_angles is None:
            pa = self.BPA * np.ones(len(X0))
        else:
            if not isinstance(position_angles, astropy.units.Quantity):
                raise ArgumentError("Position angles must be a quantity")

            if position_angles.isscalar:
                pa = position_angles * np.ones(len(X0))
            else:
                pa = position_angles
                pa[np.isnan(pa)] = self.BPA

        npix = ((nbeam / 2.0) * a / self.pixelscale).value
        xmin = np.int16(np.round(X0 - npix))
        xmax = np.int16(np.round(X0 + npix)) + 1
        ymin = np.int16(np.round(Y0 - npix))
        ymax = np.int16(np.round(Y0 + npix)) + 1
        xmin[xmin < 0] = 0
        ymin[ymin < 0] = 0
        xmax[xmax > self.data.shape[-1]] = self.data.shape[-1]
        ymax[ymax > self.data.shape[-2]] = self.data.shape[-2]

        for i in range(len(X0)):
            self._inject(
                flux[i],
                X0[i],
                Y0[i],
                xmin[i],
                xmax[i],
                ymin[i],
                ymax[i],
                a[i],
                b[i],
                pa[i],
            )

    def _numba_measure(
        self,
        X0,
        Y0,
        xmin,
        xmax,
        ymin,
        ymax,
        a,
        b,
        pa,
        allow_nan=True,
        stamps=False
    ):
        """
        flux,flux_err,chisq,DOF=_numba_measure(X0, Y0, xmin, xmax, ymin, ymax, a, b, pa, allow_nan=True, stamps=False)

        or

        flux,flux_err,chisq,DOF,data,model=_numba_measure(X0, Y0, xmin, xmax, ymin, ymax, a, b,
            pa, allow_nan=True, stamps=False)

        forced photometry for a single source
        if stamps is True, will also output data and kernel postage stamps

        :param X0: x coordinate of source to measure
        :type X0: float
        :param Y0: y coordinate of source to measure
        :type Y0: float
        :param xmin: min x coordinate of postage stamp for measuring
        :type xmin: int
        :param xmax: max x coordinate of postage stamp for measuring
        :type xmax: int
        :param ymin: min y coordinate of postage stamp for measuring
        :type ymin: int
        :param ymax: max y coordinate of postage stamp for measuring
        :type ymax: int
        :param a: fwhm along major axis in angular units
        :type a: `astropy.units.Quantity`
        :param b: fwhm along minor axis in angular units
        :type b: `astropy.units.Quantity`
        :param pa: position angle in angular units
        :type pa: `astropy.units.Quantity`
        :param allow_nan: whether or not to try to measure sources even if some RMS values
            are NaN.  Defaults to True
        :type allow_nan: bool, optional
        :param stamps: whether or not to return postage stamps of the data and model for
            a single source, defaults to False
        :type stamps: bool, optional

        :returns: flux, flux_err, chisq, DOF  or  flux, flux_err, chisq, DOF, data, model
            if stamps=True
        :rtype: float, float, float, float, optionally `np.ndarray`,`np.ndarray`
        """
        sl = tuple((slice(ymin, ymax), slice(xmin, xmax)))
        n = self.noisedata[sl]
        d = self.data[sl]
        contains_nan = np.any(np.isnan(n)) or np.any(np.isnan(d))
        
        # Handle NaN-sources immediately rather than wasting time computing the kernel
        if contains_nan:
            good = np.isfinite(n) & np.isfinite(d)
            if (not allow_nan) or (good.sum() == 0):
                if not stamps:
                    return np.nan, np.nan, np.nan, 0
                else:
                    return (
                        np.nan,
                        np.nan,
                        np.nan,
                        0,
                        d,
                        np.nan * d,
                    )

        d_orig = d
        xx, yy = _meshgrid(xmin, xmax, ymin, ymax)
        ndata = np.size(xx)
        
        kernel = get_kernel(xx,
                            yy,
                            X0,
                            Y0,
                            (a/self.pixelscale).value,
                            (b/self.pixelscale).value,
                            pa.to(u.deg).value
                            )

        # Now revisit NaN-sources as required
        if contains_nan:
            n = n[good]
            d = d[good]
            kernel = kernel[good]
            ndata = good.sum()

        flux, flux_err, chisq = _convolution(d, n, kernel)

        if not stamps:
            return flux, flux_err, chisq, ndata - 1
        else:
            return (
                flux,
                flux_err,
                chisq,
                ndata - 1,
                d_orig,
                flux * kernel,
            )

    def _measure(
        self, X0, Y0, xmin, xmax, ymin, ymax, a, b, pa, allow_nan=True, stamps=False
    ):
        """
        flux,flux_err,chisq,DOF=_measure(X0, Y0, xmin, xmax, ymin, ymax, a, b, pa, allow_nan=True, stamps=False)

        or

        flux,flux_err,chisq,DOF,data,model=_measure(X0, Y0, xmin, xmax, ymin, ymax, a, b,
            pa, allow_nan=True, stamps=False)

        forced photometry for a single source
        if stamps is True, will also output data and kernel postage stamps

        :param X0: x coordinate of source to measure
        :type X0: float
        :param Y0: y coordinate of source to measure
        :type Y0: float
        :param xmin: min x coordinate of postage stamp for measuring
        :type xmin: int
        :param xmax: max x coordinate of postage stamp for measuring
        :type xmax: int
        :param ymin: min y coordinate of postage stamp for measuring
        :type ymin: int
        :param ymax: max y coordinate of postage stamp for measuring
        :type ymax: int
        :param a: fwhm along major axis in angular units
        :type a: `astropy.units.Quantity`
        :param b: fwhm along minor axis in angular units
        :type b: `astropy.units.Quantity`
        :param pa: position angle in angular units
        :type pa: `astropy.units.Quantity`
        :param allow_nan: whether or not to try to measure sources even if some RMS values
            are NaN.  Defaults to True
        :type allow_nan: bool, optional
        :param stamps: whether or not to return postage stamps of the data and model for
            a single source, defaults to False
        :type stamps: bool, optional

        :returns: flux, flux_err, chisq, DOF  or  flux, flux_err, chisq, DOF, data, model
            if stamps=True
        :rtype: float, float, float, float, optionally `np.ndarray`,`np.ndarray`
        """
        sl = tuple((slice(ymin, ymax), slice(xmin, xmax)))
        # unfortunately we have to make a custom kernel for each object
        # since the fractional-pixel offsets change for each
        x = np.arange(xmin, xmax)
        y = np.arange(ymin, ymax)
        xx, yy = np.meshgrid(x, y)
        g = G2D(X0, Y0, (a / self.pixelscale).value, (b / self.pixelscale).value, pa)

        kernel = g(xx, yy)
        
        #print(xmin, xmax, ymin, ymax)
        #print(X0,Y0)
        #print((a / self.pixelscale).value, (b / self.pixelscale).value, pa)
        
        #xx, yy = _meshgrid(xmin, xmax, ymin, ymax)
        
        #kernel = get_kernel(xx, yy, X0, Y0, (a / self.pixelscale).value,(b / self.pixelscale).value, pa.to(u.deg).value)

        # uncertainties: see discussion in Section 3 of Condon (1997)
        # the uncertainty on the amplitude is just the noise at the position of the source
        # so do a weighted average over the beam
        n = self.noisedata[sl]
        d = self.data[sl]
        ndata = np.prod(xx.shape)

        if np.any(np.isnan(n)):
            # protect against NaNs in the data or rms map
            good = np.isfinite(n) & np.isfinite(d)
            n = n[good]
            d = d[good]
            kernel = kernel[good]
            ndata = good.sum()
            if (not allow_nan) or (good.sum() == 0):
                if not stamps:
                    return np.nan, np.nan, np.nan, 0
                else:
                    return (
                        np.nan,
                        np.nan,
                        np.nan,
                        0,
                        self.data[sl],
                        np.nan * kernel,
                    )

        flux = ((d) * kernel / n ** 2).sum() / (kernel ** 2 / n ** 2).sum()
        flux_err = ((n) * kernel / n ** 2).sum() / (kernel / n ** 2).sum()
        chisq = (((d - kernel * flux) / n) ** 2).sum()

        if not stamps:
            return flux, flux_err, chisq, ndata - 1
        else:
            return (
                flux,
                flux_err,
                chisq,
                ndata - 1,
                self.data[sl],
                flux * kernel,
            )

    def _filter_out_of_range(self, X0, Y0, nbeam, edge_buffer=1.):
        """
        X0, Y0 = _filter_out_of_range(X0, Y0, nbeam)
        Filter out sources which are beyond the image range.

        :param X0: x coordinate of source to measure
        :type X0: float
        :param Y0: y coordinate of source to measure
        :type Y0: float
        nbeam: Diameter of the square cutout for fitting in units of
            the major axis. Defaults to 3.
        """
        npix = round((nbeam / 2. * self.BMAJ.to('arcsec') / self.pixelscale).value)

        npix = int(round(npix * edge_buffer))

        X0_mask = (X0 < npix) | (X0 > self.NAXIS1 - npix)
        Y0_mask = (Y0 < npix) | (Y0 > self.NAXIS2 - npix)

        final_mask = np.logical_or(
            X0_mask, Y0_mask
        )

        logger.debug(
            "Removing %i sources that are outside of the image range",
            np.sum(final_mask)
        )
        # save the mask to reconstruct arrays
        self.idx_mask = final_mask

        return X0[~final_mask], Y0[~final_mask]

    def _inject(self, flux, X0, Y0, xmin, xmax, ymin, ymax, a, b, pa):
        """
        _inject(flux, X0, Y0, xmin, xmax, ymin, ymax, a, b, pa)
        injection for a single source

        :param flux: flux of source to inject, in the same units as the image
        :type flux: float
        :param X0: x coordinate of source to measure
        :type X0: float
        :param Y0: y coordinate of source to measure
        :type Y0: float
        :param xmin: min x coordinate of postage stamp for measuring
        :type xmin: int
        :param xmax: max x coordinate of postage stamp for measuring
        :type xmax: int
        :param ymin: min y coordinate of postage stamp for measuring
        :type ymin: int
        :param ymax: max y coordinate of postage stamp for measuring
        :type ymax: int
        :param a: fwhm along major axis in angular units
        :type a: `astropy.units.Quantity`
        :param b: fwhm along minor axis in angular units
        :type b: `astropy.units.Quantity`
        :param pa: position angle in angular units
        :type pa: `astropy.units.Quantity`
        """
        sl = tuple((slice(ymin, ymax), slice(xmin, xmax)))
        # unfortunately we have to make a custom kernel for each object
        # since the fractional-pixel offsets change for each
        x = np.arange(xmin, xmax)
        y = np.arange(ymin, ymax)
        xx, yy = np.meshgrid(x, y)
        g = G2D(X0, Y0, (a / self.pixelscale).value, (b / self.pixelscale).value, pa)

        kernel = g(xx, yy).value
        self.data[sl] += kernel * flux

    def _measure_cluster(
        self,
        X0,
        Y0,
        xmin,
        xmax,
        ymin,
        ymax,
        a,
        b,
        pa,
        allow_nan=True,
        stamps=False,
        fitter=fitting.LevMarLSQFitter(),
    ):
        """
        flux,flux_err,chisq,DOF=_measure(X0, Y0, xmin, xmax, ymin, ymax, a, b, pa, allow_nan=True,
            stamps=False, fitter = fitting.LevMarLSQFitter())
        or
        flux,flux_err,chisq,DOF,data,model=_measure(X0, Y0, xmin, xmax, ymin, ymax, a, b,
            pa, allow_nan=True, stamps=False, fitter = fitting.LevMarLSQFitter())

        forced photometry for a cluster of sources using astropy fitting

        :param X0: x coordinates of source to measure
        :type X0: `numpy.ndarray`
        :param Y0: y coordinates of source to measure
        :type Y0: `numpy.ndarray`
        :param xmin: min x coordinate of postage stamp for measuring
        :type xmin: int
        :param xmax: max x coordinate of postage stamp for measuring
        :type xmax: int
        :param ymin: min y coordinate of postage stamp for measuring
        :type ymin: int
        :param ymax: max y coordinate of postage stamp for measuring
        :type ymax: int
        :param a: fwhm of each source along major axis in angular units
        :type a: `astropy.units.Quantity`
        :param b: fwhm of each source along minor axis in angular units
        :type b: `astropy.units.Quantity`
        :param pa: position angle of each source in angular units
        :type pa: `astropy.units.Quantity`
        :param allow_nan: whether or not to try to measure sources even if some RMS values
            are NaN.  Defaults to True
        :type allow_nan: bool, optional
        :param stamps: whether or not to return postage stamps of the data and model,
            defaults to False
        :type stamps: bool, optional
        :param fitter: fitting object, defaults to `fitting.LevMarLSQFitter()`
        :type fitter: optional

        :returns: flux, flux_err, chisq, DOF  or  flux, flux_err, chisq, DOF, data, model
            if stamps=True
        :rtype: numpy.ndarray, numpy.ndarray, numpy.ndarray, numpy.ndarray, optionally
            `np.ndarray`,`np.ndarray`

        """
        x = np.arange(xmin, xmax)
        y = np.arange(ymin, ymax)
        xx, yy = np.meshgrid(x, y)

        g = None
        for k in range(len(X0)):
            if g is None:
                g = models.Gaussian2D(
                    x_mean=X0[k],
                    y_mean=Y0[k],
                    x_stddev=(a[k] / self.pixelscale).value
                    / 2
                    / np.sqrt(2 * np.log(2)),
                    y_stddev=(b[k] / self.pixelscale).value
                    / 2
                    / np.sqrt(2 * np.log(2)),
                    theta=(pa[k] - pa_offset).to(u.rad).value,
                    fixed={
                        "x_mean": True,
                        "y_mean": True,
                        "x_stddev": True,
                        "y_stddev": True,
                        "theta": True,
                    },
                )
            else:
                g = g + models.Gaussian2D(
                    x_mean=X0[k],
                    y_mean=Y0[k],
                    x_stddev=(a[k] / self.pixelscale).value
                    / 2
                    / np.sqrt(2 * np.log(2)),
                    y_stddev=(b[k] / self.pixelscale).value
                    / 2
                    / np.sqrt(2 * np.log(2)),
                    theta=(pa[k] - pa_offset).to(u.rad).value,
                    fixed={
                        "x_mean": True,
                        "y_mean": True,
                        "x_stddev": True,
                        "y_stddev": True,
                        "theta": True,
                    },
                )

        sl = tuple((slice(ymin, ymax), slice(xmin, xmax)))
        n = self.noisedata[sl]
        d = self.data[sl]
        # protect against NaNs in the data or rms map
        good = np.isfinite(n) & np.isfinite(d)

        flux = np.zeros(len(X0))
        flux_err = np.zeros(len(X0))

        if (np.any(~good) and (not allow_nan)) or (good.sum() == 0):
            # either >=1 bad point and no bad points allowed
            # OR
            # no good points left
            if stamps:
                return (
                    flux * np.nan,
                    flux_err * np.nan,
                    flux * np.nan,
                    flux,
                    d,
                    d * np.nan,
                )
            else:
                return (
                    flux * np.nan,
                    flux_err * np.nan,
                    flux * np.nan,
                    flux,
                )

        try:
            out = fitter(g, xx[good], yy[good], d[good], weights=1.0 / n[good])
        except TypeError as err:
            logger.debug("Unable to fit cluster: %s", err)
            if stamps:
                return (
                    flux * np.nan,
                    flux_err * np.nan,
                    flux * np.nan,
                    flux,
                    d,
                    d * np.nan,
                )
            else:
                return (
                    flux * np.nan,
                    flux_err * np.nan,
                    flux * np.nan,
                    flux,
                )
        model = out(xx, yy)
        chisq = np.zeros(len(X0)) + (((d[good] - model[good]) / n[good]) ** 2).sum()
        dof = np.zeros(len(X0), dtype=np.int16) + np.prod(xx[good].shape) - len(X0)

        for k in range(len(X0)):
            flux[k] = out.__getattr__("amplitude_%d" % k).value
            # a weighted average would be better for the noise here, but
            # to simplify just use the noise map at the central source position
            flux_err[k] = self.noisedata[np.int16(round(Y0[k])), np.int16(round(X0[k]))]

        if stamps:
            return flux, flux_err, chisq, dof, d, model
        else:
            return flux, flux_err, chisq, dof

    def _measure_astropy(
        self, fi, fb, fn, X0, Y0, xmin, xmax, ymin, ymax, a, b, pa, nbeam=3, stamps=False
    ):
        """
        flux, flux_err, chisq, DOF = _measure_astropy(
            X0, Y0, xmin, xmax, ymin, ymax, a, b, pa, nbeam=3, stamps=False
        )
        or
        flux, flux_err, chisq, DOF, data,model = _measure_astropy(
            X0, Y0, xmin, xmax, ymin, ymax, a, b, pa, nbeam=3, stamps=False
        )


        forced photometry for a single source using our astropy version
        X0, Y0, xmin, ymin, xmax, ymax all in pixels
        a, b, pa all Quantities
        nbeam is the size of the cutout for fitting in units of the major axis
        if stamps is True, will also output data and kernel postage stamps

        this accomplishes the same task as _measure() with astropy
        it seems very similar but a factor of 2-3 slower

        JUST FOR DEBUGGING
        """
        p = astropy.wcs.utils.pixel_to_skycoord(X0, Y0, self.w)
        im = astropy.nddata.Cutout2D(fi[0].data, p, nbeam * a, wcs=self.w)
        bg = fb[0].data[
            im.ymin_original : im.ymax_original + 1,
            im.xmin_original : im.xmax_original + 1,
        ]
        ns = fn[0].data[
            im.ymin_original : im.ymax_original + 1,
            im.xmin_original : im.xmax_original + 1,
        ]

        x = np.arange(im.data.shape[1])
        y = np.arange(im.data.shape[0])
        xx, yy = np.meshgrid(x, y)
        x0, y0 = astropy.wcs.utils.skycoord_to_pixel(p, im.wcs)
        g = G2D(x0, y0, (a / self.pixelscale).value, (b / self.pixelscale).value, pa)
        kernel = g(xx, yy)
        flux = ((im.data - bg) * kernel / ns ** 2).sum() / (kernel ** 2 / ns ** 2).sum()
        flux_err = ((ns) * kernel / ns ** 2).sum() / (kernel / ns ** 2).sum()

        chisq = (((im.data - flux * kernel) / ns.data) ** 2).sum()
        dof = np.prod(xx.shape) - 1
        if not stamps:
            return flux, flux_err, chisq, dof
        else:
            return flux, flux_err, chisq, dof, im.data, flux * kernel

    @staticmethod
    def reshape_output(inputs_list, mask):
        out = []
        for el in inputs_list:
            myarr = np.zeros(mask.shape)
            myarr[np.where(mask == False)] = el
            out.append(myarr)
        return tuple(out)
