"""Class implementation for ZTF Observing Strategy"""
from __future__ import absolute_import, print_function, division
import os
import glob
import logging
from matplotlib.patches import Polygon
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import healpy as hp
import sncosmo
from astropy.cosmology import Planck15
from astropy.time import Time, TimeDelta
from opsimsummary import (AllSkySNVisualization, split_PolygonSegments,
                          convertToCelestialCoordinates, convertToSphericalCoordinates,
                          healpix_boundaries)
__all__ = ['ztfbandcolors', 'ZTFSNViz']
ztfbandcolors = dict(g='g', r='r', i='y')
#logger = logging.getLogger(__name__)



class ZTFSNViz(AllSkySNVisualization):
    def __init__(self, bandColorDict=ztfbandcolors, radius_deg=4.,
                 showVisibleFields=False, showVarScatter=False,
                 showMW=True,
                 # ztf variables
                 depths=None, data_dir=None, offset=0.):

        AllSkySNVisualization.__init__(self, bandColorDict, radius_deg,
                                       showMW=True,
                                       showVisibleFields=showVisibleFields,
                                       showVarScatter=showVarScatter)
        #self.logger = logging.getLogger('ztf')
        self.offset = offset
        if showVisibleFields:
            if data_dir is None:
                raise ValueError(
                    'data_dir cannot be None if showVisibleFields is True')
            self._data_dir = data_dir
            self.visible_fields = glob.glob(data_dir + '/*.da')

    def filenameFromMjd(self, mjd):
        jd = np.int(mjd + 2400000.5 + self.offset)
        fname = '{:d}'.format(jd)
        return os.path.join(self._data_dir, fname + '.da')

    def _minutes_since_noon2mjd(self, expMJD, minutes=0, utc_offset=-5):
        # Datetime object corresponding to Noon at UTC
        dt = Time(np.floor(expMJD), format='mjd')
        # Use utc_offset to switch to Palomar local time
        DeltaT = TimeDelta(utc_offset * 3600 + minutes * 60, format='sec')
        xx = dt + DeltaT
        return xx.mjd

    def fieldCoords(self, mjd):
        fname = self.filenameFromMjd(mjd)
        df = pd.read_csv(fname, delim_whitespace=True, skiprows=1, names=['ind', 'ra', 'dec', 'start', 'end', 'on'],
                         index_col='ind')
        # convert ra to degrees from -180. to 180. from hours
        df.ra = df.ra * 360. / 24.0 - 180.

        # Start and end times
        df.start = list(self._minutes_since_noon2mjd(mjd, v)
                        for v in df.start.values)
        df.end = list(self._minutes_since_noon2mjd(mjd, v)
                      for v in df.end.values)
        return df

    def healpixels(self, mjd, nside=8):
        df = self.fieldCoords(mjd).query(
            'on > 0.5 and @mjd < end and @mjd > start')
        theta, phi = convertToSphericalCoordinates(
            df.ra.values, df.dec.values, unit='degrees')
        phi = phi + np.pi
        ipix = hp.ang2pix(nside=nside, theta=theta, phi=phi, nest=True)
        ipix.sort()
        return np.unique(ipix)

    def split_boundaries(self, lon, lat, step, numObjs, split_lon=180.):

        numPixels = numObjs
        numPoints = 4 * step
        assert len(lon) == len(lat)
        assert len(lon) == numPoints * numPixels

        ra = lon.reshape(numPixels, numPoints)
        dec = lat.reshape(numPixels, numPoints)

        polypts = []
        for (llo, lla) in zip(ra, dec):
            polypts += self.split_boundary(llo, lla, split_lon=split_lon)
        return polypts

    def split_boundary(self, lon, lat, split_lon=180):

        if 1 > 0:
            # if any(lon<=split_lon) and any(lon>split_lon):
            mask = lon <= split_lon + 0.000001
            lon0 = lon[mask]
            lat0 = lat[mask]

            lon1 = lon[~mask]
            lat1 = lat[~mask]
            return list(((lon0, lat0), (lon1, lat1)))
        else:
            return list((lon, lat))

    def get_visible_field_polygons(self, mjd, m, facecolor='k', alpha=1., **kwargs):
        patches = []
        nside = kwargs.get('nside', 8)
        ipix = kwargs.get('ipix', None)
        split_lon = kwargs.get('split_lon', 180)

        if ipix is None:
            ipix = self.healpixels(mjd=mjd, nside=nside)
        lon, lat = healpix_boundaries(ipix, nside=nside, units='degrees',
                                      convention='celestial', step=10,
                                      nest=True)
        polypts = self.split_boundaries(
            lon, lat, step=10, numObjs=len(ipix), split_lon=split_lon)

        for poly in polypts:
            ra, dec = poly
            if len(ra) > 0:
                patches.append(self._build_poly(
                    ra, dec, m, facecolor=facecolor, alpha=alpha, edgecolor='k'))
        return patches

    def _build_poly(self, ra, dec, m, facecolor='k', alpha=1., edgecolor='k'):
        x, y = m(ra, dec)
        xy = zip(x, y)
        p = Polygon(xy, facecolor=facecolor, fill=True,
                    alpha=alpha, edgecolor=edgecolor, lw=0)
        return p

    @staticmethod
    def maglims(band):
        """depth for bands"""
        maglim = dict(g=22, r=22, i=22)
        return maglim[band]

    def scale_mags_size(self, mags, band):
        """Size of points to be used in a scatter plot"""
        return 0.1 * (self.maglims(band) - mags)

    def generate_var_scatter(self, mjd, band, simsdf):
        """use a simulated catalog of SN to generate the data
        for a scatter plot of the SN
        
        Parameters
        ----------
        mjd : float, mandatory, units of days
            modified julian day of observation
        band : string, 
            denoting bandpass filter used for observation
        simsdf : `pd.DataFrame`
            object containing catalog of simulated SN including all
            SN to a redshift of 0.2 going off in the Northern Sky, their
            redshifts `z`, locations in terms of `ra` and `dec` and peak
            times `t_0`
        """
        regband = 'sdss' + band
        simsdf['time'] = mjd - simsdf.t0
        simsdf = simsdf.query('time > -30 and time < 50')
        model = sncosmo.Model(source='salt2')
        x0 = np.zeros(len(simsdf))
        mag = np.zeros(len(simsdf))
        i = 0
        for z, time in zip(simsdf.z.values, simsdf.time.values):
            model.set(z=z)
            model.set_source_peakabsmag(-19.3, 'bessellb', 'ab', cosmo=Planck15)
            x0[i] = model.get('x0')
            mag[i] = model.bandmag(time=time, band=regband, magsys='ab')
            i += 1
        simsdf['mag'] = mag
        simsdf['x0'] = x0
        simsdf['rad'] = self.scale_mags_size(simsdf.mag, band)
        return simsdf

    def generate_images_from(self, obsHistIDs, obsdf, snsims,
                             savefig=False, outdir='./',
                             rootname='ztf_obsHistID_'):
        """
        Conveniently generate a set of images
        """
        for obsHistID, row in obsdf.loc[obsHistIDs].iterrows():
            mjd, ra, dec, band = row[['expMJD', 'fieldRA', 'fieldDec',
                                      'filter']]
            fname = os.path.join(outdir,
                                 rootname +'{:06d}.png'.format(int(obsHistID)))
            # Keep the image generation process in a try except block so that a
            # failure does not derail the entire process.
            try:
                fig, ax = self.generate_image(ra=np.degrees(ra),
                                              dec=np.degrees(dec),
                                              radius_deg=4., mjd=mjd,
                                              band=band, sndf=snsims) 
                _ = fig.savefig(fname)
                fig.clf()
                plt.close()
                logging.info('generated image for obsHistID {}'.format(obsHistID))
            except:
                logging.error('Failed to generate image for obsHistID {}'.format(obsHistID))


        

        
                            
