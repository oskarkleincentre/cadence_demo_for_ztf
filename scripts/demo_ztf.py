#!/usr/bin/env python
import os
#from copy import deepcopy
#import datetime
#import glob
import time
import logging


import numpy as np
import pandas as pd

#import matplotlib.pyplot as plt
#from matplotlib.patches import Polygon
#from mpl_toolkits.basemap import Basemap
#import sncosmo
#import healpy as hp
#from opsimsummary import OpSimOutput
#from opsimsummary import convertToSphericalCoordinates, healpix_boundaries
from vizztf import ZTFSNViz

from sqlalchemy import create_engine
#from astropy.cosmology import Planck15
#from astropy.time import Time, TimeDelta, TimeDeltaJD

from joblib import Parallel, delayed

def read_scheduler(fname):
    """read the scheduler"""
    filename = 'sqlite:///' + fname
    print(filename)
    assert os.path.exists(fname)
    engine = create_engine(filename)

    df = pd.read_sql_table('Summary', con=engine)
    df.obsHistID = df.obsHistID.astype(np.int)
    df.set_index('obsHistID', inplace=True)
    df.expMJD = df.expMJD + 365 * 2 # ZTF is not done
    return df

def read_sim(fname='/home/rbisw/data/ZTF/sn_coord_time.dat'):
    """read the sn catalog"""
    df = pd.read_csv(fname,
                     skiprows=1, delim_whitespace=True,
                     names=('z', 'ra', 'dec', 't0'))
    #df.ra = df.ra * 360./24.0 - 180.
    df.t0 = df.t0 - 2400000.5
    return df

#class Pointings(object):
#    def __init__(self, data_dir, offset=0.):
#        self._data_dir = data_dir
#        self.visible_fields = glob.glob(data_dir + '/*.da')
#        self.dates = np.array(list(np.float(x.split('/')[-1].strip('.da')) for x in self.visible_fields))
#        self.dates.sort()
#        self.offset = offset
#        self.colordict=dict(g='g', r='r', i='y')
#    def filenameFromMjd(self, mjd):
#        jd = np.int(mjd + 2400000.5 + self.offset)
#        fname = '{:d}'.format(jd)
#        return os.path.join(self._data_dir, fname + '.da')
#    
#    def _minutes_since_noon2mjd(self, expMJD, minutes=0, utc_offset=-5):
#        # Datetime object corresponding to Noon at UTC
#        dt = Time(np.floor(expMJD), format='mjd')
#        # Use utc_offset to switch to Palomar local time
#        DeltaT = TimeDelta(utc_offset * 3600 + minutes*60, format='sec')
#        xx = dt + DeltaT
#        return xx.mjd
#    
#    def fieldCoords(self, mjd):
#        fname = self.filenameFromMjd(mjd)
#        df = pd.read_csv(fname, delim_whitespace=True, skiprows=1,
#                         names=['ind', 'ra', 'dec', 'start', 'end', 'on'],
#                         index_col='ind')
#        # convert ra to degrees from -180. to 180. from hours
#        df.ra = df.ra * 360./24.0 - 180.
#        
#        # Start and end times
#        df.start = list(self._minutes_since_noon2mjd(mjd, v) for v in df.start.values)
#        df.end = list(self._minutes_since_noon2mjd(mjd, v) for v in df.end.values)
#        
#        return df
#    def healpixels(self, mjd, nside=8):
#        df = self.fieldCoords(mjd).query('on > 0.5 and @mjd < end and @mjd > start')
#        theta, phi = convertToSphericalCoordinates(df.ra.values, df.dec.values, unit='degrees')
#        phi = phi + np.pi
#        ipix = hp.ang2pix(nside=nside, theta=theta, phi=phi, nest=True)
#        ipix.sort()
#        return np.unique(ipix)
#    def pixel_patches(self, mjd, m, nside=8, facecolor='k', alpha=1., ipix=None):
#        patches = []
#        if ipix is None:
#            ipix = self.healpixels(mjd=mjd, nside=nside)
#        for pix in ipix:
#            lon, lat = healpix_boundaries(pix, nside=nside, units='degrees',
#                                          convention='celestial', step=10,
#                                          nest=True)
#            x, y = m(lon, lat)
#            xy = zip(x, y)
#            p = Polygon(xy, facecolor=facecolor, fill=True,alpha=alpha, edgecolor='k', lw=0)
#            patches.append(p)
#        return patches
#    def generate_image(self, obsHistID, df, snsims=None):
#        
#        mjd = df.ix[obsHistID, 'expMJD']
#        band = df.ix[obsHistID, 'filter']
#        ra = np.degrees(df.ix[obsHistID, 'fieldRA'])
#        dec = np.degrees(df.ix[obsHistID, 'fieldDec'])
#        # generate fig
#        fig, ax = plt.subplots()#1, 2)
#        #ax = axx[0]
#        success = True
#        try:
#            m = Basemap(#llcrnrlat=-32., llcrnrlon=48.,
#                        #urcrnrlat=-22., urcrnrlon=58,
#                        projection='moll', lon_0=0., lat_0=0.,
#                        ax=ax, celestial=True)
#            _ = m.drawparallels(np.arange(-91.,91.,20.))
#            _ = m.drawmeridians(np.arange(-180., 181., 30.))
#            # color the sky blue
#            _ = m.drawmapboundary(color='b', fill_color='b')
#            ####_allpatches = self.pixel_patches(mjd, m=m, nside=32, facecolor='b', alpha=0.8, 
#            ####                                 ipix=np.arange(hp.nside2npix(32)))
#            ####_ = list(ax.add_patch(p) for p in _allpatches)
#            # color visible fields black
#            vispatches = self.pixel_patches(mjd=mjd, m=m, nside=8)
#            _ = list(ax.add_patch(p) for p in vispatches)
#            # put the field of view
#            m.tissot(lon_0=ra, lat_0=dec, radius_deg=4., npts=100, ax=ax,
#             **dict(fill=False, edgecolor=self.colordict[band], lw=2))
#            if snsims is not None:
#                scatter_vals = self.sn_scatter(obsHistID, df, snsims)
#                x, y = m(scatter_vals.ra.values, scatter_vals.dec.values)
#                ax.scatter(x, y, s=scatter_vals.rad.values, c='w', edgecolors='w', zorder=10)
#        except:
#            fig, ax = plt.subplots()
#            success = False
#            m = None
#        return fig, ax, m, success
#    def generate_images(self, obsHistIDs, df, snsims=None, savefig=False, outdir='./'):
#        figstuff = list(self.generate_image(obsHistID, df, snsims=snsims) for obsHistID in obsHistIDs)
#        mjd = df.ix[obsHistIDs, 'expMJD'].values
#        if savefig:
#            for obsHistID, figstuff in zip(obsHistIDs, figstuff):
#                if figstuff[-1]:
#                    fig = figstuff[0]
#                    fname = os.path.join(outdir, 'ztf_obsHistID_{:06d}.png'.format(np.int(obsHistID)))
#                    #print(type(fig), fname)
#                    fig.savefig(fname)
#        return figstuff, obsHistIDs
#    def sn_scatter(self, obsHistID, df, simsdf):
#        time = df.ix[obsHistID, 'expMJD']
#        band = 'sdss' + df.ix[obsHistID, 'filter']
#        depth = dict(sdssg=22, sdssr=22, sdssi=22)
#        simsdf = deepcopy(simsdf.query('z > 0.02'))
#        simsdf['time'] = time - simsdf.t0
#        simsdf = simsdf.query('time > -30 & time < 50')
#        model = sncosmo.Model(source='salt2')
#        x0 = np.zeros(len(simsdf))
#        mag = np.zeros(len(simsdf))
#        i = 0
#        for z, time in zip(simsdf.z.values, simsdf.time.values):
#            model.set(z=z, t0=time)
#            model.set_source_peakabsmag(-19.3, 'bessellb', 'ab', cosmo=Planck15)
#            x0[i] = model.get('x0')
#            mag[i] = model.bandmag(time=time, band=band, magsys='ab')
#            i += 1
#        simsdf['mag'] = mag
#        simsdf['x0'] = x0
#        simsdf['area'] = 0.1* (depth[band] - simsdf.mag)
#        return simsdf

if __name__=='__main__':
    
    logger = logging.getLogger('ztf')
    logger.setLevel(logging.INFO)
    logging.basicConfig(filename='multiImage.log', level=logging.INFO)

    startRow = 30000 
    numSplits = 1000
    numRows = 5000
    ztf_data = '/home/rbisw/data/ZTF'
    scheduler_fname = os.path.join(ztf_data, 'test_schedule_v3.db')
    df = read_scheduler(scheduler_fname).loc[startRow : startRow + numRows]
    dfs = np.array_split(df, numSplits)

    def make_images(i):
        ztf_data = '/home/rbisw/data/ZTF'
        sim_fname = os.path.join(ztf_data, 'sn_coord_time.dat')
        ztf_fields = '/nfs/brahe/ZTF/Scheduler/NugentSim/tst/data_year'
        ztf_figs = '/nfs/brahe/scratch/rbiswas/ztf_figs/'
        logfname = 'images_{}.log'.format(i)
        logfname = os.path.join(ztf_figs, logfname)
        tsplitstart = time.time()
        with open(logfname, 'w') as f:
            f.write('starting split {0} at time {1}\n'.format(i, tsplitstart))
        simsdf = read_sim(sim_fname)
        ztfsky = ZTFSNViz(showVisibleFields=True, showVarScatter=True, showMW=True,
                          data_dir=ztf_fields)
        ztfsky.generate_images_from(dfs[i].index.values, dfs[i], snsims=simsdf, savefig=True) 
        tsplitend = time.time()
        with open(logfname, mode='a+') as f:
            f.write('images generated at time {} \n'.format(tsplitend))
        return 0

    ndf = Parallel(n_jobs=-1)(delayed(make_images)(i) for i in range(numSplits))
