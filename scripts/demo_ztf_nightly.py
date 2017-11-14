#!/usr/bin/env python
import os
#from copy import deepcopy
#import datetime
#import glob
import time
import logging
from datetime import datetime
from astropy.time import Time, TimeDelta
import pytz
import traceback
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
    df.t0 = df.t0 - 2400000.5
    return df


if __name__=='__main__':
    
    logger = logging.getLogger('ztf_nights')
    logger.setLevel(logging.INFO)
    logging.basicConfig(filename='multiImage.log', level=logging.INFO)

    numSplits = 24
    ztf_data = '/home/rbisw/data/ZTF'
    scheduler_fname = os.path.join(ztf_data, 'test_schedule_v3.db')
    df = read_scheduler(scheduler_fname)

    dt = datetime(2018, 1, 3, 0, 0, tzinfo=pytz.timezone('US/Pacific'))
    tt = Time(dt)
    mjds = np.array(list((tt + TimeDelta(n, format='jd')).mjd for n in range(365))) 
    obsHistIDs = list((df['expMJD']-mjd).abs().argsort()[:1].values[0] for mjd in mjds)  
    dftodo = df.loc[obsHistIDs]
    dfs = np.array_split(dftodo, numSplits)

    def make_images(i):
        ztf_data = '/home/rbisw/data/ZTF'
        sim_fname = os.path.join(ztf_data, 'sn_coord_time.dat')
        ztf_fields = '/nfs/brahe/ZTF/Scheduler/NugentSim/tst/data_year'
        ztf_figs = '/nfs/brahe/scratch/rbiswas/ztf_figs_year2/'
        logfname = 'images_{}.log'.format(i)
        logfname = os.path.join(ztf_figs, logfname)
        tsplitstart = time.time()
        with open(logfname, 'w') as f:
            f.write('starting split {0} at time {1}\n'.format(i, tsplitstart))
        simsdf = read_sim(sim_fname)
        ztfsky = ZTFSNViz(showVisibleFields=True, showVarScatter=True, showMW=True,
                          data_dir=ztf_fields)
        ztfsky.generate_images_from(dfs[i].index.values, dfs[i], snsims=simsdf, savefig=True,
                                    outdir=ztf_figs, surveystart=df.expMJD.min()) 
        tsplitend = time.time()
        with open(logfname, mode='a+') as f:
            f.write('images generated at time {} \n'.format(tsplitend))
        return 0

    ndf = Parallel(n_jobs=-1)(delayed(make_images)(i) for i in range(numSplits))
