import os
import pandas as pd

from astropy.time import Time
from astropy.coordinates import SkyCoord

from ztfquery import query

import argparse

parser = argparse.ArgumentParser(description='Download ZTF image metadata')
parser.add_argument('--table', type=str, required=True,
                    help='Path of BTS table in csv format')
parser.add_argument('--outdir', type=str, default='metadata',
                    help='Directory where metadata tables are saved')
parser.add_argument('--ifirst', type=int, default=0,
                    help='Index of first object in the list')
parser.add_argument('--ilast', type=int, default=0,
                    help='Index of last object in the list. Default is the last one')
parser.add_argument('--phase-min', type=int, default=50,
                    help='Number of days before peak flux')
parser.add_argument('--phase-max', type=int, default=150,
                    help='Number of days after peak flux')
parser.add_argument('--overwrite', action="store_true",
                    help='If set, will re-download all metatables')
args = parser.parse_args()

def read_bts_table(table_name):
    ''' Reads Bright Transient Survey list of objects,
        which can be obtained with:
        from ztfquery import bts
        databts = bts.download_bts_table()
        databts.to_csv('BTS_transients.csv')
    '''
    print('Reading BTS table of transients..')
    tab = pd.read_csv(table_name)
    print(f'Found {len(tab)} objects')

    name = tab['ZTFID'].values.astype(str)
    sky_coord = SkyCoord(tab['RA'], tab['Dec'], unit=('hourangle', 'deg'))
    ra = sky_coord.ra.deg
    dec = sky_coord.dec.deg
    jd_peak = Time(tab['peakt_jd'], format='jd').jd

    return name, ra, dec, jd_peak

def get_metadata(name, ra, dec, jd_peak, 
                  phase_min=50, phase_max=150,
                  overwrite=False):
    ''' Loads metadata information about this transient
        If metadata is available locally, it reads it
        otherwise it makes a query
    '''

    zquery = query.ZTFQuery()
    meta_name = f'{args.outdir}/{name}_{phase_min}_{phase_max}.csv'
    if os.path.exists(meta_name) and overwrite==False:
        print(f'{name}: Metatable already exists')
        return 
    else:
        print(f'{name}: Querying metatable between -{phase_min} and {phase_max} days of peak')
        jd_start = jd_peak - phase_min 
        jd_end = jd_peak + phase_max 
        
        #-- Make query
        zquery.load_metadata(radec=[ra, dec], 
                            sql_query=f"obsjd>{jd_start} and obsjd<{jd_end}")
        print(f'{name}: found {len(zquery.metatable)} exposures')

        #-- Save info to file 
        zquery.metatable.to_csv(meta_name)




names, ras, decs, jd_peaks = read_bts_table(args.table)
os.makedirs(args.outdir, exist_ok=True)

ilast = len(name) if args.ilast == 0 else args.ilast

for i in range(args.ifirst, ilast):
    print(f'Progress: {args.ifirst} <= {i} < {args.ilast}')
    get_metadata(names[i], ras[i], decs[i], jd_peaks[i], 
                  phase_min=args.phase_min, phase_max=args.phase_max,
                  overwrite=args.overwrite)
        
    