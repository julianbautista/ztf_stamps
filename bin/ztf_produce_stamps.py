import os, sys
import pandas as pd
import numpy as np

from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from astropy.io import fits

from ztfquery import query
import argparse



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
    #ra = sky_coord.ra.deg
    #dec = sky_coord.dec.deg
    jd_peak = Time(tab['peakt_jd'], format='jd').jd

    return name, sky_coord, jd_peak

def get_metadata(name, meta_dir='metadata', phase_min=50, phase_max=150):
    ''' Loads metadata information about this transient
        If metadata is available locally, it reads it
        otherwise returns None
    '''
    zquery = query.ZTFQuery()
    meta_name = f'{meta_dir}/{name}_{phase_min}_{phase_max}.csv'
    #meta_name = f'{meta_dir}/{name}.csv'

    if os.path.exists(meta_name):
        zquery = query.ZTFQuery().from_metafile(meta_name)
        print(f'{name}: Loading local metatable with '
                f'{len(zquery.metatable)} epochs from {meta_name}')
        return zquery
    else:
        print(f'{name}: Metatable not found {meta_name}')
        return None

def exist(files):
    ''' Takes a list of filenames and returns an array of booleans 
        stating if the files exist
    '''
    files_exist = np.array([os.path.exists(f) for f in files])
    return files_exist 

def get_filenames_from_meta(zquery, suffix='sciimg.fits'):
    
    #-- Create filenames from metadata but do not download them
    files_sci = zquery.get_data(suffix=suffix, downloadit=False)

    #-- Check how many exist
    files_sci_exist = exist(files_sci)

    return files_sci, files_sci_exist

def get_stampnames_from_meta(zquery, name, outdir):

    mjd = Time(zquery.metatable['obsjd'], format='jd').mjd
    filt = {1:'g', 2:'r', 3:'i'}
    cam = np.array([filt[fid] for fid in zquery.metatable['fid'].values])
    stamp_names = [f'{outdir}/{name}/{name}-{cam[i]}-{mjd[i]:.5f}.fits' for i in range(mjd.size)]
    return stamp_names

def get_stamp(image_file, sky_coord, pixels=32):
    ''' This version uses astropy and saves the header from original image
    '''
    if not os.path.exists(image_file):
        print(f' Image not found: {image_file}')
        return None
    image_full = fits.open(image_file)[0]
    wcs = WCS(header=image_full.header)

    cutout = Cutout2D(  image_full.data, 
                        position = sky_coord, 
                        size = (pixels, pixels),
                        wcs = wcs,
                        mode = 'partial',
                        fill_value = np.nan,
                        copy=True)

    header = image_full.header
    header.update(cutout.wcs.to_header())
    cutout.header = header
    cutout.filter = header['FILTER'][-1:]
    cutout.mjd = header['OBSMJD']
                    
    return cutout
    
def get_stamps(image_files, sky_coord, pixels=32):

    stamps = []
    nfiles = image_files.size
    for i in range(nfiles):
        print(f' Producing stamp for image {i} out of {nfiles}')
        stamp = get_stamp(image_files[i], sky_coord, pixels=pixels)
        if not stamp is None:
            stamps.append(stamp)
    return stamps

def write_stamp(name, stamp, directory='stamps_v1', overwrite=False):
    
    filter = stamp.filter
    mjd = stamp.mjd 
    
    os.makedirs(f'{directory}/{name}', exist_ok=True)
    stamp_name = (f'{directory}/{name}/{name}-{filter}-{mjd:.5f}.fits')
    
    if os.path.exists(stamp_name) and overwrite==False:
        return

    hdu = fits.PrimaryHDU(data=stamp.data, header=stamp.header)
    hdu.writeto(stamp_name, overwrite=True)

def main():
    parser = argparse.ArgumentParser(description='Download ZTF image metadata')
    parser.add_argument('--table', type=str, required=True,
                        help='Path of BTS table in csv format')
    parser.add_argument('--meta-dir', type=str, default='metadata',
                        help='Directory containing metadata tables')
    parser.add_argument('--outdir', type=str, default='stamps',
                        help='Directory where stamps are saved')
    parser.add_argument('--pixels', type=int, default=32,
                        help='Number of pixels on a side for the stamps.')
    parser.add_argument('--ifirst', type=int, default=0,
                        help='Index of first object in the list')
    parser.add_argument('--ilast', type=int, default=0,
                        help='Index of last object in the list. Default is the last one')
    parser.add_argument('--overwrite', action="store_true", default=False,
                        help='If set, will re-write all stamps')
    parser.add_argument('--check-stamps', action="store_true", default=False,
                        help='If set, will only print the summary of how many stamps were produced')
    parser.add_argument('--suffix', type=str, default='sciimg.fits', 
                        help='''
    Here is the list of available options depending on you image kind:
        
    # Science image (kind="sci"):
    - sciimg.fits (primary science image) # (default)
    - mskimg.fits (bit-mask image)
    - psfcat.fits (PSF-fit photometry catalog)
    - sexcat.fits (nested-aperture photometry catalog)
    - sciimgdao.psf (spatially varying PSF estimate in DAOPhot's lookup table format)
    - sciimgdaopsfcent.fits (PSF estimate at science image center as a FITS image)
    - sciimlog.txt (log output from instrumental calibration pipeline)
    - scimrefdiffimg.fits.fz (difference image: science minus reference; fpack-compressed)
    - diffimgpsf.fits (PSF estimate for difference image as a FITS image)
    - diffimlog.txt (log output from image subtraction and extraction pipeline)
    - log.txt (overall system summary log from realtime pipeline)
    
    # Reference image (kind="ref"):
    -log.txt
    -refcov.fits
    -refimg.fits # (default)
    -refimlog.txt
    -refpsfcat.fits
    -refsexcat.fits
    -refunc.fits
    # Raw images (kind="raw")
    No Choice so suffix is ignored for raw data
    
    # Calibration (kind="cal")
    - None (#default) returns `caltype`.fits
    - log:            returns `caltype`log.txt
    - unc:            returns `caltype`unc.fits

''')
    args = parser.parse_args()

    names, sky_coords, jd_peaks = read_bts_table(args.table)
    os.makedirs(args.outdir, exist_ok=True)

    ilast = len(names) if args.ilast == 0 else args.ilast

    for i in range(args.ifirst, ilast):
        print(f'\nProgress: {args.ifirst} <= {i} < {args.ilast}')
        name = names[i]
        sky_coord = sky_coords[i]

        zquery = get_metadata(names[i], meta_dir=args.meta_dir)
        stamp_names = get_stampnames_from_meta(zquery, name, args.outdir)
        stamp_exist = exist(stamp_names)
        print(f'{name} : {np.sum(stamp_exist)} stamps exist out of {len(stamp_exist)}')
        if args.check_stamps:
            continue

        #-- Get image paths and check if they exist
        image_files, image_exist = get_filenames_from_meta(zquery, suffix=args.suffix)
        if np.sum(image_exist) == 0:
            print(f' ERROR: No images found for {name}')
            continue
        if np.sum(image_exist) < image_exist.size:
            print(' Warning: Not all images are downloaded, we will use those that exist')

        #-- If not overwriting, only produce stamps for those missing
        if args.overwrite == False:
            image_files = image_files[~stamp_exist]
            print(f'{name}: Producing the missing {image_files.size} stamps')

        #-- Produce stamps
        stamps = get_stamps(image_files, sky_coord, pixels=args.pixels)

        #-- Write them
        for stamp in stamps:
            write_stamp(name, stamp, directory=args.outdir)

if __name__ == '__main__' and len(sys.argv) > 1:
    main()