import argparse
import time
import os 
import numpy as np

from ztfquery import query


parser = argparse.ArgumentParser(description='Download ZTF images')
parser.add_argument('--table', type=str, required=True,
                    help='Path of metadata table in csv format')
parser.add_argument('--ncpu', type=int, default=20,
                    help='Number of machines for dask')
parser.add_argument('--overwrite', action="store_true", default=False,
                    help='If set, will re-download all images')
parser.add_argument('--check-images', action="store_true", default=False,
                    help='If set, will only check if images exist')
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

def limit_numpy(nthreads=4):
    """ This is required to avoid multithreading within a cpu"""
    import os
    threads = str(nthreads)
    #print(f"threads {threads}")
    os.environ["NUMEXPR_NUM_THREADS"] = threads
    os.environ["OMP_NUM_THREADS"] = threads
    os.environ["OPENBLAS_NUM_THREADS"] = threads
    os.environ["MKL_NUM_THREADS"] = threads
    os.environ["VECLIB_MAXIMUM_THREADS"] = threads

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
    print(f'{len(files_sci)} sci files, {np.sum(files_sci_exist)} exist')
    
    return files_sci, files_sci_exist

def main():
    limit_numpy(4)

    args = parser.parse_args()


    zquery = query.ZTFQuery().from_metafile(args.table)
    image_files, image_exist = get_filenames_from_meta(zquery, suffix=args.suffix)

    print(f'Found {len(image_files)} exposures in {args.table}, {sum(image_exist)} exist already')

    if args.check_images:
        return

    print('Creating dask client...')
    from dask_jobqueue import SGECluster
    from dask.distributed import Client
    cluster = SGECluster(name="dask-worker",  walltime="02:00:00", 
                        memory='5GB', death_timeout=120, 
                        project="P_lsst", resource_spec='sps=1', 
                        cores=1, processes=1)
    cluster.scale(args.ncpu) 
    client = Client(cluster)



    print('Downloading data...')
    t0 = time.time()
    tt = zquery.get_data(args.suffix, show_progress=False, 
                        indexes=zquery.metatable.index,
                        client=client,
                        overwrite=args.overwrite)
    t1 = time.time()
    print(f'Time elapsed in download: {(t1-t0)/60:.3f} minutes')
    print('Shutting down dask client...')
    client.shutdown()
    t2 = time.time()
    print(f'Time elapsed {(t2-t1)/60:.3f} minutes')
    
if __name__ == '__main__':
    main()