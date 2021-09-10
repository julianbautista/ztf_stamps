import argparse
import time 

parser = argparse.ArgumentParser(description='Download ZTF images')
parser.add_argument('--table', type=str, required=True,
                    help='Path of metadata table in csv format')
parser.add_argument('--ncpu', type=int, default=20,
                    help='Number of machines for dask')
parser.add_argument('--overwrite', action="store_true", default=False,
                    help='If set, will re-download all images')
args = parser.parse_args()

def limit_numpy(nthreads=4):
    """ """
    import os
    threads = str(nthreads)
    #print(f"threads {threads}")
    os.environ["NUMEXPR_NUM_THREADS"] = threads
    os.environ["OMP_NUM_THREADS"] = threads
    os.environ["OPENBLAS_NUM_THREADS"] = threads
    os.environ["MKL_NUM_THREADS"] = threads
    os.environ["VECLIB_MAXIMUM_THREADS"] = threads
limit_numpy(4)

import dask
from dask_jobqueue import SGECluster
from dask.distributed import Client, LocalCluster

cluster = SGECluster(name="dask-worker",  walltime="02:00:00", 
                     memory='5GB', death_timeout=120, 
                     project="P_lsst", resource_spec='sps=1', 
                     cores=1, processes=1)
cluster.scale(args.ncpu) 
client = Client(cluster)

from ztfquery import query
zquery = query.ZTFQuery().from_metafile(args.table)
print(f'Found {len(zquery.metatable)} exposures')

t0 = time.time()
tt = zquery.get_data("sciimg.fits", show_progress=False, 
                     indexes=zquery.metatable.index,
                     client=client,
                     overwrite=args.overwrite)
t1 = time.time()
print('Time elapsed: {(t1-t0)/60:.3f} minutes')
client.shutdown()