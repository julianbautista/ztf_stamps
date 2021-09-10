# ztf_stamps

This package is meant to download all images associated with a list of transists. 
It works at the CC, if you have an account and ztf access. 

1. First download the metadata associated with the transient list

`python bin/ztf_download_metadata.py --table BTS_transients_2021-07-05.csv --outdir metatables --phase-min 50 --phase-max 150 --ifirst 0 --ilast 130 --overwrite`

You can add `--ifirst` and `--ilast` options if you want to download only some of them. These are the indexes in the table. 

2. Merge the tables you want into a single one using 

`python bin/ztf_merge_metadata.py ...`

3. Download the data using 

`python bin/ztf_download_images.py ...`

