# ZTF Stamps

This package is meant to download all images associated with a list of transists, given a time window. 
It has only been used at the Centre de Calcul IN2P3. An account and access to ztf data are required. 

1. First download the metadata associated with the transient list. The script will read the name of the object, ra, dec, and the mjd of peak flux and write the metadata associated to all available images between `phase-min` and `phase-max` in days. 

 `python bin/ztf_download_metadata.py --table BTS_transients_2021-07-05.csv --outdir metatables --phase-min 50 --phase-max 150 [--ifirst 0] [--ilast 10] [--overwrite]`

 You can add `--ifirst` and `--ilast` if you want to download only some of them. These are the indexes of the transient in list. 

 You can merge several metadata tables into a single one with 

 `python bin/ztf_merge_metadata.py --meta-dir metatables --outname ZTF_bigtable.csv `

2. Download the data from a metadata table. First login into a computing node (not a login node)

 `qlogin  -P P_lsst -l ct=15:00:00,vmem=10000M`

 Then launch the downloading script : 

 `python bin/ztf_download_images.py --table ZTF_bigtable.csv [--suffix sciimg.fits] [--ncpu 20] [--overwrite] [--check-images]`

 The `--check-images` option allows to only print the summary of how many images are already downloaded.

 To follow the progress, login in the same machine in a new terminal with the following command, replacing XXXX with the correct number:

 `ssh -J username@cca.in2p3.fr -L 8787:localhost:8787 username@ccwigeXXXX.in2p3.fr`

 Then open the following link in a browser 

 https://localhost:8787/status

4. Produce the stamps with the following command

 `python bin/ztf_produce_stamps.py --table BTS_transients_2021-07-05.csv --meta-dir metatables --outdir stamps [--ifirst 0] [--ilast 10] [--overwrite] [--check-stamps] [--suffix sciimg.fits] `



