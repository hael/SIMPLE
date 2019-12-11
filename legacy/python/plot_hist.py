#!/usr/bin/env python3
##Plot the histogram of a MRC file
##  Created by Michael Eager, 2018
##
##  plot_hist.py is part of the Simple cryoEM software suite

import sys,os
import numpy as np
import matplotlib.pyplot as plt
import mrcfile
import argparse

def plot_histogram(img,nbins):

    hist, bin_edges = np.histogram(img, bins=nbins)
    bin_centers = 0.5*(bin_edges[:-1] + bin_edges[1:])

    plt.subplot(121)
    plt.plot(bin_centers, hist, lw=2)
#    plt.axvline(threshold  , color='r', ls='--', lw=2)
    plt.text(0.57, 0.8, 'Histogram', fontsize=20, transform = plt.gca().transAxes)

    plt.yticks([])
    plt.subplot(122)
    if img.ndim == 3 and img.shape[2] > 1:
        mid = np.floor(img.shape[2] / 2)
        plt.imshow(img[:,:,mid])
    else:
        plt.imshow(img)
    plt.show()

def check_valid(filename):
    import io
    output = io.StringIO()
    res= mrcfile.validate(filename, print_file=output)
    print(output.getvalue().strip())
    return res

if __name__ == '__main__':

    # Parse command line arguments and validate img directory

    parser = argparse.ArgumentParser(usage='''mydenoise.py -i dir ''',
                                     description='''
                                     ''')
    parser.add_argument('-i', '--inputfile',help='''Input file. Must be a
        image MRC file ''', required=True)
    parser.add_argument('-b', '--bins',help='''Number of bins in histogram ''')
    parser.add_argument(
        '-v', '--verbose', help='Verbose.', action="store_true")


    args = parser.parse_args()
    # Check input folder exists
    if not os.path.exists(args.inputfile):
        print('Error: File \'' + args.inputfile + '\' does not exist.')
        sys.exit(1)

    if not args.inputfile.endswith('mrc'):
        print('Error: File \'' + args.inputfile + '\' must be an MRC file.')
        sys.exit(1)

    if args.bins:
        nbins= int(args.bins)
    else:
        nbins= 256

    if not check_valid(args.inputfile):
        print('Error: File \'' + args.inputfile + '\' is an invalid MRC file.')
        sys.exit(1)

    mrc = mrcfile.open(args.inputfile, mode='r')
    #with mrcfile.open(args.inputfile) as mrc:
    print('Data size :', mrc.data.shape)
    # if mrc.is_single_image():
    #     img = np.squeeze(np.array(mrc.data, dtype=np.int16))
    # if mrc.is_image_stack():
    #     img = np.squeeze(np.array(mrc.data, dtype=np.int16))
    # if mrc.is_volume():
    #     print('Error: File \'' + args.inputfile + '\' does not exist.')
    #     sys.exit(1)
    if mrc.is_volume_stack():
        print('Error: Cannot process volume stack.')
        sys.exit(1)

    img = np.squeeze(np.array(mrc.data, dtype=np.int16))
    plot_histogram(img, nbins)
    mrc.close()
