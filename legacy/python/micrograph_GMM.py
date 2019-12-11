#!/usr/bin/env python
## Plot the estimated threshold, histogram and mask of an MRC file
##  Created by Michael Eager, 2018
##
##  micrograph_GMM.py is part of the Simple cryoEM software suite
import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
import mrcfile
import argparse


def plot_GMM(image, nbins):

    np.random.seed(1)
    n = 10
    l = image.shape
    im = np.zeros(l)
    points = l*np.random.random((2, n**2))
    im[(points[0]).astype(np.int), (points[1]).astype(np.int)] = 1
    im = ndimage.gaussian_filter(im, sigma=l[0]/(4.*n))

    mask = (im > im.mean()).astype(np.float)


    img = mask + 0.3*np.random.randn(*mask.shape)

    hist, bin_edges = np.histogram(img, bins=60)
    bin_centers = 0.5*(bin_edges[:-1] + bin_edges[1:])

    classif = GaussianMixture(n_components=2)
    classif.fit(img.reshape((img.size, 1)))

    threshold = np.mean(classif.means_)
    binary_img = img > threshold

    hist2, bin_edges2 = np.histogram(image, bins=nbins)
    bin_centers2 = 0.5*(bin_edges2[:-1] + bin_edges2[1:])

    classif2 = GaussianMixture(n_components=2)
    classif2.fit(image.reshape((image.size, 1)))

    threshold2 = np.mean(classif2.means_)
    binary_img2 = image > threshold2


    plt.figure(figsize=(11,8))

    plt.subplot(231)
    plt.imshow(img)
    plt.axis('off')
    plt.subplot(232)
    plt.plot(bin_centers, hist, lw=2)
    plt.axvline(threshold, color='r', ls='--', lw=2)
    plt.text(0.57, 0.8, 'histogram', fontsize=20, transform = plt.gca().transAxes)
    plt.yticks([])
    plt.subplot(233)
    plt.imshow(binary_img, cmap=plt.cm.gray, interpolation='nearest')
    plt.axis('off')


    plt.subplot(234)
    plt.imshow(image)
    plt.axis('off')
    plt.subplot(235)
    plt.plot(bin_centers2, hist2, lw=2)
    plt.axvline(threshold2, color='r', ls='--', lw=2)
    plt.text(0.57, 0.8, 'histogram', fontsize=20, transform = plt.gca().transAxes)
    plt.yticks([])
    plt.subplot(236)
    plt.imshow(binary_img2, cmap=plt.cm.gray, interpolation='nearest')
    plt.axis('off')

    plt.subplots_adjust(wspace=0.02, hspace=0.3, top=1, bottom=0.1, left=0, right=1)
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
    plot_GMM(img, nbins)
    mrc.close()
