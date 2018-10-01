#!/usr/bin/env python3
"""Denoise an MRC file

  Created by Michael Eager, 2018

  micrograph_denoise.py is part of the Simple cryoEM software suite

"""
import sys, os
import numpy as np
import skimage
import scipy
import matplotlib.pyplot as plt
import argparse
from skimage import data, img_as_float
from skimage.restoration import denoise_nl_means
from skimage.measure import compare_psnr
import mrcfile

from skimage.morphology import disk
from skimage.filters import rank
from skimage.filters.rank import entropy



def normalise(image):
    """Normalise numpy image
    Args:
        image (np.ndarray): Input image

    Returns: Normalised [0,1] numpy ndarray
    """
    nimg = (image - min(image.flatten())).astype(np.float32)
    nimg = nimg / max(nimg.flatten())
    return nimg


def means_filter_perct(image):
    selem = disk(20)
    return rank.mean_percentile(image, selem=selem, p0=.1, p1=.9)


def create_mask(image,frac):
    sz= np.array(image.shape,np.int32)
    downscaled_img = scipy.misc.imresize(image, size=sz/frac, interp='cubic')
    imstd = scipy.ndimage.generic_filter(downscaled_img, np.std, 10)
    imstd = means_filter_perct(imstd)
    mask = imstd > skimage.filters.threshold_otsu(imstd)
    return  scipy.misc.imresize(mask, size=sz, interp='cubic'), imstd

def estimate_H(image,mask):
    """Estimate the h-factor using the stdev of the background
    Args:
        image (np.ndarray): Input image
        mask  (np.ndarray): Background mask
     Returns: Float value, h
    """
    if not mask:
        mask, = create_mask(image,4)
    return ndimage.standard_deviation(image,mask=mask)

def getMRCimage(filename):
    """
    Args:
        filename (string): MRC filename
    Returns:
        Numpy ndarray image

    """

    mrc = mrcfile.open(args.filename, mode='r')
    print('Data size :', mrc.data.shape)
    if mrc.is_volume_stack():
        print('Error: Cannot process volume stack.')
        sys.exit(1)
    img = np.squeeze(np.array(mrc.data, dtype=np.float))
    mrc.close()
    return img


def mydenoise(image, sigma=0., patchsz=5, boxsz=13, mode=3, doplot=True):
    """ Non-local means denoise

    Args:
        image (np.ndarray) : Input image
        sigma (float)      : Noise estimate variance
        patchsz (int)      : Patch size  (total patch width=2*patchsz+1)
        boxsz (int)        : Box size (total box width= 2*boxsz+1)
    Returns:
        Numpy ndimage (denoised image)
    """
    noisy = (image - min(image.flatten())).astype(np.float32)
    noisy = noisy / max(noisy.flatten())
    patch_kw = dict(patch_size=patchsz,      # 5x5 patches
                patch_distance=boxsz,        # 13x13 search area
                multichannel=True)
    if sigma <= 0.0:
        sigma_est = estimate_H(noisy,[])
    else:
        sigma_est = sigma

    #denoise = denoise_nl_means(noisy,  sigma, **patch_kw)
    if mode == 1:
        # slow algorithm
        denoise = denoise_nl_means(noisy, h=1.15 * sigma_est, fast_mode=False,
                           **patch_kw)
    elif mode == 2:
        # slow algorithm, sigma provided
        denoise = denoise_nl_means(noisy, h=0.8 * sigma_est, sigma=sigma_est,
                            fast_mode=False, **patch_kw)
    elif mode == 3:
        # fast algorithm
        denoise = denoise_nl_means(noisy, h=0.8 * sigma_est, fast_mode=True,
                                **patch_kw)
    elif mode == 4:
        # fast algorithm, sigma provided
        denoise = denoise_nl_means(noisy, h=0.6 * sigma_est, sigma=sigma_est,
                                 fast_mode=True, **patch_kw)
    else:
        print(' mydenoise invalid mode ')
        sys.exit(1)

    #fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(8, 6),
    #                       sharex=True, sharey=True)
    psnr = compare_psnr(noisy,denoise)
    if doPlot:
        fig, ax = plt.subplots(ncols=2, figsize=(8, 4), sharex=True, sharey=True,
                               subplot_kw={'adjustable': 'box-forced'})
        ax[0].imshow(noisy)
        ax[0].axis('off')
        ax[0].set_title('noisy')
        ax[1].imshow(denoise)
        ax[1].axis('off')
        ax[1].set_title('non-local means')
        ax[1].text(0.57, 0.8, ['PSNR ' + str(psnr)], fontsize=20, transform = plt.gca().transAxes)
        fig.tight_layout()
        plt.show()
    return denoise

def check_valid(filename):
    """Run MRCfile's validation script on filename
    """
    import io
    output = io.StringIO()
    res= mrcfile.validate(filename, print_file=output)
    print(output.getvalue().strip())
    return res


if __name__ == '__main__':

    # Parse command line arguments and validate img directory

    parser = argparse.ArgumentParser(usage='''mydenoise.py -i <file|dir> ''',
                                     description='''
                                     ''')
    parser.add_argument('-i', '--input',help='''Input directory. Must be an
        image directory containing MRC files ''', required=True)
    parser.add_argument(
        '-v', '--verbose', help='Verbose.', action="store_true")


    args = parser.parse_args()
    # Check input folder exists
    if not os.path.exists(args.input):
        print('Error: Folder \'' + args.input + '\' does not exist.')
        sys.exit(1)


    if os.path.isdir(args.input):
        files = os.listdir(args.input)
        if args.verbose:
            print(files)
        mrcfiles = [f for f in files if f.endswith('mrc')]
    elif os.path.isfile(args.input):
        mrcfiles = [args.input ]

    for mfile in mrcfiles:
        if not check_valid(mfile):
            print('Error: File \'' + mfile + '\' is an invalid MRC file.')
            sys.exit(1)

        img = getMRCimage(mfile)
        mask,imgstd = create_mask(image,4)
        h_est= ndimage.standard_deviation(image,mask=mask)
        denoised_img = mydenoise(img) #sigma=0., patchsz=5, boxsz=13, mode=3, doplot=True):
        mrc = mrcfile.new(['denoised_'+ os.path.basename(mfile)])
        mrc.set_data(denoised_img)
        mrc.close()

"""


    imdata = skimage.data.imread('_intg00003-1.jpg', as_grey=True, flatten=True)
    noisy = (imdata - min(imdata.flatten()))
    noisy = noisy / max(noisy.flatten())
    imshow(noisy)



help(scipy.std)
scipy.std(noisy.flatten())
sqrt(0.289)

denoise = denoise_nl_means(noisy, 5, 21, 0.08, multichannel=True)

scipy.misc.imsave('denoised-1.jpg', denoise)


graddenoise = np.zeros(denoise.shape)
scipy.ndimage.gaussian_gradient_magnitude(denoise, 20, graddenoise, mode='reflect')
scipy.misc.imsave('graddenoised-1.jpg', graddenoise)


denoise = denoise_nl_means(noisy, 15, 31, 0.16, multichannel=True)

scipy.ndimage.gaussian_gradient_magnitude(denoise, 9,graddenoise, mode='reflect')
imshow(graddenoise)
scipy.ndimage.standard_deviation(denoise, 9,stddenoise, mode='reflect')
stddenoise = np.zeros(denoise.shape)
scipy.ndimage.standard_deviation(denoise, 9,stddenoise, mode='reflect')
help(scipy.ndimage.standard_deviation)
help(scipy.ndimage.grey_opening)
odenoise = scipy.ndimage.grey_opening(denoise,size=(5,5))
imshow(odenoise-denoise)
imshow(denoise-odenoise)
odenoise = scipy.ndimage.grey_opening(denoise,size=(15,15))
imshow(denoise-odenoise)
odenoise = scipy.ndimage.grey_opening(denoise,size=(25,25))
imshow(denoise-odenoise)
imshow(np.abs(denoise-odenoise))
imshow(np.abs(odenoise-denoise))
odenoise = scipy.ndimage.grey_opening(denoise,size=(51,51))
imshow(np.abs(odenoise-denoise))
imshow(odenoise)
odenoise = scipy.ndimage.grey_opening(denoise,size=(5,5))
imshow(odenoise)
odenoise = scipy.ndimage.grey_opening(denoise,size=(15,15))
imshow(odenoise)
odenoise = scipy.ndimage.grey_opening(denoise,size=(3,3))
imshow(odenoise)
imshow(odenoise-denoise)
imshow(denoise-odenoise)
imdown = scipy.misc.imresize(imdata, size=(1024,1024), interp='bilinear')
imshow(imdown)
imdown = scipy.misc.imresize(imdata, size=(1024,1024), interp='nearest')

imdown = (imdown - min(imdown.flatten()))/(max(imdown.flatten()) - min(imdown.flatten()))

imdown = scipy.misc.imresize(imdata, size=(512,512), interp='cubic')
imshow(imdown)
imdown = (imdown - min(imdown.flatten()))/(max(imdown.flatten()) - min(imdown.flatten()))
denoise = denoise_nl_means(imdown, 3,7, 0.25, multichannel=True)


    # Thresholding
    a=scipy.ndimage.gaussian_filter(1-denoise, 3)
    othresh = skimage.filters.threshold_otsu(a)
    bgmask = a < othresh

    """
