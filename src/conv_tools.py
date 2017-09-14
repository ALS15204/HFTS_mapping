import os, pdb, sys, pyfits
from math import *
from numpy import *
from scipy import interpolate

## CONVOLUTION TOOLS ##
def beamprof_fwhm(beam2dmap):
    
    X = arange(beam2dmap.shape[0])
    Y = arange(beam2dmap.shape[1])
    Xgrid, Ygrid = meshgrid(X, Y)
    
    beampeak = (beam2dmap == nanmax(beam2dmap))
    Xpeak, Ypeak = Xgrid[beampeak], Ygrid[beampeak]

    beamhalf = (beam2dmap <= 0.5*nanmax(beam2dmap))
    Xhalf, Yhalf = Xgrid[beamhalf], Ygrid[beamhalf]

    beamfwhm = 2.*min(sqrt((Xhalf-Xpeak)**2. + (Yhalf-Ypeak)**2.))
    
    return beamfwhm

def make_gaussian_psf(fwhm, pixscale = 1., xsize = 257, ysize = 257):
    
    X = arange(xsize)
    Y = arange(ysize)
    Xgrid, Ygrid = meshgrid(X, Y)
    Xctr, Yctr = xsize/2, ysize/2
    fwhm_pix = fwhm/pixscale
    sigma_gaussian = fwhm_pix/2./sqrt(2.*log(2.))
    
    gaussian_psf = zeros((xsize, ysize))
    gaussian_psf = exp(-((Xgrid-Xctr)**2. + (Ygrid-Yctr)**2.)/ \
                       2./sigma_gaussian**2.)
        
    return gaussian_psf



def make_kernel(psfmap, refpsfmap):

    nbx, nby = refpsfmap.shape
    nbx_psf, nby_psf = psfmap.shape

    if nbx != nbx_psf:
        newpsf = zeros((nbx, nby))
        halfpsf = (nbx_psf-1)/2
        center_x = (nbx-1)/2
        newpsf[center_x-halfpsf: center_x+halfpsf+1, \
               center_x-halfpsf: center_x+halfpsf+1] = psfmap
    else:
        newpsf = psfmap.copy()

    fft_beamratio = fft.fftpack.fft2(refpsfmap)/ \
                    fft.fftpack.fft2(newpsf)

    idx = prod(refpsfmap.shape)
    idx = reshape(arange(idx), refpsfmap.shape)
    xx = mod(idx,nbx)-(nbx-1)/2
    yy = (idx/nby)-(nbx-1)/2
    dist = xx**2.+yy**2.
    freq = dist/double(nbx)
    fwhm_ideal = sqrt(beamprof_fwhm(refpsfmap)**2.-beamprof_fwhm(psfmap)**2.)
    
    fwhm_test = (arange(100)+1.)/10.    
    rx_ref = 100.
    tmplist = []
    for thisfwhm in fwhm_test[::-1]:
        nyquist_freq = 1./thisfwhm
        hanning = zeros(refpsfmap.shape, dtype = 'f')
        temp = where(freq < nyquist_freq)
        if len(temp) > 0:
            hanning[temp] = 0.5*(1+cos(pi*freq[temp]/nyquist_freq))

        hanning_shift = roll(roll(hanning, (nbx+1)/2,axis = 0), \
                             (nby+1)/2, \
                             axis = 1)
        outkernel = fft.fftpack.ifft2(hanning_shift*fft_beamratio)
        outkernel=roll(roll(outkernel, nbx/2, axis=0), \
                       nby/2, \
                       axis=1)
        outkernel=outkernel.real
        
        k1dx = outkernel[nbx/2:, nby/2]
        k1dx_slope = k1dx[1:]-k1dx[:-1]
        tmp_x = k1dx_slope > 0.
        k1dy = outkernel[nbx/2, nby/2:]
        k1dy_slope = k1dy[1:]-k1dy[:-1]
        tmp_y = k1dy_slope > 0.
        rx = min([min(arange(len(k1dx))[tmp_x]), \
                  min(arange(len(k1dy))[tmp_y])])
        if len(tmplist) >= 3:
            if rx < tmplist[-1] \
               and tmplist[-1] < tmplist[-2] \
               and tmplist[-2] < tmplist[-3]:
                break
        tmplist.append(rx)            
        if rx <= rx_ref:
            rx_ref = rx
            outkernel_ref = outkernel
        if rx > rx_ref:
            break

    k1dx = outkernel[nbx/2:, nby/2]
    k1dy = outkernel[nbx/2, nby/2:]
    rx = arange(len(k1dx))
    tvx = sum(array([abs(roll(k1dx, x) - k1dx) for x in rx]), axis = 0)
    tvy = sum(array([abs(roll(k1dy, x) - k1dy) for x in rx]), axis = 0)
    rx_cut = min([rx[(tvx == tvx.min())][0], rx[(tvy == tvy.min())][0]])
    X, Y = meshgrid(arange(outkernel_ref.shape[0]), \
                    arange(outkernel_ref.shape[1]))
    outkernel_ref[(((X-nbx/2)**2.+(Y-nby/2)**2.) > rx_cut**2.)] = 0.
    outkernel_ref[(outkernel < 0.)]  = 0.
    outkernel_ref = outkernel_ref/sum(outkernel_ref)

    return outkernel

# inpixsize and outpixsize should have the same angular units.
def img_resample(inarr, inpixsize, outpixsize, \
                 xyctr = None, newdim = None, normalize = True):

    inarr = inarr[1:-1,1:-1].copy()
    nkx, nky = inarr.shape
    X, Y = arange(nkx)+0.5, arange(nky)+0.5

    # Define the desired center from the input map. 
    # Xctr and Yctr can both be non-integers.
    if xyctr == None:
        Xctr, Yctr = nkx/2., nky/2.
    else:
        Xctr, Yctr = xyctr[0], xyctr[1]

    # Define the new dimension of the regridded map if not given
    if newdim == None:
        new_nkx, new_nky = array(inarr.shape)/(outpixsize/inpixsize)
        new_nkx = int(new_nkx)
        new_nky = int(new_nky)
        if mod(new_nkx,2) == 0:
            new_nkx = new_nkx-1
            new_nky = new_nky-1
    else:
        new_nkx, new_nky = newdim[0], newdim[1]

    new_Xctr, new_Yctr = (new_nkx)/2., (new_nky)/2.
    
    halfnewpix = outpixsize/inpixsize*0.5
    new_X = (arange(new_nkx)-new_Xctr+0.5)*outpixsize/inpixsize +  \
            Xctr
    new_Y = (arange(new_nky)-new_Yctr+0.5)*outpixsize/inpixsize + \
            Yctr
    valuesum \
        = nansum(inarr[int(new_Y[0]-halfnewpix):int(new_Y[-1]+halfnewpix), \
                       int(new_X[0]-halfnewpix):int(new_X[-1]+halfnewpix)])
    f = interpolate.RectBivariateSpline(X,Y,inarr,kx=3,ky=3)

    outarr = f(new_Y,new_X)
    outarr[where(outarr == -999.)] = nan
    # Conserve the sum of the initial array
    outarr = outarr.copy()/nansum(outarr.copy())* \
             valuesum*(outpixsize/inpixsize)**(-2.)
    if normalize:
        outarr = outarr.copy()/nansum(outarr.copy())
        
    return outarr
