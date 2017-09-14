############################ REVISION NOTE ############################
## --  Although not recommended, the map-making now support
##     regridding observations processed with different versions
##     of hipe and sampling modes.
## --  For the current version, only square maps are allowed.
## --  Convolution is now added to this library.
#######################################################################
# Load python intrinsic libraries
import os, pdb, sys, pyfits
from math import *
from types import *
from numpy import *
from matplotlib.pyplot import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy import convolution

# Load my own dependencies
from fitline_leastsq import *
from conv_tools import *
from rwu_obs_record import *

############### Define some universal constants ####################

sr_to_ascsq = (3600.*180./pi)**2.
ascsq_to_sr = 1./sr_to_ascsq

######################### DEPENDENCIES #############################
def load_observations(targetcode, datapath = '', datasuffix = ''):

    exfits = pyfits.open('../data/1342256377_avgspectrum_HR_unapod_0_0.fits')
    exfits.verify('silentfix')
    bolarr = ['SSW', 'SLW']
    obj_info = target_obsinfo(targetcode)
    nobs = obj_info['nobs']
    datapath = datapath + targetcode + '/'
    
    corr_fits = pyfits.open('../data/ExtCal_corr.fits')
    
    fits_out = []
    for obs in range(nobs):
        njiggle, jigctr = fts_obsmode_info(obj_info['sampletype'][obs])
        fits_jig = []
        for i in range(len(njiggle)):
            fits_jig.append(None)
        fits_out.append(fits_jig)
    
    ractr_all = []
    decctr_all = []
    ## read-in all fits files for map-gridding
    for obs in range(nobs):
        njiggle, jigctr = fts_obsmode_info(obj_info['sampletype'][obs])
        ## read the central bolometer fits file for the RA and DEC info
        thisdatapath = datapath + obj_info['obsid'][obs] + '/level2/'
        if obj_info['sampletype'][obs] != 'sparse':
            for bol in bolarr:
                thisdatapath = datapath + obj_info['obsid'][obs] + \
                               '/level2/'+ \
                               'HR_' + bol + '_spectrum2d/'
                thispathfiles = os.listdir(thisdatapath)
                corr_factor = corr_fits[bol].data['factor']
                for ifile in thispathfiles:
                    thisfits = pyfits.open(thisdatapath + ifile)
                    all_bolnames = [thisfits[1].data[x][7] \
                                    for x in arange(len(thisfits[1].data))]
                    all_jigid = [thisfits[1].data[x][9] \
                                 for x in arange(len(thisfits[1].data))]
                    for ijig in njiggle:
                        for iext in range(len(exfits)):
                            bolname = exfits[iext].name
                            tmp = (array(all_bolnames) == bolname) & \
                                  (array(all_jigid) == ijig)
                            if not (True in tmp):
                                continue
                            idx = arange(len(thisfits[1].data))[tmp][0]
                            wave \
                                = pyfits.Column(name = 'wave', \
                                                format = 'D', \
                                                array \
                                                = thisfits[1].data[idx][0])
                            flux \
                                = pyfits.Column(name = 'flux', \
                                                format = 'D', \
                                                array \
                                                = thisfits[1].data[idx][1]* \
                                                corr_factor)
                            error \
                                = pyfits.Column(name = 'error', \
                                                format = 'D', \
                                                array \
                                                = thisfits[1].data[idx][2])
                            all_columns = pyfits.ColDefs([wave, flux, error])
                            thistablehdu \
                                = pyfits.BinTableHDU.from_columns(all_columns)

                            thistablehdu.header.set('RA', \
                                                    thisfits[1].data[idx][3])
                            thistablehdu.header.set('DEC', \
                                                    thisfits[1].data[idx][4])
                            thistablehdu.name = bolname
                            if fits_out[obs][ijig] is None:
                                fits_out[obs][ijig] \
                                    = pyfits.HDUList(thistablehdu)
                            else:
                                fits_out[obs][ijig].append(thistablehdu)

                            if bolname == 'SLWC3' and ijig == jigctr:
                                ractr_all.append(thisfits[1].data[idx][3])
                                decctr_all.append(thisfits[1].data[idx][4])
                thisfits.close()
        else:
            thisdatapath = datapath + obj_info['obsid'][obs] + \
                           '/level2/HR_spectrum_ext/'
            thispathfiles = os.listdir(thisdatapath)
            for ifile in thispathfiles:
                thisfits = pyfits.open(thisdatapath + ifile)
                all_bolnames = array([thisfits[x].name \
                                      for x in arange(len(thisfits))])
                tmp = array([('SLW' in x or 'SSW' in x) \
                             for x in all_bolnames])
                all_bolnames = all_bolnames[tmp]
                for iext in all_bolnames:
                    corr_factor = corr_fits[iext[:3]].data['factor']
                    wave \
                        = pyfits.Column(name = 'wave', \
                                        format = 'D', \
                                        array \
                                        = thisfits[iext].data['wave'])
                    flux \
                        = pyfits.Column(name = 'flux', \
                                        format = 'D', \
                                        array \
                                        = thisfits[iext].data['flux']* \
                                        corr_factor)
                    error \
                        = pyfits.Column(name = 'error', \
                                        format = 'D', \
                                        array \
                                        = thisfits[iext].data['error'])
                    all_columns = pyfits.ColDefs([wave, flux, error])
                    thistablehdu \
                        = pyfits.BinTableHDU.from_columns(all_columns)

                    thistablehdu.header.set('RA', \
                                            thisfits[iext].header['RA'])
                    thistablehdu.header.set('DEC', \
                                            thisfits[iext].header['DEC'])
                    thistablehdu.name = iext
                    if fits_out[obs][0] is None:
                        fits_out[obs][0] \
                            = pyfits.HDUList(thistablehdu)
                    else:
                        fits_out[obs][0].append(thistablehdu)

                    if iext == 'SLWC3':
                        ractr_all.append(thisfits[iext].header['RA'])
                        decctr_all.append(thisfits[iext].header['DEC'])
                thisfits.close()
                        
                    
                                
    return fits_out, array(ractr_all), array(decctr_all)

def make_dummy_wcs(ra, dec, nxy, pixscale):
    """
    Function: make the WCS for the output maps
    """
    nxdiv=2.
    nydiv=2.

    xc=floor(nxy[0]-nxy[0]/nxdiv)+1.
    yc=floor(nxy[1]-nxy[1]/nydiv)+1.
    dummy=wcs.WCS()
    dummy.wcs.crpix=array([xc,yc])
    dummy.wcs.crval=array([ra,dec])
    dummy.wcs.cdelt=array([-pixscale,pixscale])*1./3600.
    dummy.wcs.ctype=['RA---TAN','DEC--TAN']

    return dummy

def prepare_hdr(inputwcs, arraydim, zarr = None, obsname = ''):

    outhdr = inputwcs.to_header()
    outhdr.set('SIMPLE', \
               value = True, \
               comment = 'created by Astropy', \
               before = 'WCSAXES')
    outhdr.set('BITPIX', \
               value = -64, \
               comment = 'IEEE single precision floating point', \
               before = 'WCSAXES')
    outhdr.rename_keyword('WCSAXES', 'NAXIS')
    outhdr['NAXIS'] = len(arraydim)
    outhdr.set('NAXIS1', \
               value = arraydim[-2], \
               before = 'CRPIX1')
    outhdr.set('NAXIS2',\
               value = arraydim[-1], \
               before = 'CRPIX1')
    if len(arraydim) == 2:
        outhdr.set('EXTEND', \
                   value = True, \
                   after = 'NAXIS2')
        outhdr.set('BUNIT', \
                   value = 'W m-2 sr-1', \
                   comment = 'Integrated Surface Brightness', \
                   before = 'LONPOLE')
    if len(arraydim) == 3:
        if zarr is None:
            zarr = arange(len(arraydim[0]))
        outhdr.set('NAXIS3', \
                   value = arraydim[0], \
                   before='CRPIX1')
        outhdr.set('CRVAL3', \
                   value = zarr[0], \
                   after = 'CRVAL2')
        outhdr.set('CRPIX3', \
                   value = 1.0, \
                   after = 'CRPIX2')
        outhdr.set('CDELT3', \
                   value = average(zarr[1:]-zarr[:-1]), \
                   after = 'CDELT2', \
                   comment = '[GHz]')
        outhdr.set('CUNIT3', \
                   value = 'GHz', \
                   after = 'CUNIT2')
        outhdr.set('CTYPE3', \
                   value = 'Frequency', \
                   after = 'CTYPE2')
        outhdr.set('EXTEND', \
                   value = True, \
                   after = 'NAXIS3')
        outhdr.set('BUNIT', \
                   value = 'W m-2 GHz-1 sr-1', \
                   comment = 'Surface Brightness', \
                   before = 'LONPOLE')
    outhdr.set('BSCALE', \
               value = 1.0, \
               after = 'BUNIT')       
    outhdr.set('PCOUNT', \
               value = 0, \
               after = 'EXTEND')
    outhdr.set('GCOUNT', \
               value = 1, \
               after = 'EXTEND')
    if obsname != '':
        outhdr.set('OBJECT', \
                   value = obsname, \
                   before = 'LONPOLE')

    imghdr=pyfits.Header(outhdr)
    imghdr.update(outhdr.items())

    return imghdr

def get_linecomp(x, p):

    cont = continuum_level(x, p[:3])
    nlines = get_nlines(len(p))
    alllines = []
    for i in range(nlines):
        alllines.append(one_line(x, p[3+2*i:5+2*i]))

    return cont, alllines

def get_lineint(p):

    nlines = get_nlines(len(p))

    allint = []
    for i in range(nlines):
        allint.append(abs(p[3+2*i])*1e9*linewidth)

    return allint

## MAP PLOTTING TOOLS ##
def plot_fits2darray(intmapfits, maptitle):

    intwcs = wcs.WCS(intmapfits[0].header)
    nxy = intmapfits[0].data.shape

    rect = (0.1, 0.1, 0.8, 0.8)
    mapfig, mapfig_ax = plot_grid_opendev(maptitle, \
                                          rect = rect, \
                                          fignum = 21)
    ext=0.,1.,0.,1.
    mapfig_ax.imshow(intmapfits[0].data, \
                     extent = ext, \
                     origin = 'lower', \
                     cmap = get_cmap('hot'), \
                     interpolation = 'none')

    mapfig_ax = plot_grid_addcorrd(intwcs, mapfig_ax, nxy)
    colorbar_ax = mapfig.add_axes([0.92,0.1,0.02,0.8])
    clr_map = matplotlib.cm.hot
    norm \
        = matplotlib.colors.Normalize(vmin = \
                                      nanmin(intmapfits[0].data/1e-9), \
                                      vmax = \
                                      nanmax(intmapfits[0].data/1e-9))
    cb = matplotlib.colorbar.ColorbarBase(colorbar_ax, \
                                          cmap = clr_map, \
                                          norm = norm, \
                                          orientation='vertical', \
                                          format = '%4.1f')
    colorbar_ax.tick_params(labelsize=18.)
    colorbar_ax.text(-0.6, 0.55, \
                     '(10$^{-9}$ W m$^{-2}$ sr$^{-1}$)', \
                     rotation = 90.)
    return mapfig

def plot_grid_addcorrd(mapwcs, map_dev, map_dim, nbox = 4.):

    array_tics = arange(nbox+1)
    xtics = ytics = array_tics/nbox

    xytics = list()
    for i in range(len(xtics)):
        xytics.append(array([xtics[i], ytics[i]]))
        
    xytics = array(xytics)
    axxy = xytics*map_dim-array([0.5, 0.5])
    axradec = mapwcs.wcs_pix2world(axxy[::-1], 0)
    xtlbl = list()
    ytlbl = list()
    for i in range(len(axradec)):
        c = SkyCoord(axradec[i][0], axradec[i][1], \
                     unit = (u.degree, u.degree), \
                     frame = 'icrs')
        xtlbl.append(str(int(c.ra.hms[0]))+'h\n'+ \
                     str(c.ra.hms[1])+'m\n'+ \
                     '{:5.2f}'.format(c.ra.hms[2])+'s')
        ytlbl.append(str(int(c.dec.dms[0]))+'d\n'+ \
                     str(-c.dec.dms[1])+'m\n'+ \
                     '{:5.2f}'.format(-c.dec.dms[2])+'s')
        map_dev.set_xticks(xtics)
        map_dev.set_xticklabels(xtlbl[::-1], \
                                fontsize = 'large')
        map_dev.set_yticks(ytics)
        map_dev.set_yticklabels(ytlbl[::-1], \
                                fontsize = 'large', x = -0.02)
    return map_dev

def plot_grid_onegrid(plotdevice, axrect, \
                      xr = None, yr = None, \
                      gridannotate = None):

    dev_subax = plotdevice.add_axes(axrect)
    if xr != None:
        dev_subax.set_xlim(xr)
    if yr != None:
        dev_subax.set_xlim(yr)
    dev_subax.set_xticklabels([])
    dev_subax.set_yticklabels([])
    if gridannotate != None:
        dev_subax.annotate(gridannotate, [0.1, 1.0])

    return dev_subax

def plot_grid_opendev(plot_title, \
                      rect = (0.1, 0.1, 0.8, 0.8), \
                      fignum = 0):

    gridfig = figure(fignum, figsize = [12., 12.], \
                     linewidth = 3.0, frameon = True, dpi = 100)
    gridfig.text(0.45, 0.95, plot_title, fontsize = 'x-large')
    gridfig_ax = gridfig.add_axes(rect)

    return gridfig, gridfig_ax

class HFTS_grid(object):
    
    def __init__(self, targetcode, datapath, \
                 linefile = '../data/FTS_lines.dat', montecarlo = False, \
                 pointcal = True, datasuffix = '', \
                 pixscale = -1.0, ftsfov = -1.0, \
                 ractr = -1.0, decctr = -1.0, \
                 conv_res = 43.):
                 
        self.targetcode = targetcode
        self.datapath = datapath
        self.linefile = linefile
        self.montecarlo = montecarlo
        self.pointcal = pointcal
        self.datasuffix = datasuffix
        self.pixscale = pixscale
        self.ftsfov = ftsfov
        self.ractr = ractr
        self.decctr = decctr
        self.allfits, self.ractr_all, self.decctr_all \
            = load_observations(targetcode, \
                                datapath = datapath, \
                                datasuffix = datasuffix)
        self.calfits, self.allbeamfreq, self.allbeamprof \
            = load_calinfo(datapath = datapath)
        self.beampixscale = 1.
        self.conv_res = conv_res
        
    def view_bolspec(self, bolnames = [], \
                     obsid = [-1], jigid = [0], **plot_kwrds):

        obsid = array(obsid)
        obj_info = target_obsinfo(self.targetcode)
        all_obsid = array(obj_info['obsid'])
        if obsid == [-1]:
            obsid = all_obsid

        fts_bolinfo = fts_detinfo()
        if bolnames == []:
            bolnames = append(fts_bolinfo[0]['all_det'], \
                              fts_bolinfo[1]['all_det'])
        bolfig = figure(22, figsize = [15., 8.])
        bolfig_ax = subplot(111)
        for iobs in obsid:
            tmp = (all_obsid == iobs)
            idx = arange(len(self.allfits))[tmp][0]
            for ijig in jigid:
                for ibol in bolnames:
                    bolfig_ax.plot(self.allfits[idx][ijig][ibol].data['wave'], \
                                   self.allfits[idx][ijig][ibol].data['flux'], \
                                   label = iobs + ' ' + ibol, **plot_kwrds)
        
        bolfig_ax.set_xlim(array(bolfig_ax.get_xlim())*array([1., 1.2]))
        bolfig_ax.legend(loc = 0, bbox_to_anchor = (0.8, 1.1), \
                         fontsize = 8., ncol = 2, fancybox = True)

        return bolfig            
        
    def select_beamprof(self, freq):

        if freq < 1000.:
            beamfreq = self.allbeamfreq[1]
            beamprof = self.allbeamprof[1]
        else:
            beamfreq = self.allbeamfreq[0]
            beamprof = self.allbeamprof[0]

        beammrk = where(abs(freq-beamfreq) == abs(freq-beamfreq).min())
        beam2dmap = beamprof[beammrk][0]

        return beam2dmap
    
    def prepare_conv_kernels(self, targetcode, linename):
        
        obj_info = target_obsinfo(self.targetcode)
        z = obj_info['z']
        thisline = load_thisline(linename, linefile = self.linefile)
        beamprof = self.select_beamprof(thisline['freq']/(1.+z))
        refbeam = make_gaussian_psf(self.conv_res)
        kernel = make_kernel(beamprof, refbeam)
        kernel_resamp = img_resample(kernel, \
                                     self.beampixscale, \
                                     self.pixscale)
        return kernel_resamp
    
    #########################  MAP-MAKING ####################################
    ## targetcode: (str)                                                     #
    ##             A specific code assigned to the target as recorder in the #
    ##             "fts_obs_record.py".                                      #
    ## linename: (str)                                                       #
    ##           has to follow the same convention given in the <linefile>.  #
    ## use_pipe_error: (bool)                                                #
    ##                Whether to use the pipeline uncertainty for gridding   #
    ##                or not.                                                #
    ## montecarlo: (bool)                                                    #
    ##             Whether this operation is for a MC experiment or not.     #
    ## linefile: (str)                                                       #
    ##           Filename specifying the record of all lines included in     #
    ##           the range                                                   #
    ## (refer to README for the required format of the <linefile>)           #
    ##########################################################################
    def make_line_cube(self, linename):

        # PREPARE GLOBAL INFORMATION #
        thisline = load_thisline(linename, linefile = self.linefile)
        ## Read in information about the requested target
        ### [!!! CHANGE 1 NEEDED BGN !!!]
        ### This definition of data-path has forced the user to structure
        ### their data as "<targetcode>/HIPE.xx.x.xxxx/unapodize/", which 
        ### is not convenient.
        obj_info = target_obsinfo(self.targetcode)
        nobs = obj_info['nobs']
        z = obj_info['z']
        ### [!!! CHANGE 1 NEEDED END !!!]

        ## Calibration uncertainty:
        ### The calibration uncertainty for FTS is 10% (Swinyard et al. 2014)
        calerr = 0.1

        # Read bolometer info
        fts_bolinfo = fts_detinfo()

        ## Take average of <ractr_all> and <decctr_all> to be
        ## the default ceter of the map (if not defined).
        if self.ractr == -1.0 and self.decctr == -1.0:
            self.ractr = 0.5*(min(self.ractr_all) + max(self.ractr_all))
            self.decctr = 0.5*(min(self.decctr_all) + max(self.decctr_all))
            ## Define default pixel scale (15 arcsec)
        if self.pixscale == -1.0:
            self.pixscale = 15.
            ## Define total map size
        if self.ftsfov == -1.0:
            ra_size = self.ractr_all.ptp()*3600.*cos(self.decctr*pi/180.) \
                      + 105.*2.
            dec_size = self.decctr_all.ptp()*3600. + 105.*2.
            self.ftsfov \
                = self.pixscale*ceil(max([dec_size, ra_size])/self.pixscale)
            
        ### calculate the necessary pixel numbers
        nxy = int(self.ftsfov/self.pixscale)*ones(2, dtype = 'int')
        
        dwcs = make_dummy_wcs(self.ractr, self.decctr, nxy, self.pixscale)

        ## Define the corner locations of a given pixel w.r.t. its center.
        xyc = array([array([-0.5,-0.5]), \
                     array([0.5,-0.5]), \
                     array([0.5,0.5]), \
                     array([-0.5,0.5])])

        rect = 0.1, 0.1, 0.8, 0.8
        bolspec, bolspec_ax \
            = plot_grid_opendev(linename, rect = rect, \
                                fignum = 11)
        bolspec_subax_all = []
        ystackmax=0.

        for obs in range(nobs):

            if self.montecarlo:
                random.seed()
                mcpar_cal = random.randn(1)[0]
            else:
                mcpar_cal = 0.
                
            if thisline['freq'] <= 1000.:
                this_bolinfo = fts_bolinfo[1]
            if thisline['freq'] >= 1000.:
                this_bolinfo = fts_bolinfo[0]

            all_bol = this_bolinfo['all_det']
            obs_linefreq = thisline['freq']/(1.+z)
            beamprof = self.select_beamprof(obs_linefreq)
            beamarea = sum(beamprof)*ascsq_to_sr*self.beampixscale
            njiggle = len(self.allfits[obs])
            for ijig in range(njiggle):
                fts = self.allfits[obs][ijig]
                ftsnames = array([x.name for x in fts])
                for idet in range(len(this_bolinfo['all_det'])):
                    ftsmrk = (ftsnames == this_bolinfo['all_det'][idet])
                    ftsidx = arange(len(fts))
                    try:
                        data = fts[ftsidx[ftsmrk][0]].data
                    except IndexError:
                        continue
                    calfactor \
                        = \
                        self.calfits[this_bolinfo['all_det'][idet]]. \
                        data['pointConv']*1e-26/(beamarea)

                    hdr = fts[ftsidx[ftsmrk][0]].header
                    radec = array([hdr['RA'],hdr['DEC']])
                    xy = dwcs.wcs_world2pix([radec],0)
                    xi = int(round(xy[0][0]))
                    yi = int(round(xy[0][1]))
                    datar = (data['wave'] >= obs_linefreq - 15.) & \
                            (data['wave'] <= obs_linefreq + 15.)

                    if self.montecarlo:
                        random.seed()
                        mcpar_std = random.randn(len(data[datar]))
                    else:
                        mcpar_std = 0.

                    inidata = None

                    inidata = data[datar].copy()
                    inidata['flux'][:] \
                        = inidata['flux'].copy()+ \
                        (mcpar_cal*calerr*inidata['flux'].copy()) + \
                        (mcpar_std*inidata['error'].copy())

                    if self.pointcal:
                        inidata['flux'][:] = inidata['flux'].copy()* \
                                             calfactor[datar]

                    pfit = fit_one_line(inidata['wave'], \
                                        inidata['flux'], \
                                        linename, \
                                        z = z, linefile = self.linefile)
                    nlines = get_nlines(len(pfit[0]))
                    cont, lines = get_linecomp(inidata['wave'], pfit[0])

                    if self.montecarlo:
                        errmaskr=(inidata['wave'] <= obs_linefreq - 5.) | \
                            (inidata['wave'] >= obs_linefreq + 5.)
                        errmed = median((inidata['flux']-cont)[errmaskr])
                        errsig = std((inidata['flux']-cont)[errmaskr])
                        errmc = abs(errmed + \
                                    errsig*random.randn(len(inidata['flux'])))
                        inidata['error'][:] = errmc.copy()
                        inidata['error'][:] \
                            = sqrt((mcpar_cal*calerr*inidata['flux'])**2. \
                                   + (mcpar_std*inidata['error'])**2.)

                    cont_rmv_flux = inidata['flux'].copy()-cont
                    cont_slines_rmv_flux = cont_rmv_flux
                    if len(lines) > 1:
                        for iline in range(nlines-1):
                            cont_slines_rmv_flux = cont_slines_rmv_flux - \
                                                   lines[iline+1]

                    if max(cont_slines_rmv_flux) > ystackmax:
                        ystackmax = max(cont_slines_rmv_flux)

                    if obs == ijig == idet == 0:
                        xmin = inidata['wave'].min()
                        xmax = inidata['wave'].max()
                        mapcube = zeros([len(inidata['wave']), \
                                         nxy[1], nxy[0]])
                        mapcube[:] = nan
                        errorcube = zeros([len(inidata['wave']), \
                                           nxy[1], nxy[0]])
                        errorcube[:] = nan
                        area_sum = zeros([len(inidata['wave']), \
                                          nxy[1], nxy[0]])

                    for pix in range(len(xyc)):
                        # Note that there are three conventions here:
                        # one: Python array (origin = (0, 0))
                        # two: Image coordinate (origin = (0.5, 0.5))
                        # three: astropy wcs python (origin = (-0.5, -0.5))
                        # the computation here is following convention one.

                        # convert from convention three to one
                        # since xc and yc are given in convention three.
                        xc = int(round(xyc[pix][0]-0.5))
                        yc = int(round(xyc[pix][1]-0.5))

                        # upper-right corner of the imaginary pixel around
                        # the pointing (convention one)
                        xur = int(floor(xy[0][0]+1.))
                        yur = int(floor(xy[0][1]+1.))

                        # mapidx is in convention one
                        mapidx = (yur+yc,xur+xc)
                        if (mapidx[0] > nxy[0]-1) | (mapidx[1] > nxy[1]-1) | \
                           (mapidx[0] < 0) | (mapidx[1] < 0):
                            continue

                        # calculate the area in convention one
                        area = abs(xy[0][0]+1.+xyc[pix][0]-(xur))* \
                               abs(xy[0][1]+1.+xyc[pix][1]-(yur))
                        area_sum[:, mapidx[0], mapidx[1]] \
                            = area_sum[:, mapidx[0], mapidx[1]].copy()+ \
                            area**2./inidata['error']**2.

                        if isnan(mapcube[0, mapidx[0],mapidx[1]]):

                            mapcube[:, mapidx[0],mapidx[1]] = 0.
                            errorcube[:, mapidx[0],mapidx[1]] = 0.

                        mapcube[:, mapidx[0],mapidx[1]] \
                            = \
                            mapcube[:, mapidx[0],mapidx[1]].copy()+ \
                            area*cont_slines_rmv_flux/inidata['error']**2.
                        errorcube[:, mapidx[0],mapidx[1]] \
                            = \
                            errorcube[:, mapidx[0],mapidx[1]].copy()+ \
                            area/inidata['error']**2.

                    subrect = rect[0] + (xi)*(rect[2]/nxy[0]*1.), \
                              rect[1] + (yi)*(rect[3]/nxy[1]*1.), \
                              rect[2]/nxy[0]*1., \
                              rect[3]/nxy[1]*1.
                    bolspec_subax \
                        = plot_grid_onegrid(bolspec, subrect, \
                                            xr = [xmin, xmax])
                    bolspec_subax.plot(inidata['wave'], \
                                       cont_slines_rmv_flux)
                    bolspec_subax_all.append(bolspec_subax)

        for i in range(len(bolspec_subax_all)):
            bolspec_subax_all[i].set_ylim(array([-0.1, 1.1])*ystackmax)

        bolspec_ax = plot_grid_addcorrd(dwcs, bolspec_ax, nxy)
        mapcube=mapcube.copy()/errorcube.copy()
        errorcube=sqrt(area_sum.copy())/errorcube.copy()

        outhdr = prepare_hdr(dwcs, mapcube.shape, \
                             zarr = inidata['wave'], \
                             obsname = obj_info['name'])
        outhdulist = pyfits.HDUList()
        outhdulist.append(pyfits.PrimaryHDU(mapcube, header = outhdr))
        outhdulist.append(pyfits.PrimaryHDU(errorcube, header = outhdr))

        return bolspec, outhdulist

    def convolution_linecube(self, cubefits, linename):
        
        ## Prepare convolution kernels
        kernelfile = '../kernels/kernel_' + \
                     self.targetcode + '_' + \
                     linename + '_to_gaussian_' + \
                     str(self.conv_res) + '.fits'
        
        if not os.path.exists(kernelfile):
            kernel = self.prepare_conv_kernels(self.targetcode, linename)
            kernel_hdu = pyfits.PrimaryHDU(kernel)
            kernel_hdu.writeto(kernelfile, clobber = True)
        else:
            kernelfits = pyfits.open(kernelfile)
            kernel = kernelfits[0].data

        datacube = cubefits[0].data
        datacube_conv = zeros(datacube.shape)
        
        temp = where(isnan(datacube[0, :, :]) == True)
        nxy = datacube[0, :, :].shape
        mask = zeros(nxy)
        mask[-1,:] = mask[0,:] = mask[:,-1] = mask[:,0] = 1.
        for i in range(len(temp[0])):
            if temp[0][i] < nxy[0]-1 and temp[1][i] < nxy[1]-1:
                mask[temp[0][i],temp[1][i]] = 1.
                mask[temp[0][i]+1,temp[1][i]+1] = 1.
                mask[temp[0][i]+1,temp[1][i]-1] = 1.
                mask[temp[0][i]-1,temp[1][i]-1] = 1.
                mask[temp[0][i]-1,temp[1][i]+1] = 1.
                
        for i in range(datacube_conv.shape[0]):
            onelayer \
                = convolution.convolve_fft(datacube[i, :, :], \
                                           kernel, \
                                           normalize_kernel=True)
            tmp = (isnan(datacube[i, :, :]) == True)
            onelayer[tmp] = nan
            onelayer[(mask == 1.)] = nan
            datacube_conv[i, :, :] = onelayer
            
        cubefits[0].data = datacube_conv
        
        return cubefits[0]
    
    def intmeasure_linecube(self, cube_primhdu, linename):
        
        obj_info = target_obsinfo(self.targetcode)

        obswcs = wcs.WCS(cube_primhdu.header, naxis = 2)
        int_cube = cube_primhdu.data
        freq = cube_primhdu.header['CRVAL3'] + \
               cube_primhdu.header['CDELT3']* \
               (arange(cube_primhdu.header['NAXIS3'])- \
                (cube_primhdu.header['CRPIX3']-1.))

        nxy = int_cube[0, :, :].shape
        intmap = zeros(nxy, dtype = 'float')
        interr = zeros(nxy, dtype = 'float')

        rect = 0.1, 0.1, 0.8, 0.8
        gridspec, gridspec_ax = plot_grid_opendev(linename, \
                                                  rect = rect, \
                                                  fignum = 13)
        yrmax = nanmax(int_cube)
        xr = [freq.min(), freq.max()]
        yr = [-0.1*yrmax, 1.1*yrmax]
        for ix in range(nxy[0]):
            for iy in range(nxy[1]):
                if isnan(int_cube[:, iy, ix][0]):
                    continue
                else:
                    linefit = fit_one_line(freq, int_cube[:, iy, ix], \
                                           linename, \
                                           z = obj_info['z'], \
                                           linefile = self.linefile, \
                                           nosidelines = True)
                    linemod = continuum_level(freq, linefit[0][:3]) + \
                              one_line(freq, linefit[0][3:])
                    all_int = get_lineint(linefit[0])
                    intmap[iy, ix] = all_int[0]
                    interr[iy, ix] \
                        = sqrt(sum(((int_cube[:, iy, ix] - linemod)* \
                                    cube_primhdu.header['CDELT3']*1e9)**2.))
                    
                    subrect = rect[0] + (ix)*(rect[2]/nxy[0]*1.), \
                              rect[1] + (iy)*(rect[3]/nxy[1]*1.), \
                              rect[2]/nxy[0]*1., \
                              rect[3]/nxy[1]*1.
                    gridspec_subax \
                        = plot_grid_onegrid(gridspec, subrect, xr = xr)
                    gridspec_subax.plot(freq, int_cube[:, iy, ix], \
                                        linewidth = 3., color = 'red')
                    gridspec_subax.plot(freq, linemod, color = 'cyan')
                    gridspec_subax.set_ylim(yr)

        tmp = intmap == 0.
        intmap[tmp] = nan
        interr[tmp] = nan
        outhdr = prepare_hdr(obswcs, intmap.shape, \
                             obsname = obj_info['name'])
        outhdulist = pyfits.HDUList()
        outhdulist.append(pyfits.PrimaryHDU(intmap, header = outhdr))
        outhdulist.append(pyfits.PrimaryHDU(interr, header = outhdr))

        gridspec_ax = plot_grid_addcorrd(obswcs, gridspec_ax, nxy)
        return gridspec, outhdulist
        
