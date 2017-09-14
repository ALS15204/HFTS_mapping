import os, sys, pdb, pyfits
from numpy import *
from scipy.optimize import leastsq
from astropy import wcs

def read_linelist(linefile):
    
    linetbl = loadtxt(linefile, delimiter = ',', dtype = 'str')
    linehdr = map(str.strip, linetbl[0])
    linedat = linetbl[1:]

    outdict = []
    for i in range(len(linedat)):
        thisdict = dict.fromkeys(linehdr)
        for j in range(len(linehdr)):
            if linehdr[j] == 'name':
                thisdict[linehdr[j]] \
                    = (linedat[i, j].strip()).replace('J', '')
            else:
                thisdict[linehdr[j]] = float(linedat[i, j])
        outdict.append(thisdict)

    return outdict

def load_thisline(linename, linefile = None):

    if linefile == None:
        linelist = read_linelist('FTS_lines.dat')
    else:
        linelist = read_linelist(linefile)

    idx = arange(len(linelist))
    tmp = (array([xx['name'] for xx in linelist]) == linename)
    thisline = array(linelist)[tmp]

    return thisline[0]

def get_nlines(np):

    nlines = (np-3)/2

    return nlines
    
def continuum_level(x, p):

    y_cont = (p[0]*x**2. + p[1]*x + p[2])

    return y_cont

def one_line(x, p):
    
    y_oneline = abs(p[0])*sinc((x-p[1])/linewidth)
    
    return y_oneline

def minfunc(p, x ,y ,err):

    ymod = continuum_level(x, p[:3])
    nlines = get_nlines(len(p))
    for i in range(nlines):
        ymod += one_line(x, p[3+2*i:5+2*i])

    tomin = (y-ymod)/err
    
    return tomin

# Estimate the initial parameters.
def initial_values(x, y, z = None, \
                   linename = None, linefile = None, \
                   nosidelines = None):

    if linefile == None:
        linelist = read_linelist('FTS_lines.dat')
    else:
        linelist = read_linelist(linefile)
        
    all_linename = array([xx['name'] for xx in linelist])
    all_linefreq = array([xx['freq'] for xx in linelist])/(1.+z)
    
    thisline_mrk = (all_linename == linename)
    thislinefreq = all_linefreq[thisline_mrk][0]-min(x)

    sidelines_inrange = (x[0]+1. < all_linefreq) & \
                        (all_linefreq < x[-1]-1.) & \
                        (all_linename != linename)
    n_sidelines = 0
    if True in sidelines_inrange and not nosidelines:
        sidelinefreq_inrange = all_linefreq[sidelines_inrange]-min(x)
        n_sidelines = len(sidelinefreq_inrange)
        
    p0 = []
    # parabola 2nd-order coefficient
    p0.append(1e-8)

    # parabola 1st-order coefficient
    p0.append(1e-8)

    # parabola 0th-order coefficient
    p0.append(median(y))

    mrk = (thislinefreq-1. < x-min(x)) & (x-min(x) < thislinefreq+1.)
    p0.append(max(y[mrk]))
    
    p0.append(thislinefreq + \
              0.5*random.uniform(low = -1., high = 1.))

    if n_sidelines > 0:
        for iline in range(n_sidelines):
            mrk = (sidelinefreq_inrange[iline]-1. < x-min(x)) & \
                  (x-min(x) < sidelinefreq_inrange[iline]+1.) 
            p0.append(max(y[mrk]))
            p0.append(sidelinefreq_inrange[iline] + \
                      0.5*random.uniform(low = -1., high = 1.))

    return p0

# spec: spectrum trimmed to the range for the fitting.
def fit_one_line(xghz, yint, linename, z = None, \
                 error = None, linefile = None, \
                 nosidelines = None):

    x_0 = min(xghz)
    y_0 = min(abs(yint))
    ngood = 10

    x = xghz
    x_fit = x - x_0
    y_fit = yint/y_0

    if error == None:
        error = ones(len(y_fit))
    else:
        error = error/y_0
        
    fa = (x_fit, y_fit, error)
    p0 = initial_values(x, y_fit, z = z, \
                            linename = linename, linefile = linefile, \
                            nosidelines = nosidelines)
    
    if len(p0) > 5:
        ntry = 300
    else:
        ntry = 50

    itry = 1
    igood = 0
    while itry <= ntry and igood < ngood:
        ylinefit_tmp = leastsq(minfunc,
                               p0, \
                               args = fa, \
                               xtol = 1e-8, \
                               factor = 0.1, \
                               epsfcn = 1e-6
        )
        chisq_tmp = sum(minfunc(ylinefit_tmp[0], x_fit, y_fit, error)**2.) \
                    /(len(y_fit)-len(p0)-1.)
        if len(p0) > 5:
            mrk = (abs(ylinefit_tmp[0][4]-p0[4]) < 1.) & \
                  (abs(ylinefit_tmp[0][6]-p0[6]) < 1.)
        else:
            mrk = (abs(ylinefit_tmp[0][4]-p0[4]) < 1.5)
            
        p0 = initial_values(x, y_fit, z = z, \
                            linename = linename, linefile = linefile, \
                            nosidelines = nosidelines)
        itry += 1
        if mrk:
            igood += 1
            if 'minf_val' not in vars():
                minf_val = chisq_tmp
                ylinefit = copy(ylinefit_tmp)
            else:
                if chisq_tmp < minf_val:
                    minf_val = chisq_tmp
                    ylinefit = copy(ylinefit_tmp)
        else:
            if itry == ntry:
                minf_val = chisq_tmp
                ylinefit = copy(ylinefit_tmp)
            else:
                continue
        
    ylinefit[0][2] = (ylinefit[0][0]*x_0**2. - \
                      ylinefit[0][1]*x_0 + \
                      ylinefit[0][2])*y_0
    ylinefit[0][1] = (ylinefit[0][1] - 2.*ylinefit[0][0]*x_0)*y_0
    ylinefit[0][0] = (ylinefit[0][0])*y_0
    nlines = get_nlines(len(p0))
    for i in range(nlines):
        ylinefit[0][3+2*i] = ylinefit[0][3+2*i]*y_0
        ylinefit[0][4+2*i] = ylinefit[0][4+2*i]+x_0
    return ylinefit

linewidth = 0.37733*pi # GHz
