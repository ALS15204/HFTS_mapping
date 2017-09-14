from map_gridding import *
import matplotlib, pylab
from multiprocessing import Pool
from numpy import *
import resource, sys

switch_backend('Agg')
targetcode = 'Carina_Nebula'
datapath = '/home/rwu/Obs_data/Herschel/'
#datapath = '../../../data/FTS/'
outputroot = '../../../results/' + targetcode + '/MC_v4/'
if not os.path.exists(outputroot):
    os.system('mkdir ' + outputroot)
    
# [OPTIONAL INPUTS]
## -- a file that gives the emission line information
# linefile = './FTS_lines.dat'
## -- WCS information for the output maps
##### uncomment the three lines bellow to define the WCS
# pixscale = 15.
# ractr, decctr = 0.0, 0.0
# ftsfov = 300.
## uncomment montecarlo if wish to simulate map-making procedure with
## 300 random initial conditions
# montecarlo = False
## turn pointcal off if one wishes to use extended-source calibrated data
# pointcal = True
# datasuffix = ''
####################################################################
def mc_process(imc):
    
    linenames = linetab[1:, 0]
    for thisline in linenames:
        outputpath = outputroot + thisline.replace('(', '').replace(')', '') \
                     + '/'
        if not os.path.exists(outputpath):
            os.system('mkdir ' + outputpath)
            
        testfile = outputpath + \
                   targetcode + '_' + \
                   thisline + '_conv_intmap_' + str(imc) + '.png'
        if os.path.exists(testfile):
            continue
        gridfig, cubehdu = foo.make_line_cube(thisline)

        cube_fig, inthdu = foo.intmeasure_linecube(cubehdu[0], thisline)
        intfig = plot_fits2darray(inthdu, thisline)

        cube_conv = foo.convolution_linecube(cubehdu, thisline)

        cube_conv_fig, int_conv_hdu \
            = foo.intmeasure_linecube(cube_conv, thisline)
        int_conv_fig = plot_fits2darray(int_conv_hdu, thisline)

        gridfig.savefig(outputpath + \
                        targetcode + '_' + \
                        thisline + '_grid_' + str(imc) + '.png', \
                        format = 'PNG')
        cubehdu.writeto(outputpath + \
                        targetcode + '_' + \
                        thisline + '_grid_' + str(imc) + '.fits', \
                        clobber = True)
        pylab.close(gridfig)
        cubehdu.close()
        cube_fig.savefig(outputpath + \
                         targetcode + '_' + \
                         thisline + '_cube_' + str(imc) + '.png', \
                         format = 'PNG')
        pylab.close(cube_fig)
        cube_conv.writeto(outputpath + \
                          targetcode + '_' + \
                          thisline + '_grid_conv_' + str(imc) + '.fits', \
                          clobber = True)
        cube_conv_fig.savefig(outputpath + \
                       targetcode + '_' + \
                       thisline + '_cube_conv.png', format = 'PNG')
        pylab.close(cube_conv_fig)
        inthdu.writeto(outputpath + \
                       targetcode + '_' + \
                       thisline + '_intmap_' + str(imc) + '.fits', \
                       clobber = True)
        inthdu.close()
        intfig.savefig(outputpath + \
                        targetcode + '_' + \
                        thisline + '_intmap_' + str(imc) + '.png', \
                       format = 'PNG')
        pylab.close(intfig)
        int_conv_hdu.writeto(outputpath + \
                             targetcode + '_' + \
                             thisline + '_conv_intmap_' + str(imc) + '.fits', \
                             clobber = True)
        int_conv_hdu.close()
        int_conv_fig.savefig(outputpath + \
                             targetcode + '_' + \
                             thisline + '_conv_intmap_' + str(imc) + '.png', \
                             format = 'PNG')
        pylab.close(int_conv_fig)
        
    return

foo = HFTS_grid(targetcode, datapath, montecarlo = True)
linetab = loadtxt('../data/FTS_lines.dat', delimiter = ',', \
                  dtype = 'str')
### prepare to assign obsids to each node
nloop = 300
allloop = arange(nloop)

pool = Pool(processes = 6)
mcarr = []
for i in range(len(allloop)):
    testfilename = '../../../results/Carina_Nebula/MC_v4/13CO14-13/' + \
                   targetcode + '_13CO(14-13)_conv_intmap_' + \
                   str(i) + '.png'
    if not os.path.exists(testfilename):
        mcarr.append(i)
mcarr = array(mcarr)
pool.map(mc_process, mcarr)
