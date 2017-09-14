from map_gridding import *
import matplotlib

switch_backend('Agg')
targetcode = 'Horsehead'
datapath = '/home/rwu/Obs_data/Herschel/'
#datapath = '../../../data/FTS/'
outputpath = '../../../results/' + targetcode + '/'
if not os.path.exists(outputpath):
    os.system('mkdir ' + outputpath)
    
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

foo = HFTS_grid(targetcode, datapath)
linetab = loadtxt('../data/FTS_lines.dat', delimiter = ',', dtype = 'str')

linenames = linetab[1:, 0]
for thisline in linenames:
    gridfig, cubehdu = foo.make_line_cube(thisline)
    
    cube_fig, inthdu = foo.intmeasure_linecube(cubehdu[0], thisline)
    intfig = plot_fits2darray(inthdu, thisline)
    
    cube_conv = foo.convolution_linecube(cubehdu, thisline)
    
    cube_conv_fig, int_conv_hdu = foo.intmeasure_linecube(cube_conv, thisline)
    int_conv_fig = plot_fits2darray(int_conv_hdu, thisline)

    gridfig.savefig(outputpath + \
                    targetcode + '_' + \
                    thisline + '_grid.png', format = 'PNG')
    cubehdu.writeto(outputpath + \
                    targetcode + '_' + \
                    thisline + '_grid.fits', clobber = True)
    close(gridfig)
    cubehdu.close()
    cube_fig.savefig(outputpath + \
                     targetcode + '_' + \
                     thisline + '_cube.png', format = 'PNG')
    close(cube_fig)
    cube_conv.writeto(outputpath + \
                      targetcode + '_' + \
                      thisline + '_grid_conv.fits', clobber = True)
    cube_conv_fig.savefig(outputpath + \
                   targetcode + '_' + \
                   thisline + '_cube_conv.png', format = 'PNG')
    close(cube_conv_fig)
    inthdu.writeto(outputpath + \
                   targetcode + '_' + \
                   thisline + '_intmap.fits', clobber = True)
    inthdu.close()
    intfig.savefig(outputpath + \
                    targetcode + '_' + \
                    thisline + '_intmap.png', format = 'PNG')
    close(intfig)
    int_conv_hdu.writeto(outputpath + \
                         targetcode + '_' + \
                         thisline + '_conv_intmap.fits', clobber = True)
    int_conv_hdu.close()
    int_conv_fig.savefig(outputpath + \
                         targetcode + '_' + \
                         thisline + '_conv_intmap.png', format = 'PNG')
    close(int_conv_fig)
