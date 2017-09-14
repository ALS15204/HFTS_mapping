from numpy import *
import pdb, pyfits
def target_obsinfo(targetcode):

    obs_dict = array([{'code': '30Dor', \
                       'name': 'LMC-30 Doradus', \
                       'nobs': 3, \
                       'obsid': ['1342219550', \
                                 '1342257932', \
                                 '1342262908'], \
                       'sampletype': ['intermediate', \
                                      'intermediate', \
                                      'intermediate'], \
                       'z': 0.002}, \
                      {'code': 'Carina_Nebula', \
                       'name': 'Carina Nebula-Trumpler 14', \
                       'nobs': 6, \
                       'obsid': ['1342256376', \
                                 '1342256377', \
                                 '1342262913', \
                                 '1342262914', \
                                 '1342262915', \
                                 '1342262916'], \
                       'sampletype': ['intermediate', \
                                      'intermediate', \
                                      'intermediate', \
                                      'intermediate', \
                                      'intermediate', \
                                      'intermediate'], \
                       'z': 0.0}, \
                      {'code': 'Horsehead', \
                       'name': 'Horse Head', \
                       'nobs': 3, \
                       'obsid': ['1342216880', \
                                 '1342216881', \
                                 '1342216882', \
                                 '1342228732'], \
                       'sampletype': ['full', \
                                      'full', \
                                      'full', \
                                      'full'], \
                       'z': 0.0}, \
                      {'code': 'M17', \
                       'name': 'M17', \
                       'nobs': 2, \
                       'obsid': ['1342228703', \
                                 '1342231045'], \
                       'sampletype': ['full', \
                                      'sparse'], \
                       'z': 0.0}, \
                      {'code': 'M83', \
                       'name': 'M83', \
                       'nobs': 1, \
                       'obsid': ['1342212345'], \
                       'sampletype': ['full'], \
                       'z': 0.001711}, \
                      {'code': 'Cen_A', \
                       'hipevsn': ['11.0.2825'], \
                       'name': 'Centaurus A', \
                       'nobs': 1, \
                       'obsid': ['1342204036'], \
                       'sampletype': ['full'], \
                       'z': 0.001823}, \
                      {'code': 'NGC7023', \
                       'name': 'NGC 7023', \
                       'nobs': 3, \
                       'obsid': ['1342198923', \
                                 '1342201204', \
                                 '1342201205'], \
                       'sampletype': ['full', \
                                      'full', \
                                      'full'], \
                       'z': 0.0}, \
                      {'code': 'M82', \
                       'name': 'M82', \
                       'nobs': 1, \
                       'obsid': ['1342208388'], \
                       'sampletype': ['full'], \
                       'z': 0.000677}, \
                      {'code': 'NGC4418', \
                       'name': 'NGC4418', \
                       'nobs': 1, \
                       'obsid': ['1342210848'], \
                       'sampletype': ['intermediate'], \
                       'z': 0.007268}, \
                      {'code': 'He2-10', \
                       'name': 'Henize 2-10', \
                       'nobs': 1, \
                       'obsid': ['1342245083'], \
                       'sampletype': ['intermediate'], \
                       'z': 0.002912}, \
                      {'code': 'IC10', \
                       'name': 'IC 10', \
                       'nobs': 1, \
                       'obsid': ['1342246982'], \
                       'sampletype': ['sparse'], \
                       'z': -0.001161}, \
                      {'code': 'LMC-Diffuse', \
                       'name': 'LMC-Diffuse', \
                       'nobs': 1, \
                       'obsid': ['1342256081'], \
                       'sampletype': ['intermediate'], \
                       'z': 0.002}, \
                      {'code': 'LMC-N11', \
                       'name': 'LMC-N11', \
                       'nobs': 1, \
                       'obsid': ['1342257915'], \
                       'sampletype': ['intermediate'], \
                       'z': 0.002}, \
                      {'code': 'LMC-N157', \
                       'name': 'LMC-N157', \
                       'nobs': 1, \
                       'obsid': ['1342262907'], \
                       'sampletype': ['intermediate'], \
                       'z': 0.002}, \
                      {'code': 'LMC-N158', \
                       'name': 'LMC-N158', \
                       'nobs': 1, \
                       'obsid': ['1342257914'], \
                       'sampletype': ['intermediate'], \
                       'z': 0.002}, \
                      {'code': 'LMC-N159', \
                       'name': 'LMC-N159', \
                       'nobs': 1, \
                       'obsid': ['1342259066'], \
                       'sampletype': ['intermediate'], \
                       'z': 0.001}, \
                      {'code': 'LMC-N160', \
                       'name': 'LMC-N160', \
                       'nobs': 1, \
                       'obsid': ['1342262909'], \
                       'sampletype': ['intermediate'], \
                       'z': 0.002}, \
                      {'code': 'LMC-N180', \
                       'name': 'LMC-N180', \
                       'nobs': 1, \
                       'obsid': ['1342257913'], \
                       'sampletype': ['intermediate'], \
                       'z': 0.002}, \
                      {'code': 'LMC-N44', \
                       'name': 'LMC-N44', \
                       'nobs': 1, \
                       'obsid': ['1342262905'], \
                       'sampletype': ['intermediate'], \
                       'z': 0.002}, \
                      {'code': 'LMC-N79', \
                       'name': 'LMC-N79', \
                       'nobs': 2, \
                       'obsid': ['1342245112', \
                                 '1342245113'], \
                       'sampletype': ['full', \
                                      'full'], \
                       'z': 0.002}, \
                      {'code': 'SMC-N76', \
                       'name': 'SMC-N76', \
                       'nobs': 1, \
                       'obsid': ['1342270032'], \
                       'sampletype': ['intermediate'], \
                       'z': 0.001}, \
                      {'code': 'SMC-NGC249', \
                       'name': 'SMC-NGC249', \
                       'nobs': 1, \
                       'obsid': ['1342270031'], \
                       'sampletype': ['intermediate'], \
                       'z': 0.001}])
    output = None
    for i in range(len(obs_dict)):
        if obs_dict[i]['code'] == targetcode:
            output = obs_dict[i]

    return output
def fts_detinfo():

    det_dict = [{'module' : 'SSW', \
                 'center_det' : 'SSWD4', \
                 'all_det' : ['SSWB2','SSWB3','SSWB4', \
                              'SSWC2','SSWC3','SSWC4','SSWC5', \
                              'SSWD2','SSWD3','SSWD4','SSWD6', \
                              'SSWE2','SSWE3','SSWE4','SSWE5', \
                              'SSWF2','SSWF3']}, \
                {'module' : 'SLW', \
                 'center_det' : 'SLWC3', \
                 'all_det' : ['SLWB2','SLWB3', \
                              'SLWC2','SLWC3','SLWC4', \
                              'SLWD2','SLWD3']}]
    return det_dict

def fts_obsmode_info(obsmode):
    
    if obsmode == 'full':
        njiggle = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
        jigctr = 6

    if obsmode == 'intermediate':
        njiggle = [0,1,2,3]
        jigctr = 0

    if obsmode == 'sparse':
        njiggle = [0]
        jigctr = 0

    return njiggle, jigctr

# Always use the latest spire_cal files.
def load_calinfo(datapath = None):
    
    caldir = datapath + 'spire_cal_14_3/' \
             'herschel.spire.ia.dataset.SpecBeamParam/'
    calfits = pyfits.open(caldir+'SCalSpecBeamParam_HR_20120218_v4.fits')
    beamdir = datapath + 'spire_cal_14_3/' \
              'herschel.spire.ia.dataset.SpecBeamProf/'
    slwbeam = pyfits.open(beamdir+'SCalSpecBeamProf_SLW_v5.fits')
    sswbeam = pyfits.open(beamdir+'SCalSpecBeamProf_SSW_v5.fits')
    beamfreq = [
        sswbeam[1].header['CRVAL3']+ \
        arange(sswbeam[1].header['NAXIS3'])* \
        sswbeam[1].header['CDELT3'], \
        slwbeam[1].header['CRVAL3']+ \
        arange(slwbeam[1].header['NAXIS3'])* \
        slwbeam[1].header['CDELT3']
    ]
    beamprof = [sswbeam[1].data, slwbeam[1].data]

    slwbeam.close()
    sswbeam.close()
    
    return calfits, beamfreq, beamprof
