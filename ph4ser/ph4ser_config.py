"""
Configuration of different template parameters for self-calibration and imaging.
This is intended to be used as a first trial of self-calibration.
"""


visibility_info = {'path':"/media/sagauga/galnet/LIRGI_Sample/VLA-Archive/A_config/23A-324/C_band"
                          "/MCG+12-02-001/autoselfcal/",
                   'field':'MCG12-02-001',
                   'vis_name':'MCG12-02-001.calibrated.avg',
                   'savename':'_A_C_sf'
       }


FIELD = ''
SPWS = ''
ANTENNAS = ''
refantmode = 'flex' # 'strict' or 'flex' -- flex may be better in most situations.
refant = ''  # If '', will compute automatically -- this is VERY recommended. 
minblperant = 2 # Minimum number of baselines per antenna
solnorm = False # True or False (globally used). If '', will use False for phases and True for
             # amplitudes.
             
quiet = True # WSClean logging (printing stuff during deconvolution)
do_additional_images = False # Not working yet, leave to False.
run_mode = 'terminal'        # Not tested, leave to 'terminal' for now 

new_phasecentre = None
# new_phasecentre = '09:36:30.804 +61.18.30.641' #UGC5101 L band

#VLA
receiver = 'C'
instrument = 'EVLA' # 'EVLA' or 'eM'
if instrument == 'eM':
       refantmode = 'flex'

#number of channels-out for wsclean's MFS imager.
if instrument == 'EVLA':
       nc = 4 #number of bandwidth split during convolution (number of sub-band WSClean images)
if instrument == 'eM':
       nc = 3 #number of bandwidth split during convolution (number of sub-band WSClean images)

negative_arg='no-negative'
# negative_arg='negative'
steps = [
    'startup',          # create directory structure, start variables and clear visibilities.
#     'fov_image',      # create a FOV image
    'save_init_flags',  # save (or restore) the initial flags
    'statwt',           # run statwt on the initial data
    'autoflag_init',    # run rflag on the initial data -- use for L and S band data.
    'test_image',       # create a test image
    'select_refant',    # select reference antenna
    'p0',               # initial test  of selfcal, phase only (p)
    'p1',               # redo phase-only selfcal (if enough flux density); ignores p0
    'p2',               # continue phase-only selfcal (incremental)
    'ap1',              # amp-selfcal (ap); uses p0 or (p1 and p2)
    'split_trial_1',    # split the data after first trial (and run wsclean)
#     'autoflag_final',   # run rflag on the final data (use for L and S band data)
    'report_results'    # report results of first trial
]


cell_sizes_JVLA = {'L':'0.2arcsec',
                   'S':'0.1arcsec',
                   'C':'0.06arcsec',
                   'X':'0.04arcsec',
                   'Ku':'0.02arcsec',
                   'K':'0.01arcsec',
                   'Ka':'0.01arcsec'}

cell_sizes_eMERLIN = {'L':'0.03arcsec',
                      'C':'0.008arcsec'}



taper_sizes_eMERLIN = {'L':'0.3arcsec',
                       'C':'0.06arcsec'}

taper_sizes_JVLA = {'L':'2.0arcsec',
                    'S':'1.0arcsec',
                    'C':'0.6arcsec',
                    'X':'0.4arcsec',
                    'Ku':'0.2arcsec', #Ku-A
                    # 'Ku':'1.0arcsec', #Ku-C
                    'K':'0.1arcsec',
                    # 'Ka':'1.0arcsec' #Ka-C
                    'Ka':'0.1arcsec' #Ka-A
                    }


if instrument == 'eM':
    cell_size = cell_sizes_eMERLIN[receiver]
    taper_size = taper_sizes_eMERLIN[receiver]
if instrument == 'EVLA':
    cell_size = cell_sizes_JVLA[receiver]
    taper_size = taper_sizes_JVLA[receiver]

init_parameters = {'fov_image': {'imsize': 1024*8,
                                'cell': '0.25arcsec',
                                'basename': 'FOV_phasecal_image',
                                'FIELD_SHIFT': None,
                                'niter': 1000,
                                'robust': 0.5},
                  'test_image': {'imsize': int(1024*3),
                                 'imsizey': int(1024*3),
                                 'FIELD_SHIFT': None,
                                 'cell': cell_size,
                                 'prefix': 'test_image',
                                 'uvtaper': [''],
                                 'niter': 10000,
                                 'robust': 0.5 if receiver in ('K', 'Ka') or instrument == 'eM' else 0.0}
                  }

global_parameters = {'imsize': init_parameters['test_image']['imsize'],
                     'imsizey': init_parameters['test_image']['imsizey'],
                     'FIELD_SHIFT': init_parameters['test_image']['FIELD_SHIFT'],
                     'cell': init_parameters['test_image']['cell'],
                     'nsigma_automask' : '3.0',
                     'nsigma_autothreshold' : '1.5',
                     'mask_grow_iterations': 2,
                     'uvtaper' : [''],
                     'with_multiscale' : True,
                     'scales' : 'None',
                     'niter':100000}

general_settings = {'timebin_statw': '12s',#timebin for statwt
                    'statwt_statalg': 'chauvenet',
                    'calwt': False,#calibrate/update the weights during applycal?
                    'calwt_ap': True,#calibrate/update the weights during ap-applycal?
                    'applymode_p': 'calflag',
                    'applymode_ap': 'calflag',
                    'timebin':None, #not yet implemented
                    'channel_width': None,#not yet implemented
                    'do_average_time': False,#not yet implemented
                    'do_average_freq': False,#not yet implemented
                    'correlations':'RR,LL',
                    'new_phasecentre': new_phasecentre,
                    'force_combine_spw' : True, #force the combination of spectral windows for faint sources. 
                    'allow_combine_spw' : False, #allow the combination of spectral windows for bright sources.
                    'allow_tapper' : False, #allow the use of tapering.
                    }


"""
Selfcal parameters to be used for very faint sources, 
with a total integrated flux density lower than 10 mJy.
"""
params_very_faint = {'name': 'very_faint',
                     'p0': {'robust': 0.5,
                            'solint': 'inf' if instrument == 'eM' else ('240s' if receiver in ('K', 'Ka', 'Ku') else '120s'),
                            # 'solint' : '120s',
                            # 'sigma_mask': 6.0 if receiver in ('K', 'Ka', 'Ku') or instrument == 'eM' else 15.0,#C-Config
                            'sigma_mask': 6.0 if receiver in ('K', 'Ka', 'Ku') or instrument == 'eM' else 12.0,#A-Config
                            'mask_grow_iterations': 4,
                            'combine': 'spw' if general_settings['force_combine_spw'] else '',
                            # 'combine': '',# if joint array configs/instruments
                            'gaintype': 'T',
                            'calmode': 'p',
                            # 'minsnr': 0.75 if instrument == 'eM' else 1.0,
                            'minsnr': 0.1 if instrument == 'eM' else 0.1, #testing
                            'spwmap': [], #leavy empty here. It will be filled later if combine='spw'
                            'nsigma_automask' : '3.0',
                            'nsigma_autothreshold' : '1.5',
                            'uvtaper' : [''],
                            'with_multiscale' : True,
                            'scales': 'None',
                            'compare_solints' : False},
                     'ap1': {'robust': 0.5 if instrument == 'eM' else 1.0,
                            'solint': 'inf' if instrument == 'eM' else ('240s' if receiver in ('X', 'K', 'Ka', 'Ku') else '120s'),
                            # 'solint' : '120s',
                             'sigma_mask': 6.0,
                             'mask_grow_iterations': 8,
                             'combine': 'spw' if general_settings['force_combine_spw'] else '',
                            #  'combine': '',# if joint array configs/instruments
                             'gaintype': 'T',
                             'calmode': 'ap',
                            #  'minsnr': 0.75 if instrument == 'eM' else 1.0,
                             'minsnr': 0.1 if instrument == 'eM' else 0.1, #testing
                             'spwmap': [], #leavy empty here. It will be filled later if combine='spw'
                             'nsigma_automask' : '3.0',
                             'nsigma_autothreshold' : '1.0',
                             'uvtaper' : [''],
                             'with_multiscale' : True,
                             'scales': 'None',
                             'compare_solints' : False},
                     }


"""
Selfcal parameters to be used for faint sources, 
with a total integrated flux density between 10 and 20 mJy.
"""
params_faint = {'name': 'faint',
                'p0': {'robust': 0.0,
                       'solint' : '192s' if receiver in ('X', 'K', 'Ka', 'Ku') or instrument == 'eM' else '120s',
                       'sigma_mask': 12.0 if instrument == 'eM' else 20,
                       'mask_grow_iterations': 2,
                       'combine': 'spw' if general_settings['force_combine_spw'] else '',
                       'gaintype': 'T',
                       'calmode': 'p',
                       'minsnr': 0.1 if instrument == 'eM' else 0.1,
                       'spwmap': [],
                       'nsigma_automask' : '5.0',
                       'nsigma_autothreshold' : '2.5',
                       'uvtaper' : [''],
                       'with_multiscale' : False,
                       'scales' : 'None',
                       'compare_solints' : False},
                'p1': {'robust': 0.0,
                       'solint' : '192s' if receiver in ('K', 'Ka', 'Ku') or instrument == 'eM' else '96s',
                       'sigma_mask': 20,
                       'mask_grow_iterations': 4,
                       'combine':  '' if general_settings['allow_combine_spw'] == False else ('spw' if receiver in ('K', 'Ka', 'Ku') or instrument == 'eM' else ''),
                     #   'combine': '', #if VLA-C-config
                       'gaintype': 'T',
                       'calmode': 'p',
                       'minsnr': 0.1 if instrument == 'eM' else 0.1,
                       'spwmap': [],
                       'nsigma_automask' : '4.0',
                       'nsigma_autothreshold' : '1.5',
                       'uvtaper' : [''],
                     #   'uvtaper': [taper_size] if receiver in ('Ku', 'K', 'Ka') or
                     #                              instrument == 'eM' else [''],
                       'with_multiscale': False if receiver in ('K', 'Ka', 'Ku') or instrument ==
                                                   'eM' else True,
                       # 'scales' : '0,5,20',
                       'scales': 'None',
                       'compare_solints' : False},
                'p2': {'robust': 0.5 if instrument == 'eM' else 1.0,
                     #   'solint' : '240s' if receiver in ('X', 'K', 'Ka', 'Ku') or instrument == 'eM' else '60s',
                       'solint' : '96s' if instrument == 'eM' else '96s',
                       'sigma_mask': 15,
                       'mask_grow_iterations': 4,
                       'combine':  '' if general_settings['allow_combine_spw'] == False else ('spw' if receiver in ('K', 'Ka', 'Ku') or instrument == 'eM' else ''),
                     #   'combine': '', #if VLA-C-config
                       'gaintype': 'T',
                       'calmode': 'p',
                       'minsnr': 0.1 if instrument == 'eM' else 0.1,
                       'spwmap': [],
                       'nsigma_automask': '4.0',
                       'nsigma_autothreshold': '1.5',
                       # 'uvtaper' : [''], #if VLA-C-config
                       'uvtaper' : [''] if general_settings['allow_tapper'] == False else ([taper_size] if receiver in ('C','X','K', 'Ka', 'Ku') or instrument == 'eM' else ['']),
                       'with_multiscale': True,
                       'scales': 'None',
                       'compare_solints': False},
                'ap1': {'robust': 0.5 if instrument == 'eM' else 1.0,
                     #    'solint': 'inf' if instrument == 'eM' else ('240s' if receiver in ('X', 'K', 'Ka', 'Ku') else '192s'), # testing
                        'solint': '192s' if instrument == 'eM' or receiver in ('K', 'Ka', 'Ku') else '96s',
                        'sigma_mask': 7,
                        'mask_grow_iterations': 4,
                        'combine':  '' if general_settings['allow_combine_spw'] == False else ('spw' if receiver in ('K', 'Ka', 'Ku') or instrument == 'eM' else ''),
                     #    'combine': 'spw' if receiver in ('C', 'X', 'K', 'Ka', 'Ku') or instrument == 'eM' else '',#testing
                     #    'combine': '', #if VLA-C-config
                        'gaintype': 'T',
                        'calmode': 'ap',
                        'minsnr': 0.1 if instrument == 'eM' else 0.1,
                        'spwmap': [],
                        'nsigma_automask' : '4.0',
                        'nsigma_autothreshold' : '1.5',
                        # 'uvtaper' : [''], #if VLA-C-config
                        'uvtaper' : [''] if general_settings['allow_tapper'] == False else ([taper_size] if receiver in ('C','X','K', 'Ka', 'Ku') or instrument == 'eM' else ['']),
                        'with_multiscale': True,
                        # 'scales': '0,5,20,40',
                        'scales': 'None',
                        'compare_solints' : False},
                }



"""
Selfcal parameters to be used for standard sources, 
with a total integrated flux density between 20 and 50 mJy.
"""
params_standard_1 = {'name': 'standard_1',
                   'p0': {
                          'robust': -0.5 if receiver in ('K', 'Ka') or instrument == 'eM' else -0.5,
                          'solint' : '240s' if instrument == 'eM' else '96s',
                          'sigma_mask': 20.0 if receiver in ('Ku', 'K', 'Ka', 'Q') or instrument == 'eM' else 40.0,
                          'mask_grow_iterations': 2,
                     #      'combine': 'spw' if receiver in ('X', 'K', 'Ka', 'Ku', 'Q') or instrument == 'eM' else '',
                          'combine': 'spw', #testing
                          'gaintype': 'T',
                          'calmode': 'p',
                     #      'minsnr': 1.0 if receiver in ('X', 'Ku', 'K', 'Ka') or instrument == 'eM' else 1.0,
                          'minsnr': 0.1, #testing
                          'spwmap': [],
                          'nsigma_automask' : '5.0',
                          'nsigma_autothreshold' : '2.5',
                          'uvtaper' : [''],
                          'with_multiscale' : False,
                          # 'scales' : '0,5,20',
                          'scales': 'None',
                          'compare_solints' : False},
                   'p1': {'robust': 0.0 if receiver in ('K', 'Ka') or instrument == 'eM' else 0.0,
                          'solint' : '120s' if instrument == 'eM' else '96s',
                     #      'sigma_mask': 30.0 if receiver in ('C', 'X', 'K', 'Ka', 'Ku') or instrument == 'eM' else 60.0,#testing
                          'sigma_mask': 20.0 if receiver in ('K', 'Ka', 'Ku') or instrument == 'eM' else 30.0,
                          'mask_grow_iterations': 4,
                          'combine':  '' if general_settings['allow_combine_spw'] == False else ('spw' if receiver in ('K', 'Ka', 'Ku') or instrument == 'eM' else ''),
                     #      'combine': 'spw' if instrument == 'eM' else '', #testing
                     #      'combine': 'spw', #testing
                          'gaintype': 'T' if receiver in ('Ku','K', 'Ka') or instrument == 'eM' else 'G',
                          'calmode': 'p',
                     #      'minsnr': 1.0,
                          'minsnr': 0.1 if instrument == 'eM' else 0.1, #testing
                          'spwmap': [],
                          'nsigma_automask' : '3.0',
                          'nsigma_autothreshold' : '1.5',
                          'uvtaper' : [''],
                          'with_multiscale' : True,
                          # 'scales' : '0,5,20',
                          'scales': 'None',
                          'compare_solints' : False},
                   'p2': {'robust': 0.5 if instrument == 'eM' else 1.0,
                          'solint': '48s',
                          'sigma_mask': 15.0 if receiver in ('K', 'Ka', 'Ku') or instrument == 'eM' else 15.0,
                          'mask_grow_iterations': 4,
                          'combine':  '' if general_settings['allow_combine_spw'] == False else ('spw' if receiver in ('K', 'Ka', 'Ku') or instrument == 'eM' else ''),
                     #      'combine': 'spw', #testing
                          'gaintype': 'T',
                          'calmode': 'p',
                     #      'minsnr': 1.0,
                          'minsnr': 0.1 if instrument == 'eM' else 0.1, #testing
                          'spwmap': [],
                          'nsigma_automask' : '3.0',
                          'nsigma_autothreshold' : '1.0',
                          # 'uvtaper' : [''],
                          'uvtaper' : [''] if general_settings['allow_tapper'] == False else ([taper_size] if receiver in ('C','X','K', 'Ka', 'Ku') or instrument == 'eM' else ['']),
                          #testing
                     #      'uvtaper': [taper_size] if receiver in ('C', 'Ku', 'K', 'Ka') or instrument == 'eM' else [''],#testing
                     #      'uvtaper' : [taper_size] if receiver in ('X', 'Ku', 'K', 'Ka') or instrument == 'eM' else [''],
                          'with_multiscale': True,
                          # 'scales': '0,5,10,20,40',
                          'scales': 'None',
                          'compare_solints' : False},
                   'ap1': {'robust': 0.5 if instrument == 'eM' else 1.0,
                           'solint': '240s' if receiver in ('K', 'Ka', 'Ku') or instrument == 'eM' else '96s',
                           'sigma_mask': 10.0 if receiver in ('K', 'Ka', 'Ku') or instrument == 'eM' else 10.0,
                           'mask_grow_iterations': 8,
                           'combine':  '' if general_settings['allow_combine_spw'] == False else ('spw' if receiver in ('K', 'Ka', 'Ku') or instrument == 'eM' else ''),
                     #       'combine': 'spw' if instrument == 'eM' else '', #testing
                     #       'combine': 'spw', #testing
                           'gaintype': 'T' if instrument == 'eM' else 'G',
                           'calmode': 'ap',
                     #       'minsnr': 1.0,
                           'minsnr': 0.1 if instrument == 'eM' else 0.1, #testing
                           'spwmap': [],
                           'nsigma_automask' : '3.0',
                           'nsigma_autothreshold' : '1.0',
                           # 'uvtaper' : [''],
                           'uvtaper' : [''] if general_settings['allow_tapper'] == False else ([taper_size] if receiver in ('C','X','K', 'Ka', 'Ku') or instrument == 'eM' else ['']),
                     #       'uvtaper': [taper_size] if receiver in ('C', 'Ku', 'K', 'Ka') or instrument == 'eM' else [''],#testing
                     #       'uvtaper' : [taper_size] if receiver in ('Ku', 'K', 'Ka') or instrument == 'eM' else [''],
                           'with_multiscale': True,
                     #       'with_multiscale': False if receiver in ('K', 'Ka', 'Ku') or
                     #                                   instrument == 'eM' else True,
                           # 'scales': '0,5,10,20,40',
                           'scales': 'None',
                           'compare_solints' : False},
                 }


"""
Selfcal parameters to be used for standard sources, 
with a total integrated flux density between 50 and 100 mJy.
"""
"""
Selfcal parameters to be used for standard sources, 
with a total integrated flux density between 50 and 100 mJy.
Note that some values may change if using e-MERLIN or VLA.
"""
params_standard_2 = {'name': 'standard_2',
                   'p0': {'robust': -0.5 if receiver in ('K', 'Ku', 'Ka') or instrument == 'eM' else -1.0,
                          'solint' : '96s' if instrument == 'eM' else '96s',
                          'sigma_mask': 30.0 if instrument == 'eM' else 50.0,
                          'mask_grow_iterations': 2,
                          'combine': 'spw' if instrument == 'eM' else '',
                          'gaintype': 'T',
                          'calmode': 'p',
                          'minsnr': 1.0 if receiver in ('K', 'Ku', 'Ka') or instrument == 'eM' else 2.0,
                          'spwmap': [],
                          'nsigma_automask' : '6.0',
                          'nsigma_autothreshold' : '3.0',
                          'uvtaper' : [''],
                          'with_multiscale' : False,
                          # 'scales': '0,5,10',
                          'scales': 'None',
                          'compare_solints' : False},
                   'p1': {'robust': 0.0 if receiver in ('K', 'Ku', 'Ka') or instrument == 'eM' else -0.5,
                          'solint' : '60s' if instrument == 'eM' else '60s',
                          'sigma_mask': 20.0 if receiver in ('C', 'X', 'K', 'Ka', 'Ku') or instrument == 'eM' else 20.0,
                          'mask_grow_iterations': 4,
                          'combine':  '' if general_settings['allow_combine_spw'] == False else ('spw' if receiver in ('K', 'Ka', 'Ku') or instrument == 'eM' else ''),
                     #      'combine': 'spw' if instrument == 'eM' else '', #needs to be tested
                          'gaintype': 'T' if instrument == 'eM' else 'G',
                          'calmode': 'p',
                          'minsnr': 1.0 if instrument == 'eM' else 2.0,
                          'spwmap': [],
                          'nsigma_automask' : '3.0',
                          'nsigma_autothreshold' : '1.5',
                          'uvtaper' : [''],
                     #      'uvtaper' : [taper_size] if instrument == 'eM' else [''],#testing
                          'with_multiscale' : True if instrument == 'eM' else True,
                          # 'scales': '0,5,10,20',
                          'scales': 'None',
                          'compare_solints' : False},
                   'p2': {'robust': 0.5 if instrument == 'eM' else 0.5,
                          'solint': '60s' if instrument == 'eM' else '36s',
                          'sigma_mask': 15.0 if receiver in ('C', 'X', 'K', 'Ka', 'Ku') or instrument == 'eM' else 15.0,
                          'mask_grow_iterations': 4,
                          'combine':  '' if general_settings['allow_combine_spw'] == False else ('spw' if receiver in ('K', 'Ka', 'Ku') or instrument == 'eM' else ''),
                     #      'combine': 'spw' if instrument == 'eM' else '',
                          'gaintype': 'T',
                          'calmode': 'p',
                          'minsnr': 1.0 if instrument == 'eM' else 2.0,
                          'spwmap': [],
                          'nsigma_automask' : '3.0',
                          'nsigma_autothreshold' : '1.0',
                          'uvtaper' : [''],
                     #      'uvtaper': [taper_size] if receiver in ('Ku', 'K', 'Ka') or instrument == 'eM' else [''],
                     #      'uvtaper': [taper_size] if receiver in ('S', 'C', 'Ku', 'K', 'Ka') or instrument == 'eM' else [''],#testing
                          'with_multiscale' : True,
                          # 'scales': '0,5,10,20,40',
                          'scales': 'None',
                          'compare_solints' : False},
                   'ap1': {'robust': 0.5,
                           'solint': '120s' if instrument == 'eM' else '60s',
                           'sigma_mask': 10.0 if receiver in ('K', 'Ka', 'Ku') or instrument == 'eM' else 10.0,
                           'mask_grow_iterations': 8,
                           'combine':  '' if general_settings['allow_combine_spw'] == False else ('spw' if receiver in ('K', 'Ka', 'Ku') or instrument == 'eM' else ''),
                     #       'combine': 'spw' if instrument == 'eM' else '',
                           'gaintype': 'G',
                           'calmode': 'ap',
                           'minsnr': 1.0 if instrument == 'eM' else 2.0,
                           'spwmap': [],
                           'nsigma_automask' : '3.0',
                           'nsigma_autothreshold' : '1.0',
                           'uvtaper' : [''],
                     #       'uvtaper': [taper_size] if receiver in ('S', 'C', 'Ku', 'K', 'Ka') or instrument == 'eM' else [''],#testing
                           'with_multiscale' : True if instrument == 'eM' else True,
                           # 'scales': '0,5,10,20,40',
                           'scales': 'None',
                           'compare_solints' : False},
                 }

"""
Selfcal parameters to be used for bright sources, 
with a total integrated flux density above 0.1 Jy.
"""
params_bright = {'name': 'bright',
                 'p0': {'robust': -0.5 if receiver in ('K', 'Ka') else -1.0,
                        'solint': '96s' if instrument == 'eM' else '48s',
                        'sigma_mask': 60,
                        'mask_grow_iterations': 2,
                        'combine': 'spw' if instrument == 'eM' else '',
                        'gaintype': 'G',
                        'calmode': 'p',
                        'minsnr': 1.0 if instrument == 'eM' else 2.0,
                        'spwmap': [],
                        'nsigma_automask': '6.0',
                        'nsigma_autothreshold': '3.0',
                        'uvtaper' : [''],
                        'with_multiscale' : False,
                        # 'scales': '0,5,10',
                        'scales': 'None',
                        'compare_solints' : False},
                 'p1': {'robust': 0.0 if receiver in ('K', 'Ka') or instrument == 'eM' else -0.5,
                        'solint' : '96s' if instrument == 'eM' else '24s',
                        'sigma_mask': 30.0 if receiver in ('X', 'K', 'Ka', 'Ku') or instrument == 'eM' else 60.0,
                        'mask_grow_iterations': 4,
                     #    'combine': 'spw' if instrument == 'eM' else '', #testing
                        'combine': '',
                        'gaintype': 'T' if instrument == 'eM' else 'G',
                        'calmode': 'p',
                        'minsnr': 1.0 if instrument == 'eM' else 2.0,
                        'spwmap': [],
                        'nsigma_automask': '4.0' if instrument == 'eM' else '6.0',
                        'nsigma_autothreshold': '1.5' if instrument == 'eM' else '3.0',
                        'uvtaper' : [''],
                     #    'uvtaper' : [taper_size] if instrument == 'eM' else [''], # testing
                        'with_multiscale': True,
                        # 'scales': '0,5,10,20',
                        'scales': 'None',
                        'compare_solints': False},
                 'p2': {'robust': 0.0,
                        'solint': '60s' if instrument == 'eM' else '24s',
                        'sigma_mask': 12.0 if receiver in ('X', 'K', 'Ka', 'Ku') or instrument == 'eM' else 30.0,
                        'mask_grow_iterations': 4,
                        'combine': '',
                        'gaintype': 'G',
                        'calmode': 'p',
                        'minsnr': 1.0 if instrument == 'eM' else 2.0,
                        'spwmap': [],
                        'uvtaper' : [''],
                        'nsigma_automask': '4.0',
                        'nsigma_autothreshold': '1.0',
                        'with_multiscale': True,
                        # 'scales': '0,5,10,20,40',
                        'scales': 'None',
                        'compare_solints': False},
                 'ap1': {'robust': 0.5,
                         'solint': '96s' if instrument == 'eM' else '36s',
                         'sigma_mask': 8.0 if receiver in ('X', 'K', 'Ka', 'Ku') or instrument == 'eM' else 20.0,
                         'mask_grow_iterations': 8,
                         'combine': '',
                     #     'combine': 'spw' if instrument == 'eM' else '', #testing
                         'gaintype': 'T' if instrument == 'eM' else 'G',
                         'calmode': 'ap',
                         'minsnr': 1.0 if instrument == 'eM' else 2.0,
                         'spwmap': [],
                         'uvtaper' : [''],
                         'nsigma_automask': '3.0',
                         'nsigma_autothreshold': '1.0',
                         'with_multiscale': True if instrument == 'eM' else True,
                         # 'scales': '0,5,10,20,40',
                         'scales': 'None',
                         'compare_solints': False},
                 }


params_trial_2 = None