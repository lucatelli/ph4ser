"""
Configuration of different template parameters for self-calibration and imaging.
This is intended to be used as a first trial of self-calibration.
"""
# visibility_info = {'path':"/mnt/scratch/lucatelli/astronomical-data/M82_v2/eM_C/sc_v12/CY2204/standard_manchester_test_2/",
#                    'vis_name':'M82_CY2204_eM_C.avg4s',
#                    'field':'M82',
#                    'savename':'_eM_C_CY2204_sc_v12'
#        }

visibility_info = {'path':"/media/sagauga/void/astronomical-data/LIRGI_Sample_v2/VV705/VLA_A_C/sc_v13/standard/",
                   'vis_name':'VV705_SDSSJ1518.calibrated',
                   'field':'VV705',
                   'savename':'_eM_C_CY2204_sc_v13'
       }

# path = ""
# vis_name = ""
# field = ""        # does not need to match field name in the visibility (this is just a name to be attatched to some output files).
# savename = ""     # final name will be field + savename + '.ms'

# visibility_info = {'path': path,
#                    'field': field,
#                    'vis_name': vis_name,
#                    'savename': savename
#        }

FIELD = ''
SPWS = ''
ANTENNAS = ''
refantmode = 'flex'          # 'strict' or 'flex' -- flex may be better in most situations. But use with caution.
refant = ''                  # If '', will compute automatically -- this is VERY recommended, unless you know what you are doing.
minblperant = 2              # Minimum number of baselines per antenna. For safety, set to 3.
solnorm = False              # True or False (globally used). If '', will use False for phases and True for
                             # amplitudes. As of now, solnorm=False is safer and better (why? coming soon....)
default_combine = 'scan'

quiet = True                 # WSClean logging (print or not stuff during deconvolution)
plotting_verbosity = 1       # for CASA plots (controls how many kinds of plots are generated.)
show_figures = True          # Leave True for interactive mode (e.g. jupyter notebook); set to False for terminal mode.
do_additional_images = False # Not working yet, leave to False for now.
run_mode = 'terminal'        # Not tested, leave to 'terminal' for now. Use the jupyter notebook for interactive mode (within ../examples).

new_phasecentre = None
multi_config = False # True if using multiple configurations or arrays.
receiver = 'C'
instrument = 'EVLA' # 'EVLA' or 'eM'


if instrument == 'eM':
       refantmode = 'flex'

#number of channels-out for wsclean's MFS imager.
if instrument == 'EVLA':
    if receiver in ('L', 'S', 'C'):
        nc = 6 #number of bandwidth split during convolution (number of sub-band WSClean images)
    else:
        if multi_config:
            nc = 3 #number of bandwidth split during convolution (number of sub-band WSClean images)
        else:
            nc = 6 #number of bandwidth split during convolution (number of sub-band WSClean images)
if instrument == 'eM':
       if multi_config:
              nc = 3 #number of bandwidth split during convolution (number of sub-band WSClean images)
       else:
              if receiver == 'L':
                     nc = 4
              else:
                     nc = 3 #number of bandwidth split during convolution (number of sub-band WSClean images)
nc = 6 #overwrite for testing

negative_arg='no-negative'  #dont allow negative components during WSClean cleaning.
# negative_arg='negative'     #allow negative components during WSClean cleaning.
steps = [
    'startup',          # create directory structure, visibility preparation (check, avg, ...); start variables and clear visibilities.
#     'fov_image',      # create a FOV image
    'save_init_flags',  # save (or restore) the initial flags
    # 'statwt',           # run statwt on the initial data
    'autoflag_init',    # run rflag on the initial data -- use for L and S band data.
    'test_image',       # create a test image
    'select_refant',    # select reference antenna
    # 'delay_K',
    'p0',               # initial test  of selfcal, phase only (p)
    'p1',               # redo phase-only selfcal (if enough flux density); ignores p0
    'p2',               # continue phase-only selfcal (incremental)
    'ap1',              # amp-selfcal (ap); uses p0 or (p1 and p2)
    'split_trial_1',    # split the data after first trial (and run wsclean)
    # 'autoflag_final',   # run rflag on the final data (use for L and S band data)
    'report_results'    # report results of first trial
]


cell_sizes_JVLA = {'L':'0.17arcsec',
                   'S':'0.09arcsec',
                   'C':'0.05arcsec',
                   'X':'0.035arcsec',
                   'Ku':'0.02arcsec',
                   'K':'0.012arcsec',
                   'Ka':'0.009arcsec',
                   'Q':'0.007arcsec'
                   }

cell_sizes_eMERLIN = {'L':'0.03arcsec',
                      'C':'0.008arcsec'}
#
taper_sizes_eMERLIN = {'L':'0.3arcsec',
                       'C':'0.05arcsec'}

taper_sizes_JVLA = {'L':'2.0arcsec',
                    'S':'1.0arcsec',
                    'C':'0.4arcsec',
                    'X':'0.5arcsec',
                    'Ku':'0.2arcsec', #Ku-A
                    # 'Ku':'1.0arcsec', #Ku-C
                    'K':'0.1arcsec',
                    # 'Ka':'0.8arcsec', #Ka-C
                    'Ka':'0.08arcsec', #Ka-A
                    'Q':'0.1arcsec' #Q-A
                    }


if instrument == 'eM':
    cell_size = cell_sizes_eMERLIN[receiver]
    taper_size = taper_sizes_eMERLIN[receiver]
if instrument == 'EVLA':
    cell_size = cell_sizes_JVLA[receiver]
    taper_size = taper_sizes_JVLA[receiver]

init_parameters = {'fov_image': {'imsize': 1024*8,
                                'cell': '0.4arcsec',
                                'basename': 'FOV_phasecal_image',
                                'FIELD_SHIFT': None,
                            #     'FIELD_SHIFT':"'14:57:43.145 +24.35.10.257'", #VV340a >1Jy outlier
                            #     'FIELD_SHIFT':"'13:37:29.220001 +48.18.20.59999'", #NGC5256
                                'niter': 1000,
                                'robust': -0.25 if multi_config else (0.5 if receiver in ('Ku', 'K', 'Ka', 'Q') or instrument == 'eM' else 0.0),
                                },
                  'test_image': {'imsize': int(1024*2),
                                 'imsizey': int(1024*2),
                                 'FIELD_SHIFT': None,
                            #      'FIELD_SHIFT':"'13:37:29.220001 +48.18.20.59999'", #NGC5256
                                 'cell': cell_size,
                                 'prefix': 'test_image',
                                 'uvtaper': [''],
                                 'niter': 10000,
                                 'robust': 0.0 if multi_config else (0.5 if receiver in ('Ku', 'K', 'Ka', 'Q') or instrument == 'eM' else 0.0)
                                #  'robust': 0.0,
                                 }
                  }

global_parameters = {'imsize': init_parameters['test_image']['imsize'],
                     'imsizey': init_parameters['test_image']['imsizey'],
                     'FIELD_SHIFT': init_parameters['test_image']['FIELD_SHIFT'],
                     'cell': init_parameters['test_image']['cell'],
                     'nsigma_automask' : '5.0' if multi_config else '4.0',
                     'nsigma_autothreshold' : '2.5' if multi_config else '2.0',
                     'mask_grow_iterations': 2,
                     'uvtaper' : [''],
                     'with_multiscale' : True,
                     'use_mask': True,
                     'custom_mask' : None,
                     'scales' : 'None',
                     'maxmscales' : '6',
                     'niter':100000}

general_settings = {
   'timebin_statw': '12s',#timebin for statwt
   'statwt_statalg': 'chauvenet',
   'calwt': False,#calibrate/update the weights during applycal (phases)?
   'calwt_ap': False,#calibrate/update the weights during applycal (phases+amplitudes)?
   'applymode_p': 'calflag',
   'applymode_ap': 'calflag',
   'timebin': None,          # time averaging bin (e.g. '6s', '10s'); used if do_average_time=True
   'channel_width': None,    # frequency averaging width in channels (e.g. [2], [4]); used if do_average_freq=True
   'do_average_time': False, # average in time during _prepare_visibility(); set timebin above
   'do_average_freq': False, # average in frequency during _prepare_visibility(); set channel_width above
   'correlations':'RR,LL',
   'new_phasecentre': new_phasecentre,
   'force_combine_spw' : False, #force the combination of spectral windows.
   'allow_combine_spw' : False, #allow the combination of spectral windows.
   'allow_tapper' : False, #allow the use of tapering.
   'keep_p0' : False, #keep p0 as the root of the chain (p0 > p1 > p2 > ap1); if False, p0 is discarded once p1 runs.
   }



"""
Selfcal parameters to be used for very faint sources, 
with a total integrated flux density lower than 10 mJy.
"""
# params_very_faint = {'name': 'very_faint',
#                      'p0': {
#                             # 'robust': 0.0 if multi_config else 0.5,
#                             'robust': 0.5 if multi_config else (0.5 if general_settings['allow_tapper'] else 0.75), #testing
#                             # 'robust': 0.0,
#                             # 'solint': '120s' if general_settings['force_combine_spw'] else ('240s' if receiver in ('K', 'Ka') or instrument == 'eM' else '120s'),
#                             'solint': '120s' if general_settings['force_combine_spw'] else ('120s' if general_settings['allow_combine_spw'] else ('240s' if instrument == 'eM' else '120s')),
#                             # 'solint' : '120s',
#                             # 'sigma_mask': 6.0 if receiver in ('K', 'Ka', 'Ku') or instrument == 'eM' else 15.0,#C-Config
#                          #    'sigma_mask': 6.0 if multi_config else (8.0 if receiver in ('K', 'Ka', 'Ku') or instrument == 'eM' else 12.0),#A-Config
#                             'sigma_mask': 12.0 if multi_config else (12.0 if receiver in ('K', 'Ka', 'Ku') or instrument == 'eM' else 15.0), #testing
#                             # 'sigma_mask': 12.0, #test
#                             'mask_grow_iterations': 3 if multi_config else 3,
#                             'combine': 'scan,spw' if general_settings['force_combine_spw'] else ('scan,spw' if general_settings['allow_combine_spw'] else 'scan'),
#                             # 'combine': 'spw' if general_settings['force_combine_spw'] else ('spw' if general_settings['allow_combine_spw'] else ''),
#                             # 'combine': 'scan', # testing july 2025
#                             # 'combine': '',# if joint array configs/instruments
#                             'gaintype': 'T',
#                             'calmode': 'p',
#                             # 'minsnr': 0.75 if instrument == 'eM' else 1.0,
#                             'minsnr': 0.1 if instrument == 'eM' else 0.1, #testing
#                             'spwmap': [], #leavy empty here. It will be filled later if combine='spw'
#                             'nsigma_automask' : '4.0' if multi_config or instrument == 'eM' else '4.0',
#                             'nsigma_autothreshold' : '2.0' if multi_config or instrument == 'eM' else '2.0',
#                             # 'uvtaper' : [''], #if VLA-C-config
#                             # 'uvtaper': [taper_size] if receiver in ('X', 'Ku', 'K', 'Ka') or
#                             #                             instrument == 'eM' else [''],#testing
#                             # 'uvtaper': [taper_size] if receiver in ('Ku', 'K', 'Ka') or
#                             #                             instrument == 'eM' else [''],
#                             'uvtaper' : [taper_size] if general_settings['allow_tapper'] else [''],
#                             'with_multiscale' : True,
#                             # 'with_multiscale' : False if multi_config else True,
#                             # 'with_multiscale' : False,
#                             'scales': 'None',
#                             'maxmscales': '3',
#                             # 'maxmscales': '8', #testing eM + VLA
#                             'compare_solints' : False},
#                      'ap1': {
#                             # 'robust': 0.5 if multi_config else (0.5 if instrument == 'eM' else 1.0),
#                              'robust': 0.25 if multi_config else (0.75 if instrument == 'eM' else 1.0), #testing
#                             # 'robust': 0.5,
#                             # 'solint': '240s' if instrument == 'eM' else ('240s' if receiver in ('Ku', 'K', 'Ka') else '120s'),
#                             #  'solint': '120s' if general_settings['force_combine_spw'] else ('120s' if general_settings['allow_combine_spw'] else '240s'),
#                              'solint': '120s' if general_settings['force_combine_spw'] else ('120s' if general_settings['allow_combine_spw'] else ('240s' if instrument == 'eM' else '120s')),
#                             # 'solint' : '120s',
#                              'sigma_mask': 10.0 if multi_config else (10.0 if receiver in ('K', 'Ka', 'Ku') or instrument == 'eM' else 12.0),
#                             #  'sigma_mask': 12.0, #testing
#                              'mask_grow_iterations': 4 if multi_config else 4,
#                              'combine': 'scan,spw' if general_settings['force_combine_spw'] else ('scan,spw' if general_settings['allow_combine_spw'] else 'scan'),
#                             #  'combine': 'spw' if general_settings['force_combine_spw'] else ('spw' if general_settings['allow_combine_spw'] else ''),
#                             #  'combine': 'scan', # testing july 2025
#                             #  'combine': '',# if joint array configs/instruments
#                              'gaintype': 'T',
#                              'calmode': 'ap',
#                             #  'minsnr': 0.75 if instrument == 'eM' else 1.0,
#                              'minsnr': 0.1 if instrument == 'eM' else 0.1, #testing
#                              'spwmap': [], #leavy empty here. It will be filled later if combine='spw'
#                              'nsigma_automask' : '4.0' if multi_config or instrument == 'eM' else '4.0',
#                              'nsigma_autothreshold' : '2.0' if multi_config or instrument == 'eM' else '2.0',
#                             #  'uvtaper' : [''], #if VLA-C-config
#                             #  'uvtaper': [taper_size] if receiver in ('X', 'Ku', 'K', 'Ka') or
#                             #                             instrument == 'eM' else [''],#testing
#                             #  'uvtaper': [taper_size] if receiver in ('Ku', 'K', 'Ka') or
#                             #                             instrument == 'eM' else [''],
#                             #  'with_multiscale' : False if receiver in ('K', 'Ka', 'Ku') or
#                             #                               instrument == 'eM' else True,
#                              'uvtaper' : [taper_size] if general_settings['allow_tapper'] else [''],
#                              'with_multiscale' : True,
#                             #  'with_multiscale' : False if multi_config else True,
#                             #  'with_multiscale' : False,
#                              'scales': 'None',
#                              'maxmscales': '4',
#                             #  'maxmscales': '8', #testing eM + VLA
#                              'compare_solints' : False},
#                      }

params_very_faint = {'name': 'very_faint', #global very_faint template for e-MERLIN
                     'p0': {
                            'robust': 0.5, #testing
                            # 'solint': '90s' if general_settings['allow_combine_spw'] else ('240s' if instrument == 'eM' else '120s'),
                            'solint': '120s',
                            'sigma_mask': 12.0 if instrument == 'eM' else 12.0, #testing
                            'mask_grow_iterations': 3,
                            'combine': 'scan,spw' if general_settings['force_combine_spw'] else ('scan,spw' if general_settings['allow_combine_spw'] else 'scan'),
                            'gaintype': 'G',
                            'calmode': 'p',
                            'minsnr': 0.1 if instrument == 'eM' else 0.1, #testing
                            'spwmap': [], #leavy empty here. It will be filled later if combine='spw'
                            'nsigma_automask' : '4.0',
                            'nsigma_autothreshold' : '2.0',
                            'uvtaper' : [''],
                            # 'uvtaper' : [taper_size] if general_settings['allow_tapper'] else [''],
                            'with_multiscale' : True,
                            'scales': 'None',
                            'maxmscales': '3',
                            # 'maxmscales': '8', #testing eM + VLA
                            'compare_solints' : False},
                     'ap1': {
                            #  'robust': 1.0 if general_settings['allow_tapper'] else 1.25,
                             'robust': 0.5,
                            #  'solint': '90s' if general_settings['allow_combine_spw'] else ('240s' if instrument == 'eM' else '120s'),
                             'solint': '120s',
                             'sigma_mask': 12.0 if instrument == 'eM' else 10.0, #testing
                             'mask_grow_iterations': 3,
                             'combine': 'scan,spw' if general_settings['force_combine_spw'] else ('scan,spw' if general_settings['allow_combine_spw'] else 'scan'),
                             'gaintype': 'T',
                             'calmode': 'ap',
                             'minsnr': 0.1 if instrument == 'eM' else 0.1, #testing
                             'spwmap': [], #leavy empty here. It will be filled later if combine='spw'
                             'nsigma_automask' : '4.0',
                             'nsigma_autothreshold' : '2.0',
                            #  'uvtaper' : [''],
                             'uvtaper' : [taper_size] if general_settings['allow_tapper'] else [''],
                             'with_multiscale' : True,
                             'scales': 'None',
                             'maxmscales': '6',
                            #  'maxmscales': '8', #testing eM + VLA
                             'compare_solints' : False},
                     }

# params_very_faint = {'name': 'very_faint', #global very_faint template for VLA-Ka band
#                      'p0': {
#                             'robust': 1.0, #testing
#                             # 'solint': '60s' if general_settings['allow_combine_spw'] else '90s',
#                             'solint': '90s' if general_settings['allow_combine_spw'] else '120s',
#                             'sigma_mask': 12.0 if instrument == 'eM' else 10.0, #testing
#                             'mask_grow_iterations': 3,
#                             'combine': 'scan,spw' if general_settings['force_combine_spw'] else ('scan,spw' if general_settings['allow_combine_spw'] else 'scan'),
#                             'gaintype': 'G',
#                             'calmode': 'p',
#                             'minsnr': 0.1, #testing
#                             'spwmap': [], #leavy empty here. It will be filled later if combine='spw'
#                             'nsigma_automask' : '4.0',
#                             'nsigma_autothreshold' : '2.0',
#                             'uvtaper' : [''],
#                             # 'uvtaper' : [taper_size] if general_settings['allow_tapper'] else [''],
#                             'with_multiscale' : True,
#                             'scales': 'None',
#                             'maxmscales': '4',
#                             # 'maxmscales': '8', #testing eM + VLA
#                             'compare_solints' : False},
#                      'ap1': {
#                              'robust': 1.5,
#                             #  'solint': '60s' if general_settings['allow_combine_spw'] else '90s',
#                              'solint': '90s' if general_settings['allow_combine_spw'] else '120s',
#                              'sigma_mask': 12.0 if instrument == 'eM' else 8.0, #testing
#                              'mask_grow_iterations': 4,
#                              'combine': 'scan,spw' if general_settings['force_combine_spw'] else ('scan,spw' if general_settings['allow_combine_spw'] else 'scan'),
#                              'gaintype': 'T',
#                              'calmode': 'ap',
#                              'minsnr': 0.1, #testing
#                              'spwmap': [], #leavy empty here. It will be filled later if combine='spw'
#                              'nsigma_automask' : '4.0',
#                              'nsigma_autothreshold' : '2.0',
#                             #  'uvtaper' : [''],
#                              'uvtaper' : [taper_size] if general_settings['allow_tapper'] else [''],
#                              'with_multiscale' : True,
#                              'scales': 'None',
#                              'maxmscales': '6',
#                             #  'maxmscales': '8', #testing eM + VLA
#                              'compare_solints' : False},
#                      }


# params_very_faint = {'name': 'very_faint', #global very_faint template for inf - scan based solints
#                      'p0': {
#                             'robust': 0.75, 
#                             'solint': 'inf',
#                             'sigma_mask': 12.0 if instrument == 'eM' else 12.0,
#                             'mask_grow_iterations': 3,
#                             'combine': 'spw' if general_settings['force_combine_spw'] else ('spw' if general_settings['allow_combine_spw'] else ''),
#                             'gaintype': 'G',
#                             'calmode': 'p',
#                             'minsnr': 0.1, 
#                             'spwmap': [], 
#                             'nsigma_automask' : '4.0',
#                             'nsigma_autothreshold' : '2.0',
#                             # 'uvtaper' : [''],
#                             'uvtaper' : [taper_size] if general_settings['allow_tapper'] else [''],
#                             'with_multiscale' : True,
#                             'scales': 'None',
#                             'maxmscales': '3',
#                             'compare_solints' : False},
#                      'ap1': {
#                              'robust': 0.75,
#                              'solint': 'inf',
#                              'sigma_mask': 12.0 if instrument == 'eM' else 10.0,
#                              'mask_grow_iterations': 4,
#                              'combine': 'spw' if general_settings['force_combine_spw'] else ('spw' if general_settings['allow_combine_spw'] else ''),
#                              'gaintype': 'T',
#                              'calmode': 'ap',
#                              'minsnr': 0.1,
#                              'spwmap': [], 
#                              'nsigma_automask' : '4.0',
#                              'nsigma_autothreshold' : '2.0',
#                             #  'uvtaper' : [''],
#                              'uvtaper' : [taper_size] if general_settings['allow_tapper'] else [''],
#                              'with_multiscale' : True,
#                              'scales': 'None',
#                              'maxmscales': '4',
#                              'compare_solints' : False},
#                      }


# params_very_faint = {'name': 'very_faint', #global very_faint template: test for relaively faint but diffuse sources at C band
#                      'p0': {
#                             'robust': 0.75, 
#                             'solint': '60s',
#                             'sigma_mask': 12.0 if instrument == 'eM' else 12.0,
#                             'mask_grow_iterations': 3,
#                             'combine': 'scan,spw' if general_settings['force_combine_spw'] else ('scan,spw' if general_settings['allow_combine_spw'] else 'scan'),
#                             'gaintype': 'G',
#                             'calmode': 'p',
#                             'minsnr': 0.1, 
#                             'spwmap': [], 
#                             'nsigma_automask' : '5.0',
#                             'nsigma_autothreshold' : '2.0',
#                             'uvtaper' : [''],
#                             # 'uvtaper' : [taper_size] if general_settings['allow_tapper'] else [''],
#                             'with_multiscale' : True,
#                             'scales': 'None',
#                             'maxmscales': '3',
#                             'compare_solints' : False},
#                      'ap1': {
#                              'robust': 1.0,
#                              'solint': '120s',
#                              'sigma_mask': 12.0 if instrument == 'eM' else 12.0,
#                              'mask_grow_iterations': 4,
#                              'combine': 'scan,spw' if general_settings['force_combine_spw'] else ('scan,spw' if general_settings['allow_combine_spw'] else 'scan'),
#                              'gaintype': 'T',
#                              'calmode': 'ap',
#                              'minsnr': 0.1,
#                              'spwmap': [], 
#                              'nsigma_automask' : '5.0',
#                              'nsigma_autothreshold' : '2.0',
#                             #  'uvtaper' : [''],
#                              'uvtaper' : [taper_size] if general_settings['allow_tapper'] else [''],
#                              'with_multiscale' : True,
#                              'scales': 'None',
#                              'maxmscales': '6',
#                              'compare_solints' : False},
#                      }


# params_very_faint = {'name': 'very_faint', #global very_faint template; test for Q-band
#                      'p0': {
#                             'robust': 0.5, #testing
#                             # 'solint': '60s' if general_settings['allow_combine_spw'] else '90s',
#                             'solint': '120s',
#                             'sigma_mask': 12.0 if instrument == 'eM' else 15.0, #testing
#                             'mask_grow_iterations': 3,
#                             'combine': 'scan,spw' if general_settings['force_combine_spw'] else ('scan,spw' if general_settings['allow_combine_spw'] else 'scan'),
#                             'gaintype': 'G',
#                             'calmode': 'p',
#                             'minsnr': 0.33, #testing
#                             'spwmap': [], #leavy empty here. It will be filled later if combine='spw'
#                             'nsigma_automask' : '4.0',
#                             'nsigma_autothreshold' : '2.0',
#                             # 'uvtaper' : [''],
#                             'uvtaper' : [taper_size] if general_settings['allow_tapper'] else [''],
#                             'with_multiscale' : False,
#                             'scales': 'None',
#                             'maxmscales': '3',
#                             # 'maxmscales': '8', #testing eM + VLA
#                             'compare_solints' : False},
#                      'ap1': {
#                              'robust': 0.5,
#                             #  'solint': '60s' if general_settings['allow_combine_spw'] else '90s',
#                              'solint': '120s',
#                              'sigma_mask': 12.0 if instrument == 'eM' else 15.0, #testing
#                              'mask_grow_iterations': 4,
#                              'combine': 'scan,spw' if general_settings['force_combine_spw'] else ('scan,spw' if general_settings['allow_combine_spw'] else 'scan'),
#                              'gaintype': 'T',
#                              'calmode': 'ap',
#                              'minsnr': 0.33, #testing
#                              'spwmap': [], #leavy empty here. It will be filled later if combine='spw'
#                              'nsigma_automask' : '4.0',
#                              'nsigma_autothreshold' : '2.0',
#                             #  'uvtaper' : [''],
#                              'uvtaper' : [taper_size] if general_settings['allow_tapper'] else [''],
#                              'with_multiscale' : False,
#                              'scales': 'None',
#                              'maxmscales': '4',
#                             #  'maxmscales': '8', #testing eM + VLA
#                              'compare_solints' : False},
#                      }

# params_very_faint = {'name': 'very_faint', #global very_faint template
#                      'p0': {
#                             'robust': 1.5, #testing
#                             'solint': '960s' if instrument == 'eM' else '120s',
#                             'sigma_mask': 12.0 if instrument == 'eM' else 12.0, #testing
#                             'mask_grow_iterations': 3,
#                             'combine': 'scan,spw' if general_settings['force_combine_spw'] else ('scan,spw' if general_settings['allow_combine_spw'] else 'scan'),
#                             'gaintype': 'T',
#                             'calmode': 'p',
#                             'minsnr': 0.1 if instrument == 'eM' else 0.1, #testing
#                             'spwmap': [], #leavy empty here. It will be filled later if combine='spw'
#                             'nsigma_automask' : '4.0',
#                             'nsigma_autothreshold' : '2.0',
#                             'uvtaper' : [''],
#                             # 'uvtaper' : [taper_size] if general_settings['allow_tapper'] else [''],
#                             'with_multiscale' : True,
#                             'scales': 'None',
#                             'maxmscales': '3',
#                             # 'maxmscales': '8', #testing eM + VLA
#                             'compare_solints' : False},
#                      'ap1': {
#                             'robust': 1.5,
#                              'solint': '960s' if instrument == 'eM' else '120s',
#                              'sigma_mask': 12.0 if instrument == 'eM' else 12.0, #testing
#                              'mask_grow_iterations': 6,
#                              'combine': 'scan,spw' if general_settings['force_combine_spw'] else ('scan,spw' if general_settings['allow_combine_spw'] else 'scan'),
#                              'gaintype': 'T',
#                              'calmode': 'ap',
#                              'minsnr': 0.1 if instrument == 'eM' else 0.1, #testing
#                              'spwmap': [], #leavy empty here. It will be filled later if combine='spw'
#                              'nsigma_automask' : '4.0',
#                              'nsigma_autothreshold' : '2.0',
#                             #  'uvtaper' : [''],
#                              'uvtaper' : [taper_size] if general_settings['allow_tapper'] else [''],
#                              'with_multiscale' : True,
#                              'scales': 'None',
#                              'maxmscales': '6',
#                             #  'maxmscales': '8', #testing eM + VLA
#                              'compare_solints' : False},
#                      }

"""
Selfcal parameters to be used for faint sources, 
with a total integrated flux density between 10 and 20 mJy.
"""
params_faint = {'name': 'faint',
                'p0': {
                    #    'robust': -0.25 if multi_config else 0.5,
                       'robust': 0.5 if multi_config else 0.5,
                    #    'robust': 0.5,
                    #    'solint' : '190s' if receiver in ('X', 'K', 'Ka', 'Ku', 'Q') or instrument == 'eM' else '120s',
                       'solint' : '240s' if instrument == 'eM' else ('90s' if general_settings['force_combine_spw'] else ('90s' if general_settings['allow_combine_spw'] else '120s')),
                       'sigma_mask': 10.0 if multi_config else (15.0 if receiver in ('X', 'K', 'Ka', 'Ku', 'Q') or instrument == 'eM' else 20.0),
                    #    'sigma_mask': 12.0, #test
                       'mask_grow_iterations': 4 if multi_config else 2,
                     #   'combine': 'scan,spw' if general_settings['force_combine_spw'] else 'scan',
                       'combine':  'scan' if general_settings['allow_combine_spw'] == False else ('scan,spw' if receiver in ('X', 'K', 'Ka', 'Ku', 'Q') or instrument == 'eM' else 'scan'),
                    #    'combine':  '' if general_settings['allow_combine_spw'] == False else ('spw' if receiver in ('X', 'K', 'Ka', 'Ku', 'Q') or instrument == 'eM' else ''),
                     #   'combine': 'spw' if general_settings['force_combine_spw'] else '',
                     #   'combine': 'scan', # testing july 2025
                       'gaintype': 'G',
                       'calmode': 'p',
                       'minsnr': 0.1 if instrument == 'eM' else 0.1,
                       'spwmap': [],
                       'nsigma_automask' : '5.0',
                       'nsigma_autothreshold' : '2.5',
                       'uvtaper' : [''],
                       'with_multiscale' : True,
                    #    'with_multiscale' : False if multi_config else True,
                    #    'with_multiscale' : False,
                       'scales' : 'None',
                       'maxmscales': '3',
                    #    'maxmscales': '8', #testing eM + VLA
                       'compare_solints' : False},
                'p1': {
                       'robust': 0.75 if instrument == 'eM' else 0.5,
                    #    'robust': 1.0 if instrument == 'eM' else 0.5, #testing
                    #    'solint' : '120s' if receiver in ('K', 'Ka', 'Ku', 'Q') or instrument == 'eM' else '90s',
                       'solint' : '120s' if general_settings['force_combine_spw'] else ('120s' if general_settings['allow_combine_spw'] else ('120s' if receiver in ('K', 'Ka', 'Ku', 'Q') or instrument == 'eM' else '90s')),
                       'sigma_mask': 12.0 if receiver in ('K', 'Ka', 'Ku', 'Q') or instrument == 'eM' else 20.0,
                       'mask_grow_iterations': 2 if multi_config else 4,
                    #    'combine':  'scan' if general_settings['allow_combine_spw'] == False else ('scan,spw' if receiver in ('X', 'K', 'Ka', 'Ku', 'Q') or instrument == 'eM' else 'scan'),
                       'combine':  'scan' if general_settings['allow_combine_spw'] == False else ('scan,spw' if receiver in ('K', 'Ka', 'Ku', 'Q') or instrument == 'eM' else 'scan'), #testing
                     #   'combine':  '' if general_settings['allow_combine_spw'] == False else ('spw' if receiver in ('K', 'Ka', 'Ku') else ''), #testing
                     #   'combine': 'scan', # testing july 2025
                     #   'combine': 'scan,spw' if receiver in ('C', 'X', 'K', 'Ka', 'Ku') or instrument == 'eM' else '',#testing
                     #   'combine': 'scan,', #if VLA-C-config
                       'gaintype': 'G',
                       'calmode': 'p',
                       'minsnr': 0.1 if instrument == 'eM' else 0.1,
                       'spwmap': [],
                       'nsigma_automask' : '5.0',
                       'nsigma_autothreshold' : '2.5',
                       'uvtaper' : [''],
                     #   'uvtaper': [taper_size] if receiver in ('Ku', 'K', 'Ka', 'Q') or
                     #                              instrument == 'eM' else [''],
                       'with_multiscale': False if multi_config else (True if receiver in ('K', 'Ka', 'Ku', 'Q') or instrument == 'eM' else True),
                       # 'scales' : '0,5,20',
                       'scales': 'None',
                       'maxmscales': '3',
                       'compare_solints' : False},
                'p2': {
                    #    'robust': 0.5 if instrument == 'eM' else 0.75,
                    #    'robust': 0.5 if general_settings['allow_tapper'] else (0.75 if instrument == 'eM' else 1.0),
                       'robust': 0.5 if general_settings['allow_tapper'] else 1.0, #testing
                     #   'solint' : '240s' if receiver in ('X', 'K', 'Ka', 'Ku', 'Q') or instrument == 'eM' else '60s',
                    #    'solint' : '90s' if receiver in ('K', 'Ka', 'Ku', 'Q') or instrument == 'eM' else '60s',
                       'solint' : '90s' if general_settings['force_combine_spw'] else ('90s' if general_settings['allow_combine_spw'] else ('120s' if receiver in ('K', 'Ka', 'Ku', 'Q') or instrument == 'eM' else '60s')),
                       'sigma_mask': 10.0,
                       'mask_grow_iterations': 2 if multi_config else 4,
                    #    'combine':  'scan' if general_settings['allow_combine_spw'] == False else ('scan,spw' if receiver in ('X', 'K', 'Ka', 'Ku') or instrument == 'eM' else 'scan'),
                       'combine':  'scan' if general_settings['allow_combine_spw'] == False else ('scan,spw' if receiver in ('K', 'Ka', 'Ku', 'Q') or instrument == 'eM' else 'scan'), #testing
                     #   'combine':  '' if general_settings['allow_combine_spw'] == False else ('spw' if receiver in ('K', 'Ka', 'Ku') else ''), #testing
                     #   'combine': 'scan', # testing july 2025
                     #   'combine': 'spw' if receiver in ('C', 'X', 'K', 'Ka', 'Ku', 'Q') or instrument == 'eM' else '',#testing
                     #   'combine': '', #if VLA-C-config
                       'gaintype': 'T',
                       'calmode': 'p',
                       'minsnr': 0.1 if instrument == 'eM' else 0.1,
                       'spwmap': [],
                       'nsigma_automask': '4.0',
                       'nsigma_autothreshold': '2.0',
                       # 'uvtaper' : [''], #if VLA-C-config
                     #   'uvtaper': [taper_size] if receiver in ('S', 'C', 'X', 'Ku', 'K', 'Ka', 'Q') or
                     #                              instrument == 'eM' else [''],#testing
                       'uvtaper' : [''] if general_settings['allow_tapper'] == False else ([taper_size] if receiver in ('C','X','K', 'Ka', 'Ku', 'Q') or instrument == 'eM' else ['']),
                       'with_multiscale': False if multi_config else True,
                     #   'with_multiscale': False if receiver in ('K', 'Ka', 'Ku', 'Q') or instrument ==
                     #                               'eM' else True,
                       # 'scales': '0,5,10,20,40',
                       'scales': 'None',
                       'maxmscales': '4',
                       'compare_solints': False},
                'ap1': {
                    #    'robust': -0.25 if multi_config else (0.5 if instrument == 'eM' else 1.0),
                        # 'robust': 0.5 if multi_config else (0.5 if general_settings['allow_tapper'] else 0.75), #testing
                       'robust': 0.5 if general_settings['allow_tapper'] else 1.0, #testing
                     #    'solint': 'inf' if instrument == 'eM' else ('240s' if receiver in ('X', 'K', 'Ka', 'Ku') else '192s'), # testing
                        # 'solint': '240s' if instrument == 'eM' or receiver in ('K', 'Ka', 'Ku', 'Q') else '96s',
                        'solint': '120s' if general_settings['force_combine_spw'] else ('120s' if general_settings['allow_combine_spw'] else ('240s' if instrument == 'eM' else '90s')),
                        # 'sigma_mask': 10.0 if multi_config else (12.0 if receiver in ('K', 'Ka', 'Ku', 'Q') or instrument == 'eM' else 15.0),
                        'sigma_mask': 10.0, #test
                        'mask_grow_iterations': 4 if multi_config else 4,
                        'combine':  'scan' if general_settings['allow_combine_spw'] == False else ('scan,spw' if receiver in ('K', 'Ka', 'Ku', 'Q') or instrument == 'eM' else 'scan'), #testing
                        #   'combine':  'scan', #testing
                        # 'combine':  '' if general_settings['allow_combine_spw'] == False else ('spw' if receiver in ('K', 'Ka', 'Ku', 'Q') or instrument == 'eM' else ''), #testing
                     #    'combine': 'scan', # testing july 2025
                    #     'combine': 'spw' if receiver in ('C', 'X', 'K', 'Ka', 'Ku', 'Q') or instrument == 'eM' else '',#testing
                     #    'combine': '', #if VLA-C-config
                        'gaintype': 'T',
                        'calmode': 'ap',
                        'minsnr': 0.1 if instrument == 'eM' else 0.1,
                        'spwmap': [],
                        'nsigma_automask' : '4.0' if multi_config else '4.0',
                        'nsigma_autothreshold' : '2.0' if multi_config else '2.0',
                        # 'uvtaper' : [''], #if VLA-C-config
                     #   'uvtaper': [taper_size] if receiver in ('S', 'C', 'X', 'Ku', 'K', 'Ka') or
                     #                              instrument == 'eM' else [''],#testing
                        'uvtaper' : [''] if general_settings['allow_tapper'] == False else ([taper_size] if receiver in ('C','X','K', 'Ka', 'Ku', 'Q') or instrument == 'eM' else ['']),
                        'with_multiscale': True,
                        # 'with_multiscale': False if multi_config else True,
                        # 'with_multiscale' : False,
                     #    'with_multiscale': False if receiver in ('K', 'Ka', 'Ku') or instrument ==
                     #                                'eM' else True,
                        # 'scales': '0,5,20,40',
                        'scales': 'None',
                        'maxmscales': '4',
                        # 'maxmscales': '8', #testing eM + VLA
                        'compare_solints' : False},
                }


# params_faint = {'name': 'faint', #test e-MERLIN obss
#                  'p0': {
#                         'robust': -0.5,
#                         'solint': '480s',
#                         'sigma_mask': 40.0,
#                         'mask_grow_iterations': 1,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'T' if instrument == 'eM' else 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1,
#                         'spwmap': [],
#                         'nsigma_automask': '5.0',
#                         'nsigma_autothreshold': '2.5',
#                         'uvtaper' : [''],
#                         'with_multiscale': False,
#                         'scales': 'None',
#                         'maxmscales': '3',
#                         'compare_solints' : False},
#                  'p1': {
#                         'robust': -0.25,
#                         'solint': '480s',
#                         'sigma_mask': 20.0,
#                         'mask_grow_iterations': 2,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1,
#                         'spwmap': [],
#                         'nsigma_automask': '5.0',
#                         'nsigma_autothreshold': '2.5',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '3',
#                         'compare_solints' : False},
#                  'p2': {
#                         'robust': 0.25 if general_settings['allow_tapper'] else 0.5,
#                         'solint': '240s',
#                         'sigma_mask': 15.0,
#                         'mask_grow_iterations': 3,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'T' if instrument == 'eM' else 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1,
#                         'spwmap': [],
#                         'nsigma_automask': '4.0',
#                         'nsigma_autothreshold': '2.0',
#                         'uvtaper' : [taper_size] if general_settings['allow_tapper'] else [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '4',
#                         'compare_solints' : False},
#                  'ap1': {
#                         'robust': 0.5 if general_settings['allow_tapper'] else 0.5,
#                         'solint': '240s',
#                         'sigma_mask': 12.0,
#                         'mask_grow_iterations': 3,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'T' if instrument == 'eM' else 'G',
#                         'calmode': 'ap',
#                         'minsnr': 0.1,
#                         'spwmap': [],
#                         'nsigma_automask': '4.0',
#                         'nsigma_autothreshold': '2.0',
#                         'uvtaper' : [taper_size] if general_settings['allow_tapper'] else [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '4' if instrument == 'eM' else '6',
#                         'compare_solints' : False},
#                  }


"""
Selfcal parameters to be used for standard sources, 
with a total integrated flux density between 20 and 50 mJy.
"""
params_standard_1 = {'name': 'standard_1',
                   'p0': {
                        #   'robust': -1.0 if multi_config else (-0.25 if receiver in ('K', 'Ka') or instrument == 'eM' else 0.0),
                        #   'robust': -0.25 if multi_config else (-0.25 if receiver in ('K', 'Ka') or instrument == 'eM' else 0.0),
                        #   'robust': -0.25 if multi_config else 0.0,
                          'robust': -0.25 if multi_config else (0.75 if receiver in ('Q') else (-0.5 if receiver in ('C', 'S', 'L') else 0.25)),
                        #   'robust': 0.5,
                        #   'solint': '196s' if instrument == 'eM' or receiver in ('K', 'Ka', 'Ku', 'Q') else '96s',
                          'solint' : '120s' if general_settings['force_combine_spw'] else ('120s' if general_settings['allow_combine_spw'] else ('240s' if receiver in ('K', 'Ka', 'Ku', 'Q') or instrument == 'eM' else '90s')),
                          'sigma_mask': 8.0 if multi_config else (20.0 if receiver in ('Ku', 'K', 'Ka') or instrument == 'eM' else (8.0 if receiver in ('Q') else 50.0)),
                          'mask_grow_iterations': 2 if multi_config else 3,
                          'combine':  'scan' if general_settings['allow_combine_spw'] == False else ('scan,spw' if receiver in ('X', 'K', 'Ka', 'Ku', 'Q') or instrument == 'eM' else 'scan'),
                     #      'combine': 'scan', # testing july 2025
                     #      'combine': 'spw', #testing
                          'gaintype': 'G',
                          'calmode': 'p',
                     #      'minsnr': 1.0 if receiver in ('X', 'Ku', 'K', 'Ka') or instrument == 'eM' else 1.0,
                          'minsnr': 0.1, #testing
                          'spwmap': [],
                        #   'nsigma_automask' : '5.0',
                        #   'nsigma_autothreshold' : '2.5',
                          'nsigma_automask' : '2.0' if receiver in ('Q') else '5.0',       #test M82 A-Ka and A-Q band 
                          'nsigma_autothreshold' : '1.1' if receiver in ('Q')  else '2.0', #test M82 A-Ka and A-Q band 
                          'uvtaper' : [''],
                          'with_multiscale' : True,
                          # 'scales' : '0,5,20',
                          'scales': 'None',
                          'maxmscales': '3',
                          'compare_solints' : False},
                   'p1': {
                        #   'robust': 0.5,
                         #  'robust': 0.0 if receiver in ('K', 'Ka') or instrument == 'eM' else 0.0,
                          'robust': -0.25 if multi_config else (0.25 if receiver in ('K', 'Ka', 'Q') or instrument == 'eM' else (-0.25 if receiver in ('C', 'S', 'L') else 0.0)),
                          'solint' : '120s' if receiver in ('K', 'Ka', 'Ku', 'Q') or instrument == 'eM' else '96s',
                     #      'sigma_mask': 30.0 if receiver in ('C', 'X', 'K', 'Ka', 'Ku') or instrument == 'eM' else 60.0,#testing
                          'sigma_mask': 20.0 if receiver in ('K', 'Ka', 'Ku', 'Q') or instrument == 'eM' else 30.0,
                          'mask_grow_iterations': 4,
                          'combine':  'scan' if general_settings['allow_combine_spw'] == False else ('scan,spw' if receiver in ('K', 'Ka', 'Ku', 'Q') or instrument == 'eM' else 'scan'),
                     #      'combine': 'scan', # testing july 2025
                     #      'combine': 'spw' if instrument == 'eM' else '', #testing
                     #      'combine': 'spw', #testing
                          'gaintype': 'G' if receiver in ('Ku','K', 'Ka', 'Q') or instrument == 'eM' else 'G',
                          'calmode': 'p',
                     #      'minsnr': 1.0,
                          'minsnr': 0.1 if instrument == 'eM' else 0.1, #testing
                          'spwmap': [],
                          'nsigma_automask' : '5.0',
                          'nsigma_autothreshold' : '2.5',
                          'uvtaper' : [''],
                          # 'uvtaper' : [taper_size] if general_settings['allow_tapper'] else [''], #testing
                          'with_multiscale': False if multi_config else True,
                          # 'scales' : '0,5,20',
                          'scales': 'None',
                          'maxmscales': '4',
                          'compare_solints' : False},
                   'p2': {
                     #      'robust': 0.5 if instrument == 'eM' else 0.5,
                          'robust': 0.75 if receiver in ('Q') else 0.5,
                          'solint': '60s',
                     #      'sigma_mask': 20.0 if receiver in ('C', 'X', 'K', 'Ka', 'Ku') or instrument == 'eM' else 30.0,#testing
                          'sigma_mask': 15.0 if receiver in ('K', 'Ka', 'Ku', 'Q') or instrument == 'eM' else 20.0,
                          'mask_grow_iterations': 5,
                          'combine':  'scan' if general_settings['allow_combine_spw'] == False else ('scan,spw' if receiver in ('K', 'Ka', 'Ku', 'Q') or instrument == 'eM' else 'scan'),
                     #      'combine': 'scan', # testing july 2025
                     #      'combine': 'spw' if receiver in ('K', 'Ka', 'Ku', 'Q') or instrument == 'eM' else '',
                     #      'combine': 'spw' if instrument == 'eM' else '', #testing
                     #      'combine': 'spw', #testing
                          'gaintype': 'T',
                          'calmode': 'p',
                     #      'minsnr': 1.0,
                          'minsnr': 0.1 if instrument == 'eM' else 0.1, #testing
                          'spwmap': [],
                          'nsigma_automask' : '4.0',
                          'nsigma_autothreshold' : '2.0',
                          # 'uvtaper' : [''],
                          'uvtaper' : [''] if general_settings['allow_tapper'] == False else ([taper_size] if receiver in ('C','X','K', 'Ka', 'Ku', 'Q') or instrument == 'eM' else ['']),
                          #testing
                     #      'uvtaper': [taper_size] if receiver in ('C', 'Ku', 'K', 'Ka') or instrument == 'eM' else [''],#testing
                     #      'uvtaper' : [taper_size] if receiver in ('X', 'Ku', 'K', 'Ka') or instrument == 'eM' else [''],
                          'with_multiscale': False if multi_config else True,
                          # 'scales': '0,5,10,20,40',
                          'scales': 'None',
                          'maxmscales': '6',
                          'compare_solints' : False},
                   'ap1': {
                        #    'robust': -1.0 if multi_config else (0.5 if instrument == 'eM' else 1.0),
                           'robust': 0.25 if multi_config else (0.5 if instrument == 'eM' else (0.75  if receiver in ('K','Ka','Q') else 0.5)),
                        #    'robust': 0.5,
                        #    'solint': '196s' if receiver in ('K', 'Ka', 'Ku', 'Q') or instrument == 'eM' else '96s',
                           'solint' : '120s' if general_settings['force_combine_spw'] else ('120s' if general_settings['allow_combine_spw'] else ('240s' if receiver in ('K', 'Ka', 'Ku', 'Q') or instrument == 'eM' else '90s')),
                           'sigma_mask': 6.0 if multi_config else (12.0 if receiver in ('K', 'Ka', 'Ku') or instrument == 'eM' else (6.0 if receiver in ('Q') else 12.0)),
                           'mask_grow_iterations': 2 if multi_config else 5,
                           'combine':  'scan' if general_settings['allow_combine_spw'] == False else ('scan,spw' if receiver in ('X', 'K', 'Ka', 'Ku', 'Q') or instrument == 'eM' else 'scan'),
                     #       'combine': 'scan', # testing july 2025
                     #       'combine': 'spw' if instrument == 'eM' else '', #testing
                     #       'combine': 'spw', #testing
                           'gaintype': 'T' if instrument == 'eM' else 'G',
                           'calmode': 'ap',
                     #       'minsnr': 1.0,
                           'minsnr': 0.1 if instrument == 'eM' else 0.1, #testing
                           'spwmap': [],
                        #    'nsigma_automask' : '4.0',
                        #    'nsigma_autothreshold' : '2.0',
                           'nsigma_automask' : '2.0' if receiver in ('Q') else '4.0',       #test M82 A-Ka and A-Q band 
                           'nsigma_autothreshold' : '1.1' if receiver in ('Q')  else '2.0', #test M82 A-Ka and A-Q band 
                           # 'uvtaper' : [''],
                           'uvtaper' : [''] if general_settings['allow_tapper'] == False else ([taper_size] if receiver in ('C','X','K', 'Ka', 'Ku', 'Q') or instrument == 'eM' else ['']),
                     #       'uvtaper': [taper_size] if receiver in ('C', 'Ku', 'K', 'Ka') or instrument == 'eM' else [''],#testing
                     #       'uvtaper' : [taper_size] if receiver in ('Ku', 'K', 'Ka') or instrument == 'eM' else [''],
                           'with_multiscale': False if multi_config else True,
                     #       'with_multiscale': False if receiver in ('K', 'Ka', 'Ku') or
                     #                                   instrument == 'eM' else True,
                           # 'scales': '0,5,10,20,40',
                           'scales': 'None',
                           'maxmscales': '4' if receiver in ('K', 'Ka', 'Q') or instrument == 'eM' else '6',
                           'compare_solints' : False},
                 }


# params_standard_1 = {'name': 'standard_1', #test e-MERLIN obss
#                  'p0': {
#                         'robust': -0.5,
#                         'solint': '120s',
#                         'sigma_mask': 50.0,
#                         'mask_grow_iterations': 1,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1,
#                         'spwmap': [],
#                         'nsigma_automask': '5.0',
#                         'nsigma_autothreshold': '2.5',
#                         'uvtaper' : [''],
#                         'with_multiscale': False,
#                         'scales': 'None',
#                         'maxmscales': '3',
#                         'compare_solints' : False},
#                  'p1': {
#                         'robust': 0.0,
#                         'solint': '90s',
#                         'sigma_mask': 40.0,
#                         'mask_grow_iterations': 2,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1,
#                         'spwmap': [],
#                         'nsigma_automask': '5.0',
#                         'nsigma_autothreshold': '2.5',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '3',
#                         'compare_solints' : False},
#                  'p2': {
#                         'robust': 0.5,
#                         'solint': '60s',
#                         'sigma_mask': 25.0,
#                         'mask_grow_iterations': 3,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'T',
#                         'calmode': 'p',
#                         'minsnr': 0.1,
#                         'spwmap': [],
#                         'nsigma_automask': '4.0',
#                         'nsigma_autothreshold': '2.0',
#                         'uvtaper' : [taper_size] if general_settings['allow_tapper'] else [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '4',
#                         'compare_solints' : False},
#                  'ap1': {
#                         'robust': 0.75,
#                         'solint': '120s',
#                         'sigma_mask': 15.0,
#                         'mask_grow_iterations': 3,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'T',
#                         'calmode': 'ap',
#                         'minsnr': 0.1,
#                         'spwmap': [],
#                         'nsigma_automask': '4.0',
#                         'nsigma_autothreshold': '2.0',
#                         'uvtaper' : [taper_size] if general_settings['allow_tapper'] else [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '4' if instrument == 'eM' else '6',
#                         'compare_solints' : False},
#                  }

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
                   'p0': {
                        #  'robust': -1.0 if multi_config else (-0.5 if receiver in ('Ku', 'K', 'Ka') or instrument == 'eM' else -1.0),
                         'robust': -0.25 if multi_config else (-0.5 if receiver in ('Ku', 'K', 'Ka') or instrument == 'eM' else (0.5 if receiver in ('Q') else -1.0)),
                        #  'robust': -0.25 if multi_config else (0.5 if receiver in ('Ku', 'K', 'Ka') or instrument == 'eM' else -1.0), #test Arp299 eM-C
                        #   'robust': -0.25 if multi_config else (-0.25 if receiver in ('Ku', 'K', 'Ka') or instrument == 'eM' else -1.0), #test2 Arp299 eM-C
                        #   'solint' : '120s' if instrument == 'eM' else '90s',
                          'solint' : '90s' if general_settings['force_combine_spw'] else ('90s' if general_settings['allow_combine_spw'] else ('120s' if receiver in ('K', 'Ka', 'Ku', 'Q') or instrument == 'eM' else '90s')),
                        #   'sigma_mask': 10.0 if multi_config else (15.0 if instrument == 'eM' else 60.0),
                          'sigma_mask': 15.0 if multi_config else (20.0 if instrument == 'eM' else (30.0 if receiver in ('X', 'K', 'Ka', 'Ku') else (12.0 if receiver in ('Q') else 60.0))),
                        #   'sigma_mask': 40.0, #testing
                        #   'mask_grow_iterations': 2 if multi_config else 4,
                          'mask_grow_iterations': 4 if multi_config else 4, #testing
                        #   'combine': 'spw' if instrument == 'eM' else '',
                          'combine': 'scan' if general_settings['allow_combine_spw'] == False else ('scan,spw' if receiver in ('K', 'Ka', 'Ku', 'Q') or instrument == 'eM' else 'scan'),
                        #   'combine':  'scan' if general_settings['allow_combine_spw'] == False
                        #   'combine': 'spw', # testing july 2025
                          'gaintype': 'T',
                          'calmode': 'p',
                          'minsnr': 0.1 if receiver in ('K', 'Ku', 'Ka', 'Q') or instrument == 'eM' else 2.0,
                          'spwmap': [],
                          'nsigma_automask' : '4.0' if instrument == 'eM' else ('2.0' if receiver in ('Ka', 'Q') else '6.0'),
                          'nsigma_autothreshold' : '2.0' if instrument == 'eM' else ('1.5' if receiver in ('Ka', 'Q') else '3.0'),
                          'uvtaper' : [''],
                          'with_multiscale': True if multi_config else (False if receiver in ('K', 'Ka', 'Ku', 'Q') or instrument == 'eM' else True),
                          # 'scales': '0,5,10',
                          'scales': 'None',
                          'maxmscales': '3',
                          'compare_solints' : False},
                   'p1': {
                        #   'robust': -0.25 if multi_config else (0.0 if receiver in ('K', 'Ku', 'Ka', 'Q') or instrument == 'eM' else -0.5),
                        #   'robust': -0.25 if multi_config else (0.25 if receiver in ('K', 'Ku', 'Ka', 'Q') or instrument == 'eM' else -0.5), #testing
                          'robust': -0.25 if multi_config else (0.0 if receiver in ('Ku', 'K') or instrument == 'eM' else (0.25 if receiver in ('Ka', 'Q') else -0.5)),
                        #   'robust': -0.25 if multi_config else (0.5 if receiver in ('K', 'Ku', 'Ka', 'Q') or instrument == 'eM' else -0.5), #test Arp299 eM-C
                        #   'robust': -0.25 if multi_config else (-0.25 if receiver in ('K', 'Ku', 'Ka', 'Q') or instrument == 'eM' else -0.5), #test2 Arp299 eM-C
                          'solint' : '90s' if instrument == 'eM' else '60s',
                        #   'sigma_mask': 20.0 if receiver in ('C', 'X', 'K', 'Ka', 'Ku') or instrument == 'eM' else 30.0,
                          'sigma_mask': 15.0 if multi_config else (15.0 if instrument == 'eM' else (30.0 if receiver in ('X', 'K', 'Ka', 'Ku', 'Q') else 40.0)),
                        #   'mask_grow_iterations': 2 if multi_config else 4,
                          'mask_grow_iterations': 4 if multi_config else 4, #testing
                          'combine':  'scan' if general_settings['allow_combine_spw'] == False else ('scan,spw' if receiver in ('K', 'Ka', 'Ku', 'Q') or instrument == 'eM' else 'scan'),
                         #  'combine': 'spw', # testing july 2025
                     #      'combine': 'spw' if instrument == 'eM' else '', #needs to be tested
                          'gaintype': 'T' if instrument == 'eM' else 'G',
                          'calmode': 'p',
                          'minsnr': 0.1 if receiver in ('K', 'Ku', 'Ka', 'Q') or instrument == 'eM' else 2.0,
                          'spwmap': [],
                        #   'nsigma_automask' : '5.0',
                          'nsigma_automask' : '3.5' if instrument == 'eM' else '5.0',
                        #   'nsigma_autothreshold' : '2.5',
                          'nsigma_autothreshold' : '1.8' if instrument == 'eM' else '2.5',
                          'uvtaper' : [''],
                     #      'uvtaper' : [taper_size] if instrument == 'eM' else [''],#testing
                          'with_multiscale': False if multi_config else True,
                          # 'scales': '0,5,10,20',
                          'scales': 'None',
                          'maxmscales': '3' if instrument == 'eM' or receiver in ('K', 'Ka', 'Ku', 'Q') else '4',
                          'compare_solints' : False},
                   'p2': {
                        #   'robust': 0.0 if multi_config else (0.25 if instrument == 'eM' else (0.75 if receiver in ('K', 'Ku', 'Ka', 'Q') else 0.5)),
                          'robust': 0.0 if multi_config else (0.5 if receiver in ('Ku', 'K') or instrument == 'eM' else (0.75 if receiver in ('Ka', 'Q') else 0.5)),
                        #   'robust': 0.0 if multi_config else (0.5 if instrument == 'eM' else 0.5), #test Arp299 eM-C
                        #   'robust': 0.0 if multi_config else (-0.25 if instrument == 'eM' else 0.5), #test2 Arp299 eM-C
                          'solint': '60s' if instrument == 'eM' else '30s',
                        #   'sigma_mask': 15.0 if receiver in ('C', 'X', 'K', 'Ka', 'Ku') or instrument == 'eM' else 20.0,
                          'sigma_mask': 10.0 if multi_config else (12.0 if instrument == 'eM' else (20.0 if receiver in ('X', 'K', 'Ka', 'Ku', 'Q') else 30.0)),
                        #   'mask_grow_iterations': 3 if multi_config else 4,
                          'mask_grow_iterations': 6 if multi_config else 4, #testing
                          'combine':  'scan' if general_settings['allow_combine_spw'] == False else ('scan,spw' if receiver in ('K', 'Ka', 'Ku', 'Q') or instrument == 'eM' else 'scan'),
                         #  'combine': 'spw', # testing july 2025
                     #      'combine': 'spw' if instrument == 'eM' else '',
                          'gaintype': 'T',
                          'calmode': 'p',
                          'minsnr': 0.1 if receiver in ('K', 'Ku', 'Ka') or instrument == 'eM' else 2.0,
                          'spwmap': [],
                        #   'nsigma_automask' : '5.0',
                          'nsigma_automask' : '4.0',
                        #   'nsigma_autothreshold' : '2.5',
                          'nsigma_autothreshold' : '2.0',
                        #   'uvtaper' : [''],
                          'uvtaper' : [''] if general_settings['allow_tapper'] == False else ([taper_size] if receiver in ('C','X','K', 'Ka', 'Ku', 'Q') or instrument == 'eM' else ['']),
                     #      'uvtaper': [taper_size] if receiver in ('Ku', 'K', 'Ka') or instrument == 'eM' else [''],
                     #      'uvtaper': [taper_size] if receiver in ('S', 'C', 'Ku', 'K', 'Ka') or instrument == 'eM' else [''],#testing
                          'with_multiscale': False if multi_config else True,
                          # 'scales': '0,5,10,20,40',
                          'scales': 'None',
                          'maxmscales': '4' if instrument == 'eM' or receiver in ('K', 'Ka', 'Ku', 'Q') else '6',
                          'compare_solints' : False},
                   'ap1': {
                        #    'robust': -0.25 if multi_config else 0.75,
                        #    'robust': 0.25 if multi_config else 0.5,
                        #    'robust': 0.25 if multi_config else (0.5 if instrument == 'eM' else (1.0 if receiver in ('K', 'Ka', 'Q') else 0.5)), #testing
                           'robust': 0.25 if multi_config else (0.5 if receiver in ('Ku', 'K') or instrument == 'eM' else (0.75 if receiver in ('Ka', 'Q') else 0.5)),
                        #    'robust': 0.25 if multi_config else 0.5, #test Arp299 eM-C
                        #    'robust': 0.25 if multi_config else -0.25, #test2 Arp299 eM-C
                        #    'solint': '120s' if instrument == 'eM' else '60s',
                           'solint' : '60s' if general_settings['force_combine_spw'] else ('60s' if general_settings['allow_combine_spw'] else ('120s' if receiver in ('K', 'Ka', 'Ku', 'Q') or instrument == 'eM' else '60s')),
                        #    'sigma_mask': 6.0 if multi_config else (10.0 if receiver in ('K', 'Ka', 'Ku') or instrument == 'eM' else 12.0),
                           'sigma_mask': 8.0 if multi_config else (10.0 if instrument == 'eM' else (15.0 if receiver in ('X', 'K', 'Ka', 'Ku') else (10.0 if receiver in ('Q') else 20.0))),
                        #    'mask_grow_iterations': 3 if multi_config else 5,
                           'mask_grow_iterations': 6 if multi_config else 5, #testing
                           'combine':  'scan' if general_settings['allow_combine_spw'] == False else ('scan,spw' if receiver in ('K', 'Ka', 'Ku', 'Q') or instrument == 'eM' else 'scan'),
                         #   'combine': 'spw', # testing july 2025
                     #       'combine': 'spw' if instrument == 'eM' else '',
                           'gaintype': 'G',
                           'calmode': 'ap',
                           'minsnr': 0.1 if receiver in ('K', 'Ku', 'Ka', 'Q') or instrument == 'eM' else 2.0,
                           'spwmap': [],
                           'nsigma_automask' : '5.0' if multi_config else ('2.0' if receiver in ('Ka', 'Q') else '4.0'),
                           'nsigma_autothreshold' : '2.5' if multi_config else ('1.5' if receiver in ('Ka', 'Q') else '2.0'),
                           # 'uvtaper' : [''],
                           'uvtaper' : [''] if general_settings['allow_tapper'] == False else ([taper_size] if receiver in ('C','X','K', 'Ka', 'Ku', 'Q') or instrument == 'eM' else ['']),
                     #       'uvtaper': [taper_size] if receiver in ('S', 'C', 'Ku', 'K', 'Ka') or instrument == 'eM' else [''],#testing
                           'with_multiscale': False if multi_config else True,
                           # 'scales': '0,5,10,20,40',
                           'scales': 'None',
                           'maxmscales': '5' if instrument == 'eM' or receiver in ('K', 'Ka', 'Ku', 'Q') else '8',
                           'compare_solints' : False},
                 }

"""
Selfcal parameters to be used for bright sources, 
with a total integrated flux density above 0.1 Jy.
"""
# params_bright = {'name': 'bright',
#                  'p0': {
#                         # 'robust': -0.5 if receiver in ('K', 'Ku', 'Ka') or instrument == 'eM' else -1.5,
#                         'robust': -0.5 if multi_config else (-0.75 if receiver in ('Ku', 'K', 'Ka') or instrument == 'eM' else -1.5),
#                         'solint': '60s' if general_settings['allow_combine_spw'] and instrument == 'eM' else ('90s' if instrument == 'eM' else '60s'),
#                      #    'sigma_mask': 60,
#                         # 'sigma_mask': 50.0 if multi_config else 80.0,
#                         # 'sigma_mask': 50.0 if multi_config else 100.0, #testing
#                         'sigma_mask': 30.0 if multi_config else (40.0 if instrument == 'eM' else (50.0 if receiver in ('X', 'K', 'Ka', 'Ku') else 120.0)),  #needs condition for Q band also
#                         # 'mask_grow_iterations': 2,
#                         'mask_grow_iterations': 3 if multi_config else 3, #testing
#                         # 'combine': 'scan,spw' if instrument == 'eM' else '',
#                         'combine': 'scan' if general_settings['allow_combine_spw'] == False else ('scan,spw' if receiver in ('K', 'Ka', 'Ku') or instrument == 'eM' else 'scan'),
#                         'gaintype': 'T' if instrument == 'eM' else 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' or receiver in ('K', 'Ka', 'Ku', 'Q') else 1.0,
#                         'spwmap': [],
#                         'nsigma_automask': '6.0',
#                         'nsigma_autothreshold': '3.0',
#                         'uvtaper' : [''],
#                         # 'with_multiscale' : False,
#                         'with_multiscale': True if multi_config else (False if receiver in ('K', 'Ka', 'Ku', 'Q') or instrument == 'eM' else True),
#                         # 'scales': '0,5,10',
#                         'scales': 'None',
#                         'maxmscales': '3',
#                         'compare_solints' : False},
#                  'p1': {
#                         # 'robust': -0.5 if receiver in ('K', 'Ka') or instrument == 'eM' else -0.75, #needs condition for other setups
#                         'robust': -0.5 if receiver in ('K', 'Ka') or instrument == 'eM' else -1.0, #testing 
#                         'solint' : '60s' if general_settings['allow_combine_spw'] and instrument == 'eM' else ('90s' if instrument == 'eM' else '36s'),
#                      #    'sigma_mask': 30.0 if receiver in ('X', 'K', 'Ka', 'Ku') or instrument == 'eM' else 60.0,
#                         # 'sigma_mask': 50.0 if receiver in ('X', 'K', 'Ka', 'Ku') or instrument == 'eM' else 60.0, #testing
#                         # 'sigma_mask': 50.0 if receiver in ('X', 'K', 'Ka', 'Ku') or instrument == 'eM' else 80.0, #testing2
#                         'sigma_mask': 20.0 if multi_config else (30.0 if instrument == 'eM' else (50.0 if receiver in ('X', 'K', 'Ka', 'Ku') else 80.0)),
#                         # 'mask_grow_iterations': 4,
#                         'mask_grow_iterations': 4 if multi_config else 3, #testing
#                      #    'combine': 'spw' if instrument == 'eM' else '', #testing
#                         'combine': 'scan,spw' if instrument == 'eM' and general_settings['allow_combine_spw'] else 'scan', #testing
#                         # 'gaintype': 'T' if instrument == 'eM' else 'G',
#                         'gaintype': 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' or receiver in ('K', 'Ka', 'Ku', 'Q') else 1.0,
#                         'spwmap': [],
#                         'nsigma_automask': '5.0' if instrument == 'eM' else '6.0',
#                         'nsigma_autothreshold': '3.0' if instrument == 'eM' else '3.0',
#                         'uvtaper' : [''],
#                      #    'uvtaper' : [taper_size] if instrument == 'eM' else [''], # testing
#                         # 'with_multiscale': False if multi_config else True,
#                         'with_multiscale': True,
#                         # 'with_multiscale' : False, #testing M82
#                         # 'scales': '0,5,10,20',
#                         'scales': 'None',
#                         'maxmscales': '4',
#                         'compare_solints': False},
#                  'p2': {
#                     #    'robust': -0.25 if multi_config else (0.25 if receiver in ('K', 'Ka') or instrument == 'eM' else 0.0), #testing M82
#                        'robust': -0.25 if multi_config else (0.25 if receiver in ('K', 'Ka') or instrument == 'eM' else -0.25),
#                         # 'robust': -0.25 if multi_config else 0.0,
#                         'solint': '30s' if general_settings['allow_combine_spw'] and instrument == 'eM' else ('60s' if instrument == 'eM' else '18s'),
#                      #    'sigma_mask': 12.0 if receiver in ('X', 'K', 'Ka', 'Ku') or instrument == 'eM' else 30.0,
#                         # 'sigma_mask': 30.0 if receiver in ('X', 'K', 'Ka', 'Ku') or instrument == 'eM' else 40.0, #testing
#                         # 'sigma_mask': 30.0 if receiver in ('X', 'K', 'Ka', 'Ku') or instrument == 'eM' else 60.0, #testing2
#                         'sigma_mask': 15.0 if multi_config else (15.0 if instrument == 'eM' else (30.0 if receiver in ('X', 'K', 'Ka', 'Ku') else 60.0)), #testing2
#                         # 'mask_grow_iterations': 6,
#                         'mask_grow_iterations': 6 if multi_config else 6, #testing
#                         'combine': 'scan,spw' if instrument == 'eM' and general_settings['allow_combine_spw'] else 'scan', #testing
#                         'gaintype': 'T' if instrument == 'eM' else 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' or receiver in ('K', 'Ka', 'Ku', 'Q') else 1.0,
#                         'spwmap': [],
#                         'uvtaper' : [taper_size] if general_settings['allow_tapper'] else [''], #testing
#                         'nsigma_automask': '4.0' if instrument == 'eM' or receiver in ('K', 'Ka', 'Ku', 'Q') else '5.0',
#                         'nsigma_autothreshold': '2.0' if instrument == 'eM' or receiver in ('K', 'Ka', 'Ku', 'Q') else '2.5',
#                         # 'with_multiscale': False if multi_config else True,
#                         'with_multiscale': True,
#                         # 'with_multiscale' : False, #testing M82
#                         # 'scales': '0,5,10,20,40',
#                         'scales': 'None',
#                         'maxmscales': '5' if instrument == 'eM' or receiver in ('K', 'Ka', 'Ku', 'Q') else '5',
#                         'compare_solints': False},
#                  'ap1': {
#                          'robust': 0.25 if multi_config else (0.5 if receiver in ('K', 'Ka') or instrument == 'eM' else 0.25), #default
#                         #  'robust': -0.25 if multi_config else (0.5 if receiver in ('K', 'Ka') or instrument == 'eM' else 1.0), #2026: testing Arp299 (VLA-A Cband)
#                         #  'robust': 0.25 if multi_config else 0.0,
#                          'solint': '60s' if general_settings['allow_combine_spw'] and instrument == 'eM' else ('90s' if instrument == 'eM' else '36s'),
#                      #     'sigma_mask': 8.0 if receiver in ('X', 'K', 'Ka', 'Ku') or instrument == 'eM' else 20.0,
#                         #  'sigma_mask': 15.0 if receiver in ('X', 'K', 'Ka', 'Ku') or instrument == 'eM' else 30.0, #testing
#                          'sigma_mask': 10.0 if multi_config else (15.0 if instrument == 'eM' else (15.0 if receiver in ('X', 'K', 'Ka', 'Ku') else 15.0)), #testing2
#                         #  'mask_grow_iterations': 2 if multi_config else 6,
#                          'mask_grow_iterations': 6 if multi_config else 4, #testing
#                         #  'combine': 'scan,spw' if instrument == 'eM' and general_settings['allow_combine_spw'] else 'scan', #testing
#                         'combine': 'scan' if general_settings['allow_combine_spw'] == False else ('scan,spw' if receiver in ('K', 'Ka', 'Ku') or instrument == 'eM' else 'scan'),
#                      #     'combine': 'spw' if instrument == 'eM' else '', #testing
#                         #  'gaintype': 'T' if instrument == 'eM' else 'G',
#                          'gaintype': 'T' if instrument == 'eM' else 'G',
#                          'calmode': 'ap',
#                          'minsnr': 0.1 if instrument == 'eM' or receiver in ('K', 'Ka', 'Ku', 'Q') else 1.0,
#                          'spwmap': [],
#                          'uvtaper' : [taper_size] if general_settings['allow_tapper'] else [''], #testing
#                          'nsigma_automask': '4.0' if instrument == 'eM' or receiver in ('K', 'Ka', 'Ku', 'Q') else '5.0',
#                          'nsigma_autothreshold': '2.0' if instrument == 'eM' or receiver in ('K', 'Ka', 'Ku', 'Q') else '2.5',
#                         #  'with_multiscale': False if multi_config else True,
#                          'with_multiscale': True,
#                         #  'with_multiscale' : False, #testing M82
#                          # 'scales': '0,5,10,20,40',
#                          'scales': 'None',
#                          'maxmscales': '6' if instrument == 'eM' or receiver in ('K', 'Ka', 'Ku', 'Q') else '6',
#                          'compare_solints': False},
#                  }



# params_bright = {'name': 'bright', #test Arp220 (C-band); M87 (C; X-band); M82 (C-band); Mrk231 (C-band)
#                  'p0': {
#                         'robust': -2.0,
#                         'solint': '72s',
#                         'sigma_mask': 150.0,
#                         'mask_grow_iterations': 1,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'T' if instrument == 'eM' else 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.2,
#                         'spwmap': [],
#                         'nsigma_automask': '8.0',
#                         'nsigma_autothreshold': '4.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': False,
#                         'scales': 'None',
#                         'maxmscales': '3',
#                         'compare_solints' : False},
#                  'p1': {
#                         'robust': -1.5,
#                         'solint': '36s',
#                         'sigma_mask': 120.0,
#                         'mask_grow_iterations': 3,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.2,
#                         'spwmap': [],
#                         'nsigma_automask': '6.0',
#                         'nsigma_autothreshold': '4.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '3',
#                         'compare_solints' : False},
#                  'p2': {
#                         'robust': -1.0,
#                         'solint': 'int',
#                         'sigma_mask': 80.0,
#                         'mask_grow_iterations': 4,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'T' if instrument == 'eM' else 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.2,
#                         'spwmap': [],
#                         'nsigma_automask': '6.0',
#                         'nsigma_autothreshold': '3.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '4',
#                         'compare_solints' : False},
#                  'ap1': {
#                         'robust': -1.0,
#                         'solint': 'int',
#                         'sigma_mask': 40.0,
#                         'mask_grow_iterations': 4,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'T' if instrument == 'eM' else 'G',
#                         'calmode': 'ap',
#                         'minsnr': 0.2,
#                         'spwmap': [],
#                         'nsigma_automask': '4.0',
#                         'nsigma_autothreshold': '2.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '6',
#                         'compare_solints' : False},
#                  }


params_bright = {'name': 'bright', #test Arp220 (C-band); M87 (C; X-band); M82 (C-band); Mrk231 (C-band)
                 'p0': {
                        'robust': -1.0,
                        'solint': '120s',
                        'sigma_mask': 120.0,
                        'mask_grow_iterations': 1,
                        'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
                        'gaintype': 'T' if instrument == 'eM' else 'G',
                        'calmode': 'p',
                        'minsnr': 0.1 if instrument == 'eM' else 0.1,
                        'spwmap': [],
                        'nsigma_automask': '8.0',
                        'nsigma_autothreshold': '4.0',
                        'uvtaper' : [''],
                        'with_multiscale': False,
                        'scales': 'None',
                        'maxmscales': '3',
                        'compare_solints' : False},
                 'p1': {
                        'robust': -0.5,
                        'solint': '60s',
                        'sigma_mask': 80.0,
                        'mask_grow_iterations': 3,
                        'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
                        'gaintype': 'G',
                        'calmode': 'p',
                        'minsnr': 0.1 if instrument == 'eM' else 0.1,
                        'spwmap': [],
                        'nsigma_automask': '6.0',
                        'nsigma_autothreshold': '4.0',
                        'uvtaper' : [''],
                        'with_multiscale': True,
                        'scales': 'None',
                        'maxmscales': '4',
                        'compare_solints' : False},
                 'p2': {
                        'robust': 0.0,
                        'solint': '16s' if instrument == 'eM' else '30s',
                        'sigma_mask': 40.0,
                        'mask_grow_iterations': 4,
                        'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
                        'gaintype': 'T' if instrument == 'eM' else 'G',
                        'calmode': 'p',
                        'minsnr': 0.1 if instrument == 'eM' else 0.1,
                        'spwmap': [],
                        'nsigma_automask': '4.0',
                        'nsigma_autothreshold': '2.0',
                        'uvtaper' : [''],
                        'with_multiscale': True,
                        'scales': 'None',
                        'maxmscales': '5',
                        'compare_solints' : False},
                 'ap1': {
                        'robust': 0.5,
                        'solint': '32s' if instrument == 'eM' else '30s',
                        'sigma_mask': 15.0,
                        'mask_grow_iterations': 4,
                        'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
                        'gaintype': 'G',
                        'calmode': 'ap',
                        'minsnr': 0.1 if instrument == 'eM' else 0.1,
                        'spwmap': [],
                        'nsigma_automask': '4.0',
                        'nsigma_autothreshold': '2.0',
                        'uvtaper' : [''],
                        'with_multiscale': True,
                        'scales': 'None',
                        'maxmscales': '6',
                        'compare_solints' : False},
                 }

# params_bright = {'name': 'bright', #test Arp220 (C-band); M87 (C; X-band); M82 (C-band); Mrk231 (C-band)
#                  'p0': {
#                         'robust': -1.0,
#                         'solint': '64s',
#                         'sigma_mask': 120.0,
#                         'mask_grow_iterations': 1,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'T' if instrument == 'eM' else 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 3.0,
#                         'spwmap': [],
#                         'nsigma_automask': '8.0',
#                         'nsigma_autothreshold': '4.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': False,
#                         'scales': 'None',
#                         'maxmscales': '3',
#                         'compare_solints' : False},
#                  'p1': {
#                         'robust': -0.5,
#                         'solint': '32s',
#                         'sigma_mask': 80.0,
#                         'mask_grow_iterations': 3,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 3.0,
#                         'spwmap': [],
#                         'nsigma_automask': '6.0',
#                         'nsigma_autothreshold': '4.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '4',
#                         'compare_solints' : False},
#                  'p2': {
#                         'robust': 0.0,
#                         'solint': '16s' if instrument == 'eM' else 'int',
#                         'sigma_mask': 40.0,
#                         'mask_grow_iterations': 4,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'T' if instrument == 'eM' else 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 3.0,
#                         'spwmap': [],
#                         'nsigma_automask': '4.0',
#                         'nsigma_autothreshold': '2.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '6',
#                         'compare_solints' : False},
#                  'ap1': {
#                         'robust': 0.5,
#                         'solint': '32s' if instrument == 'eM' else 'int',
#                         'sigma_mask': 20.0,
#                         'mask_grow_iterations': 4,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'G',
#                         'calmode': 'ap',
#                         'minsnr': 0.1 if instrument == 'eM' else 3.0,
#                         'spwmap': [],
#                         'nsigma_automask': '4.0',
#                         'nsigma_autothreshold': '2.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '8',
#                         'compare_solints' : False},
#                  }


# params_bright = {'name': 'bright', #test Arp220 (C-band); M87 (C; X-band); M82 (C-band); Mrk231 (C-band)
#                  'p0': {
#                         'robust': -0.5,
#                         'solint': '64s',
#                         'sigma_mask': 120.0,
#                         'mask_grow_iterations': 1,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'T' if instrument == 'eM' else 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 3.0,
#                         'spwmap': [],
#                         'nsigma_automask': '8.0',
#                         'nsigma_autothreshold': '4.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': False,
#                         'scales': 'None',
#                         'maxmscales': '3',
#                         'compare_solints' : False},
#                  'p1': {
#                         'robust': 0.0,
#                         'solint': '32s',
#                         'sigma_mask': 80.0,
#                         'mask_grow_iterations': 3,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 3.0,
#                         'spwmap': [],
#                         'nsigma_automask': '6.0',
#                         'nsigma_autothreshold': '4.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '4',
#                         'compare_solints' : False},
#                  'p2': {
#                         'robust': 0.75,
#                         'solint': '16s' if instrument == 'eM' else 'int',
#                         'sigma_mask': 40.0,
#                         'mask_grow_iterations': 4,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'T' if instrument == 'eM' else 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 3.0,
#                         'spwmap': [],
#                         'nsigma_automask': '4.0',
#                         'nsigma_autothreshold': '2.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '6',
#                         'compare_solints' : False},
#                  'ap1': {
#                         'robust': 1.25,
#                         'solint': '32s' if instrument == 'eM' else 'int',
#                         'sigma_mask': 20.0,
#                         'mask_grow_iterations': 4,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'G',
#                         'calmode': 'ap',
#                         'minsnr': 0.1 if instrument == 'eM' else 3.0,
#                         'spwmap': [],
#                         'nsigma_automask': '4.0',
#                         'nsigma_autothreshold': '2.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '8',
#                         'compare_solints' : False},
#                  }

# params_bright = {'name': 'bright', #test Arp220 (C-band); M87 (C; X-band); M82 (C-band); Mrk231 (C-band)
#                  'p0': {
#                         'robust': 0.0,
#                         'solint': '64s',
#                         'sigma_mask': 15.0,
#                         'mask_grow_iterations': 1,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'T' if instrument == 'eM' else 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 3.0,
#                         'spwmap': [],
#                         'nsigma_automask': '8.0',
#                         'nsigma_autothreshold': '4.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': False,
#                         'scales': 'None',
#                         'maxmscales': '3',
#                         'compare_solints' : False},
#                  'p1': {
#                         'robust': 0.5,
#                         'solint': '32s',
#                         'sigma_mask': 15.0,
#                         'mask_grow_iterations': 3,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 3.0,
#                         'spwmap': [],
#                         'nsigma_automask': '6.0',
#                         'nsigma_autothreshold': '4.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '4',
#                         'compare_solints' : False},
#                  'p2': {
#                         'robust': 1.0,
#                         'solint': '16s' if instrument == 'eM' else 'int',
#                         'sigma_mask': 15.0,
#                         'mask_grow_iterations': 4,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'T' if instrument == 'eM' else 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 3.0,
#                         'spwmap': [],
#                         'nsigma_automask': '4.0',
#                         'nsigma_autothreshold': '2.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '6',
#                         'compare_solints' : False},
#                  'ap1': {
#                         'robust': 1.5,
#                         'solint': '32s' if instrument == 'eM' else 'int',
#                         'sigma_mask': 12.0,
#                         'mask_grow_iterations': 4,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'G',
#                         'calmode': 'ap',
#                         'minsnr': 0.1 if instrument == 'eM' else 3.0,
#                         'spwmap': [],
#                         'nsigma_automask': '4.0',
#                         'nsigma_autothreshold': '2.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '8',
#                         'compare_solints' : False},
#                  }


# params_bright = {'name': 'bright', #second pass of sc (e.g. selfcal from a previous sc run)
#                  'p0': {
#                         'robust': 0.5,
#                         'solint': 'inf',
#                         'sigma_mask': 30.0,
#                         'mask_grow_iterations': 3,
#                         'combine': '',
#                         'gaintype': 'T' if instrument == 'eM' else 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 2.0,
#                         'spwmap': [],
#                         'nsigma_automask': '4.0',
#                         'nsigma_autothreshold': '2.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '3',
#                         'compare_solints' : False},
#                  'p1': {'robust': 0.5,
#                         'solint' : '60s',
#                         'sigma_mask': 20.0,
#                         'mask_grow_iterations': 4,
#                         'combine': 'scan,spw' if instrument == 'eM' and general_settings['allow_combine_spw'] else 'scan', #testing
#                         'gaintype': 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 2.0,
#                         'spwmap': [],
#                         'nsigma_automask': '4.0',
#                         'nsigma_autothreshold': '2.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '4',
#                         'compare_solints': False},
#                  'p2': {'robust': 0.75,
#                         'solint': '30s',
#                         'sigma_mask': 15.0,
#                         'mask_grow_iterations': 6,
#                         'combine': 'scan,spw' if instrument == 'eM' and general_settings['allow_combine_spw'] else 'scan', #testing
#                         'gaintype': 'T' if instrument == 'eM' else 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 2.0,
#                         'spwmap': [],
#                         'uvtaper' : [taper_size] if general_settings['allow_tapper'] else [''], #testing
#                         'nsigma_automask': '4.0',
#                         'nsigma_autothreshold': '2.0',
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '4' if instrument == 'eM' else '5',
#                         'compare_solints': False},
#                  'ap1': {
#                          'robust': 0.75,
#                          'solint': '60s',
#                          'sigma_mask': 10.0,
#                          'mask_grow_iterations': 6,
#                         'combine': 'scan' if general_settings['allow_combine_spw'] == False else ('scan,spw' if receiver in ('K', 'Ka', 'Ku') or instrument == 'eM' else 'scan'),
#                          'gaintype': 'T' if instrument == 'eM' else 'G',
#                          'calmode': 'ap',
#                          'minsnr': 0.1 if instrument == 'eM' else 2.0,
#                          'spwmap': [],
#                          'uvtaper' : [taper_size] if general_settings['allow_tapper'] else [''], #testing
#                          'nsigma_automask': '4.0',
#                          'nsigma_autothreshold': '2.0',
#                          'with_multiscale': True,
#                          'scales': 'None',
#                          'maxmscales': '4' if instrument == 'eM' else '6',
#                          'compare_solints': False},
#                  }

# params_bright = {'name': 'bright', #snapshot observations (this is a test for VLA for short observations (e.g. 2min on source).)
#                  'p0': {
#                         # 'robust': -0.5 if receiver in ('K', 'Ku', 'Ka') or instrument == 'eM' else -1.5,
#                         'robust': -0.75,
#                         'solint': '90s' if instrument == 'eM' else '60s',
#                      #    'sigma_mask': 60,
#                         # 'sigma_mask': 50.0 if multi_config else 80.0,
#                         # 'sigma_mask': 50.0 if multi_config else 100.0, #testing
#                         'sigma_mask': 40.0 if instrument == 'eM' else (25.0 if receiver in ('X', 'K', 'Ka', 'Ku') else 50.0),
#                         # 'mask_grow_iterations': 2,
#                         'mask_grow_iterations': 2,
#                         # 'combine': 'scan,spw' if instrument == 'eM' else '',
#                         'combine': 'scan' if general_settings['allow_combine_spw'] == False else ('scan,spw' if receiver in ('K', 'Ka', 'Ku') or instrument == 'eM' else 'scan'),
#                         'gaintype': 'T' if instrument == 'eM' else 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 2.0,
#                         'spwmap': [],
#                         'nsigma_automask': '6.0',
#                         'nsigma_autothreshold': '3.0',
#                         'uvtaper' : [''],
#                         # 'with_multiscale' : False,
#                         'with_multiscale': True if multi_config else True,
#                         # 'scales': '0,5,10',
#                         'scales': 'None',
#                         'maxmscales': '3',
#                         'compare_solints' : False},
#                  'p1': {'robust': -0.25,
#                         'solint' : '90s' if instrument == 'eM' else '30s',
#                      #    'sigma_mask': 30.0 if receiver in ('X', 'K', 'Ka', 'Ku') or instrument == 'eM' else 60.0,
#                         # 'sigma_mask': 50.0 if receiver in ('X', 'K', 'Ka', 'Ku') or instrument == 'eM' else 60.0, #testing
#                         # 'sigma_mask': 50.0 if receiver in ('X', 'K', 'Ka', 'Ku') or instrument == 'eM' else 80.0, #testing2
#                         'sigma_mask': 40.0 if instrument == 'eM' else (25.0 if receiver in ('X', 'K', 'Ka', 'Ku') else 40.0),
#                         # 'mask_grow_iterations': 4,
#                         'mask_grow_iterations': 3,
#                      #    'combine': 'spw' if instrument == 'eM' else '', #testing
#                         'combine': 'scan,spw' if instrument == 'eM' and general_settings['allow_combine_spw'] else 'scan', #testing
#                         'gaintype': 'T' if instrument == 'eM' else 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 2.0,
#                         'spwmap': [],
#                         'nsigma_automask': '5.0',
#                         'nsigma_autothreshold': '2.5',
#                         'uvtaper' : [''],
#                      #    'uvtaper' : [taper_size] if instrument == 'eM' else [''], # testing
#                         # 'with_multiscale': False if multi_config else True,
#                         'with_multiscale': True,
#                         # 'with_multiscale' : False, #testing M82
#                         # 'scales': '0,5,10,20',
#                         'scales': 'None',
#                         'maxmscales': '4',
#                         'compare_solints': False},
#                  'p2': {'robust': 0.25,
#                         # 'robust': -0.25 if multi_config else 0.0,
#                         'solint': '60s' if instrument == 'eM' else '15s',
#                      #    'sigma_mask': 12.0 if receiver in ('X', 'K', 'Ka', 'Ku') or instrument == 'eM' else 30.0,
#                         # 'sigma_mask': 30.0 if receiver in ('X', 'K', 'Ka', 'Ku') or instrument == 'eM' else 40.0, #testing
#                         # 'sigma_mask': 30.0 if receiver in ('X', 'K', 'Ka', 'Ku') or instrument == 'eM' else 60.0, #testing2
#                         'sigma_mask': 15.0 if instrument == 'eM' else (15.0 if receiver in ('X', 'K', 'Ka', 'Ku') else 20.0), #testing2
#                         # 'mask_grow_iterations': 6,
#                         'mask_grow_iterations': 4, #testing
#                         'combine': 'scan,spw' if instrument == 'eM' and general_settings['allow_combine_spw'] else 'scan', #testing
#                         'gaintype': 'T' if instrument == 'eM' else 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 2.0,
#                         'spwmap': [],
#                         'uvtaper' : [''],
#                         'nsigma_automask': '4.0',
#                         'nsigma_autothreshold': '2.0',
#                         # 'with_multiscale': False if multi_config else True,
#                         'with_multiscale': True,
#                         # 'with_multiscale' : False, #testing M82
#                         # 'scales': '0,5,10,20,40',
#                         'scales': 'None',
#                         'maxmscales': '5',
#                         'compare_solints': False},
#                  'ap1': {
#                          'robust': 0.5, #default
#                         #  'robust': -0.25 if multi_config else (0.5 if receiver in ('K', 'Ka') or instrument == 'eM' else 1.0), #2026: testing Arp299 (VLA-A Cband)
#                         #  'robust': 0.25 if multi_config else 0.0,
#                          'solint': '90s' if instrument == 'eM' else '30s',
#                      #     'sigma_mask': 8.0 if receiver in ('X', 'K', 'Ka', 'Ku') or instrument == 'eM' else 20.0,
#                         #  'sigma_mask': 15.0 if receiver in ('X', 'K', 'Ka', 'Ku') or instrument == 'eM' else 30.0, #testing
#                          'sigma_mask': 10.0, #testing2
#                         #  'mask_grow_iterations': 2 if multi_config else 6,
#                          'mask_grow_iterations': 5, #testing
#                         #  'combine': 'scan,spw' if instrument == 'eM' and general_settings['allow_combine_spw'] else 'scan', #testing
#                         'combine': 'scan' if general_settings['allow_combine_spw'] == False else ('scan,spw' if receiver in ('K', 'Ka', 'Ku') or instrument == 'eM' else 'scan'),
#                      #     'combine': 'spw' if instrument == 'eM' else '', #testing
#                         #  'gaintype': 'T' if instrument == 'eM' else 'G',
#                          'gaintype': 'G',
#                          'calmode': 'ap',
#                          'minsnr': 0.1 if instrument == 'eM' else 2.0,
#                          'spwmap': [],
#                          'uvtaper' : [''],
#                          'nsigma_automask': '3.0',
#                          'nsigma_autothreshold': '1.5',
#                         #  'with_multiscale': False if multi_config else True,
#                          'with_multiscale': True,
#                         #  'with_multiscale' : False, #testing M82
#                          # 'scales': '0,5,10,20,40',
#                          'scales': 'None',
#                          'maxmscales': '6',
#                          'compare_solints': False},
#                  }


# params_bright = {'name': 'bright', #test M82
#                  'p0': {
#                         # 'robust': -0.5 if receiver in ('K', 'Ku', 'Ka') or instrument == 'eM' else -1.5,
#                         'robust': -0.5 if multi_config else (-0.5 if receiver in ('Ku', 'K', 'Ka') or instrument == 'eM' else -1.5),
#                         'solint': '90s' if instrument == 'eM' else '60s',
#                      #    'sigma_mask': 60,
#                         # 'sigma_mask': 50.0 if multi_config else 80.0,
#                         # 'sigma_mask': 50.0 if multi_config else 100.0, #testing
#                         'sigma_mask': 30.0 if multi_config else (25.0 if instrument == 'eM' else (50.0 if receiver in ('X', 'K', 'Ka', 'Ku') else 120.0)),
#                         # 'mask_grow_iterations': 2,
#                         'mask_grow_iterations': 3 if multi_config else 3, #testing
#                         # 'combine': 'scan,spw' if instrument == 'eM' else '',
#                         'combine': '' if general_settings['allow_combine_spw'] == False else ('spw' if receiver in ('K', 'Ka', 'Ku') or instrument == 'eM' else ''),
#                         'gaintype': 'T' if instrument == 'eM' else 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 2.0,
#                         'spwmap': [],
#                         'nsigma_automask': '6.0',
#                         'nsigma_autothreshold': '3.0',
#                         'uvtaper' : [''],
#                         'with_multiscale' : True,
#                         # 'with_multiscale': True if multi_config else False,
#                         # 'scales': '0,5,10',
#                         'scales': 'None',
#                         'maxmscales': '3',
#                         'compare_solints' : False},
#                  'p1': {'robust': -0.5 if receiver in ('K', 'Ka') or instrument == 'eM' else -1.0,
#                         'solint' : '90s' if instrument == 'eM' else '30s',
#                      #    'sigma_mask': 30.0 if receiver in ('X', 'K', 'Ka', 'Ku') or instrument == 'eM' else 60.0,
#                         # 'sigma_mask': 50.0 if receiver in ('X', 'K', 'Ka', 'Ku') or instrument == 'eM' else 60.0, #testing
#                         # 'sigma_mask': 50.0 if receiver in ('X', 'K', 'Ka', 'Ku') or instrument == 'eM' else 80.0, #testing2
#                         'sigma_mask': 20.0 if multi_config else (20.0 if instrument == 'eM' else (50.0 if receiver in ('X', 'K', 'Ka', 'Ku') else 80.0)),
#                         # 'mask_grow_iterations': 4,
#                         'mask_grow_iterations': 4 if multi_config else 4, #testing
#                      #    'combine': 'spw' if instrument == 'eM' else '', #testing
#                         'combine': 'spw' if instrument == 'eM' and general_settings['allow_combine_spw'] else '', #testing
#                         'gaintype': 'T' if instrument == 'eM' else 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 2.0,
#                         'spwmap': [],
#                         'nsigma_automask': '5.0' if instrument == 'eM' else '6.0',
#                         'nsigma_autothreshold': '3.0' if instrument == 'eM' else '3.0',
#                         'uvtaper' : [''],
#                      #    'uvtaper' : [taper_size] if instrument == 'eM' else [''], # testing
#                         # 'with_multiscale': False if multi_config else True,
#                         'with_multiscale': True,
#                         'maxmscales': '4',
#                         # 'with_multiscale' : False, #testing M82
#                         # 'scales': '0,5,10,20',
#                         'scales': 'None',
#                         'compare_solints': False},
#                  'p2': {'robust': -0.25 if multi_config else (0.0 if receiver in ('K', 'Ka') or instrument == 'eM' else 0.0), #testing M82
#                         # 'robust': -0.25 if multi_config else 0.0,
#                         'solint': '60s' if instrument == 'eM' else '15s',
#                      #    'sigma_mask': 12.0 if receiver in ('X', 'K', 'Ka', 'Ku') or instrument == 'eM' else 30.0,
#                         # 'sigma_mask': 30.0 if receiver in ('X', 'K', 'Ka', 'Ku') or instrument == 'eM' else 40.0, #testing
#                         # 'sigma_mask': 30.0 if receiver in ('X', 'K', 'Ka', 'Ku') or instrument == 'eM' else 60.0, #testing2
#                         'sigma_mask': 15.0 if multi_config else (15.0 if instrument == 'eM' else (30.0 if receiver in ('X', 'K', 'Ka', 'Ku') else 40.0)), #testing2
#                         # 'mask_grow_iterations': 6,
#                         'mask_grow_iterations': 6 if multi_config else 6, #testing
#                         'combine': 'spw' if instrument == 'eM' and general_settings['allow_combine_spw'] else '', #testing
#                         'gaintype': 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 2.0,
#                         'spwmap': [],
#                         'uvtaper' : [''],
#                         'nsigma_automask': '5.0',
#                         'nsigma_autothreshold': '2.5',
#                         # 'with_multiscale': False if multi_config else True,
#                         'with_multiscale': True,
#                         # 'with_multiscale' : False, #testing M82
#                         # 'scales': '0,5,10,20,40',
#                         'scales': 'None',
#                         'maxmscales': '8',
#                         'compare_solints': False},
#                  'ap1': {
#                          'robust': 0.25 if multi_config else (0.5 if receiver in ('K', 'Ka') or instrument == 'eM' else 0.0), #default
#                         #  'robust': -0.25 if multi_config else (0.5 if receiver in ('K', 'Ka') or instrument == 'eM' else 1.0), #2026: testing Arp299 (VLA-A Cband)
#                         #  'robust': 0.25 if multi_config else 0.0,
#                          'solint': '90s' if instrument == 'eM' else '30s',
#                      #     'sigma_mask': 8.0 if receiver in ('X', 'K', 'Ka', 'Ku') or instrument == 'eM' else 20.0,
#                         #  'sigma_mask': 15.0 if receiver in ('X', 'K', 'Ka', 'Ku') or instrument == 'eM' else 30.0, #testing
#                          'sigma_mask': 15.0 if multi_config else (15.0 if instrument == 'eM' else (15.0 if receiver in ('X', 'K', 'Ka', 'Ku') else 20.0)), #testing2
#                         #  'mask_grow_iterations': 2 if multi_config else 6,
#                          'mask_grow_iterations': 6 if multi_config else 6, #testing
#                         #  'combine': 'scan,spw' if instrument == 'eM' and general_settings['allow_combine_spw'] else 'scan', #testing
#                         'combine': '' if general_settings['allow_combine_spw'] == False else ('spw' if receiver in ('K', 'Ka', 'Ku') or instrument == 'eM' else ''),
#                      #     'combine': 'spw' if instrument == 'eM' else '', #testing
#                         #  'gaintype': 'T' if instrument == 'eM' else 'G',
#                          'gaintype': 'G',
#                          'calmode': 'ap',
#                          'minsnr': 0.1 if instrument == 'eM' else 2.0,
#                          'spwmap': [],
#                          'uvtaper' : [''],
#                          'nsigma_automask': '5.0',
#                          'nsigma_autothreshold': '2.5',
#                         #  'with_multiscale': False if multi_config else True,
#                          'with_multiscale': True,
#                         #  'with_multiscale' : False, #testing M82
#                          # 'scales': '0,5,10,20,40',
#                          'scales': 'None',
#                          'maxmscales': '8',
#                          'compare_solints': False},
#                  }

# params_bright = {'name': 'bright', #test Arp220 (C-band); M87 (C; X-band); M82 (C-band); Mrk231 (C-band)
#                  'p0': {
#                         'robust': -1.0,
#                         'solint': '60s',
#                         'sigma_mask': 120.0,
#                         'mask_grow_iterations': 1,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'T' if instrument == 'eM' else 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.2,
#                         'spwmap': [],
#                         'nsigma_automask': '8.0',
#                         'nsigma_autothreshold': '4.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': False,
#                         'scales': 'None',
#                         'maxmscales': '3',
#                         'compare_solints' : False},
#                  'p1': {
#                         'robust': -1.0,
#                         'solint': '30s',
#                         'sigma_mask': 80.0,
#                         'mask_grow_iterations': 3,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.2,
#                         'spwmap': [],
#                         'nsigma_automask': '6.0',
#                         'nsigma_autothreshold': '4.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '4',
#                         'compare_solints' : False},
#                  'p2': {
#                         'robust': -1.0,
#                         'solint': '10s',
#                         'sigma_mask': 40.0,
#                         'mask_grow_iterations': 4,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'T' if instrument == 'eM' else 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.2,
#                         'spwmap': [],
#                         'nsigma_automask': '4.0',
#                         'nsigma_autothreshold': '2.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '4',
#                         'compare_solints' : False},
#                  'ap1': {
#                         'robust': -1.0,
#                         'solint': '10s',
#                         'sigma_mask': 20.0,
#                         'mask_grow_iterations': 4,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'T' if instrument == 'eM' else 'G',
#                         'calmode': 'ap',
#                         'minsnr': 0.2,
#                         'spwmap': [],
#                         'nsigma_automask': '4.0',
#                         'nsigma_autothreshold': '2.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '4',
#                         'compare_solints' : False},
#                  }

# params_very_faint = {'name': 'very_faint', #global for everything; first loop / faint sources or with diffuse/complex emission
#                  'p0': {
#                         'robust': 0.75 if receiver in ('Q') else 1.0,
#                         'solint': '240s' if instrument == 'eM' else '120s',
#                         'sigma_mask': 10.0 if instrument == 'eM' else (12.0 if receiver in ('Q') else 15.0),
#                         'mask_grow_iterations': 3,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'T' if instrument == 'eM' else 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 0.1,
#                         'spwmap': [],
#                         'nsigma_automask': '3.0' if receiver in ('Q') else '4.0',
#                         'nsigma_autothreshold': '2.0',
#                         'uvtaper' : [taper_size] if general_settings['allow_tapper'] else [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '3',
#                         'compare_solints' : False},
#                  'ap1': {
#                         # 'robust': 0.75 if general_settings['allow_tapper'] else 1.0,
#                         'robust': 1.5 if receiver in ('Q') else 1.0,
#                         'solint': '240s' if instrument == 'eM' else '120s',
#                         'sigma_mask': 8.0 if instrument == 'eM' else (12.0 if receiver in ('Q') else 15.0),
#                         'mask_grow_iterations': 4,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'T' if instrument == 'eM' else 'G',
#                         'calmode': 'ap',
#                         'minsnr': 0.1 if instrument == 'eM' else 0.1,
#                         'spwmap': [],
#                         'nsigma_automask': '3.0' if receiver in ('Q') else '4.0',
#                         'nsigma_autothreshold': '2.0',
#                         'uvtaper' : [taper_size] if general_settings['allow_tapper'] else [''],
#                         # 'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '6',
#                         'compare_solints' : False},
#                  }



params_global = {'name': 'global', #global for everything; first loop / faint sources or with diffuse/complex emission
                 'p0': {
                        'robust': 0.25 if receiver in ('Q') else  0.0,
                        'solint': '120s',
                        'sigma_mask': 18 if receiver in ('Q') else 40.0,
                        'mask_grow_iterations': 3,
                        'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
                        'gaintype': 'T' if instrument == 'eM' else 'G',
                        'calmode': 'p',
                        'minsnr': 0.1 if instrument == 'eM' else 0.1,
                        'spwmap': [],
                        'nsigma_automask': '4.0' if receiver in ('Q') else '6.0',
                        'nsigma_autothreshold': '2.0' if receiver in ('Q') else '2.5',
                        'uvtaper' : [''],
                        'with_multiscale': True,
                        'scales': 'None',
                        'maxmscales': '3',
                        'compare_solints' : False},
                 'p1': {
                        'robust': 0.5 if receiver in ('Q') else 0.25,
                        'solint': '90s',
                        'sigma_mask': 14 if receiver in ('Q') else 30.0,
                        'mask_grow_iterations': 3,
                        'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
                        'gaintype': 'G',
                        'calmode': 'p',
                        'minsnr': 0.1 if instrument == 'eM' else 0.1,
                        'spwmap': [],
                        'nsigma_automask': '4.0' if receiver in ('Q') else '5.0',
                        'nsigma_autothreshold': '1.5' if receiver in ('Q') else '2.0',
                        'uvtaper' : [''],
                        'with_multiscale': True,
                        'scales': 'None',
                        'maxmscales': '4',
                        'compare_solints' : False},
                 'p2': {
                        'robust': 0.25 if general_settings['allow_tapper'] else (1.0 if receiver in ('Q') else 0.75),
                        'solint': '60s',
                        'sigma_mask': 12 if receiver in ('Q') else 15.0,
                        'mask_grow_iterations': 4,
                        'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
                        'gaintype': 'T',
                        'calmode': 'p',
                        'minsnr': 0.1 if instrument == 'eM' else 0.1,
                        'spwmap': [],
                        'nsigma_automask': '3.0' if receiver in ('Q') else '4.0',
                        'nsigma_autothreshold': '1.5' if receiver in ('Q') else '2.0',
                        'uvtaper' : [taper_size] if general_settings['allow_tapper'] else [''],
                        # 'uvtaper' : [''],
                        'with_multiscale': True,
                        'scales': 'None',
                        'maxmscales': '5',
                        'compare_solints' : False},
                 'ap1': {
                        # 'robust': 0.75 if general_settings['allow_tapper'] else 1.0,
                        'robust': 0.5 if general_settings['allow_tapper'] else (1.5 if receiver in ('Q') else 0.75),
                        'solint': '90s',
                        'sigma_mask': 10 if receiver in ('Q') else 12.0,
                        'mask_grow_iterations': 4,
                        'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
                        'gaintype': 'G',
                        'calmode': 'ap',
                        'minsnr': 0.1 if instrument == 'eM' else 0.1,
                        'spwmap': [],
                        'nsigma_automask': '3.0' if receiver in ('Q') else '4.0',
                        'nsigma_autothreshold': '1.5' if receiver in ('Q') else '2.0',
                        'uvtaper' : [taper_size] if general_settings['allow_tapper'] else [''],
                        # 'uvtaper' : [''],
                        'with_multiscale': True,
                        'scales': 'None',
                        'maxmscales': '6',
                        'compare_solints' : False},
                 }


# params_global = {'name': 'global', #global for everything; initial loop / very bright sources / lots of initial artefacts
#                  'p0': {
#                         'robust': -2.0,
#                         'solint': '60s',
#                         'sigma_mask': 150.0,
#                         'mask_grow_iterations': 1,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'T' if instrument == 'eM' else 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 0.5,
#                         'spwmap': [],
#                         'nsigma_automask': '8.0',
#                         'nsigma_autothreshold': '4.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': False,
#                         'scales': 'None',
#                         'maxmscales': '3',
#                         'compare_solints' : False},
#                  'p1': {
#                         'robust': -1.0,
#                         'solint': '60s',
#                         'sigma_mask': 120.0,
#                         'mask_grow_iterations': 3,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 0.5,
#                         'spwmap': [],
#                         'nsigma_automask': '8.0',
#                         'nsigma_autothreshold': '3.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '4',
#                         'compare_solints' : False},
#                  'p2': {
#                         'robust': -0.5,
#                         'solint': '30s',
#                         'sigma_mask': 60.0,
#                         'mask_grow_iterations': 4,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'T' if instrument == 'eM' else 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 0.5,
#                         'spwmap': [],
#                         'nsigma_automask': '6.0',
#                         'nsigma_autothreshold': '2.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '5',
#                         'compare_solints' : False},
#                  'ap1': {
#                         'robust': 0.0,
#                         'solint': '30s',
#                         'sigma_mask': 40.0,
#                         'mask_grow_iterations': 4,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'G',
#                         'calmode': 'ap',
#                         'minsnr': 0.1 if instrument == 'eM' else 0.5,
#                         'spwmap': [],
#                         'nsigma_automask': '6.0',
#                         'nsigma_autothreshold': '2.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '6',
#                         'compare_solints' : False},
#                  }

# params_global = {'name': 'global', #global for everything; initial loop / same as before, but for high-frequency
#                  'p0': {
#                         'robust': -2.0,
#                         'solint': '120s',
#                         'sigma_mask': 80.0,
#                         'mask_grow_iterations': 1,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'T' if instrument == 'eM' else 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 0.5,
#                         'spwmap': [],
#                         'nsigma_automask': '6.0',
#                         'nsigma_autothreshold': '3.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': False,
#                         'scales': 'None',
#                         'maxmscales': '3',
#                         'compare_solints' : False},
#                  'p1': {
#                         'robust': -1.0,
#                         'solint': '60s',
#                         'sigma_mask': 50.0,
#                         'mask_grow_iterations': 2,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 0.5,
#                         'spwmap': [],
#                         'nsigma_automask': '6.0',
#                         'nsigma_autothreshold': '3.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '4',
#                         'compare_solints' : False},
#                  'p2': {
#                         'robust': -1.0,
#                         'solint': '30s',
#                         'sigma_mask': 25.0,
#                         'mask_grow_iterations': 4,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'T' if instrument == 'eM' else 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 0.5,
#                         'spwmap': [],
#                         'nsigma_automask': '6.0',
#                         'nsigma_autothreshold': '2.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '5',
#                         'compare_solints' : False},
#                  'ap1': {
#                         'robust': -1.0,
#                         'solint': '30s',
#                         'sigma_mask': 20.0,
#                         'mask_grow_iterations': 4,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'G',
#                         'calmode': 'ap',
#                         'minsnr': 0.1 if instrument == 'eM' else 0.5,
#                         'spwmap': [],
#                         'nsigma_automask': '6.0',
#                         'nsigma_autothreshold': '2.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '6',
#                         'compare_solints' : False},
#                  }



# params_global = {'name': 'global', #global for everything; second loop / very bright sources / moderate initial artefacts leftovers from a previous run
#                  'p0': {
#                         'robust': 0.0,
#                         'solint': '30s',
#                         'sigma_mask': 80.0,
#                         'mask_grow_iterations': 2,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'T' if instrument == 'eM' else 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 0.5,
#                         'spwmap': [],
#                         'nsigma_automask': '8.0',
#                         'nsigma_autothreshold': '3.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '4',
#                         'compare_solints' : False},
#                  'p1': {
#                         'robust': 0.0,
#                         'solint': '18s',
#                         'sigma_mask': 40.0,
#                         'mask_grow_iterations': 3,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 0.5,
#                         'spwmap': [],
#                         'nsigma_automask': '8.0',
#                         'nsigma_autothreshold': '3.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '6',
#                         'compare_solints' : False},
#                  'p2': {
#                         'robust': 0.5,
#                         'solint': 'int',
#                         'sigma_mask': 30.0,
#                         'mask_grow_iterations': 4,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'T' if instrument == 'eM' else 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 0.5,
#                         'spwmap': [],
#                         'nsigma_automask': '8.0',
#                         'nsigma_autothreshold': '2.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '8',
#                         'compare_solints' : False},
#                  'ap1': {
#                         'robust': 0.5,
#                         'solint': '60s',
#                         'sigma_mask': 20.0,
#                         'mask_grow_iterations': 4,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'G',
#                         'calmode': 'ap',
#                         'minsnr': 0.1 if instrument == 'eM' else 0.5,
#                         'spwmap': [],
#                         'nsigma_automask': '8.0',
#                         'nsigma_autothreshold': '2.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '8',
#                         'compare_solints' : False},
#                  }



# params_global = {'name': 'global', #global for everything; first loop / moderate bright sources / moderate initial artefacts / significant extended emission
#                  'p0': {
#                         'robust': -0.5,
#                         'solint': '120s',
#                         'sigma_mask': 80.0,
#                         'mask_grow_iterations': 1,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'T' if instrument == 'eM' else 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 0.1,
#                         'spwmap': [],
#                         'nsigma_automask': '8.0',
#                         'nsigma_autothreshold': '3.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '4',
#                         'compare_solints' : False},
#                  'p1': {
#                         'robust': -0.5,
#                         'solint': '60s',
#                         'sigma_mask': 40.0,
#                         'mask_grow_iterations': 3,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 0.1,
#                         'spwmap': [],
#                         'nsigma_automask': '8.0',
#                         'nsigma_autothreshold': '3.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '4',
#                         'compare_solints' : False},
#                  'p2': {
#                         'robust': 0.5,
#                         'solint': '30s',
#                         'sigma_mask': 20.0,
#                         'mask_grow_iterations': 4,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'T' if instrument == 'eM' else 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 0.1,
#                         'spwmap': [],
#                         'nsigma_automask': '6.0',
#                         'nsigma_autothreshold': '2.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '6',
#                         'compare_solints' : False},
#                  'ap1': {
#                         'robust': 1.0,
#                         'solint': '30s',
#                         'sigma_mask': 15.0,
#                         'mask_grow_iterations': 4,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'G',
#                         'calmode': 'ap',
#                         'minsnr': 0.1 if instrument == 'eM' else 0.1,
#                         'spwmap': [],
#                         'nsigma_automask': '4.0',
#                         'nsigma_autothreshold': '2.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '6',
#                         'compare_solints' : False},
#                  }


# params_global = {'name': 'global', #global for everything; second loop / moderate bright sources / significant extended emission
#                  'p0': {
#                         'robust': 0.0,
#                         'solint': '60s',
#                         'sigma_mask': 30.0,
#                         'mask_grow_iterations': 1,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'T' if instrument == 'eM' else 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 0.1,
#                         'spwmap': [],
#                         'nsigma_automask': '8.0',
#                         'nsigma_autothreshold': '3.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '4',
#                         'compare_solints' : False},
#                  'p1': {
#                         'robust': 0.5,
#                         'solint': '60s',
#                         'sigma_mask': 20.0,
#                         'mask_grow_iterations': 3,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 0.1,
#                         'spwmap': [],
#                         'nsigma_automask': '8.0',
#                         'nsigma_autothreshold': '3.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '4',
#                         'compare_solints' : False},
#                  'p2': {
#                         'robust': 1.5,
#                         'solint': '18s',
#                         'sigma_mask': 15.0,
#                         'mask_grow_iterations': 4,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'T' if instrument == 'eM' else 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 0.1,
#                         'spwmap': [],
#                         'nsigma_automask': '6.0',
#                         'nsigma_autothreshold': '2.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '6',
#                         'compare_solints' : False},
#                  'ap1': {
#                         'robust': 1.5,
#                         'solint': '18s',
#                         'sigma_mask': 12.0,
#                         'mask_grow_iterations': 4,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'G',
#                         'calmode': 'ap',
#                         'minsnr': 0.1 if instrument == 'eM' else 0.1,
#                         'spwmap': [],
#                         'nsigma_automask': '4.0',
#                         'nsigma_autothreshold': '2.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '6',
#                         'compare_solints' : False},
#                  }

# params_global = {'name': 'global', #global for everything - working with M82; first loop / e-MERLIN
#                  'p0': {
#                         'robust': -0.5,
#                         'solint': '240s',
#                         'sigma_mask': 80.0,
#                         'mask_grow_iterations': 1,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'T' if instrument == 'eM' else 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 0.5,
#                         'spwmap': [],
#                         'nsigma_automask': '6.0',
#                         'nsigma_autothreshold': '3.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '3',
#                         'compare_solints' : False},
#                  'p1': {
#                         'robust': 0.0,
#                         'solint': '120s',
#                         'sigma_mask': 40.0,
#                         'mask_grow_iterations': 3,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 0.5,
#                         'spwmap': [],
#                         'nsigma_automask': '4.0',
#                         'nsigma_autothreshold': '1.5',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '4',
#                         'compare_solints' : False},
#                  'p2': {
#                         'robust': 0.0,
#                         'solint': '40s' if instrument == 'eM' else 'int',
#                         'sigma_mask': 20.0,
#                         'mask_grow_iterations': 4,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'T' if instrument == 'eM' else 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 0.5,
#                         'spwmap': [],
#                         'nsigma_automask': '3.0',
#                         'nsigma_autothreshold': '1.5',
#                         # 'uvtaper' : [''],
#                         'uvtaper' : [taper_size] if general_settings['allow_tapper'] else [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '5',
#                         'compare_solints' : False},
#                  'ap1': {
#                         'robust': 0.5,
#                         'solint': '80s' if instrument == 'eM' else 'int',
#                         'sigma_mask': 12.0,
#                         'mask_grow_iterations': 4,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'G',
#                         'calmode': 'ap',
#                         'minsnr': 0.1 if instrument == 'eM' else 0.5,
#                         'spwmap': [],
#                         'nsigma_automask': '3.0',
#                         'nsigma_autothreshold': '1.5',
#                         # 'uvtaper' : [''],
#                         'uvtaper' : [taper_size] if general_settings['allow_tapper'] else [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '6',
#                         'compare_solints' : False},
#                  }

# params_global = {'name': 'global', #global for everything; tests with moderately bright and complex eM emission
#                  'p0': {
#                         'robust': -0.5,
#                         'solint': '120s',
#                         'sigma_mask':40.0,
#                         'mask_grow_iterations': 4,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'T' if instrument == 'eM' else 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 0.1,
#                         'spwmap': [],
#                         'nsigma_automask': '6.0',
#                         'nsigma_autothreshold': '3.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '3',
#                         'compare_solints' : False},
#                  'p1': {
#                         'robust': -0.25,
#                         'solint': '80s',
#                         'sigma_mask': 25.0,
#                         'mask_grow_iterations': 4,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 0.1,
#                         'spwmap': [],
#                         'nsigma_automask': '6.0',
#                         'nsigma_autothreshold': '2.5',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '4',
#                         'compare_solints' : False},
#                  'p2': {
#                         'robust': 0.25,
#                         'solint': '40s',
#                         'sigma_mask': 15.0,
#                         'mask_grow_iterations': 4,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'T',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 0.1,
#                         'spwmap': [],
#                         'nsigma_automask': '4.0',
#                         'nsigma_autothreshold': '2.0',
#                         'uvtaper' : [taper_size] if general_settings['allow_tapper'] else [''],
#                         # 'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '6',
#                         'compare_solints' : False},
#                  'ap1': {
#                         # 'robust': 0.75 if general_settings['allow_tapper'] else 1.0,
#                         'robust': 0.75,
#                         'solint': '80s',
#                         'sigma_mask': 12.0,
#                         'mask_grow_iterations': 6,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'G',
#                         'calmode': 'ap',
#                         'minsnr': 0.1 if instrument == 'eM' else 0.1,
#                         'spwmap': [],
#                         'nsigma_automask': '4.0',
#                         'nsigma_autothreshold': '2.0',
#                         'uvtaper' : [taper_size] if general_settings['allow_tapper'] else [''],
#                         # 'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '6',
#                         'compare_solints' : False},
#                  }


# params_global = {'name': 'global', #global for everything; tests with moderately bright and complex eM emission; second sc run
#                  'p0': {
#                         'robust': 0.0,
#                         'solint': '240s',
#                         'sigma_mask':40.0,
#                         'mask_grow_iterations': 4,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'G' if instrument == 'eM' else 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 0.1,
#                         'spwmap': [],
#                         'nsigma_automask': '5.0',
#                         'nsigma_autothreshold': '2.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '3',
#                         'compare_solints' : False},
#                  'p1': {
#                         'robust': 0.5,
#                         'solint': '120s',
#                         'sigma_mask': 25.0,
#                         'mask_grow_iterations': 4,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 0.1,
#                         'spwmap': [],
#                         'nsigma_automask': '5.0',
#                         'nsigma_autothreshold': '2.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '4',
#                         'compare_solints' : False},
#                  'p2': {
#                         'robust': 0.75,
#                         'solint': '80s',
#                         'sigma_mask': 15.0,
#                         'mask_grow_iterations': 4,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 0.1,
#                         'spwmap': [],
#                         'nsigma_automask': '4.0',
#                         'nsigma_autothreshold': '1.5',
#                         'uvtaper' : [taper_size] if general_settings['allow_tapper'] else [''],
#                         # 'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '6',
#                         'compare_solints' : False},
#                  'ap1': {
#                         # 'robust': 0.75 if general_settings['allow_tapper'] else 1.0,
#                         'robust': 1.0,
#                         'solint': '240s',
#                         'sigma_mask': 12.0,
#                         'mask_grow_iterations': 6,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'G',
#                         'calmode': 'ap',
#                         'minsnr': 0.1 if instrument == 'eM' else 0.1,
#                         'spwmap': [],
#                         'nsigma_automask': '4.0',
#                         'nsigma_autothreshold': '1.5',
#                         'uvtaper' : [taper_size] if general_settings['allow_tapper'] else [''],
#                         # 'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '6',
#                         'compare_solints' : False},
#                  }



# params_very_faint = {'name': 'very_faint', #very faint global for everything; combined selfcal for faint sources or with diffuse/complex emission
#                  'p0': {
#                         'robust': 0.75,
#                         'solint': '120s' if instrument == 'eM' else '120s',
#                         'sigma_mask': 12.0 if instrument == 'eM' else (12.0 if receiver in ('Q') else 15.0),
#                         'mask_grow_iterations': 3,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'T' if instrument == 'eM' else 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 0.1,
#                         'spwmap': [],
#                         'nsigma_automask': '3.0' if receiver in ('Q') else '4.0',
#                         'nsigma_autothreshold': '2.0',
#                         'uvtaper' : [taper_size] if general_settings['allow_tapper'] else [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '6',
#                         'compare_solints' : False},
#                  'ap1': {
#                         'robust': 0.75,
#                         'solint': '120s' if instrument == 'eM' else '120s',
#                         'sigma_mask': 12.0 if instrument == 'eM' else (12.0 if receiver in ('Q') else 15.0),
#                         'mask_grow_iterations': 4,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'T' if instrument == 'eM' else 'G',
#                         'calmode': 'ap',
#                         'minsnr': 0.1 if instrument == 'eM' else 0.1,
#                         'spwmap': [],
#                         'nsigma_automask': '3.0' if receiver in ('Q') else '4.0',
#                         'nsigma_autothreshold': '2.0',
#                         'uvtaper' : [taper_size] if general_settings['allow_tapper'] else [''],
#                         # 'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '6',
#                         'compare_solints' : False},
#                  }


# params_global = {'name': 'global', #global for everything; combined selfcal for faint sources or with diffuse/complex emission
#                  'p0': {
#                         'robust': 0.75,
#                         'solint': '120s',
#                         'sigma_mask': 12.0 if instrument == 'eM' else (12.0 if receiver in ('Q') else 15.0),
#                         'mask_grow_iterations': 3,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 0.1,
#                         'spwmap': [],
#                         'nsigma_automask': '4.0',
#                         'nsigma_autothreshold': '2.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '6',
#                         'compare_solints' : False},
#                  'p1': {
#                         'robust': 0.5,
#                         'solint': '90s',
#                         'sigma_mask': 14 if receiver in ('Q') else 30.0,
#                         'mask_grow_iterations': 3,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 0.1,
#                         'spwmap': [],
#                         'nsigma_automask': '4.0',
#                         'nsigma_autothreshold': '2.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '4',
#                         'compare_solints' : False},
#                  'p2': {
#                         'robust': 0.5,
#                         'solint': '60s',
#                         'sigma_mask': 12 if receiver in ('Q') else 15.0,
#                         'mask_grow_iterations': 4,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'T',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 0.1,
#                         'spwmap': [],
#                         'nsigma_automask': '4.0',
#                         'nsigma_autothreshold': '2.0',
#                         'uvtaper' : [taper_size] if general_settings['allow_tapper'] else [''],
#                         # 'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '5',
#                         'compare_solints' : False},
#                  'ap1': {
#                         # 'robust': 0.75 if general_settings['allow_tapper'] else 1.0,
#                         'robust': 0.75,
#                         'solint': '120s',
#                         'sigma_mask': 12.0 if instrument == 'eM' else (12.0 if receiver in ('Q') else 15.0),
#                         'mask_grow_iterations': 3,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'G',
#                         'calmode': 'ap',
#                         'minsnr': 0.1 if instrument == 'eM' else 0.1,
#                         'spwmap': [],
#                         'nsigma_automask': '4.0',
#                         'nsigma_autothreshold': '2.0',
#                         'uvtaper' : [taper_size] if general_settings['allow_tapper'] else [''],
#                         # 'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '6',
#                         'compare_solints' : False},
#                  }


# params_global = {'name': 'global', #global for everything; combined selfcal for moderate bright sources or with diffuse/complex emission
#                  'p0': {
#                         'robust': -0.5 if instrument == 'eM' else 0.0,
#                         'solint': '120s',
#                         'sigma_mask': 12.0 if instrument == 'eM' else (12.0 if receiver in ('Q') else 15.0),
#                         'mask_grow_iterations': 3,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 0.1,
#                         'spwmap': [],
#                         'nsigma_automask': '4.0',
#                         'nsigma_autothreshold': '2.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '6',
#                         'compare_solints' : False},
#                  'p1': {
#                         'robust': -0.5 if instrument == 'eM' else 0.0,
#                         'solint': '90s',
#                         'sigma_mask': 14 if receiver in ('Q') else 30.0,
#                         'mask_grow_iterations': 3,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'G',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 0.1,
#                         'spwmap': [],
#                         'nsigma_automask': '4.0',
#                         'nsigma_autothreshold': '2.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '4',
#                         'compare_solints' : False},
#                  'p2': {
#                         'robust': -0.5 if instrument == 'eM' else 0.5,
#                         'solint': '60s',
#                         'sigma_mask': 12 if receiver in ('Q') else 15.0,
#                         'mask_grow_iterations': 4,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'T',
#                         'calmode': 'p',
#                         'minsnr': 0.1 if instrument == 'eM' else 0.1,
#                         'spwmap': [],
#                         'nsigma_automask': '4.0',
#                         'nsigma_autothreshold': '2.0',
#                         'uvtaper' : [taper_size] if general_settings['allow_tapper'] else [''],
#                         # 'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '5',
#                         'compare_solints' : False},
#                  'ap1': {
#                         # 'robust': 0.75 if general_settings['allow_tapper'] else 1.0,
#                         'robust': -0.5 if instrument == 'eM' else 0.5,
#                         'solint': '120s',
#                         'sigma_mask': 12.0 if instrument == 'eM' else (12.0 if receiver in ('Q') else 15.0),
#                         'mask_grow_iterations': 3,
#                         'combine': 'scan,spw' if general_settings['allow_combine_spw'] else 'scan',
#                         'gaintype': 'G',
#                         'calmode': 'ap',
#                         'minsnr': 0.1 if instrument == 'eM' else 0.1,
#                         'spwmap': [],
#                         'nsigma_automask': '4.0',
#                         'nsigma_autothreshold': '2.0',
#                         'uvtaper' : [taper_size] if general_settings['allow_tapper'] else [''],
#                         # 'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': 'None',
#                         'maxmscales': '6',
#                         'compare_solints' : False},
#                  }

params_very_faint = params_very_faint.copy()
params_very_faint['name'] = 'very_faint'
params_faint = params_global.copy()
params_faint['name'] = 'faint'
params_standard_1 = params_global.copy()
params_standard_1['name'] = 'standard_1'
params_standard_2 = params_global.copy()
params_standard_2['name'] = 'standard_2'
params_bright = params_global.copy()
params_bright['name'] = 'bright'


params_trial_2 = None # comment this and uncomment the following lines
                      # if this is the second pass of self-calibration.


# params_trial_2 = {'name': 'trial_2',
#                  'p0': {'robust': 0.0,
#                         'solint' : '36s',
#                         'sigma_mask': 12,
#                         'combine': '',
#                         'gaintype': 'G',
#                         'calmode': 'p',
#                         'minsnr': 3.0,
#                         'spwmap': [],
#                         'nsigma_automask': '6.0',
#                         'nsigma_autothreshold': '3.0',
#                         'uvtaper' : [''],
#                         'with_multiscale' : True,
#                         'scales': '0,5,20,50',
#                         'compare_solints' : False},
#                  'p1': {'robust': 0.5,
#                         'solint' : '12s',
#                         'sigma_mask': 8,#set to 15 if e-MERLIN
#                         'combine': '',
#                         'gaintype': 'G',
#                         'calmode': 'p',
#                         'minsnr': 3.0,
#                         'spwmap': [],
#                         'nsigma_automask': '6.0',
#                         'nsigma_autothreshold': '3.0',
#                         'uvtaper' : [''],
#                         'with_multiscale': True,
#                         'scales': '0,5,20,50',
#                         'compare_solints': False},
#                  'ap1': {'robust': 0.5,
#                          'solint': '60s',
#                          'sigma_mask': 8,
#                          'combine': '',
#                          'gaintype': 'G',
#                          'calmode': 'ap',
#                          'minsnr': 3.0,
#                          'spwmap': [],
#                          'uvtaper' : [''],
#                          'nsigma_automask': '3.0',
#                          'nsigma_autothreshold': '1.5',
#                          'with_multiscale': True,
#                          'scales': '0,5,20,50',
#                          'compare_solints': False},
#                  }