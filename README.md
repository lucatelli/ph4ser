# ph4ser
The module `ph4ser` is a self-calibration pipeline that uses `WSClean` as imager. 
This is an improved version of the self-calibration module contained in the `morphen_dev` repository.
`ph4ser` is still in development, so expect bugs and issues. Feel free to send any feedback and idea
that you may have. 


## Install dependencies
To install the dependencies, just follow the instructions found in 
[https://github.com/lucatelli/morphen](https://github.com/lucatelli/morphen). 

Then, clone the `ph4ser` repository. 
```
git clone https://github.com/lucatelli/ph4ser.git
```


## Running the code

The module needs the configuration file `ph4ser_config.py`. In there, 
you can manually input the information needed for `ph4ser`. 

### Set visibility name/path.
The first step is to set the path to the visibility data:
```
visibility_info = {'path':"/data/C_band/MCG+12-02-001/autoselfcal/",
                   'field':'MCG12-02-001',
                   'vis_name':'MCG12-02-001.calibrated.avg',
                   'savename':'_A_C_sf'
       }
```
1. Copy the splited calibrated visibility into `path`. 
2. Provide the visibility filename to `vis_name` (without the `.ms` extensions)
3. `savename` is just a convenient name so that at the end self-calibrated files
will be saved with that prefix.


### Set additional info
The second step is to tell `ph4ser` the data's instrument/receiver information.
```
receiver = 'C'
instrument = 'EVLA' # 'EVLA' or 'eM'
```

Complementary information, such as the image size can be specified within the dictionary
`init_parameters`. To set the image size in x and y direction to 2048 pixels, you can use 
the keys `imsize` and `imsizey`, respectively. For example:
```
init_parameters = {'fov_image': {'imsize': 1024*8,
                                'cell': '0.4arcsec',
                                'basename': 'FOV_phasecal_image',
                                'FIELD_SHIFT': None,
                                'niter': 1000,
                                'robust': 0.5},
                  'test_image': {'imsize': int(1024*2),
                                 'imsizey': int(1024*2),
                                 'FIELD_SHIFT': None,
                                 'cell': cell_size,
                                 'prefix': 'test_image',
                                 'uvtaper': [''],
                                 'niter': 10000,
                                 'robust': 0.5 if receiver in ('K', 'Ka') or instrument == 'eM' else 0.0}
                  }
```

### Selecting steps that will be run.
The number of self-calibration steps that are going to be run depends on the collective 
flux density of all sources within the region to be imaged. If the total flux density of 
all sources is high, more steps of self-calibration will be performed 
(e.g. p -> p -> p -> ap). If the value is low, we must avoid 
multiple self-calibration steps, so usually a single phase step will be 
executed (and an amplitude will be attempted).

The general configuration is to set all steps, and the exact number will be determined 
during runtime, when a test image will be generated.

The default configuration contained within `ph4ser_config` is:
```
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
```

### Set the path to the WSClean singularity container. 
The imaging routines with be run with the help of the script `imaging_with_wsclean.py`. The only change in 
that file is to set where your singularity container with `WSClean` is saved. Check line `219`, just edit the 
variable `wsclean_dir`. In my case, I have:
```
    if running_container == 'singularity':
        mount_dir = root_dir_sys + ':/mnt'
        root_dir = '/mnt/'
        wsclean_dir = '/media/sagauga/xfs_evo/wsclean3.4-idg-everybeam-eMERLIN_portable.sif'
```


### Running the code


The code `ph4ser` can be run on a terminal or on a Jupyter notebook. If run in a 
notebook, additional plots generated during runtime will stay there. This is useful 
to easily check how self-calibration improves the visibilities. 

If you prefer to run `ph4ser` on the terminal, after making the changes to the configuration file, you can 
simply do 
```
python ph4ser.py
```
Remember to activate your conda environment in order to run `ph4ser`. 


An alternative, is to run the entire pipeline within a Jupyter notebook. See below. 


### Example notebook. 
The Jupyter notebook [examples/ph4ser_example.ipynb](examples/ph4ser_example.ipynb) 
is an example containing the required set of commands to run self-calibration on a visibility. 






