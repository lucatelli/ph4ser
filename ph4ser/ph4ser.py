"""
                                                          ..___|**_
                                                  .|||||||||*+@+*__*++.
                                              _||||.           .*+;].,#_
          Morphen                        _|||*_                _    .@@@#@.
          ph4ser                   _|||||_               .@##@#| _||_
   Radio Self-Calibration     |****_                   .@.,/\..@_.
          Module             #///#+++*|    .       .@@@;#.,.\@.
                              .||__|**|||||*||*+@#];_.  ;,;_
     Geferson Lucatelli                        +\*_.__|**#
                                              |..      .]]
                                               ;@       @.*.
                                                #|       _;]];|.
                                                 ]_          _+;]@.
                                                 _/_             |]\|    .  _
                                              ...._@* __ .....     ]]+ ..   _
                                                  .. .       . .. .|.|_ ..


This module consists of performing interferometric imaging with wsclean and
running CASA's task gaincal for self-calibration.
It was tested for VLA (L,S,C,X,Ku) and eMERLIN (C band) observations.

Faint sources or higher-frequency observations (e.g. < 10 mJy)
may not work well. So, more experiments are required for
K and Ka VLA bands and eMERLIN fainter sources.

The user is advised to run the code in an interactive session (ipython),
step-by-step, and check the results of each step.
Check the config.py file at:
https://github.com/lucatelli/morphen/blob/main/selfcal/config.py

Note that the pure automated self-calibration is still experimental,
but showed to be good in most cases.

Check https://github.com/lucatelli/morphen/blob/main/selfcal/README.md for more information.

"""
__versions__ = ('0.3.1alpha-1', '0.4.0alpha-1', '0.5.0alpha-1')
__codenames__ = ('Pelicoto', 'Saurinho', '')
__dates__ =  ('2024 03 25','2024 11 13','2024 12 18')
__version__ = '0.4.0alpha-1'
__codename__ = 'Saurinho'
__author__ = 'Geferson Lucatelli'
__coauthors__ = ('')
__email__ = 'geferson.lucatelli@postgrad.manchester.ac.uk'
__date__ = '2024 12 18'
# print(__doc__)

import os
import sys

sys.path.append('./')
import mlibs as mlibs
import glob
import pandas as pd
from casatasks import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit

try:
    # import casatools
    # from casatasks import *
    import casatasks, casatools, casaplotms
    from casaplotms import plotms
    from casatasks import gaincal, phaseshift, applycal, bandpass, split, importasdm, statwt
    from casatasks import flagdata, flagmanager, plotants, plotweather, hanningsmooth
    from casatasks import clearcal, delmod, setjy, fluxscale, tclean, imhead, listobs
except:
    print('Not importing casatools. '
          'Maybe you are inside inside CASA? '
          'Or check your modular installation?')
    pass

from casaplotms import plotms
from casaviewer.imview import imview


def reload_libs():
    """
    This allows you to reload the library file without restarting the kernel.

    Usage:
    In [1]: reload_libs()
    """
    import mlibs as mlibs
    from importlib import reload
    reload(mlibs)


# reload_libs()


msmd = casatools.msmetadata()
ms = casatools.ms()
tb = casatools.table()

# mlibs.reset_rc_params()
import matplotlib as mpl
import os
import tableprint
from importlib import reload


import ph4ser_config as cf
# import ph4ser_config_combined as cf



class Configuration:
    """
    Configuration Class to specify basic parameters for data processing and visualization.
    """

    def __init__(self):
        """
        Initialize configuration with default values and load external settings.
        """
        # Load external configuration
        reload(cf)

        # Visibility and processing parameters
        self.field = cf.FIELD
        self.antennas = cf.ANTENNAS
        self.refantmode = cf.refantmode
        self.spws = cf.SPWS
        self.minblperant = cf.minblperant

        # Instrument-specific parameters
        self.cell_sizes_jvla = cf.cell_sizes_JVLA
        self.cell_sizes_emerlin = cf.cell_sizes_eMERLIN
        self.taper_sizes_emerlin = cf.taper_sizes_eMERLIN
        self.taper_sizes_jvla = cf.taper_sizes_JVLA




        # Visibility information
        # self.visibility_info = cf.run_mode
        self.visibility_info = cf.visibility_info
        self.path = cf.visibility_info['path']
        self.field = cf.visibility_info['field']
        self.vis_name = cf.visibility_info['vis_name']
        self.savename = cf.visibility_info['savename']

        # Processing steps and settings
        self.steps = cf.steps
        self.refant = cf.refant
        self.instrument = cf.instrument

        # Parameter sets for different scenarios
        self.init_parameters = cf.init_parameters
        self.global_parameters = cf.global_parameters
        self.general_settings = cf.general_settings
        self.params_very_faint = cf.params_very_faint
        self.params_faint = cf.params_faint
        self.params_standard_1 = cf.params_standard_1
        self.params_standard_2 = cf.params_standard_2
        self.params_bright = cf.params_bright
        self.params_trial_2 = cf.params_trial_2
        self.do_additional_images = cf.do_additional_images

        # Imaging parameters
        self.receiver = cf.receiver
        self.cell_size = cf.cell_size
        self.taper_size = cf.taper_size
        self.nc = cf.nc
        self.negative_arg = cf.negative_arg
        self.solnorm = cf.solnorm
        self.quiet = cf.quiet
        self.run_mode = cf.run_mode
        self.imsize = cf.global_parameters['imsize']
        self.imsizey = cf.global_parameters['imsizey']

        # Set default matplotlib parameters
        self.reset_rc_params()

    @staticmethod
    def reset_rc_params():
        """
        Set global configuration for matplotlib.pyplot
        """
        mpl.rcParams.update({
            'font.size': 16,
            'text.usetex': False,
            'font.family': 'STIXGeneral',
            'mathtext.fontset': 'stix',
            'font.weight': 'medium',
            'xtick.labelsize': 16,
            'ytick.labelsize': 16,
            'figure.figsize': (6, 4),
            'axes.labelsize': 16,
            'xtick.major.width': 1,
            'ytick.major.width': 1,
            'axes.linewidth': 1.5,
            'axes.edgecolor': 'orange',
            'lines.linewidth': 2,
            'legend.fontsize': 14,
            'grid.linestyle': '--',
            'axes.grid.which': 'major',
            'axes.grid.axis': 'both',
            'axes.spines.right': True,
            'axes.grid': True,
            'axes.titlesize': 16,
            'legend.framealpha': 1.0
        })

    def update_config(self, **kwargs):
        """
        Update configuration parameters dynamically.

        Args:
            **kwargs: Key-value pairs of configuration parameters to update
        """
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                raise AttributeError(f"Configuration has no attribute '{key}'")


class Pipeline:
    def __init__(self, configuration=None):
        """Initialize pipeline with configuration"""
        self.config = configuration if configuration else Configuration()
        self.steps_performed = []
        self.initialize_storage()
        self.freq_ranges = {
            (1, 2): "L",
            (2, 4): "S",
            (4, 8): "C",
            (8, 12): "X",
            (12, 18): "Ku",
            (18, 26.5): "K",
            (26.5, 40): "Ka",
            (40, 50): "Q"
        }

    def select_parameters(self, total_flux, snr=None):
        """Select appropriate parameters based on total flux"""
        if total_flux < 10:
            params = self.config.params_very_faint.copy()
        elif 10 <= total_flux < 20:
            params = self.config.params_faint.copy()
        elif 20 <= total_flux < 50:
            params = self.config.params_standard_1.copy()
        elif 50 <= total_flux < 100:
            params = self.config.params_standard_2.copy()
        else:
            params = self.config.params_bright.copy()
        return params

    def get_spwids(self, vis):
        """Get spectral window IDs from visibility data"""
        lobs = listobs(vis=vis)
        extract_spwids = {key: lobs[key] for key in lobs if 'scan_' in key}

        unique_spwids = set()
        for key in extract_spwids:
            nested_dict = extract_spwids[key]
            for inner_key in nested_dict:
                spwids = nested_dict[inner_key]['SpwIds']
                unique_spwids.add(tuple(sorted(spwids)))

        unique_spwids_lists = sorted([list(t) for t in unique_spwids])
        unique_elements = set(element for sublist in unique_spwids_lists
                              for element in sublist)
        return sorted(list(unique_elements))

    # def get_spwmap(self, vis):
    #     """Generate spectral window mapping"""
    #     unique_spwids_lists = self._get_unique_spwids(vis)
    #     counts = self._count_spw_occurrences(unique_spwids_lists)
    #     spwmap_i = [item for item, count in counts.items()
    #                for _ in range(count)]
    #     return [spwmap_i[:len(self.get_spwids(vis))]]

    def get_spwmap(self, vis):
        lobs = listobs(vis=vis)
        extract_spwids = {key: lobs[key] for key in lobs if 'scan_' in key}

        unique_spwids = set()

        for key in extract_spwids:
            nested_dict = extract_spwids[key]
            for inner_key in nested_dict:
                spwids = nested_dict[inner_key]['SpwIds']
                # Convert the array to a sorted tuple and add to the set
                unique_spwids.add(tuple(sorted(spwids)))

        # Convert the set of tuples back to a sorted list of lists
        unique_spwids_lists = sorted([list(t) for t in unique_spwids])

        counts = {}
        for lst in unique_spwids_lists:
            if lst[0] not in counts:
                counts[lst[0]] = len(lst)

        # Construct the spwmap
        spwmap_i = [item for item, count in counts.items() for _ in range(count)]
        spwmap = [spwmap_i[:len(self.get_spwids(vis))]]
        return spwmap

    def get_all_chan_freqs(self, vis):
        """
        Retrieve channel frequencies for all spectral windows.

        Args:
            vis: Input visibility file

        Returns:
            tuple: (Flattened channel frequencies array, Channel frequencies per SPW)
        """
        spw_list = self.get_spwids(vis)
        chan_freqs = np.empty(len(spw_list), dtype=object)

        msmd.open(vis)
        for i in range(len(spw_list)):
            chan_freq = msmd.chanfreqs(spw_list[i])
            chan_freqs[i] = chan_freq
        msmd.done()

        chan_freqs_flat = np.hstack(chan_freqs)
        return chan_freqs_flat, chan_freqs

    def get_spwids(self, vis):
        """
        Get unique spectral window IDs from visibility data.

        Args:
            vis: Input visibility file

        Returns:
            numpy.ndarray: Array of unique spectral window IDs
        """
        lobs = listobs(vis=vis)
        extract_spwids = {key: lobs[key] for key in lobs if 'scan_' in key}

        unique_spwids = set()
        for key in extract_spwids:
            nested_dict = extract_spwids[key]
            for inner_key in nested_dict:
                spwids = nested_dict[inner_key]['SpwIds']
                unique_spwids.add(tuple(sorted(spwids)))

        unique_spwids_lists = sorted([list(t) for t in unique_spwids])
        unique_elements = set(element for sublist in unique_spwids_lists
                              for element in sublist)

        return np.asarray(sorted(list(unique_elements))).astype(int)

    def get_cell_size(self, vis):
        """Calculate optimal cell size if not defined in config"""
        if not self.config.cell_size:
            ms.open(vis)
            ms.selectinit(datadescid=0)
            uvw = ms.getdata('uvw')['uvw']
            ms.close()

            uvdist_meters = np.sqrt(uvw[0] ** 2 + uvw[1] ** 2)
            longest_baseline_meters = np.nanmax(uvdist_meters)

            spw_list = self.get_spwids(vis)
            chan_freqs = np.empty(len(spw_list), dtype=object)

            msmd.open(vis)
            for i in range(len(spw_list)):
                chan_freqs[i] = msmd.chanfreqs(spw_list[i])
            msmd.done()

            max_freq = np.max(np.hstack(chan_freqs))
            wavelength_meters = 3e8 / max_freq
            longest_baseline_lambda = longest_baseline_meters / wavelength_meters

            cell_float = (180.0 * 3600 / (np.pi * 5)) * (1.0 / longest_baseline_lambda)
            return f'{cell_float:.4f}arcsec'

        return self.config.cell_size

    def check_band(self, vis):
        """
        Determine the frequency band of the observations and calculate frequency statistics.

        Args:
            vis: Input visibility file

        Returns:
            tuple: (band_name, mean_frequency, max_frequency, min_frequency)
        """
        if not self.config.receiver:
            msmd.open(vis)
            bandwidth = msmd.bandwidths()
            nspw = len(bandwidth)

            chan_freqs_all = np.empty(nspw, dtype=object)
            spws_freq = np.zeros(nspw)

            for nch in range(nspw):
                chan_freqs_all[nch] = msmd.chanfreqs(nch)
                spws_freq[nch] = np.nanmean(chan_freqs_all[nch])

            msmd.done()

            mean_freq = np.nanmean(spws_freq) * 1e-9  # Convert to GHz
            max_freq = np.nanmax(spws_freq) * 1e-9
            min_freq = np.nanmin(spws_freq) * 1e-9

            band_name = None
            for freq_range, band in self.freq_ranges.items():
                if freq_range[0] <= mean_freq <= freq_range[1]:
                    band_name = band
                    break
            return band_name
        return self.config.receiver

    def print_table(self, data):
        """Print dictionary as formatted table"""
        rows = []
        for key, value in data.items():
            if isinstance(value, list):
                value = ', '.join(map(str, value))
            rows.append((key, value))
        headers = ["Parameter", "Value"]
        tableprint.table(rows, headers)

    def compute_image_stats(self, path, image_list, image_statistics, prefix='', sigma=None,
                            selfcal_step=None):
        """
        Compute statistics of cleaned images from a wsclean run at a given self-calibration step.
        """
        if sigma is None:
            if (selfcal_step == 'test_image') or (selfcal_step == 'p0'):
                """
                We must be more conservative when creating masks to compute the total flux density 
                before self-calibration. The image may contain artifacts above the default sigma 
                threshold of 6.0, and may lead to overestimation of the total flux density.
                An alternative sigma is 8. Note that the mask dilation is a very powerful approach and 
                very sensitive to the sigma threshold. A sigma of 8 results large differences in 
                relation to a sigma of 6. 
                """
                sigma = 10.0
            else:
                sigma = 6.0

        file_list = glob.glob(f"{path}*{prefix}*MFS-image.fits")
        file_list.sort(key=os.path.getmtime, reverse=False)

        try:
            image_list[prefix] = file_list[-1]
        except:
            image_list[prefix] = file_list
        image_list[prefix + '_residual'] = image_list[prefix].replace(
            'MFS-image.fits', 'MFS-residual.fits')
        image_list[prefix + '_model'] = image_list[prefix].replace(
            'MFS-image.fits', 'MFS-model.fits')

        if os.path.exists(image_list[prefix].replace('MFS-image.fits', 'MFS-image-pb.fits')):
            has_pb_image = True
            image_to_analyse = image_list[prefix].replace('MFS-image.fits', 'MFS-image-pb.fits')
            residual_to_analyse = image_list[prefix + '_residual'].replace('MFS-residual.fits',
                                                                           'MFS-residual-pb.fits')
        else:
            has_pb_image = False
            image_to_analyse = image_list[prefix]
            residual_to_analyse = image_list[prefix + '_residual']

        try:
            # mask should be created from normal -MFS-image, and not from -MFS-image-pb!!!!
            _, mask_dilated = mlibs.mask_dilation(image_list[prefix],
                                                  sigma=sigma, PLOT=False)
            level_stats = mlibs.level_statistics(image_to_analyse, sigma=sigma,
                                                 mask=mask_dilated)
        except:
            try:
                _, mask_dilated = mlibs.mask_dilation(
                    image_list[prefix].replace('-MFS-image', '-MFS-dirty'),
                    sigma=sigma, PLOT=False)
                level_stats = mlibs.level_statistics(image_to_analyse, sigma=sigma,
                                                     mask=mask_dilated)
            except:
                _, mask_dilated = mlibs.mask_dilation(
                    image_list[prefix].replace('-MFS-image', '-MFS-dirty'),
                    sigma=5.0, PLOT=False)
                level_stats = mlibs.level_statistics(image_to_analyse, sigma=sigma,
                                                     mask=mask_dilated)

        image_stats = mlibs.get_image_statistics(imagename=image_to_analyse,
                                                 dic_data=level_stats,
                                                 sigma_mask=sigma,
                                                 mask=mask_dilated,
                                                 )
        img_props = mlibs.compute_image_properties(image_to_analyse,
                                                   residual_to_analyse,
                                                   results=image_stats,
                                                   sigma_mask=sigma,
                                                   mask=mask_dilated,
                                                   do_fit_ellipse=False,
                                                   show_figure=False)[-1]

        image_statistics[prefix] = img_props

        try:
            centre = mlibs.nd.maximum_position(
                np.nan_to_num(mlibs.ctn(image_to_analyse) * mask_dilated, nan=0))[::-1]
            # centre = (int(image_statistics[prefix]['y0']),int(image_statistics[prefix]['x0']))
            fig = plt.figure(figsize=(16, 8))
            ax0 = fig.add_subplot(1, 2, 2)
            ax0 = mlibs.eimshow(imagename=image_to_analyse,
                                center=centre,
                                projection='offset', plot_colorbar=True,
                                rms=mlibs.mad_std(
                                    np.nan_to_num(mlibs.ctn(residual_to_analyse), nan=0)),
                                # crop=True,box_size=300,
                                figsize=(8, 8), ax=ax0, fig=fig,
                                crop=True, box_size=int(4 * image_statistics[prefix]['C95radii']),
                                # save_name=image_list[prefix].replace('.fits', '_map'),
                                # CM='magma',
                                add_beam=True)
            ax0.set_title(f'Radio Map')
            ax1 = fig.add_subplot(1, 2, 1)
            ax1 = mlibs.eimshow(imagename=residual_to_analyse,
                                center=centre,
                                projection='offset',
                                vmin_factor=-3.0, vmax_factor=0.99,
                                add_contours=False,
                                figsize=(8, 8), ax=ax1, fig=fig,
                                crop=True, box_size=int(4 * image_statistics[prefix]['C95radii']),
                                save_name=image_to_analyse.replace('.fits', '_map'),
                                plot_title=f'Residual Map', plot_colorbar=False,
                                # CM='magma',
                                add_beam=False)

        except:
            print('--==>> Error on plotting radio map.')
            pass
        try:
            if os.path.exists(image_list[prefix].replace('MFS-image.fits', 'MFS-image-pb.fits')):
                sub_band_images = glob.glob(
                    image_list[prefix].replace('-MFS-image.fits', '') + '-????-image-pb.fits')
                sub_band_residuals = glob.glob(
                    image_list[prefix + '_residual'].replace('-MFS-residual.fits',
                                                             '') + '-????-residual-pb.fits')
            else:
                sub_band_images = glob.glob(
                    image_list[prefix].replace('-MFS-image.fits', '') + '-????-image.fits')
                sub_band_residuals = glob.glob(
                    image_list[prefix + '_residual'].replace('-MFS-residual.fits',
                                                             '') + '-????-residual.fits')

            _FLUXES = []
            _FLUXES_err = []
            for i in range(len(sub_band_images)):
                # flux_density, flux_density_err = mlibs.compute_flux_density(sub_band_images[i],
                #                                                       sub_band_residuals[i],
                #                                                       mask=None)

                img_props = mlibs.compute_image_properties(sub_band_images[i],
                                                           sub_band_residuals[i],
                                                           sigma_mask=sigma,
                                                           last_level=1.5,
                                                           #    mask=None,
                                                           mask=mask_dilated,
                                                           verbose=1,
                                                           save_csv=True,
                                                           do_fit_ellipse=False,
                                                           show_figure=False)[-1]

                flux_density, flux_density_err = img_props['total_flux_mask'], img_props[
                    'flux_error_res_3']

                print('Flux density = ', flux_density)
                _FLUXES.append(flux_density)
                _FLUXES_err.append(flux_density_err)
            FLUXES = mlibs.np.asarray(_FLUXES)
            FLUXES_err = mlibs.np.asarray(_FLUXES_err)
            freqlist = mlibs.getfreqs(sub_band_images)

            mini, result_1, param_dict = (
                mlibs.do_fit_spec_RC_linear(freqlist,
                                            FLUXES * 1000,
                                            FLUXES_err * 1000,
                                            basename_save=image_to_analyse,
                                            title_text=r'WSClean Sub-band Images',
                                            plot_errors_shade=True, do_mcmc_fit=True,
                                            verbose=2))

            plt.close()
        except:
            print('--==>> Some error found in the sub-bands images.')
            pass
        return (image_statistics, image_list)

    def create_mask(self, imagename, rms_mask, sigma_mask, mask_grow_iterations, PLOT=False):

        valid_sigma_mask = sigma_mask
        while True:
            mask_valid = mlibs.mask_dilation(imagename,
                                             PLOT=PLOT,
                                             rms=rms_mask,
                                             dilation_size=None,
                                             sigma=valid_sigma_mask,
                                             iterations=mask_grow_iterations)[1]
            if mask_valid.sum() > 0:
                break
            print(' ++>> No mask found with sigma_mask:', valid_sigma_mask)
            print(' ++>> Reducing sigma_mask by 2 until valid mask is found...')
            valid_sigma_mask = valid_sigma_mask - 2.0

            if valid_sigma_mask <= 6:
                print("Reached minimum sigma threshold without finding a valid mask.")
                break

        mask_valid = mlibs.mask_dilation(imagename,
                                         PLOT=PLOT,
                                         rms=rms_mask,
                                         dilation_size=None,
                                         sigma=valid_sigma_mask - 1,
                                         iterations=mask_grow_iterations)[1]

        mask = mask_valid
        mask_wslclean = mask * 1.0  # mask in wsclean is inverted
        mask_name = imagename.replace('.fits', '') + '_mask.fits'
        mlibs.pf.writeto(mask_name, mask_wslclean, overwrite=True)
        return (mask_name)

    def plot_visibilities(self, g_vis, name, with_DATA=True, with_MODEL=False,
                          with_CORRECTED=False, with_RESIDUAL=False):

        if with_DATA == True:
            plotfile = os.path.dirname(g_vis) + '/selfcal/plots/' + name + '_uvwave_amp_data.jpg'
            if not os.path.isfile(plotfile):
                plotms(vis=g_vis, xaxis='UVwave', yaxis='amp', avgantenna=False, avgscan=False,
                       # antenna=self.config.antennas, spw=self.config.spws,
                       coloraxis='baseline',
                       ydatacolumn='data', avgchannel='9999', avgtime='9999',
                       correlation='RR', plotrange=[0, 0, 0, 0],
                       width=2000, height=800, showgui=False, overwrite=True, dpi=1200,
                       highres=True,
                       customsymbol=True, symbolsize=4, symbolshape='diamond',
                       # plotrange=[-1,-1,-1,0.3],
                       plotfile=plotfile)
            else:
                pass

            plotfile = os.path.dirname(g_vis) + '/selfcal/plots/' + name + '_freq_amp_data.jpg'
            if not os.path.isfile(plotfile):
                plotms(vis=g_vis, xaxis='freq', yaxis='amp', avgantenna=False, avgscan=False,
                       # antenna=self.config.antennas, spw=self.config.spws,
                       coloraxis='baseline',
                       ydatacolumn='data',
                       #    avgchannel='8',
                       avgtime='9999',
                       correlation='RR', plotrange=[0, 0, 0, 0],
                       width=2000, height=800, showgui=False, overwrite=True, dpi=1200,
                       highres=True,
                       customsymbol=True, symbolsize=4, symbolshape='diamond',
                       plotfile=plotfile)
            else:
                pass

        if with_CORRECTED == True:
            plotfile = os.path.dirname(
                g_vis) + '/selfcal/plots/' + name + '_uvwave_amp_corrected_div_model.jpg'
            if not os.path.isfile(plotfile):
                plotms(vis=g_vis, xaxis='UVwave', yaxis='amp',
                       # antenna=self.config.antennas, spw=self.config.spws,
                       coloraxis='baseline', avgantenna=False, avgscan=False,
                       ydatacolumn='corrected/model', avgchannel='9999', avgtime='9999',
                       correlation='RR',
                       width=2000, height=800, showgui=False, overwrite=True, dpi=1200,
                       highres=True,
                       customsymbol=True, symbolsize=4, symbolshape='diamond',
                       plotrange=[0, 0, 0, 10],
                       plotfile=plotfile)
            else:
                pass

            plotfile = os.path.dirname(
                g_vis) + '/selfcal/plots/' + name + '_uvwave_amp_corrected.jpg'
            if not os.path.isfile(plotfile):
                plotms(vis=g_vis, xaxis='UVwave', yaxis='amp', avgantenna=False, avgscan=False,
                       # antenna=self.config.antennas, spw=self.config.spws,
                       coloraxis='baseline',
                       # plotrange=[-1,-1,0,0.3],
                       ydatacolumn='corrected', avgchannel='9999', avgtime='9999',
                       correlation='RR', plotrange=[0, 0, 0, 0],
                       width=2000, height=800, showgui=False, overwrite=True, dpi=1200,
                       highres=True,
                       customsymbol=True, symbolsize=4, symbolshape='diamond',
                       plotfile=plotfile)
            else:
                pass

            plotfile = os.path.dirname(g_vis) + '/selfcal/plots/' + name + '_freq_amp_corrected.jpg'
            if not os.path.isfile(plotfile):
                plotms(vis=g_vis, xaxis='freq', yaxis='amp', avgantenna=False, avgscan=False,
                       # antenna=self.config.antennas, spw=self.config.spws,
                       coloraxis='baseline',
                       ydatacolumn='corrected',
                       #    avgchannel='8',
                       avgtime='9999',
                       correlation='RR', plotrange=[0, 0, 0, 0],
                       width=2000, height=800, showgui=False, overwrite=True, dpi=1200,
                       highres=True,
                       customsymbol=True, symbolsize=4, symbolshape='diamond',
                       plotfile=plotfile)
            else:
                pass

        if with_MODEL == True:
            plotfile = os.path.dirname(g_vis) + '/selfcal/plots/' + name + '_uvwave_amp_model.jpg'
            if not os.path.isfile(plotfile):
                plotms(vis=g_vis, xaxis='UVwave', yaxis='amp', avgantenna=False, avgscan=False,
                       # antenna=self.config.antennas, spw=self.config.spws,
                       coloraxis='baseline',
                       ydatacolumn='model', avgchannel='9999', avgtime='9999',
                       correlation='RR', plotrange=[0, 0, 0, 0],
                       width=2000, height=800, showgui=False, overwrite=True, dpi=1200,
                       highres=True,
                       customsymbol=True, symbolsize=4, symbolshape='diamond',
                       plotfile=plotfile)
            else:
                pass

            plotfile = os.path.dirname(g_vis) + '/selfcal/plots/' + name + '_freq_amp_model.jpg'
            if not os.path.isfile(plotfile):
                plotms(vis=g_vis, xaxis='freq', yaxis='amp', avgantenna=False, avgscan=False,
                       # antenna=self.config.antennas, spw=self.config.spws,
                       coloraxis='baseline',
                       ydatacolumn='model',
                       #    avgchannel='',
                       avgtime='9999',
                       correlation='RR', plotrange=[0, 0, 0, 0],
                       width=2000, height=800, showgui=False, overwrite=True, dpi=1200,
                       highres=True,
                       customsymbol=True, symbolsize=4, symbolshape='diamond',
                       plotfile=plotfile)
            else:
                pass

        if with_RESIDUAL == True:
            plotfile = os.path.dirname(
                g_vis) + '/selfcal/plots/' + name + '_uvwave_amp_corrected-model.jpg'
            if not os.path.isfile(plotfile):
                plotms(vis=g_vis, xaxis='UVwave', yaxis='amp',
                       # antenna=self.config.antennas, spw=self.config.spws,
                       coloraxis='baseline', avgantenna=False, avgscan=False,
                       ydatacolumn='corrected-model', avgchannel='9999', avgtime='9999',
                       correlation='RR', plotrange=[0, 0, 0, 0],
                       width=2000, height=800, showgui=False, overwrite=True, dpi=1200,
                       highres=True,
                       customsymbol=True, symbolsize=4, symbolshape='diamond',
                       plotfile=plotfile)
            else:
                pass

        pass

    def plot_uvwave(self, g_vis, name):

        plotfile = os.path.dirname(g_vis) + '/selfcal/plots/' + name + '_uvwave.jpg'
        if not os.path.isfile(plotfile):
            plotms(vis=g_vis, xaxis='uwave', yaxis='vwave', avgantenna=False,
                   # antenna=self.config.antennas, spw=self.config.spws,
                   coloraxis='spw',
                   ydatacolumn='data', avgchannel='16', avgtime='12',
                   correlation='RR,LL',
                   width=800, height=800, showgui=False, overwrite=True, dpi=1200, highres=True,
                   # customsymbol=True,symbolsize=1,symbolshape='diamond',
                   plotfile=plotfile)
        plotfile = os.path.dirname(g_vis) + '/selfcal/plots/' + name + '_uv.jpg'
        if not os.path.isfile(plotfile):
            plotms(vis=g_vis, xaxis='u', yaxis='v', avgantenna=False,
                   # antenna=self.config.antennas, spw=self.config.spws,
                   coloraxis='spw',
                   ydatacolumn='data', avgchannel='9999', avgtime='2',
                   correlation='RR,LL',
                   width=800, height=800, showgui=False, overwrite=True, dpi=1200, highres=True,
                   # customsymbol=True,symbolsize=1,symbolshape='diamond',
                   plotfile=plotfile)

        pass

    def get_tb_data(self, table, param):
        tb.open(table)
        param_data = tb.getcol(param).ravel()
        tb.close()
        return (param_data)

    def make_plot_snr(self, caltable, cut_off, plot_snr=True, bins=50, density=True,
                      save_fig=False):
        import numpy as np
        import matplotlib.pyplot as plt
        from scipy import stats
        snr = self.get_tb_data(caltable, 'SNR')
        plt.figure(figsize=(3, 3))
        if plot_snr:
            plt.hist(snr, bins=bins, density=density, histtype='step')
            # plt.legend( loc='upper right' )
            plt.xlabel('SNR')
            # plt.semilogy()
            # plt.semilogx()
            plt.axvline(x=3, color='k', linestyle='--')
            plt.axvline(x=cut_off, color='r', linestyle='--')
            plt.grid()
            if save_fig == True:
                try:
                    plt.savefig(caltable.replace('.tb', '.jpg'), dpi=300, bbox_inches='tight')
                except:
                    plt.savefig(caltable + '.jpg', dpi=300, bbox_inches='tight')
            # plt.show()
            plt.clf()
            plt.close()

        fraction_flagged_solutions = stats.percentileofscore(snr, cut_off)

        print(f" ++==>> Fraction of flagged solutions with SNR < {cut_off} is"
              f" {fraction_flagged_solutions:.2f}%")

        # print('P(<=' + str(cut_off) + ') = {0}  ({1})'.format(
        #     stats.percentileofscore(snr, cut_off), ''))
        pass

    def calibration_table_plot(self, table, stage='selfcal',
                               table_type='gain_phase', kind='',
                               xaxis='time', yaxis='phase',
                               fields=[''], showgui=True):
        if not os.path.exists(os.path.dirname(table) + '/plots/' + stage):
            os.makedirs(os.path.dirname(table) + '/plots/' + stage)

        if yaxis == 'phase':
            plotrange = [-1, -1, -180, 180]
        else:
            plotrange = [-1, -1, -1, -1]

        if fields == '':

            plotms(vis=table, xaxis=xaxis, yaxis=yaxis, field='',
                   gridcols=1, gridrows=1, coloraxis='spw', antenna='', plotrange=plotrange,
                   width=1000, height=400, dpi=600, overwrite=True, showgui=showgui,
                   # correlation='LL,RR',
                   plotfile=os.path.dirname(
                       table) + '/plots/' + stage + '/' + table_type + '_' + xaxis + '_' + yaxis + '_field_' + str(
                       'all') + '.jpg')

            plotms(vis=table, xaxis=xaxis, yaxis=yaxis, field='', avgbaseline=True,
                   gridcols=1, gridrows=1, coloraxis='spw', antenna='', plotrange=plotrange,
                   width=1000, height=400, dpi=600, overwrite=True, showgui=showgui,
                   # correlation='LL,RR',
                   plotfile=os.path.dirname(
                       table) + '/plots/' + stage + '/' + table_type + '_' + xaxis + '_' +
                            yaxis + '_avgbaseline_field_' + str(
                       'all') + '.jpg')

        else:

            for FIELD in fields:
                plotms(vis=table, xaxis=xaxis, yaxis=yaxis, field=FIELD,
                       # gridcols=4,gridrows=4,coloraxis='spw',antenna='',iteraxis='antenna',
                       # width=2048,height=1280,dpi=256,overwrite=True,showgui=False,
                       gridcols=1, gridrows=1, coloraxis='spw', antenna='',
                       plotrange=plotrange,
                       width=1000, height=400, dpi=600, overwrite=True, showgui=showgui,
                       # correlation='LL,RR',
                       plotfile=os.path.dirname(
                           table) + '/plots/' + stage + '/' + table_type + '_' + xaxis + '_' + yaxis + '_field_' + str(
                           FIELD) + '.jpg')

        pass

    def check_solutions(self, g_name, cut_off=2.0, minsnr=2.0, n_interaction=0, uvrange='',
                        solnorm=False, combine='', calmode='p', gaintype='G', solint_factor=1.0,
                        interp='', spwmap=[],
                        gain_tables_selfcal=[''], special_name='', refant='', minblperant=4,
                        return_solution_stats=False):
        g_vis = g_name + '.ms'
        minsnr = minsnr
        solint_template = np.asarray([24, 48, 96, 192, 384])
        solints = solint_template * solint_factor

        caltable_int = (os.path.dirname(g_name) + '/selfcal/selfcal_test_' + str(
            n_interaction) + '_' + os.path.basename(g_name) + '_solint_int_minsnr_' + str(
            minsnr) + '_calmode' + calmode + '_combine' + combine + '_gtype_' +
                        gaintype + special_name + '.tb')

        caltable_1 = (os.path.dirname(g_name) + '/selfcal/selfcal_test_' + str(
            n_interaction) + '_' + os.path.basename(g_name) + '_solint_' +
                      str(int(solints[0])) + '_minsnr_' + str(
                    minsnr) + '_calmode' + calmode + '_combine' + combine + '_gtype_' +
                      gaintype + special_name + '.tb')

        caltable_2 = (os.path.dirname(g_name) + '/selfcal/selfcal_test_' + str(
            n_interaction) + '_' + os.path.basename(g_name) + '_solint_' +
                      str(int(solints[1])) + '_minsnr_' + str(
                    minsnr) + '_calmode' + calmode + '_combine' + combine + '_gtype_' +
                      gaintype + special_name + '.tb')

        caltable_3 = (os.path.dirname(g_name) + '/selfcal/selfcal_test_' + str(
            n_interaction) + '_' + os.path.basename(g_name) + '_solint_' +
                      str(int(solints[2])) + '_minsnr_' + str(
                    minsnr) + '_calmode' + calmode + '_combine' + combine + '_gtype_' +
                      gaintype + special_name + '.tb')

        caltable_4 = (os.path.dirname(g_name) + '/selfcal/selfcal_test_' + str(
            n_interaction) + '_' + os.path.basename(g_name) + '_solint_' +
                      str(int(solints[3])) + '_minsnr_' + str(
                    minsnr) + '_calmode' + calmode + '_combine' + combine + '_gtype_' +
                      gaintype + special_name + '.tb')

        caltable_5 = (os.path.dirname(g_name) + '/selfcal/selfcal_test_' + str(
            n_interaction) + '_' + os.path.basename(g_name) + '_solint_' +
                      str(int(solints[4])) + '_minsnr_' + str(
                    minsnr) + '_calmode' + calmode + '_combine' + combine + '_gtype_' +
                      gaintype + special_name + '.tb')

        caltable_inf = (os.path.dirname(g_name) + '/selfcal/selfcal_test_' + str(
            n_interaction) + '_' + os.path.basename(g_name) + '_solint_inf_minsnr_' + str(
            minsnr) + '_calmode' + calmode + '_combine' + combine + '_gtype_' +
                        gaintype + special_name + '.tb')

        if not os.path.exists(caltable_int):
            print('>> Performing test-gaincal for solint=int...')
            gaincal(vis=g_vis, caltable=caltable_int, solint='int', refant=refant,
                    interp=interp, spwmap=spwmap,
                    solnorm=solnorm, combine=combine, minblperant=minblperant,
                    calmode=calmode, gaintype=gaintype, minsnr=minsnr, uvrange=uvrange,
                    gaintable=gain_tables_selfcal)
        if not os.path.exists(caltable_1):
            print('>> Performing test-gaincal for solint=' + str(solints[0]) + 's...')
            gaincal(vis=g_vis, caltable=caltable_1, solint=str(solints[0]) + 's',
                    refant=refant, interp=interp, spwmap=spwmap,
                    solnorm=solnorm, combine=combine, minblperant=minblperant,
                    calmode=calmode, gaintype=gaintype, minsnr=minsnr, uvrange=uvrange,
                    gaintable=gain_tables_selfcal)
        if not os.path.exists(caltable_2):
            print('>> Performing test-gaincal for solint=' + str(solints[1]) + 's...')
            gaincal(vis=g_vis, caltable=caltable_2, solint=str(solints[1]) + 's',
                    refant=refant, interp=interp, spwmap=spwmap,
                    solnorm=solnorm, combine=combine, minblperant=minblperant,
                    calmode=calmode, gaintype=gaintype, minsnr=minsnr, uvrange=uvrange,
                    gaintable=gain_tables_selfcal)
        if not os.path.exists(caltable_3):
            print('>> Performing test-gaincal for solint=' + str(solints[2]) + 's...')
            gaincal(vis=g_vis, caltable=caltable_3, solint=str(solints[2]) + 's',
                    refant=refant, interp=interp, spwmap=spwmap,
                    solnorm=solnorm, combine=combine, minblperant=minblperant,
                    calmode=calmode, gaintype=gaintype, minsnr=minsnr, uvrange=uvrange,
                    gaintable=gain_tables_selfcal)
        if not os.path.exists(caltable_4):
            print('>> Performing test-gaincal for solint=' + str(solints[3]) + 's...')
            gaincal(vis=g_vis, caltable=caltable_4, solint=str(solints[3]) + 's',
                    refant=refant, interp=interp, spwmap=spwmap,
                    solnorm=solnorm, combine=combine, minblperant=minblperant,
                    calmode=calmode, gaintype=gaintype, minsnr=minsnr, uvrange=uvrange,
                    gaintable=gain_tables_selfcal)
        if not os.path.exists(caltable_5):
            print('>> Performing test-gaincal for solint=' + str(solints[4]) + 's...')
            gaincal(vis=g_vis, caltable=caltable_5, solint=str(solints[4]) + 's',
                    refant=refant, interp=interp, spwmap=spwmap,
                    solnorm=solnorm, combine=combine, minblperant=minblperant,
                    calmode=calmode, gaintype=gaintype, minsnr=minsnr, uvrange=uvrange,
                    gaintable=gain_tables_selfcal)
        if not os.path.exists(caltable_inf):
            print('>> Performing test-gaincal for solint=inf...')
            gaincal(vis=g_vis, caltable=caltable_inf, solint='inf', refant=refant,
                    interp=interp, spwmap=spwmap,
                    solnorm=solnorm, combine=combine, minblperant=minblperant,
                    calmode=calmode, gaintype=gaintype, minsnr=minsnr, uvrange=uvrange,
                    gaintable=gain_tables_selfcal)

        def make_plot_check(cut_off=cut_off, return_solution_stats=False):
            import numpy as np
            import matplotlib.pyplot as plt
            from scipy import stats
            snr_int = get_tb_data(caltable_int, 'SNR')
            # snr_5 = get_tb_data(caltable_5,'SNR')
            snr_1 = self.get_tb_data(caltable_1, 'SNR')
            snr_2 = self.get_tb_data(caltable_2, 'SNR')
            snr_3 = self.get_tb_data(caltable_3, 'SNR')
            snr_4 = self.get_tb_data(caltable_4, 'SNR')
            snr_5 = self.get_tb_data(caltable_5, 'SNR')
            snr_inf = self.get_tb_data(caltable_inf, 'SNR')

            plt.figure()
            plt.hist(snr_int, bins=50, density=True, histtype='step',
                     label='int')
            plt.hist(snr_1, bins=50, density=True, histtype='step',
                     label=str(solints[0]) + ' seconds')
            plt.hist(snr_2, bins=50, density=True, histtype='step',
                     label=str(solints[1]) + ' seconds')
            plt.hist(snr_3, bins=50, density=True, histtype='step',
                     label=str(solints[2]) + ' seconds')
            plt.hist(snr_4, bins=50, density=True, histtype='step',
                     label=str(solints[3]) + ' seconds')
            plt.hist(snr_5, bins=50, density=True, histtype='step',
                     label=str(solints[4]) + ' seconds')
            plt.hist(snr_inf, bins=50, density=True, histtype='step',
                     label='inf')
            plt.legend(loc='upper right')
            plt.xlabel('SNR')
            # plt.semilogx()
            plt.savefig(os.path.dirname(g_name) + '/selfcal/plots/' + str(n_interaction) +
                        '_' + os.path.basename(
                g_name) + '_calmode' + calmode + '_combine' + combine + '_gtype_' + gaintype
                        + special_name + '_gain_solutions_comparisons_norm.pdf')
            # plt.clf()
            # plt.close()
            # plt.show()
            plt.figure()
            plt.hist(snr_int, bins=50, density=False, histtype='step',
                     label='int')
            plt.hist(snr_1, bins=50, density=False, histtype='step',
                     label=str(solints[0]) + ' seconds')
            plt.hist(snr_2, bins=50, density=False, histtype='step',
                     label=str(solints[1]) + ' seconds')
            plt.hist(snr_3, bins=50, density=False, histtype='step',
                     label=str(solints[2]) + ' seconds')
            plt.hist(snr_4, bins=50, density=False, histtype='step',
                     label=str(solints[3]) + ' seconds')
            plt.hist(snr_5, bins=50, density=False, histtype='step',
                     label=str(solints[4]) + ' seconds')
            plt.hist(snr_inf, bins=50, density=False, histtype='step',
                     label='inf')
            plt.legend(loc='upper right')
            plt.xlabel('SNR')
            # plt.semilogx()
            plt.savefig(os.path.dirname(g_name) + '/selfcal/plots/' + str(n_interaction) +
                        '_' + os.path.basename(
                g_name) + '_calmode' + calmode + '_combine' + combine +
                        '_gtype_' + gaintype + special_name +
                        '_gain_solutions_comparisons.pdf')

            print('P(<=' + str(cut_off) + ') = {0}  ({1})'.format(
                stats.percentileofscore(snr_int, cut_off), 'int'))
            print('P(<=' + str(cut_off) + ') = {0}  ({1})'.format(
                stats.percentileofscore(snr_1, cut_off), str(solints[0]) + ' s'))
            print('P(<=' + str(cut_off) + ') = {0}  ({1})'.format(
                stats.percentileofscore(snr_2, cut_off), str(solints[1]) + ' s'))
            print('P(<=' + str(cut_off) + ') = {0}  ({1})'.format(
                stats.percentileofscore(snr_3, cut_off), str(solints[2]) + ' s'))
            print('P(<=' + str(cut_off) + ') = {0}  ({1})'.format(
                stats.percentileofscore(snr_4, cut_off), str(solints[3]) + ' s'))
            print('P(<=' + str(cut_off) + ') = {0}  ({1})'.format(
                stats.percentileofscore(snr_5, cut_off), str(solints[4]) + ' s'))
            print('P(<=' + str(cut_off) + ') = {0}  ({1})'.format(
                stats.percentileofscore(snr_inf, cut_off), 'inf'))

            # plt.show()
            # print('################################')
            # print(np.mean(snr_int))
            # print('################################')
            # print(stats.percentileofscore(snr_int, cut_off))
            SNRs = [
                np.array(snr_int),
                np.array(snr_1),
                np.array(snr_2),
                np.array(snr_3),
                np.array(snr_4),
                np.array(snr_5),
                np.array(snr_inf)]
            percentiles_SNRs = np.asarray([
                stats.percentileofscore(snr_int, cut_off),
                stats.percentileofscore(snr_1, cut_off),
                stats.percentileofscore(snr_2, cut_off),
                stats.percentileofscore(snr_3, cut_off),
                stats.percentileofscore(snr_4, cut_off),
                stats.percentileofscore(snr_5, cut_off),
                stats.percentileofscore(snr_inf, cut_off)])

            snr_data = {
                'int': SNRs[0],
                '24s': SNRs[1],
                '48s': SNRs[2],
                '96s': SNRs[3],
                '192s': SNRs[4],
                '384s': SNRs[5],
                'inf': SNRs[6]
            }

            if return_solution_stats:
                return (snr_data, percentiles_SNRs)
            else:
                pass
            plt.clf()
            plt.close()

        def compare_phase_variation():
            plotms(caltable_1, antenna='', scan='', yaxis='phase', avgbaseline=True)

            plotms(caltable_3, antenna='', scan='', yaxis='phase', plotindex=1,
                   clearplots=False, customsymbol=True, symbolsize=20, avgbaseline=True,
                   symbolcolor='ff0000', symbolshape='circle')

            plotms(caltable_2, antenna='', scan='', yaxis='phase', plotindex=2,
                   clearplots=False, customsymbol=True, symbolsize=12, avgbaseline=True,
                   symbolcolor='green', symbolshape='square')

            plotms(caltable_inf, antenna='', scan='', yaxis='phase', plotindex=3,
                   clearplots=False, customsymbol=True, symbolsize=8, avgbaseline=True,
                   symbolcolor='yellow', symbolshape='square')

            plotms(caltable_4, antenna='', scan='', yaxis='phase', plotindex=4,
                   clearplots=False, customsymbol=True, symbolsize=4, avgbaseline=True,
                   symbolcolor='purple', symbolshape='square',
                   width=1300, height=400, showgui=True, overwrite=True,
                   plotfile=os.path.dirname(g_name) + '/selfcal/plots/' + str(
                       n_interaction) +
                            '_' + os.path.basename(
                       g_name) + '_combine' + '_calmode' + calmode + combine +
                            '_gtype_' + gaintype + special_name +
                            '_phase_variation_intervals.jpg')

        def compare_amp_variation():
            plotms(caltable_1, antenna='', scan='', yaxis='amp',
                   plotrange=[0, 0, 0, 0],
                   avgbaseline=True)

            plotms(caltable_3, antenna='', scan='', yaxis='amp', plotindex=1,
                   plotrange=[0, 0, 0, 0],
                   clearplots=False, customsymbol=True, symbolsize=20, avgbaseline=True,
                   symbolcolor='ff0000', symbolshape='circle')

            plotms(caltable_2, antenna='', scan='', yaxis='amp', plotindex=2,
                   plotrange=[0, 0, 0, 0],
                   clearplots=False, customsymbol=True, symbolsize=12, avgbaseline=True,
                   symbolcolor='green', symbolshape='square')

            plotms(caltable_inf, antenna='', scan='', yaxis='amp', plotindex=3,
                   plotrange=[0, 0, 0, 0],
                   clearplots=False, customsymbol=True, symbolsize=8, avgbaseline=True,
                   symbolcolor='yellow', symbolshape='square')

            plotms(caltable_4, antenna='', scan='', yaxis='amp', plotindex=4,
                   plotrange=[0, 0, 0, 0],
                   clearplots=False, customsymbol=True, symbolsize=4, avgbaseline=True,
                   symbolcolor='purple', symbolshape='square',
                   width=1300, height=400, showgui=True, overwrite=True,
                   plotfile=os.path.dirname(g_name) + '/selfcal/plots/' + str(
                       n_interaction) +
                            '_' + os.path.basename(g_name) +
                            '_combine' + '_calmode' + calmode + combine +
                            '_gtype_' + gaintype + special_name +
                            '_amp_variation_intervals.jpg')

        #
        # def plot_gains():
        #     plotms(caltable_int,antenna='ea01',scan='',yaxis='phase',
        #         gridrows=5,gridcols=5,iteraxis='antenna',coloraxis='spw')

        if return_solution_stats == True:
            SNRs, percentiles_SNRs = \
                make_plot_check(cut_off=cut_off,
                                return_solution_stats=return_solution_stats)
        else:
            make_plot_check(cut_off=cut_off)
        compare_phase_variation()
        if calmode == 'ap':
            compare_amp_variation()

        if return_solution_stats == True:
            return (SNRs, percentiles_SNRs, caltable_int, caltable_3, caltable_inf)
        else:
            pass

    def start_image(self, g_name, n_interaction, imsize='2048', imsizey=None, cell='0.05asec',
                    robust=0.0,
                    base_name=None,
                    nsigma_automask='7.0', nsigma_autothreshold='0.1',
                    delmodel=False, niter=600,
                    opt_args='', quiet=True, shift=None,
                    nc=4, negative_arg='negative',
                    PLOT=False, datacolumn='DATA', mask=None,
                    savemodel=True, uvtaper=[""]):
        '''
        Wsclean wrapper function. It calls wslcean from the command line with some
        predifined arguments. This initial step runs on the DATA column and creates
        the initial model which is used to calculate the initial complex self-gains.
        '''
        g_vis = g_name + '.ms'
        if imsizey is None:
            imsizey = imsize
        if base_name is None:
            base_name = str(n_interaction) + '_start_image_'
        else:
            base_name = base_name

        os.system("export OPENBLAS_NUM_THREADS=1 && python imaging_with_wsclean.py --f " +
                  g_name + " --sx "
                  + str(imsize) + " --sy " + str(imsizey) + " --niter "
                  + str(niter) + " --data " + datacolumn + " --cellsize " + cell
                  + ' --nsigma_automask ' + nsigma_automask + ' --mask ' + str(mask)
                  + ' --nsigma_autothreshold ' + nsigma_autothreshold
                  # +' --opt_args '+ opt_args
                  + ' --quiet ' + str(quiet)
                  + ' --nc ' + str(nc) + ' --negative_arg ' + negative_arg
                  + ' --shift ' + str(shift)
                  + " --r " + str(robust) + " --t " + str(uvtaper)
                  + " --update_model " + str(savemodel) + " --save_basename " + base_name)

        if PLOT == True:
            self.plot_visibilities(g_vis=g_vis, name=base_name,
                                   with_MODEL=True, with_DATA=False, with_CORRECTED=True)

        pass

    def run_wsclean(self, g_name, n_interaction, imsize='2048', imsizey=None, cell='0.05asec',
                    robust=0.5, base_name=None,
                    savemodel=True, shift=None,
                    nsigma_automask='4.0', nsigma_autothreshold='2.0',
                    datacolumn='CORRECTED', mask=None,
                    niter=10000, quiet=True,
                    nc=4, negative_arg='negative',
                    with_multiscale=False, scales="'0,5,20,40'",
                    uvtaper=[], PLOT=False, with_DATA=True, with_CORRECTED=True, with_MODEL=True):

        g_vis = g_name + '.ms'
        if imsizey is None:
            imsizey = imsize
        if base_name is None:
            base_name = str(n_interaction) + '_update_model_image_'
        else:
            base_name = base_name

        os.system("export OPENBLAS_NUM_THREADS=1 && python imaging_with_wsclean.py --f " +
                  g_name + " --sx "
                  + str(imsize) + " --sy " + str(imsizey) + " --niter "
                  + str(niter) + " --data " + datacolumn + " --cellsize " + cell
                  + ' --nsigma_automask ' + nsigma_automask + ' --mask ' + str(mask)
                  + ' --nsigma_autothreshold ' + nsigma_autothreshold
                  + ' --scales ' + scales
                  + ' --nc ' + str(nc) + ' --negative_arg ' + negative_arg
                  # +' --opt_args '+ opt_args
                  + ' --quiet ' + str(quiet) + ' --with_multiscale ' + str(with_multiscale)
                  + ' --shift ' + str(shift)
                  + " --r " + str(robust) + " --t " + str(uvtaper)
                  + " --update_model " + str(savemodel) + " --save_basename " + base_name)

        if PLOT == True:
            self.plot_visibilities(g_vis=g_vis, name=base_name, with_DATA=with_DATA,
                                   with_MODEL=with_MODEL, with_CORRECTED=with_CORRECTED)

        pass

    def self_gain_cal(self, g_name, n_interaction, field='*',
                      gain_tables=[], spwmaps=[],
                      combine='', solnorm=False, normtype='median',
                      spw='*', refantmode='strict',
                      spwmap=[], uvrange='', append=False, solmode='',  # L1R
                      minsnr=5.0, solint='inf', gaintype='G', calmode='p',
                      applymode='calflag',
                      interp='', refant='', minblperant=4,
                      action='apply', flagbackup=True, calwt=False,
                      PLOT=False, with_CORRECTED=True, with_MODEL=True, with_DATA=False,
                      overwrite_gaintable=False,
                      special_name=''):
        g_vis = g_name + '.ms'
        # refantmode = 'flex' if refantmode == 'flex' else 'strict'
        cal_basename = '_selfcal_'
        base_name = str(n_interaction) + '_update_model_image_' + cal_basename
        # base_name =  str(n_interaction)+'_loop_correct_'+cal_basename

        if calmode == 'p':
            cal_basename = cal_basename + 'phase_'
            base_name = base_name + 'phase_'
        if calmode == 'ap' or calmode == 'a':
            cal_basename = cal_basename + 'ampphase_'
            base_name = base_name + 'ampphase_'
        if gain_tables != []:
            cal_basename = cal_basename + 'incremental_'

        caltable = (os.path.dirname(g_name) + '/selfcal/' + str(n_interaction) \
                    + cal_basename + os.path.basename(g_name) \
                    + '_' + '_solint_' + solint + '_minsnr_' + str(minsnr) +
                    '_combine' + combine + '_gtype_' + gaintype + special_name + '.tb')

        if solnorm == '':
            if calmode == 'ap' or calmode == 'a':
                print(' ++==> Using normalised solutions for amplitude self-calibration.')
                solnorm = True
            else:
                solnorm = False
        else:
            solnorm = solnorm

        if not os.path.exists(caltable):
            # overwrite_gaintable
            gaincal(vis=g_vis, field=field, caltable=caltable, spwmap=spwmap,
                    spw=spw,
                    solint=solint, gaintable=gain_tables, combine=combine,
                    refant=refant, calmode=calmode, gaintype=gaintype,
                    refantmode=refantmode,
                    uvrange=uvrange, append=append, solmode=solmode,
                    interp=interp,
                    minsnr=minsnr, solnorm=solnorm, normtype=normtype,
                    minblperant=minblperant)
        else:
            print(' => Using existing caltable with same parameters asked.')
            print(' => Not computing again...')

        self.calibration_table_plot(table=caltable,
                                    fields='', yaxis='phase',
                                    table_type=str(
                                        n_interaction) + '_selfcal_phase_' + os.path.basename(
                                        g_name) +
                                               '_solint_' + solint + '_minsnr_' + str(
                                        minsnr) + '_combine' + combine +
                                               '_gtype_' + gaintype + special_name)

        if calmode == 'ap' or calmode == 'a':
            self.calibration_table_plot(table=caltable,
                                        fields='', yaxis='amp',
                                        table_type=str(n_interaction) + '_selfcal_ampphase_' +
                                                   os.path.basename(g_name) + '_solint_' + solint +
                                                   '_minsnr_' + str(
                                            minsnr) + '_combine' + combine +
                                                   '_gtype_' + gaintype + special_name)

        self.make_plot_snr(caltable=caltable, cut_off=minsnr,
                           plot_snr=True, bins=50, density=True, save_fig=True)

        if action == 'apply':
            if flagbackup == True:
                print('     => Creating new flagbackup file before mode ',
                      calmode, ' selfcal ...')
                flagmanager(vis=g_vis, mode='save',
                            versionname='before_selfcal_mode_' + calmode,
                            comment='Before selfcal apply.')

            gain_tables.append(caltable)
            # if spwmap != []:
            #     spwmaps.append(spwmap[-1])
            # else:
            #     spwmaps.append(spwmap)
            print('     => Reporting data flagged before selfcal '
                  'apply interaction', n_interaction, '...')
            summary_bef = flagdata(vis=g_vis, field='', mode='summary')
            self.report_flag(summary_bef, 'field')
            self.report_flag(summary_bef, 'observation')

            # if calmode == 'ap' or calmode == 'a':

            applycal(vis=g_vis, gaintable=gain_tables, spwmap=spwmap,
                     interp=interp,
                     applymode=applymode,
                     flagbackup=False, calwt=calwt)

            print('     => Reporting data flagged after selfcal '
                  'apply interaction', n_interaction, '...')
            summary_aft = flagdata(vis=g_vis, field='', mode='summary')
            self.report_flag(summary_aft, 'field')
            self.report_flag(summary_aft, 'observation')
            # flag_data_steps[f'step_{n_interaction}'] = summary_aft.copy()

            if PLOT == True:
                self.plot_visibilities(g_vis=g_vis, name=base_name,
                                       with_CORRECTED=with_CORRECTED,
                                       with_MODEL=with_MODEL,
                                       with_DATA=with_DATA)

        return (gain_tables, spwmap)

    def run_autoflag(self, g_vis, display='report', action='calculate',
                     timedevscale=3.5, freqdevscale=3.5,
                     timecutoff=3.0, freqcutoff=3.0,
                     maxnpieces=5,
                     mode='tfcrop', ntime='120s',
                     winsize=7, datacolumn='corrected'):
        if action == 'apply':
            print(' ++==>> Flag statistics before auto-flagging')
            summary_before = flagdata(vis=g_vis, field='', mode='summary')
            self.report_flag(summary_before, 'field')
            self.report_flag(summary_before, 'observation')
            flagmanager(vis=g_vis, mode='save', versionname='selfcal_before_' + mode,
                        comment='Before ' + mode + ' at selfcal step.')
        if mode == 'clip':
            print(' ++==>> Using clip mode for flagging...')
            flagdata(vis=g_vis, mode='clip', field='', spw='',
                     datacolumn=datacolumn, clipzeros=True, clipoutside=True,
                     extendflags=False,
                     clipminmax=[0, 1.0],
                     # channelavg=True, chanbin=1, timeavg=True, timebin='24s',
                     # timedevscale=timedevscale, freqdevscale=freqdevscale,
                     action=action, flagbackup=False, savepars=False)
        if mode == 'tfcrop':
            print(' ++==>> Using tfcrop mode for flagging...')
            flagdata(vis=g_vis, mode='tfcrop', field='', spw='',
                     datacolumn=datacolumn, ntime=ntime, combinescans=False,
                     extendflags=False, winsize=winsize, maxnpieces=maxnpieces,
                     flagnearfreq=False,
                     flagneartime=False, growaround=True,
                     #  usewindowstats='both', halfwin=2,
                     usewindowstats='sum', halfwin=2,
                     #  channelavg=True,chanbin=4,
                     #  timeavg=True, timebin='36s',
                     timecutoff=timecutoff, freqcutoff=freqcutoff,
                     freqfit='line',
                     action=action, flagbackup=False, savepars=False,
                     )
        if mode == 'rflag':
            print(' ++==>> Using rflag mode for flagging...')
            flagdata(vis=g_vis, mode='rflag', field='', spw='', display=display,
                     datacolumn=datacolumn, ntime=ntime, combinescans=False,
                     extendflags=False,
                     winsize=winsize,
                     # channelavg=True,chanbin=1,
                     # timeavg=True, timebin='24s',
                     timedevscale=timedevscale, freqdevscale=freqdevscale,
                     flagnearfreq=False, flagneartime=False, growaround=True,
                     action=action, flagbackup=False, savepars=True
                     )
        if action == 'apply':
            flagdata(vis=g_vis, field='', spw='',
                     datacolumn=datacolumn,
                     mode='extend', action=action, display='report',
                     flagbackup=False, growtime=80.0, growaround=True,
                     growfreq=80.0, extendpols=False)
            # flagdata(vis=g_vis, field='', spw='',
            #          datacolumn=datacolumn,
            #          mode='extend', action=action, display='report',
            #          flagbackup=False, growtime=75.0,
            #          growfreq=75.0, extendpols=False)

        if action == 'apply':
            flagmanager(vis=g_vis, mode='save', versionname='selfcal_after_rflag',
                        comment='After rflag at selfcal step.')
            # try:
            #     print(' ++==> Running statwt...')
            #     statwt(vis=g_vis, statalg='chauvenet', timebin='24s',
            #            datacolumn='corrected',minsamp = 3)
            # except:
            #     print(' ++==> Running statwt...')
            #     statwt(vis=g_vis, statalg='chauvenet', timebin='24s',
            #            datacolumn='data', minsamp = 3)

            print(' ++==> Flag statistics after rflag:')
            summary_after = flagdata(vis=g_vis, field='', mode='summary')
            self.report_flag(summary_after, 'field')
            self.report_flag(summary_after, 'observation')

    def find_refant(self, msfile, field, tablename, combine='', minsnr=1.0):
        """
        This function comes from the e-MERLIN CASA Pipeline.
        https://github.com/e-merlin/eMERLIN_CASA_pipeline/blob/master/functions/eMCP_functions.py#L1501
        """
        # Find phase solutions per scan:
        # tablename = calib_dir +
        # if not os.path.exists(tablename):
        gaincal(vis=msfile,
                caltable=tablename,
                field=field,
                refantmode='flex',
                combine=combine,
                solint='inf',
                minblperant=3,
                gaintype='G',
                minsnr=minsnr,
                # combine='spw',
                calmode='p')
        # find_casa_problems()
        # Read solutions (phases):
        tb.open(tablename + '/ANTENNA')
        antenna_names = tb.getcol('NAME')
        tb.close()
        tb.open(tablename)
        antenna_ids = tb.getcol('ANTENNA1')
        # times  = tb.getcol('TIME')
        flags = tb.getcol('FLAG')
        phases = np.angle(tb.getcol('CPARAM'))
        snrs = tb.getcol('SNR')
        tb.close()

        self.make_plot_snr(caltable=tablename, cut_off=1.0,
                           plot_snr=True, bins=50, density=True, save_fig=True)

        self.calibration_table_plot(table=tablename,
                                    fields='', yaxis='phase',
                                    table_type='select_refant',
                                    showgui=False)

        # Analyse number of good solutions:
        good_frac = []
        good_snrs = []
        for i, ant_id in enumerate(np.unique(antenna_ids)):
            cond = antenna_ids == ant_id
            # t = times[cond]
            f = flags[0, 0, :][cond]
            p = phases[0, 0, :][cond]
            snr = snrs[0, 0, :][cond]
            frac = 1.0 * np.count_nonzero(~f) / len(f) * 100.
            snr_mean = np.nanmean(snr[~f])
            good_frac.append(frac)
            good_snrs.append(snr_mean)
        sort_idx = np.argsort(good_frac)[::-1]
        # sort_idx = np.argsort(good_snrs)[::-1]
        print('Antennas sorted by % of good solutions:')
        for i in sort_idx:
            print('{0:3}: {1:4.1f}, <SNR> = {2:4.1f}'.format(antenna_names[i],
                                                             good_frac[i],
                                                             good_snrs[i]))
        if good_frac[sort_idx[0]] < 90:
            print('Small fraction of good solutions with selected refant!')
            print('Please inspect antennas to select optimal refant')
            print('You may want to use refantmode= flex" in default_params')
        pref_ant = antenna_names[sort_idx]
        # if 'Lo' in antenna_names:
        #     priorities = ['Pi','Da','Kn','De','Cm']
        # else:
        #     priorities = ['Mk2','Pi','Da','Kn', 'Cm', 'De']
        # refant = ','.join([a for a in pref_ant if a in priorities])
        pref_ant_list = ','.join(list(pref_ant))
        return pref_ant_list

    def get_phase_centre(self, vis):
        from astropy.coordinates import SkyCoord
        import astropy.units as u

        msmd.open(vis)
        ra_radians = msmd.phasecenter()['m0']['value']
        dec_radians = msmd.phasecenter()['m1']['value']
        msmd.close()
        # Convert to SkyCoord object
        coord = SkyCoord(ra=ra_radians * u.radian, dec=dec_radians * u.radian, frame='icrs')

        # Format the output using 'hmsdms'
        formatted_coord = coord.to_string('hmsdms')
        formatted_ra, formatted_dec = formatted_coord.split()

        formatted_ra_hms = formatted_ra.replace('h', ':').replace('m', ':').replace('s', '')
        formatted_dec_dms = formatted_dec.replace('d', '.').replace('m', '.').replace('s', '')

        formatted_output = "J2000 {} {}".format(formatted_ra_hms, formatted_dec_dms)
        frame = "J2000"
        coordinates = "{} {}".format(formatted_ra_hms, formatted_dec_dms)
        # print(formatted_output)
        return (formatted_output, frame, coordinates)

    # Function to calculate percentages
    def calculate_percentages(self, data):
        percentages = {}
        for category, subcategories in data.items():
            if isinstance(subcategories, dict):
                percentages[category] = {}
                for subcategory, values in subcategories.items():
                    if 'flagged' in values and 'total' in values:
                        percentages[category][subcategory] = (values['flagged'] / values[
                            'total']) * 100
        return percentages

    # Function to plot data for a given category
    def plot_category_data(self, percentages_over_steps, category):
        run_steps = list(percentages_over_steps.keys())
        plt.figure(figsize=(10, 10))
        # Ensure category exists in all steps
        if not all(category in step for step in percentages_over_steps.values()):
            print(f"Category '{category}' not found in all data steps.")
            return

        subcategories = set()
        for step in percentages_over_steps.values():
            subcategories.update(step.get(category, {}).keys())

        for subcategory in subcategories:
            flagged_percentages = [
                percentages_over_steps[step].get(category, {}).get(subcategory, 0) for step
                in run_steps
            ]
            plt.plot(run_steps, flagged_percentages, label=subcategory, marker='o')

        plt.xlabel('Processing Step')
        plt.ylabel('Flagged Data (%)')
        plt.title(f'Flagged Data Over Processing Steps ({category.capitalize()})')
        plt.legend()
        plt.ylim(0, 100)
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(f'{os.path.dirname(self.g_name)}/selfcal/flag_stats'
                    f'_{category}_flagged_data.jpg', dpi=300,
                    bbox_inches='tight')
        # plt.show()
        plt.clf()
        plt.close()

    def report_flag(self, summary, axis):
        """Report flagged data percentage"""
        for id, stats in summary[axis].items():
            print(f'{axis} {id}: {100. * stats["flagged"] / stats["total"]:.1f} percent flagged')

    def initialize_storage(self):
        """Initialize storage containers for pipeline data"""
        self.image_list = {}
        self.residual_list = {}
        self.model_list = {}
        self.image_statistics = {}
        self.gain_tables_applied = {}
        self.spwmaps_applied = {}
        self.flag_data_steps = {}
        self.parameter_selection = {}
        self.trial_gain_tables = []
        self.final_gain_tables = []

    def _init_names(self):
        # import casalogger.__main__
        self.g_name_ = f"{self.config.path}{self.config.vis_name}"
        self.g_vis_ = f"{self.g_name_}.ms"
        if self.config.general_settings['new_phasecentre'] is None:
            self.g_name = self.g_name_
            self.g_vis = self.g_vis_
        else:
            print(" ++==>> A shift for the phase centre was provided.")
            self.g_name = self.g_name_ + '_phs'
            self.g_vis = self.g_name + '.ms'
            _, self.rame, self.or_phc = self.get_phase_centre(self.g_vis_)
            if not os.path.exists(self.g_vis):
                print(f" ++==>> going to run phaseshift on input MS.")
                print(f" ++==>> Shifting from ")
                print(f"          {self.or_phc}   ")
                print(f" ++==>> to ")
                print(f"          {self.config.general_settings['new_phasecentre']}   ")
                phaseshift(vis=self.g_vis_,
                           outputvis=self.g_vis,
                           phasecenter=f"{self.rame} {self.config.general_settings['new_phasecentre']}")
            else:
                print('--==>> Phase-shifted MS already exist in current directory.')
                print('--==>> Going to skip.')
                pass

    def _create_directories(self):
        # print(f"++==> Preparing to selfcalibrate {self.g_vis_}.")
        print('++==> Creating basic directory structure.')
        if not os.path.exists(self.config.path + 'selfcal/'):
            os.makedirs(self.config.path + 'selfcal/')
        if not os.path.exists(self.config.path + 'selfcal/plots'):
            os.makedirs(self.config.path + 'selfcal/plots')

    def _generate_listobs(self):
        print('++==> Generating listobs file.')
        listobs(vis=self.g_vis, listfile=self.g_name + '.listobs', overwrite=True)
        if self.g_name != self.g_name_:
            listobs(vis=self.g_vis_, listfile=self.g_name_ + '.listobs', overwrite=True)

    def _clear_model(self):
        print('--==> Clearing model column...')
        delmod(vis=self.g_vis, otf=True, scr=False)
        clearcal(vis=self.g_vis)

    def run_step(self, step_name):
        """Execute a specific pipeline step"""
        if step_name not in self.steps_performed:
            if hasattr(self, f'_run_{step_name}'):
                getattr(self, f'_run_{step_name}')()
                self.steps_performed.append(step_name)

    def _run_startup(self):
        """Execute startup step"""
        print(f"++==> Preparing to selfcalibrate {self.config.path}{self.config.vis_name}")
        self._init_names()
        self._create_directories()
        self._generate_listobs()
        self._clear_model()

    def _run_save_init_flags(self):
        """Execute initial flag saving step"""
        if not os.path.exists(f"{self.g_name}.ms.flagversions/flags.Original/"):
            print("     ==> Creating backup flags file 'Original'...")
            flagmanager(vis=f"{self.g_name}.ms", mode='save',
                        versionname='Original', comment='Original flags.')
        else:
            print("     --==> Skipping flagging backup init (exists).")
            print("     --==> Restoring flags to original...")
            flagmanager(vis=self.g_name + '.ms',
                        mode='restore',
                        versionname='Original')
        print(" ++==> Amount of data flagged at the start of selfcal.")
        summary = flagdata(vis=self.g_vis, field='', mode='summary')
        self.flag_data_steps['original'] = summary.copy()
        self.report_flag(summary, 'field')
        self.report_flag(summary, 'observation')
        # steps_performed.append('save_init_flags')

    def _run_statwt(self):
        if not os.path.exists(self.g_name + '.ms.flagversions/flags.statwt_1/'):
            print("     ==> Running statwt.")
            statwt(vis=self.g_vis,
                   statalg=self.config.general_settings['statwt_statalg'],
                   timebin=self.config.general_settings['timebin_statw'],
                   datacolumn='data')
        else:
            print("     --==> Skipping statwt (flag file exists).")

    def _run_initweights(self):
        print("     ++==> Running initweights.")
        initweights(vis=g_vis,dowtsp=False,wtmode='nyq')

    def _run_autoflag_init(self):
        """
        Run automatic rflag on the data before selfcalibration.
        """
        self.plot_visibilities(g_vis=self.g_vis, name='before_initial_flag',
                               with_DATA=True,
                               with_MODEL=False, with_CORRECTED=False)

        # run_autoflag(g_vis, display='both', action='apply', mode = 'tfcrop',
        #              timecutoff=3.5, freqcutoff=3.5,
        #              winsize=5, datacolumn='data')
        self.run_autoflag(self.g_vis, display='report', action='apply', mode='rflag',
                          timedevscale=4.0, freqdevscale=4.0,
                          winsize=7, datacolumn='data')
        # run_autoflag(g_vis, action='apply', mode = 'clip',datacolumn='corrected')

        summary = flagdata(vis=self.g_vis, field='', mode='summary')
        self.flag_data_steps['autoflag_init'] = summary.copy()

        self.plot_visibilities(g_vis=self.g_vis, name='after_initial_flag',
                               with_DATA=True,
                               with_MODEL=False, with_CORRECTED=False)

    def _run_select_refant(self):
        if self.config.refant == '':
            # if 'select_refant' in steps and 'select_refant' not in steps_performed:
            print(' ++==> Estimating order of best referent antennas...')
            tablename_refant = os.path.dirname(self.g_name) + '/selfcal/find_refant.phase'
            refant = self.find_refant(msfile=self.g_vis, field='',
                                      #   combine='spw',
                                      tablename=tablename_refant)
            print(' ++==> Preferential reference antenna order = ', refant)
            self.refant = refant
            # print(f"The referent antenna is {refant}")
        else:
            print(' ++==> Using provided referent antenna order.')
            print(f"     => {self.config.refant}")
            self.refant = self.config.refant

    def _run_test_image(self):
        niter_test = self.config.init_parameters['test_image']['niter']
        robust = self.config.init_parameters['test_image']['robust']
        prefix = self.config.init_parameters['test_image']['prefix']
        self.run_wsclean(self.g_name,
                         imsize=self.config.imsize,
                         imsizey=self.config.imsizey,
                         cell=self.config.cell_size,
                         robust=robust, base_name=prefix,
                         nsigma_automask='4.0', nsigma_autothreshold='2.0',
                         n_interaction='0', savemodel=False, quiet=self.config.quiet,
                         datacolumn='DATA', shift=self.config.global_parameters['FIELD_SHIFT'],
                         with_multiscale=False, scales='0,5,20,40',
                         nc=self.config.nc, negative_arg=self.config.negative_arg,
                         uvtaper=self.config.init_parameters['test_image']['uvtaper'],
                         niter=niter_test,
                         PLOT=True, with_DATA=True, with_CORRECTED=False, with_MODEL=False)

        self.image_statistics, self.image_list = self.compute_image_stats(path=self.config.path,
                                                                          image_list=self.image_list,
                                                                          image_statistics=self.image_statistics,
                                                                          prefix=prefix,
                                                                          selfcal_step='test_image')

        current_total_flux = self.image_statistics['test_image']['total_flux_mask'] * 1000
        if current_total_flux > 50.0:
            self.image_statistics, self.image_list = self.compute_image_stats(path=self.config.path,
                                                                              image_list=self.image_list,
                                                                              image_statistics=self.image_statistics,
                                                                              prefix=prefix,
                                                                              sigma=25,
                                                                              selfcal_step='test_image')

        self.modified_robust = None
        if current_total_flux < 5.0:
            """
            Sometimes, a lower robust parameter (e.g. 0.0) may result in an image 
            with a lower flux density in relation to an image recovered with a
            higher robust parameter (e.g. 0.5 or 1.0), depending of the structure 
            of the source. In such cases, we attempt an image with a higher value, 
            and check if that is actually true.
            """
            if self.config.init_parameters['test_image']['uvtaper'] != ['']:
                self.run_wsclean(self.g_name,
                                 imsize=self.config.imsize,
                                 imsizey=self.config.imsizey,
                                 cell=self.config.cell_size,
                                 robust=robust, base_name=prefix,
                                 nsigma_automask=self.config.global_parameters['nsigma_automask'],
                                 nsigma_autothreshold=self.config.global_parameters[
                                     'nsigma_autothreshold'],
                                 n_interaction='0', savemodel=False, quiet=self.config.quiet,
                                 datacolumn='DATA',
                                 shift=self.config.global_parameters['FIELD_SHIFT'],
                                 with_multiscale=False, scales='0,5,20,40',
                                 nc=self.config.nc, negative_arg=self.config.negative_arg,
                                 uvtaper=self.config.init_parameters['test_image']['uvtaper'],
                                 niter=self.config.init_parameters['test_image']['niter'],
                                 PLOT=False)
            else:
                self.modified_robust = robust + 0.5
                self.run_wsclean(self.g_name,
                                 imsize=self.config.imsize,
                                 imsizey=self.config.imsizey,
                                 cell=self.config.cell_size,
                                 robust=self.modified_robust, base_name=prefix,
                                 nsigma_automask=self.config.global_parameters['nsigma_automask'],
                                 nsigma_autothreshold=self.config.global_parameters[
                                     'nsigma_autothreshold'],
                                 n_interaction='0', savemodel=False, quiet=self.config.quiet,
                                 datacolumn='DATA',
                                 shift=self.config.global_parameters['FIELD_SHIFT'],
                                 with_multiscale=False, scales='0,5,20,40',
                                 nc=self.config.nc, negative_arg=self.config.negative_arg,
                                 uvtaper=[''],
                                 niter=self.config.init_parameters['test_image']['niter'],
                                 PLOT=False)

            self.image_statistics, self.image_list = self.compute_image_stats(path=self.config.path,
                                                                              image_list=self.image_list,
                                                                              image_statistics=self.image_statistics,
                                                                              prefix=prefix,
                                                                              selfcal_step='test_image')

        self.plot_uvwave(self.g_vis, 'vis_plot_init')

    def check_init_parameters(self):
        if self.config.params_trial_2 is not None:
            self.p0_params = self.config.params_trial_2['p0']
            self.parameter_selection['test_image'] = self.config.params_trial_2
            self.print_table(self.p0_params)
        else:
            try:
                # current_total_flux = image_statistics['test_image']['total_flux_mask'] * 1000
                self.selfcal_params = self.select_parameters(
                    self.image_statistics['test_image']['total_flux_mask'] * 1000)
                self.parameter_selection['test_image'] = self.selfcal_params
                print('Initial Template of Parameters:',
                      self.parameter_selection['test_image']['name'])
                self.p0_params = self.parameter_selection['test_image']['p0']
                self.print_table(self.p0_params)
            except:
                print('No test image found. Have you run the test_image step?')

    def check_p0_parameters(self):
        if self.config.params_trial_2 is not None:
            self.parameter_selection['p0_pos'] = self.config.params_trial_2
        else:
            try:
                self.selfcal_params = self.select_parameters(
                    self.image_statistics['selfcal_test_0']['total_flux_mask'] * 1000)
                self.parameter_selection['p0_pos'] = self.selfcal_params
                print(' ++++>> Template of Parameters to be used from now on:',
                      self.parameter_selection['p0_pos']['name'])
                if self.parameter_selection['p0_pos']['p0']['combine'] == 'spw':
                    self.parameter_selection['p0_pos']['p0']['spwmap'] = self.get_spwmap(self.g_vis)
            except:
                pass

    def _run_p0(self):
        """
        Run the test image step and check the initial parameters.
        """
        # self.run_step('test_image')
        # self.check_init_parameters()
        iteration = '0'
        ############################################################################
        #### 0. Zero interaction. Use a small/negative robust parameter,        ####
        ####    to only the bright/compact emission components.                 ####
        ############################################################################

        if self.config.params_trial_2 is None and self.modified_robust is not None:
            self.p0_params['robust'] = self.modified_robust

        # minblperant = 3
        # combine='spw'
        if (self.p0_params['combine'] == 'spw' or
                self.p0_params['combine'] == 'scan,spw' or
                self.p0_params['combine'] == 'spw,scan'):
            self.p0_params['spwmap'] = self.get_spwmap(self.g_vis)

        print('Params that are currently being used:',
              self.parameter_selection['test_image']['name'])
        self.print_table(self.p0_params)

        if 'start_image' not in self.steps_performed:

            rms_mask = None  # 1 * image_statistics['test_image']['rms_box']

            # if image_statistics['test_image']['total_flux'] * 1000 > 100.0:
            """
            If the source is too bright, it may contains lots of artifacts for a robust 
            r = 0.5 (the initial test image), and those artifacts can be printed in the mask below. 
            So, we create a new test image with a lower robust parameter. 
            """
            prefix = 'test_image_0'
            self.run_wsclean(self.g_name,
                             imsize=self.config.imsize,
                             imsizey=self.config.imsizey,
                             cell=self.config.cell_size,
                             robust=self.p0_params['robust'],
                             base_name=prefix,
                             nsigma_automask=self.p0_params['nsigma_automask'],
                             nsigma_autothreshold=self.p0_params['nsigma_autothreshold'],
                             n_interaction=iteration, savemodel=False, quiet=self.config.quiet,
                             datacolumn='DATA',
                             shift=self.config.global_parameters['FIELD_SHIFT'],
                             uvtaper=self.p0_params['uvtaper'],
                             scales=self.p0_params['scales'],
                             nc=self.config.nc, negative_arg=self.config.negative_arg,
                             niter=self.config.global_parameters['niter'],
                             PLOT=False)

            self.image_statistics, self.image_list = self.compute_image_stats(path=self.config.path,
                                                                              image_list=self.image_list,
                                                                              image_statistics=self.image_statistics,
                                                                              prefix=prefix,
                                                                              sigma=self.p0_params[
                                                                                  'sigma_mask'],
                                                                              selfcal_step='p0')

            # if image_statistics['test_image_0']['inner_flux_f'] > 0.5:
            #     mask_grow_iterations = 6
            # if image_statistics['test_image_0']['inner_flux_f'] < 0.5:
            #     mask_grow_iterations = 6
            # mask_grow_iterations = self.global_parameters['mask_grow_iterations']

            mask_name = self.create_mask(self.image_list['test_image_0'],
                                         rms_mask=rms_mask,
                                         sigma_mask=self.p0_params['sigma_mask'],
                                         mask_grow_iterations=self.config.global_parameters[
                                             'mask_grow_iterations'])

            self.start_image(self.g_name, n_interaction=iteration,
                             imsize=self.config.imsize,
                             imsizey=self.config.imsizey,
                             cell=self.config.cell_size,
                             # uvtaper=['0.05arcsec'],
                             delmodel=True,
                             # opt_args=' -multiscale -multiscale-scales 0 ',
                             nsigma_automask=self.p0_params['nsigma_automask'],
                             nsigma_autothreshold=self.p0_params['nsigma_autothreshold'],
                             # next time probably needs to use 7.0 instead of 3.0
                             niter=self.config.global_parameters['niter'],
                             shift=self.config.global_parameters['FIELD_SHIFT'],
                             quiet=self.config.quiet,
                             uvtaper=self.p0_params['uvtaper'],
                             nc=self.config.nc, negative_arg=self.config.negative_arg,
                             savemodel=True, mask=mask_name, PLOT=True,
                             robust=self.p0_params['robust'], datacolumn='DATA')

            self.image_statistics, self.image_list = self.compute_image_stats(path=self.config.path,
                                                                              image_list=self.image_list,
                                                                              image_statistics=self.image_statistics,
                                                                              prefix='start_image',
                                                                              sigma=self.p0_params[
                                                                                  'sigma_mask'],
                                                                              selfcal_step='test_image')

            if 'start_image' not in self.steps_performed:
                self.steps_performed.append('start_image')

        # SNRs, percentiles_SNRs, caltable_int, caltable_3, caltable_inf = \
        #     check_solutions(g_name,
        #                     field, cut_off=p0_params['minsnr'],
        #                     minsnr=p0_params['minsnr'],
        #                     n_interaction=iteration,
        #                     solnorm=solnorm,
        #                     combine=p0_params['combine'], spwmap=p0_params['spwmap'],
        #                     calmode=p0_params['calmode'], refant=refant,
        #                     gaintype=p0_params['gaintype'],
        #                     # interp='cubic,cubic',
        #                     gain_tables_selfcal=[],
        #                     return_solution_stats=True)

        if 'select_refant' in self.config.steps:
            # #     if refant == '':
            # # if 'select_refant' in steps and 'select_refant' not in steps_performed:
            # print(' ++==> Estimating order of best referent antennas...')
            # self.tablename_refant = os.path.dirname(self.g_name) + '/selfcal/find_refant.phase'
            # self.refant = find_refant(msfile=self.g_vis, field='',
            #                         combine=self.p0_params['combine'],
            #                         tablename=tablename_refant)
            # print(' ++==> Preferential reference antenna order = ', refant)
            # # steps_performed.append('select_refant')
            # # print(f"The referent antenna is {refant}")
            # #     else:
            # #         refant = refant
            self._run_select_refant()

        if 'p0' not in self.steps_performed:
            # p0_params['solint'] = 'inf'
            # p0_params['combine'] = 'spw'
            self.gain_tables_selfcal_temp, self.spwmaps_selfcal_temp = (
                self.self_gain_cal(self.g_name,
                                   n_interaction=iteration,
                                   minsnr=self.p0_params['minsnr'],
                                   solint=self.p0_params['solint'],
                                   flagbackup=True,
                                   gaintype=self.p0_params['gaintype'],
                                   combine=self.p0_params['combine'],
                                   refant=self.refant,
                                   refantmode=self.config.refantmode,
                                   minblperant=self.config.minblperant,
                                   calmode=self.p0_params['calmode'],
                                   spwmap=self.p0_params['spwmap'],
                                   solnorm=self.config.solnorm,
                                   calwt=self.config.general_settings['calwt'],
                                   applymode=self.config.general_settings['applymode_p'],
                                   #  interp = 'cubicPD,'
                                   #           'cubicPD',
                                   #  interp='cubic,cubic',
                                   # interp='linearPD,'
                                   #        'linearflagrel',
                                   action='apply',
                                   PLOT=True,
                                   gain_tables=[]
                                   )
            )

            summary = flagdata(vis=self.g_vis, field='', mode='summary')
            self.flag_data_steps['selfcal_p0'] = summary.copy()

            self.run_wsclean(self.g_name, robust=self.p0_params['robust'],
                             imsize=self.config.imsize,
                             imsizey=self.config.imsizey,
                             cell=self.config.cell_size,
                             base_name='selfcal_test_0',
                             nsigma_automask=self.p0_params['nsigma_automask'],
                             nsigma_autothreshold=self.p0_params['nsigma_autothreshold'],
                             n_interaction='', savemodel=False, quiet=self.config.quiet,
                             with_multiscale=self.p0_params['with_multiscale'],
                             datacolumn='CORRECTED_DATA',
                             uvtaper=self.p0_params['uvtaper'],
                             scales=self.p0_params['scales'],
                             nc=self.config.nc,
                             # negative_arg=negative_arg,
                             niter=self.config.global_parameters['niter'],
                             shift=self.config.global_parameters['FIELD_SHIFT'],
                             PLOT=False)
            self.image_statistics, self.image_list = self.compute_image_stats(path=self.config.path,
                                                                              image_list=self.image_list,
                                                                              image_statistics=self.image_statistics,
                                                                              prefix='selfcal_test_0',
                                                                              sigma=8,
                                                                              selfcal_step='p0')

            if self.config.params_trial_2 is None:
                if self.image_statistics['selfcal_test_0']['total_flux_mask'] * 1000 < 10.0:
                    """
                    After the first pass of selfcal [phase], the data quality may have
                    improved. If originaly, the total flux density was below 10 mJy,
                    now, after improved corrected phases, the flux may be above 10 mJy.

                    We can check if using a taper or a higher robust will increase the
                    flux density above 10 mJy. If yes, the source will not considered as `very faint`,
                    and we may attempt a second phase-selfcal run with the template `faint` (
                    i.e. `p1` will be executed), and then p2 and ap1 later.
                    If not, we will continue with the `very faint` template and will proceed to
                    `ap1`, e.i., `p1` and `p2` are not going to be executed.
                    """
                    # modified_robust = robust + 0.5

                    print('Deconvolving image with a taper.')
                    self.run_wsclean(self.g_name,
                                     imsize=self.config.imsize,
                                     imsizey=self.config.imsizey,
                                     cell=self.config.cell_size,
                                     robust=1.0, base_name='selfcal_test_0',
                                     nsigma_automask=self.p0_params['nsigma_automask'],
                                     nsigma_autothreshold=self.p0_params['nsigma_autothreshold'],
                                     n_interaction='0', savemodel=False, quiet=self.config.quiet,
                                     datacolumn='CORRECTED_DATA',
                                     with_multiscale=self.p0_params['with_multiscale'],
                                     scales=self.p0_params['scales'],
                                     uvtaper=self.p0_params['uvtaper'],
                                     nc=self.config.nc, negative_arg=self.config.negative_arg,
                                     niter=self.config.global_parameters['niter'],
                                     shift=self.config.global_parameters['FIELD_SHIFT'],
                                     PLOT=False)
                    self.image_statistics, self.image_list = self.compute_image_stats(
                        path=self.config.path,
                        image_list=self.image_list,
                        image_statistics=self.image_statistics,
                        prefix='selfcal_test_0')

            # parameter_selection['p0_pos']['p0']['spwmap'] = p0_params['spwmap']

            self.trial_gain_tables.append(self.gain_tables_selfcal_temp)
            self.gain_tables_applied['p0'] = self.gain_tables_selfcal_temp
            self.spwmaps_applied['p0'] = self.spwmaps_selfcal_temp
            # self.steps_performed.append('p0')

    def _run_p1(self):
        iteration = '1'
        # current_total_flux = image_statistics['selfcal_test_0']['total_flux_mask'] * 1000

        ############################################################################
        #### 1. First interaction. Increase a little the robust parameter,      ####
        ####    start to consider more extended emission.                       ####
        ############################################################################
        self.p1_params = self.parameter_selection['p0_pos']['p1']
        print('Params that are currently being used:', self.parameter_selection['p0_pos']['name'])
        self.print_table(self.p1_params)
        if 'update_model_1' not in self.steps_performed:
            mask_name = self.create_mask(self.image_list['selfcal_test_0'],
                                         rms_mask=None,
                                         sigma_mask=self.p1_params['sigma_mask'],
                                         mask_grow_iterations=self.p1_params[
                                             'mask_grow_iterations'])

            self.run_wsclean(self.g_name, robust=self.p1_params['robust'],
                             imsize=self.config.imsize,
                             imsizey=self.config.imsizey,
                             cell=self.config.cell_size,
                             nsigma_automask=self.p1_params['nsigma_automask'],
                             nsigma_autothreshold=self.p1_params['nsigma_autothreshold'],
                             n_interaction=iteration, savemodel=True, quiet=self.config.quiet,
                             with_multiscale=self.p1_params['with_multiscale'],
                             scales=self.p1_params['scales'],
                             datacolumn='CORRECTED_DATA', mask=mask_name,
                             shift=self.config.global_parameters['FIELD_SHIFT'],
                             uvtaper=self.p1_params['uvtaper'],
                             nc=self.config.nc, negative_arg=self.config.negative_arg,
                             niter=self.config.global_parameters['niter'],
                             PLOT=False)

            self.image_statistics, self.image_list = self.compute_image_stats(path=self.config.path,
                                                                              image_list=self.image_list,
                                                                              image_statistics=self.image_statistics,
                                                                              sigma=self.p1_params[
                                                                                  'sigma_mask'],
                                                                              prefix='1_update_model_image')

            self.steps_performed.append('update_model_1')

        self.phase_tables = []
        self.spwmaps = []

        if self.config.params_trial_2 is not None:
            if (self.p1_params['combine'] == 'spw' or
                    self.p1_params['combine'] == 'scan,spw' or
                    self.p1_params['combine'] == 'spw,scan'):
                self.p1_spwmap = get_spwmap(self.g_vis)[0]
                self.spwmaps.append(self.p1_spwmap)
                self.p1_params['spwmap'] = self.spwmaps
            else:
                self.p1_spwmap = []
                self.spwmaps.append(self.p1_spwmap)
                self.p1_params['spwmap'] = self.spwmaps
        else:
            if self.p1_params['combine'] == 'spw':
                self.p1_params['spwmap'] = self.get_spwmap(self.g_vis)

        if 'p1' not in self.steps_performed:
            self.gain_tables_selfcal_p1, self.spwmaps_selfcal_p1 = (
                self.self_gain_cal(self.g_name,
                                   n_interaction=iteration,
                                   minsnr=self.p1_params['minsnr'],
                                   solint=self.p1_params['solint'],
                                   flagbackup=True,
                                   gaintype=self.p1_params['gaintype'],
                                   combine=self.p1_params['combine'],
                                   refant=self.refant,
                                   minblperant=self.config.minblperant,
                                   refantmode=self.config.refantmode,
                                   solnorm=self.config.solnorm,
                                   calwt=self.config.general_settings['calwt'],
                                   applymode=self.config.general_settings['applymode_p'],
                                   spwmap=self.p1_params['spwmap'],
                                   calmode=self.p1_params['calmode'],
                                   # interp='cubic,cubic',
                                   action='apply',
                                   PLOT=False,
                                   gain_tables=self.phase_tables,
                                   # spwmaps=spwmaps,
                                   # gain_tables=gain_tables_applied[
                                   #     'p0'].copy()
                                   ))
            summary = flagdata(vis=self.g_vis, field='', mode='summary')
            self.flag_data_steps['selfcal_p1'] = summary.copy()

            self.trial_gain_tables.append(self.gain_tables_selfcal_p1)
            self.gain_tables_applied['p1'] = self.gain_tables_selfcal_p1
            self.spwmaps_applied['p1'] = self.spwmaps_selfcal_p1
            # self.steps_performed.append('p1')

    def _run_p2(self):
        iteration = '2'
        ############################################################################
        #### 2. Second interaction. Increase more the robust parameter, or use  ####
        ####    uvtapering. Consider even more extended emission (if there is). ####
        ############################################################################
        # current_total_flux = image_statistics['1_update_model_image']['total_flux_mask'] * 1000
        # selfcal_params = select_parameters(current_total_flux)
        self.selfcal_params = self.parameter_selection['p0_pos']
        self.p2_params = self.selfcal_params['p2']

        self.spwmaps = self.spwmaps_applied['p1'].copy()
        self.phase_tables = self.gain_tables_applied['p1'].copy()

        if self.p2_params['combine'] == 'spw':
            self.p2_spwmap = self.get_spwmap(self.g_vis)[0]
            self.spwmaps.append(self.p2_spwmap)
            self.p2_params['spwmap'] = self.spwmaps
        else:
            self.p2_spwmap = []
            self.spwmaps.append(self.p2_spwmap)
            self.p2_params['spwmap'] = self.spwmaps

        print('Params that are currently being used:', self.parameter_selection['p0_pos']['name'])
        self.print_table(self.p2_params)

        if 'update_model_2' not in self.steps_performed:
            self.run_wsclean(self.g_name, robust=self.p2_params['robust'],
                             imsize=self.config.imsize,
                             imsizey=self.config.imsizey,
                             cell=self.config.cell_size,
                             base_name='selfcal_test_1',
                             nsigma_automask=self.p2_params['nsigma_automask'],
                             nsigma_autothreshold=self.p2_params['nsigma_autothreshold'],
                             n_interaction='', savemodel=False, quiet=self.config.quiet,
                             with_multiscale=self.p2_params['with_multiscale'],
                             scales=self.p2_params['scales'],
                             datacolumn='CORRECTED_DATA',
                             uvtaper=self.p2_params['uvtaper'],
                             nc=self.config.nc,
                             negative_arg=self.config.negative_arg,
                             shift=self.config.global_parameters['FIELD_SHIFT'],
                             niter=self.config.global_parameters['niter'],
                             PLOT=False)

            self.image_statistics, self.image_list = self.compute_image_stats(path=self.config.path,
                                                                              image_list=self.image_list,
                                                                              image_statistics=self.image_statistics,
                                                                              sigma=self.p2_params[
                                                                                  'sigma_mask'],
                                                                              prefix='selfcal_test_1')

            mask_name = self.create_mask(self.image_list['selfcal_test_1'],
                                         rms_mask=None,
                                         sigma_mask=self.p2_params['sigma_mask'],
                                         mask_grow_iterations=self.p2_params[
                                             'mask_grow_iterations'])

            self.run_wsclean(self.g_name, robust=self.p2_params['robust'],
                             imsize=self.config.imsize,
                             imsizey=self.config.imsizey,
                             cell=self.config.cell_size,
                             nsigma_automask=self.p2_params['nsigma_automask'],
                             nsigma_autothreshold=self.p2_params['nsigma_autothreshold'],
                             n_interaction=iteration, savemodel=True, quiet=self.config.quiet,
                             with_multiscale=self.p2_params['with_multiscale'],
                             scales=self.p2_params['scales'],
                             datacolumn='CORRECTED_DATA', mask=mask_name,
                             shift=self.config.global_parameters['FIELD_SHIFT'],
                             uvtaper=self.p2_params['uvtaper'],
                             nc=self.config.nc, negative_arg=self.config.negative_arg,
                             niter=self.config.global_parameters['niter'],
                             PLOT=True, with_DATA=False, with_CORRECTED=False, with_MODEL=True)

            self.image_statistics, self.image_list = self.compute_image_stats(path=self.config.path,
                                                                              image_list=self.image_list,
                                                                              image_statistics=self.image_statistics,
                                                                              sigma=self.p2_params[
                                                                                  'sigma_mask'],
                                                                              prefix='2_update_model_image')

            self.steps_performed.append('update_model_2')

        if 'p2' not in self.steps_performed:
            self.gain_tables_selfcal_p2, self.spwmaps_selfcal_p2 = (
                self.self_gain_cal(self.g_name,
                                   n_interaction=iteration,
                                   minsnr=self.p2_params['minsnr'],
                                   solint=self.p2_params['solint'],
                                   flagbackup=True,
                                   gaintype=self.p2_params['gaintype'],
                                   combine=self.p2_params['combine'],
                                   spwmap=self.p2_params['spwmap'],
                                   refant=self.refant,
                                   refantmode=self.config.refantmode,
                                   solnorm=self.config.solnorm,
                                   calwt=self.config.general_settings['calwt'],
                                   applymode=self.config.general_settings['applymode_p'],
                                   minblperant=self.config.minblperant,
                                   # interp = 'cubic,cubic',
                                   calmode=self.p2_params['calmode'],
                                   action='apply',
                                   PLOT=True,
                                   gain_tables=self.phase_tables.copy(),
                                   spwmaps=self.spwmaps.copy()
                                   )
            )
            summary = flagdata(vis=self.g_vis, field='', mode='summary')
            self.flag_data_steps['selfcal_p2'] = summary.copy()

            self.trial_gain_tables.append(self.gain_tables_selfcal_p2)
            self.gain_tables_applied['p2'] = self.gain_tables_selfcal_p2
            self.spwmaps_applied['p2'] = self.spwmaps_selfcal_p2
            # self.steps_performed.append('p2')

    def _run_ap1(self):
        iteration = '3'
        ############################################################################
        #### 3. Third interaction. Increase more the robust parameter, or use  ####
        ####    uvtapering. Consider even more extended emission (if there is). ####
        ############################################################################
        self.vis_split_name_p_statwt = f"{os.path.dirname(self.g_name)}/{self.config.field}{self.config.savename}_p.ms"

        if not os.path.exists(self.vis_split_name_p_statwt):
            print(f' ++==> Splitting phase-only self-calibrated visibility.')
            print(f'       Filename is: {self.vis_split_name_p_statwt}')
            if self.config.general_settings['new_phasecentre'] is None:
                split(vis=self.g_name + '.ms',
                      outputvis=self.vis_split_name_p_statwt,
                      datacolumn='corrected', keepflags=True)
            else:
                split(vis=self.g_name + '.ms',
                      outputvis=self.vis_split_name_p_statwt.replace('.ms', '_temp.ms'),
                      datacolumn='corrected', keepflags=True)
                print(' ++==> Phase-shiftting to original phasecentre...')
                phaseshift(vis=self.vis_split_name_p_statwt.replace('.ms', '_temp.ms'),
                           outputvis=self.vis_split_name_p_statwt,
                           phasecenter=f"J2000 {self.or_phc}"
                           )
                os.system(f"rm -r {self.vis_split_name_p_statwt.replace('.ms', '_temp.ms')}")

        if self.config.instrument == 'eM':
            minblperant = 3
        else:
            minblperant = self.config.minblperant

        if self.config.params_trial_2 is not None:
            if 'p1' not in self.gain_tables_applied:
                self.phase_tables = self.gain_tables_applied['p0'].copy()
                self.spwmaps = self.spwmaps_applied['p0'].copy()
            else:
                self.phase_tables = self.gain_tables_applied['p1'].copy()
                self.spwmaps = self.spwmaps_applied['p1'].copy()
        else:
            if 'p2' not in self.gain_tables_applied:
                if 'p1' not in self.gain_tables_applied:
                    self.phase_tables = self.gain_tables_applied['p0'].copy()
                    self.spwmaps = self.spwmaps_applied['p0'].copy()
                else:
                    self.phase_tables = self.gain_tables_applied['p1'].copy()
                    self.spwmaps = self.spwmaps_applied['p1'].copy()
            else:
                self.phase_tables = self.gain_tables_applied['p2'].copy()
                self.spwmaps = self.spwmaps_applied['p2'].copy()

        # ap1_params['spwmap'] = get_spwmap(g_vis)
        # ap1_params['spwmap'].append(get_spwmap(g_vis)[0])

        self.selfcal_params = self.parameter_selection['p0_pos']
        self.ap1_params = self.selfcal_params['ap1']
        if (self.ap1_params['combine'] == 'spw' or
                self.ap1_params['combine'] == 'scan,spw' or
                self.ap1_params['combine'] == 'spw,scan'):
            self.ap1_spwmap = self.get_spwmap(self.g_vis)[0]
            self.spwmaps.append(self.ap1_spwmap)
            self.ap1_params['spwmap'] = self.spwmaps
        else:
            self.ap1_spwmap = []
            self.spwmaps.append(self.ap1_spwmap)
            self.ap1_params['spwmap'] = self.spwmaps

        print('Params that are currently being used:', self.parameter_selection['p0_pos']['name'])
        self.print_table(self.ap1_params)

        if 'update_model_3' not in self.steps_performed:
            self.run_wsclean(self.g_name, robust=self.ap1_params['robust'],
                             imsize=self.config.imsize,
                             imsizey=self.config.imsizey,
                             cell=self.config.cell_size,
                             base_name='selfcal_test_2',
                             nsigma_automask=self.ap1_params['nsigma_automask'],
                             nsigma_autothreshold=self.ap1_params['nsigma_autothreshold'],
                             n_interaction='', savemodel=False, quiet=self.config.quiet,
                             with_multiscale=self.ap1_params['with_multiscale'],
                             scales=self.ap1_params['scales'],
                             datacolumn='CORRECTED_DATA',
                             uvtaper=self.ap1_params['uvtaper'],
                             shift=self.config.global_parameters['FIELD_SHIFT'],
                             nc=self.config.nc,
                             # negative_arg=self.config.negative_arg,
                             niter=self.config.global_parameters['niter'],
                             PLOT=False)

            self.image_statistics, self.image_list = self.compute_image_stats(path=self.config.path,
                                                                              image_list=self.image_list,
                                                                              image_statistics=self.image_statistics,
                                                                              sigma=self.ap1_params[
                                                                                  'sigma_mask'],
                                                                              prefix='selfcal_test_2')

            # mask_grow_iterations = self.ap1_params['mask_grow_iterations']

            mask_name = self.create_mask(self.image_list['selfcal_test_2'],
                                         rms_mask=None,
                                         sigma_mask=self.ap1_params['sigma_mask'],
                                         mask_grow_iterations=self.ap1_params[
                                             'mask_grow_iterations'])

            self.run_wsclean(self.g_name, robust=self.ap1_params['robust'],
                             imsize=self.config.imsize,
                             imsizey=self.config.imsizey,
                             cell=self.config.cell_size,
                             nsigma_automask=self.ap1_params['nsigma_automask'],
                             nsigma_autothreshold=self.ap1_params['nsigma_autothreshold'],
                             n_interaction=iteration, savemodel=True, quiet=self.config.quiet,
                             with_multiscale=self.ap1_params['with_multiscale'],
                             scales=self.ap1_params['scales'],
                             datacolumn='CORRECTED_DATA', mask=mask_name,
                             shift=self.config.global_parameters['FIELD_SHIFT'],
                             uvtaper=self.ap1_params['uvtaper'],
                             nc=self.config.nc, negative_arg=self.config.negative_arg,
                             niter=self.config.global_parameters['niter'],
                             PLOT=False)

            self.image_statistics, self.image_list = self.compute_image_stats(path=self.config.path,
                                                                              image_list=self.image_list,
                                                                              image_statistics=self.image_statistics,
                                                                              sigma=self.ap1_params[
                                                                                  'sigma_mask'],
                                                                              prefix='3_update_model_image')

            self.steps_performed.append('update_model_3')

        if 'ap1' not in self.steps_performed:
            self.gain_tables_selfcal_ap1, self.spwmaps_selfcal_ap1 = (
                self.self_gain_cal(self.g_name,
                                   n_interaction=iteration,
                                   minsnr=self.ap1_params['minsnr'],
                                   solint=self.ap1_params['solint'],
                                   flagbackup=True,
                                   gaintype=self.ap1_params['gaintype'],
                                   combine=self.ap1_params['combine'],
                                   refant=self.refant,
                                   minblperant=minblperant,
                                   solnorm=self.config.solnorm,
                                   calwt=self.config.general_settings['calwt_ap'],
                                   applymode=self.config.general_settings['applymode_ap'],
                                   refantmode=self.config.refantmode,
                                   spwmap=self.ap1_params['spwmap'],
                                   # interp='cubicPD,'
                                   #        'cubicPD',
                                   # interp = 'cubic,cubic',
                                   calmode=self.ap1_params['calmode'],
                                   action='apply',
                                   PLOT=True,
                                   gain_tables=self.phase_tables.copy(),
                                   spwmaps=self.spwmaps.copy()
                                   )
            )
            summary = flagdata(vis=self.g_vis, field='', mode='summary')
            self.flag_data_steps['selfcal_ap1'] = summary.copy()

            self.trial_gain_tables.append(self.gain_tables_selfcal_ap1)
            self.gain_tables_applied['ap1'] = self.gain_tables_selfcal_ap1
            self.spwmaps_applied['ap1'] = self.spwmaps_selfcal_ap1
            # self.steps_performed.append('ap1')

    def _run_split_trial_1(self):
        self.vis_split_name_1_statwt = f"{os.path.dirname(self.g_name)}/{self.config.field}{self.config.savename}.ms"

        if not os.path.exists(self.vis_split_name_1_statwt):
            print(f' ++==> Splitting final self-calibrated (p+ap) visibility.')
            print(f'       Filename is: {self.vis_split_name_1_statwt}')
            if self.config.general_settings['new_phasecentre'] is None:
                split(vis=self.g_name + '.ms',
                      outputvis=self.vis_split_name_1_statwt,
                      datacolumn='corrected', keepflags=True)
                print(' ++==> Running statw on split data...')
                statwt(vis=self.vis_split_name_1_statwt,
                       statalg=self.config.general_settings['statwt_statalg'],
                       timebin=self.config.general_settings['timebin_statw'],
                       datacolumn='data')
            else:
                split(vis=self.g_name + '.ms',
                      outputvis=self.vis_split_name_1_statwt.replace('.ms', '_temp.ms'),
                      datacolumn='corrected', keepflags=True)
                print(' ++==> Phase-shiftting to original phasecentre...')
                phaseshift(vis=self.vis_split_name_1_statwt.replace('.ms', '_temp.ms'),
                           outputvis=self.vis_split_name_1_statwt,
                           phasecenter=f"J2000 {self.or_phc}"
                           )
                print(' ++==> Running statw on split data...')
                statwt(vis=self.vis_split_name_1_statwt,
                       statalg=self.config.general_settings['statwt_statalg'],
                       timebin=self.config.general_settings['timebin_statw'],
                       datacolumn='data')
                os.system(f"rm -r {self.vis_split_name_1_statwt.replace('.ms', '_temp.ms')}")

        if os.path.exists(self.vis_split_name_p_statwt):
            print(' ++==> Running statw on split phase-only self-calibrated visibility...')
            statwt(vis=self.vis_split_name_p_statwt,
                   statalg=self.config.general_settings['statwt_statalg'],
                   timebin=self.config.general_settings['timebin_statw'],
                   datacolumn='data')

        print(' ++==> Imaging visibilities after self-calibration...')

        niter = self.config.global_parameters['niter']

        if self.config.params_trial_2 is not None:
            ROBUSTS = [-0.75, -0.5, -0.25, 0.5]
        else:
            ROBUSTS = [0.5]

        # ROBUSTS = [0.5]
        for robust in ROBUSTS:
            # split_list = [vis_split_name_1, vis_split_name_1_statwt]
            self.split_list = [self.vis_split_name_1_statwt]
            for vis_split in self.split_list:
                try:
                    self.run_wsclean(vis_split, robust=robust,
                                     imsize=self.config.imsize,
                                     imsizey=self.config.imsizey,
                                     cell=self.config.cell_size,
                                     base_name='selfcal_image',
                                     nsigma_automask=self.config.global_parameters[
                                         'nsigma_automask'],
                                     nsigma_autothreshold=self.config.global_parameters[
                                         'nsigma_autothreshold'],
                                     n_interaction='', savemodel=False, quiet=self.config.quiet,
                                     with_multiscale=self.config.global_parameters[
                                         'with_multiscale'],
                                     scales=self.config.global_parameters['scales'],
                                     datacolumn='DATA',
                                     shift=self.config.global_parameters['FIELD_SHIFT'],
                                     uvtaper=self.config.global_parameters['uvtaper'],
                                     nc=self.config.nc,
                                     negative_arg='negative',
                                     niter=niter,
                                     PLOT=False)

                    self.image_statistics, self.image_list = self.compute_image_stats(
                        path=self.config.path,
                        image_list=self.image_list,
                        image_statistics=self.image_statistics,
                        prefix='selfcal_image')
                except:
                    pass

        # self.plot_uvwave(self.vis_split_name_1_statwt,'vis_plot_final')

    def _run_autoflag_final(self):
        # run_autoflag(g_vis, display='report', action='apply',mode='rflag',
        #           timedevscale=2.5, freqdevscale=2.5, winsize=5,
        #           datacolumn='corrected')

        self.vis_split_name_flag = self.vis_split_name_1_statwt.replace('.ms', '_post_flag.ms')

        # self.plot_visibilities(g_vis=self.g_vis,
        #                     name='selfcal_final_before_post_flag',
        #                     with_MODEL=False, with_CORRECTED=False,
        #                     with_DATA=False,with_RESIDUAL=True)

        self.run_autoflag(self.g_vis,
                          display='report', action='apply', mode='rflag',
                          timedevscale=3.5, freqdevscale=3.5, winsize=5,
                          datacolumn='residual')

        summary = flagdata(vis=self.g_vis, field='', mode='summary')
        self.flag_data_steps['autoflag_final'] = summary.copy()

        self.plot_visibilities(g_vis=self.g_vis,
                               name='selfcal_final_after_post_flag',
                               with_MODEL=False, with_CORRECTED=False,
                               with_DATA=False, with_RESIDUAL=True)

        if not os.path.exists(self.vis_split_name_flag):
            print(' ++==> Splitting data for final auto-flagging...')
            if self.config.general_settings['new_phasecentre'] is None:
                split(vis=self.g_vis,
                      outputvis=self.vis_split_name_flag,
                      datacolumn='corrected',
                      keepflags=True)
            else:
                split(vis=self.g_vis,
                      outputvis=self.vis_split_name_flag.replace('.ms', '_temp.ms'),
                      datacolumn='corrected', keepflags=True)
                print(
                    ' ++==> Phase-shiftting final (flagged) visibility to original phasecentre...')
                phaseshift(vis=self.vis_split_name_flag.replace('.ms', '_temp.ms'),
                           outputvis=self.vis_split_name_1_statwt,
                           phasecenter=f"J2000 {self.or_phc}"
                           )
                os.system(f"rm -r {self.vis_split_name_flag.replace('.ms', '_temp.ms')}")

            statwt(vis=self.vis_split_name_flag,
                   statalg=self.config.general_settings['statwt_statalg'],
                   timebin=self.config.general_settings['timebin_statw'],
                   datacolumn='data')

        self.plot_visibilities(g_vis=self.vis_split_name_flag,
                               name='selfcal_final_post_flag',
                               with_MODEL=False, with_CORRECTED=True,
                               with_DATA=False, with_RESIDUAL=False)

        niter = self.config.global_parameters['niter']
        robust = 0.5
        self.run_wsclean(self.vis_split_name_flag, robust=robust,
                         imsize=self.config.imsize,
                         imsizey=self.config.imsizey,
                         cell=self.config.cell_size,
                         base_name='selfcal_image_post_flag',
                         nsigma_automask='3.0',
                         nsigma_autothreshold='1.0',
                         n_interaction='', savemodel=False, quiet=True,
                         with_multiscale=True,
                         scales='None',
                         datacolumn='DATA',
                         shift=self.config.global_parameters['FIELD_SHIFT'],
                         uvtaper=self.config.global_parameters['uvtaper'],
                         nc=self.config.nc,
                         negative_arg='negative',
                         niter=niter,
                         PLOT=False)

        self.image_statistics, self.image_list = self.compute_image_stats(path=self.config.path,
                                                                          image_list=self.image_list,
                                                                          image_statistics=self.image_statistics,
                                                                          prefix='selfcal_image_post_flag')

        # steps_performed.append('run_autoflag_final')

    def _run_report_results(self):
        """
        To do: save and plot the results of the selfcalibration.
        """
        self.snr = [self.image_statistics[image]['snr'] for image in self.image_statistics]

        self.df = pd.DataFrame.from_dict(self.image_statistics, orient='index')

        self.df.to_csv(self.g_name + '_selfcal_statistics.csv',
                       header=True,
                       index=False)
        self.df_gt = pd.DataFrame.from_dict(self.gain_tables_applied, orient='index')
        self.df_gt.to_csv(self.g_name + '_tables_applied.csv',
                          header=True,
                          index=False)

        # try:
        # df_selfcal[['DR_SNR_E', 'DR_pk_rmsbox', 'snr', 'peak_of_flux', 'total_flux']]
        try:
            self.df_new = self.df[
                ['snr', 'peak_of_flux', 'total_flux', 'rms_box', 'DR_pk_rmsim']]
            self.df_new = self.df_new.T
            self.df_new = self.df_new[
                ['start_image', '1_update_model_image', '3_update_model_image', 'selfcal_image']]
            for col in self.df_new.columns:
                self.df_new[col] = self.df_new[col].map("{:.2e}".format)
        except:
            self.df_new = self.df[
                ['snr', 'peak_of_flux', 'total_flux', 'rms_box', 'DR_pk_rmsim']]
            self.df_new = self.df_new.T
            self.df_new = self.df_new[
                ['start_image', 'selfcal_test_0', '3_update_model_image',
                 'selfcal_image']]
            for col in self.df_new.columns:
                self.df_new[col] = self.df_new[col].map("{:.2e}".format)
        fig = plt.figure(figsize=(12, 10))
        gs = gridspec.GridSpec(5, 6)

        ax_snr = plt.subplot(gs[0, 0:2])
        ax_snr.scatter(self.df['snr'], self.df['peak_of_flux'], label='SNR vs Sp')
        ax_snr.set_title('SNR vs Sp')

        ax_flux = plt.subplot(gs[1, 0:2])
        ax_flux.scatter(self.df['total_flux'] * 1000, self.df['peak_of_flux'], label='Flux vs Sp')
        ax_flux.set_title('Flux vs Sp')

        ax_dr_rms = plt.subplot(gs[2, 0:2])
        ax_dr_rms.scatter(self.df['DR_pk_rmsbox'], self.df['rms_box'] * 1000, label='DR vs RMS')
        ax_dr_rms.set_title('DR vs RMS')

        ax_table = plt.subplot(gs[0:3, 2:6])
        ax_table.axis('off')  # Hide the axis
        # table_data = df_.set_index('Statistics').T

        table = ax_table.table(
            cellText=self.df_new.values,
            colLabels=self.df_new.columns,
            rowLabels=self.df_new.index,
            loc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(8)
        table.scale(0.5, 0.9)

        def plot_stage_images(image, model, residual, axes_index, fig, gs, phase, rms,
                              centre,
                              crop=True):
            ax_image = plt.subplot(gs[3:5, axes_index:axes_index + 2])
            # phase = 'AP'
            ax_image = mlibs.eimshow(image,
                                     rms=None, add_contours=True, center=centre,
                                     ax=ax_image, fig=fig,
                                     crop=crop, box_size=400)
            ax_image = ax_image
            # ax_image.imshow(image_placeholder)
            ax_image.set_title(f'{phase}_IMAGE')
            # ax_image.axis('off')

        self.rms = mlibs.mad_std(mlibs.ctn(self.image_list['selfcal_image_residual']))
        self.centre = mlibs.nd.maximum_position(mlibs.ctn(self.image_list['selfcal_image']))[::-1]

        plot_stage_images(image=self.image_list['selfcal_test_0'],
                          model=self.image_list['selfcal_test_0_model'],
                          residual=self.image_list['selfcal_test_0_residual'],
                          axes_index=0,
                          fig=fig,
                          gs=gs,
                          phase='O',
                          rms=self.rms,
                          centre=self.centre
                          )
        try:
            plot_stage_images(image=self.image_list['2_update_model_image'],
                              model=self.image_list['2_update_model_image_model'],
                              residual=self.image_list['2_update_model_image_residual'],
                              axes_index=2,
                              fig=fig,
                              gs=gs,
                              phase='P',
                              rms=self.rms,
                              centre=self.centre
                              )
        except:
            plot_stage_images(image=self.image_list['3_update_model_image'],
                              model=self.image_list['3_update_model_image_model'],
                              residual=self.image_list['3_update_model_image_residual'],
                              axes_index=2,
                              fig=fig,
                              gs=gs,
                              phase='P',
                              rms=self.rms,
                              centre=self.centre
                              )
        plot_stage_images(image=self.image_list['selfcal_image'],
                          model=self.image_list['selfcal_image_model'],
                          residual=self.image_list['selfcal_image_residual'],
                          axes_index=4,
                          fig=fig,
                          gs=gs,
                          phase='AP',
                          rms=self.rms,
                          centre=self.centre
                          )

        plt.subplots_adjust(wspace=0.1, hspace=0.1)
        plt.tight_layout()
        plt.savefig(self.g_name + '_selfcal_results.pdf', dpi=300, bbox_inches='tight')
        plt.clf()
        plt.close()
        # plt.show()
        # except:
        #     pass

        """
        Report flagging statitstics
        """

        # Accumulate percentages with labels
        percentages_over_steps = {}
        for label, step in self.flag_data_steps.items():
            percentages_over_steps[label] = self.calculate_percentages(step)
        # Plotting data for each category
        # categories = ['antenna', 'array', 'correlation', 'field', 'observation', 'scan', 'spw']
        categories = ['field', 'observation', 'spw', 'correlation', 'antenna']
        for category in categories:
            self.plot_category_data(percentages_over_steps, category)

        pass

    # def run_pipeline(self):
    #     """Execute all configured pipeline steps"""
    #     # for step in self.config.steps:
    #     config = Configuration()
    #     if 'startup' in self.config.steps:
    #         self.run_step('startup')
    #     if config.cell_size is None:
    #         config.cell_size = self.get_cell_size()

    #         # self.run_step(step)


if __name__ == '__main__':
    config = Configuration()
    pipeline = Pipeline(config)
    if 'startup' in config.steps:
        pipeline.run_step('startup')

    config.cell_size = None
    config.receiver = None

    if config.cell_size is None:
        config.cell_size = pipeline.get_cell_size(pipeline.g_vis)
    if config.receiver is None:
        config.receiver = pipeline.check_band(pipeline.g_vis)
    # if config.instrument is None:
    #     config.instrument = pipeline.get_instrument(pipeline.g_vis)
    if 'save_init_flags' in config.steps:
        pipeline.run_step('save_init_flags')

    if 'statwt' in config.steps:
        pipeline.run_step('statwt')

    if 'autoflag_init' in config.steps:
        pipeline.run_step('autoflag_init')

    if 'test_image' in config.steps:
        pipeline.run_step('test_image')
    pipeline.check_init_parameters()

    if 'select_refant' in config.steps:
        pipeline.run_step('select_refant')

    if 'p0' in config.steps:
        pipeline.run_step('p0')
    pipeline.check_p0_parameters()

    if (('p1' in config.steps) and ('p1' not in pipeline.steps_performed) and
            ('p1' in pipeline.parameter_selection['p0_pos'])):
        pipeline.run_step('p1')

    if (('p2' in config.steps)
            and ('p2' not in pipeline.steps_performed)
            and ('p2' in pipeline.parameter_selection['p0_pos'])):
        pipeline.run_step('p2')

    if (('ap1' in config.steps) and
            ('ap1' not in pipeline.steps_performed) and
            ('ap1' in pipeline.parameter_selection['p0_pos'].keys())):
        pipeline.run_step('ap1')

    if 'split_trial_1' in config.steps:
        pipeline.run_step('split_trial_1')

    if 'autoflag_final' in config.steps:
        pipeline.run_step('autoflag_final')

    if 'report_results' in config.steps:
        pipeline.run_step('report_results')

    # exit()