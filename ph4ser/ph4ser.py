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
import shutil
sys.path.append('./')
import mlibs as mlibs
import glob
import pandas as pd
from casatasks import *
import numpy as np
np.set_printoptions(legacy='1.25')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit
try:
    # import casatools
    # from casatasks import *
    import casatasks, casatools, casaplotms
    from casaplotms import plotms
    from casatasks import gaincal, phaseshift, applycal, bandpass, split, importasdm, statwt, mstransform
    from casatasks import flagdata, flagmanager, plotants, plotweather, hanningsmooth
    from casatasks import clearcal, delmod, setjy, fluxscale, tclean, imhead, listobs
except:
    print('Not importing casatools. '
          'Maybe you are inside inside CASA? '
          'Or check your modular installation?')
    pass
from casaplotms import plotms
from casaviewer.imview import imview

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
from plot_vis_python import plot_visibilities, plot_uvwave
# from vis_data_plots import plot_uvwave


from matplotlib import use as mpluse
mpluse('Agg')

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


def get_nspw(vis):
    """
    Total number of rows in the SPECTRAL_WINDOW table of `vis`.

    This counts *every* SPW described by the metadata, including leftovers
    that no longer carry any data (see get_science_spwids()).

    Args:
        vis: Input visibility file (.ms).

    Returns:
        int: Number of SPWs in the SPECTRAL_WINDOW sub-table.
    """
    msmd.open(vis)
    try:
        return int(msmd.nspw())
    finally:
        msmd.done()


def get_science_spwids(vis, drop_empty=True, verbose=False):
    """
    SPW ids that actually carry data in `vis`.

    A measurement set can describe more SPWs in its SPECTRAL_WINDOW table than
    it really contains data for: leftovers from a previous split/mstransform,
    ALMA WVR and channel-averaged windows, or zero-channel rows. Those must be
    excluded, otherwise any per-SPW quantity (e.g. the `chanbin` map used for
    frequency averaging) ends up misaligned with the SPW selection that the
    pipeline actually passes to CASA.

    The selection is the union of the SPWs referenced by the scans of the MS,
    minus ALMA WVR/channel-average windows and (optionally) empty SPWs. If the
    metadata exposes no scans at all, every SPW is kept so that this never
    returns an empty selection.

    Args:
        vis: Input visibility file (.ms).
        drop_empty: If True, discard SPWs with zero channels.
        verbose: If True, report which SPWs were discarded and why.

    Returns:
        list: Sorted list of int SPW ids that carry data.
    """
    msmd.open(vis)
    try:
        all_ids = list(range(int(msmd.nspw())))

        nchan = {}
        for spwid in all_ids:
            try:
                nchan[spwid] = int(len(msmd.chanfreqs(spwid)))
            except Exception:
                # SPW described in the metadata but not readable -> treat as empty
                nchan[spwid] = 0

        # SPWs referenced by at least one scan (same criterion as listobs)
        used = set()
        try:
            for scan in msmd.scannumbers():
                try:
                    used.update(int(spwid) for spwid in msmd.spwsforscan(int(scan)))
                except Exception:
                    continue
        except Exception:
            used = set()

        # ALMA-specific auxiliary windows (WVR / channel-averaged)
        aux = set()
        try:
            aux = set(int(spwid) for spwid in msmd.almaspws(wvr=True, chavg=True))
        except Exception:
            aux = set()
    finally:
        msmd.done()

    if not used:
        # metadata does not expose scans (or the query failed): keep everything
        used = set(all_ids)

    selected = sorted(spwid for spwid in used
                      if spwid in nchan
                      and spwid not in aux
                      and (nchan[spwid] > 0 or not drop_empty))

    if not selected:
        # never hand back an empty selection; fall back to all described SPWs
        selected = sorted(all_ids)

    if verbose:
        discarded = sorted(set(all_ids) - set(selected))
        print(f'++==> SPW selection for {os.path.basename(vis)}: '
              f'{len(selected)} of {len(all_ids)} SPWs carry data.')
        if discarded:
            reasons = []
            for spwid in discarded:
                if spwid in aux:
                    reasons.append(f'{spwid} (ALMA WVR/chavg)')
                elif nchan.get(spwid, 0) == 0:
                    reasons.append(f'{spwid} (no channels)')
                else:
                    reasons.append(f'{spwid} (no scans)')
            print(f'     --==> Discarded SPWs: {", ".join(reasons)}')

    return selected


def get_spw_groups(vis, spwids=None, verbose=False):
    """
    Groups of spectral windows that are observed together.

    Each scan of an MS records the set of SPWs it covers; the distinct sets are
    the groups over which solutions can be combined (`combine='spw'`). The
    information is taken from listobs, with the msmd metadata as a fallback if
    listobs fails or reports no scans, mirroring Pipeline.get_spwids().

    Groups are intersected with `spwids`, so SPWs that carry no data (leftovers
    in the SPECTRAL_WINDOW table, ALMA WVR windows, ...) never appear.

    Args:
        vis: Input visibility file (.ms).
        spwids: SPWs to keep. None (default) uses get_science_spwids().
        verbose: If True, print the groups found.

    Returns:
        list: List of sorted lists of int SPW ids; never empty.
    """
    if spwids is None:
        spwids = get_science_spwids(vis)
    spwids = sorted(set(int(spwid) for spwid in np.atleast_1d(spwids)))

    groups = set()
    try:
        lobs = listobs(vis=vis)
        for key in (k for k in lobs if 'scan_' in k):
            for inner_key in lobs[key]:
                groups.add(tuple(sorted(int(spwid) for spwid
                                        in lobs[key][inner_key]['SpwIds'])))
    except Exception as exc:
        print(f'[get_spw_groups] listobs failed on {vis} ({exc}); '
              'falling back to msmd metadata.')

    if not groups:
        msmd.open(vis)
        try:
            for scan in msmd.scannumbers():
                try:
                    groups.add(tuple(sorted(int(spwid) for spwid
                                            in msmd.spwsforscan(int(scan)))))
                except Exception:
                    continue
        except Exception:
            groups = set()
        finally:
            msmd.done()

    # drop SPWs that carry no data, then drop groups left empty
    keep = set(spwids)
    grouped = sorted({tuple(spwid for spwid in group if spwid in keep)
                      for group in groups} - {()})

    if not grouped:
        # no usable scan information: treat the whole selection as one group
        grouped = [tuple(spwids)]

    grouped = [list(group) for group in grouped]

    if verbose:
        print(f'++==> SPW groups for {os.path.basename(vis)}: {grouped}')

    return grouped


def get_spwmap(vis, spwids=None, verbose=False):
    """
    Build the `spwmap` list used by gaincal/applycal when solutions are
    combined across spectral windows (``combine='spw'``).

    CASA indexes spwmap by SPW *id*: element i names the SPW whose solutions
    are applied to SPW i. The list must therefore be long enough to cover the
    highest SPW id in use -- it is *not* one entry per selected SPW. Every SPW
    in a group is mapped onto the lowest id of that group (where the combined
    solution is stored); SPWs outside the selection are mapped onto themselves,
    which is a no-op, so leftover metadata cannot shift the mapping.

    Args:
        vis: Input visibility file (.ms).
        spwids: SPWs that carry data. None (default) uses get_science_spwids().
        verbose: If True, print the groups and the resulting map.

    Returns:
        list: ``[spwmap_i]`` -- a single-element list holding the map, matching
        the per-gaintable nesting expected by gaincal/applycal.
    """
    if spwids is None:
        spwids = get_science_spwids(vis)
    spwids = sorted(set(int(spwid) for spwid in np.atleast_1d(spwids)))

    if not spwids:
        raise ValueError(f'get_spwmap: empty SPW selection for {vis}.')

    groups = get_spw_groups(vis, spwids=spwids, verbose=verbose)

    reference = {}
    for group in groups:
        ref = min(group)
        for spwid in group:
            # deterministic if groups happen to overlap: lowest reference wins
            reference[spwid] = min(reference.get(spwid, ref), ref)

    # identity for every id up to the highest one in use, then apply the groups
    spwmap_i = list(range(max(spwids) + 1))
    for spwid, ref in reference.items():
        spwmap_i[spwid] = ref

    if verbose:
        print(f'     ==> spwmap = {spwmap_i}')

    return [spwmap_i]


def get_chan_avg_map(vis, chan_out_avg=64, spwids=None, verbose=True,
                     return_spwids=False):
    """
    Compute per-SPW channel binning factors which average each spectral window
    down to (at most) `chan_out_avg` output channels.

    The returned list is meant to be assigned to
    ``general_settings['channel_width']``, which is forwarded to the `chanbin`
    parameter of mstransform inside Pipeline._prepare_visibility(). It must
    therefore have exactly one entry per SPW *selected* for the transform, in
    ascending SPW-id order -- not one entry per row of the SPECTRAL_WINDOW
    table, which may contain leftovers that carry no data. By default the SPW
    selection is taken from get_science_spwids(), which applies the same
    criterion (SPWs referenced by scans) used by Pipeline.get_spwids().

    Args:
        vis: Input visibility file (.ms).
        chan_out_avg: Desired number of output channels per SPW.
        spwids: SPWs to map. None (default) selects the SPWs that carry data;
            'all' forces every SPW in the SPECTRAL_WINDOW table; a list/array
            of ids restricts the map to those SPWs (e.g. the output of
            ``pipeline.get_spwids(vis)``).
        verbose: If True, print a short summary of the averaging map.
        return_spwids: If True, return ``(chan_width_avg, spwids)`` so the map
            can be checked against the SPW selection used downstream.

    Returns:
        list: One integer binning factor (>= 1) per selected SPW, ordered by
        SPW id; or the tuple ``(chan_width_avg, spwids)`` if `return_spwids`.

    Usage:
        >>> config.general_settings['channel_width'] = get_chan_avg_map(vis, chan_out_avg=128)
        >>> # equivalently, from an instantiated pipeline (uses pipeline.get_spwids):
        >>> config.general_settings['channel_width'] = pipeline.get_chan_avg_map(vis, chan_out_avg=128)
    """
    if spwids is None:
        spwids = get_science_spwids(vis, verbose=verbose)
    elif isinstance(spwids, str) and spwids.lower() == 'all':
        spwids = list(range(get_nspw(vis)))
    elif isinstance(spwids, str):
        spwids = [int(spwid) for spwid in spwids.split(',') if spwid.strip() != '']
    else:
        spwids = [int(spwid) for spwid in np.atleast_1d(spwids)]
    spwids = sorted(set(int(spwid) for spwid in spwids))

    if len(spwids) == 0:
        raise ValueError(f'get_chan_avg_map: empty SPW selection for {vis}.')

    msmd.open(vis)
    try:
        nspw_total = int(msmd.nspw())
        bad = [spwid for spwid in spwids if spwid < 0 or spwid >= nspw_total]
        if bad:
            raise ValueError(f'get_chan_avg_map: SPW id(s) {bad} not present in '
                             f'{vis} (it describes {nspw_total} SPWs).')

        chan_freqs_all = np.empty(len(spwids), dtype=object)
        spws_freq = np.zeros(len(spwids))
        for i, spwid in enumerate(spwids):
            chan_freqs_all[i] = msmd.chanfreqs(spwid)
            spws_freq[i] = np.nanmean(chan_freqs_all[i])
    finally:
        msmd.done()

    mean_freq = np.nanmean(spws_freq) * 1e-9
    nchan_per_spw = np.asarray([arr.shape[0] for arr in chan_freqs_all])
    chan_width_avg = [max(1, int(nchan / chan_out_avg)) for nchan in nchan_per_spw]

    if verbose:
        nchan_out = [int(nchan // width)
                     for nchan, width in zip(nchan_per_spw, chan_width_avg)]
        print(f'++==> Channel averaging map for {os.path.basename(vis)}:')
        print(f'     ==> {len(spwids)} SPWs selected {spwids}, '
              f'mean frequency = {mean_freq:.3f} GHz.')
        print(f'     ==> Input channels per SPW  = {list(map(int, nchan_per_spw))}')
        print(f'     ==> Channel bin widths      = {chan_width_avg}')
        print(f'     ==> Output channels per SPW = {nchan_out}')

    if return_spwids:
        return chan_width_avg, spwids
    return chan_width_avg


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
        self.plotting_verbosity = cf.plotting_verbosity
        self.show_figures = cf.show_figures

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
        self.multi_config = cf.multi_config
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
            # 'text.usetex' : True,
            # 'font.family' : 'serif',
            # 'font.serif' : ['Garamond Libre', 'EB Garamond', 'Cormorant Garamond', 'serif'],
            # 'text.latex.preamble': r'''
            # \usepackage{ebgaramond-maths}
            # \usepackage{garamondlibre}
            # \usepackage{amsmath}
            # \usepackage{amssymb}
            # ''',
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

    plot_visibilities = plot_visibilities
    plot_uvwave       = plot_uvwave
    
    def select_parameters(self, total_flux, snr=None):
        """
        Select appropriate parameters based on total flux
        These are empirical references values based on extensive testing 
        with VLA (L,S,C,X,Ku,K,Ka and Q bands) and e-MERLIN (L and C bands).
        """
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

    # def select_parameters(self, total_flux, snr=None):
    #     """Select appropriate parameters based on total flux"""
    #     if total_flux < 10:
    #         params = self.config.params_very_faint.copy()
    #     elif 10 <= total_flux < 20:
    #         params = self.config.params_faint.copy()
    #     elif 20 <= total_flux < 50:
    #         params = self.config.params_standard_1.copy()
    #     elif 50 <= total_flux < 100:
    #         params = self.config.params_standard_2.copy()
    #     else:
    #         params = self.config.params_bright.copy()
    #     return params

    # def get_spwids(self, vis, return_string=False):
    #     """Get spectral window IDs from visibility data"""
    #     lobs = listobs(vis=vis)
    #     extract_spwids = {key: lobs[key] for key in lobs if 'scan_' in key}

    #     unique_spwids = set()
    #     for key in extract_spwids:
    #         nested_dict = extract_spwids[key]
    #         for inner_key in nested_dict:
    #             spwids = nested_dict[inner_key]['SpwIds']
    #             unique_spwids.add(tuple(sorted(spwids)))

    #     unique_spwids_lists = sorted([list(t) for t in unique_spwids])
    #     unique_elements = set(element for sublist in unique_spwids_lists
    #                           for element in sublist)
    #     if return_string:
    #         return ','.join(map(str, sorted(list(unique_elements))))
    #     else:
    #         return sorted(list(unique_elements))

    # def get_spwmap(self, vis):
    #     """Generate spectral window mapping"""
    #     unique_spwids_lists = self._get_unique_spwids(vis)
    #     counts = self._count_spw_occurrences(unique_spwids_lists)
    #     spwmap_i = [item for item, count in counts.items()
    #                for _ in range(count)]
    #     return [spwmap_i[:len(self.get_spwids(vis))]]

    def get_spwmap(self, vis, spwids=None, verbose=False):
        """
        spwmap for gaincal/applycal when solutions are combined across
        spectral windows; see the module-level get_spwmap().

        The SPW selection defaults to ``self.get_spwids(vis)``, the same
        selection used everywhere else in the pipeline, so SPWs described by
        the metadata but carrying no data cannot shift the mapping.

        Args:
            vis: Input visibility file (.ms).
            spwids: SPWs that carry data; None (default) uses self.get_spwids(vis).
            verbose: If True, print the groups and the resulting map.

        Returns:
            list: ``[spwmap_i]``, indexed by SPW id.
        """
        if spwids is None:
            spwids = self.get_spwids(vis)
        return get_spwmap(vis, spwids=spwids, verbose=verbose)

    def get_spw_groups(self, vis, spwids=None, verbose=False):
        """
        Groups of SPWs observed together; see the module-level
        get_spw_groups(). Defaults to the ``self.get_spwids(vis)`` selection.

        Args:
            vis: Input visibility file (.ms).
            spwids: SPWs to keep; None (default) uses self.get_spwids(vis).
            verbose: If True, print the groups found.

        Returns:
            list: List of sorted lists of int SPW ids.
        """
        if spwids is None:
            spwids = self.get_spwids(vis)
        return get_spw_groups(vis, spwids=spwids, verbose=verbose)

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

    def get_chan_avg_map(self, vis, chan_out_avg=64, spwids=None, verbose=True,
                         return_spwids=False):
        """
        Per-SPW channel binning factors to average down to `chan_out_avg`
        output channels; see the module-level get_chan_avg_map().

        The SPW selection defaults to ``self.get_spwids(vis)`` -- exactly the
        selection that _prepare_visibility() passes to mstransform -- so the
        returned map is guaranteed to line up with `chanbin` even when the MS
        metadata still describes leftover SPWs that carry no data.

        This is a thin wrapper so that the same helper is reachable both as
        ``ph4ser.get_chan_avg_map(vis, ...)`` and as
        ``pipeline.get_chan_avg_map(vis, ...)``.

        Args:
            vis: Input visibility file (.ms).
            chan_out_avg: Desired number of output channels per SPW.
            spwids: SPWs to map; None (default) uses self.get_spwids(vis).
            verbose: If True, print a short summary of the averaging map.
            return_spwids: If True, return ``(chan_width_avg, spwids)``.

        Returns:
            list: One integer binning factor (>= 1) per selected SPW, ordered
            by SPW id; or ``(chan_width_avg, spwids)`` if `return_spwids`.
        """
        if spwids is None:
            spwids = self.get_spwids(vis)
        return get_chan_avg_map(vis, chan_out_avg=chan_out_avg, spwids=spwids,
                                verbose=verbose, return_spwids=return_spwids)

    def get_spwids(self, vis, return_string=False):
        """
        Get unique spectral window IDs from visibility data.

        Args:
            vis: Input visibility file

        Returns:
            numpy.ndarray: Array of unique spectral window IDs
        """
        unique_elements = set()
        try:
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
        except Exception as exc:
            print(f'[get_spwids] listobs failed on {vis} ({exc}); '
                  'falling back to msmd metadata.')

        if not unique_elements:
            # listobs reported no scans (or failed): use the msmd-based selection
            unique_elements = set(get_science_spwids(vis))

        if return_string:
            return ','.join(map(str, sorted(list(unique_elements))))
        else:
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

            sub_sampling_factor = 5
            # sub_sampling_factor = 10 #test for Arp200
            cell_float = (180.0 * 3600 / (np.pi * sub_sampling_factor)) * (1.0 / longest_baseline_lambda)
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
        
        last_level = 3
        
        if sigma is None:
            if (selfcal_step == 'test_image' or selfcal_step == 'p0') and (self.config.multi_config == False):
                """
                We must be more conservative when creating masks to compute the total flux density 
                before self-calibration. The image may contain artifacts above the default sigma 
                threshold of 6.0, and may lead to overestimation of the total flux density.
                An alternative sigma is 10. Note that the mask dilation is a very powerful approach and 
                very sensitive to the sigma threshold. A sigma of 10 results large differences in 
                relation to a sigma of 6. 
                """
                sigma = 12.0
            else:
                sigma = 12.0 #the previous value of 10 was unstable. 
                

        if (selfcal_step == 'test_image' or selfcal_step == 'p0') and (self.config.multi_config == False):
            last_level = 6
        else:
            last_level = 3

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
            # residual_to_analyse = image_list[prefix + '_residual'].replace('MFS-residual.fits',
            #                                                                'MFS-residual-pb.fits')
            """#testing; use normal residual instead of pb residual, 
                         because the latter may be ``too noisy'' at the boundaries.
            """
            residual_to_analyse = image_list[prefix + '_residual'].replace('MFS-residual.fits',
                                                                           'MFS-residual.fits')
        else:
            has_pb_image = False
            image_to_analyse = image_list[prefix]
            residual_to_analyse = image_list[prefix + '_residual']

        try:
            # mask should be created from normal -MFS-image, and not from -MFS-image-pb!!!!
            _, mask_dilated = mlibs.mask_dilation(image_list[prefix],
                                                  sigma=sigma, PLOT=False,show_figure=False)
            # print('mask sum', np.nansum(mask_dilated))
            level_stats = mlibs.level_statistics(image_to_analyse, sigma=sigma,
                                                 mask=mask_dilated)
        except:
            try:
                _, mask_dilated = mlibs.mask_dilation(
                    image_list[prefix].replace('-MFS-image', '-MFS-dirty'),
                    sigma=sigma, PLOT=False,show_figure=False)
                # print('mask sum', np.nansum(mask_dilated))
                level_stats = mlibs.level_statistics(image_to_analyse, sigma=sigma,
                                                     mask=mask_dilated)
            except:
                sigma = 8.0
                _, mask_dilated = mlibs.mask_dilation(
                    image_list[prefix].replace('-MFS-image', '-MFS-dirty'),
                    sigma=sigma, PLOT=False,show_figure=False)
                if np.nansum(mask_dilated) < 1:
                    sigma = 6.0
                    _, mask_dilated = mlibs.mask_dilation(
                        image_list[prefix].replace('-MFS-image', '-MFS-dirty'),
                        sigma=sigma, PLOT=False,show_figure=False)
                
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
                                                   last_level = last_level,
                                                   do_fit_ellipse=False,
                                                   show_figure=self.config.show_figures)[-1]

        image_statistics[prefix] = img_props

        try:
            centre = mlibs.nd.maximum_position(
                np.nan_to_num(mlibs.load_fits_data(image_to_analyse) * mask_dilated, nan=0))[::-1]
            # centre = (int(image_statistics[prefix]['y0']),int(image_statistics[prefix]['x0']))
            fig = plt.figure(figsize=(16, 8))
            ax0 = fig.add_subplot(1, 2, 2)
            ax0 = mlibs.eimshow(imagename=image_to_analyse,
                                center=centre,
                                projection='offset', 
                                rms=mlibs.mad_std(
                                    np.nan_to_num(mlibs.load_fits_data(residual_to_analyse), nan=0)),
                                # crop=True,box_size=300,
                                figsize=(8, 8), ax=ax0, fig=fig,
                                plot_colorbar=True,vmin_factor=1.0, vmax_factor=1.0,
                                plot_rms=True,
                                crop=True, box_size=int(4 * image_statistics[prefix]['C95radii']),
                                # save_name=image_list[prefix].replace('.fits', '_map'),
                                # CM='magma',
                                add_beam=True)
            ax0.set_title(f'Radio Map')
            ax1 = fig.add_subplot(1, 2, 1)
            ax1 = mlibs.eimshow(imagename=residual_to_analyse,
                                center=centre,
                                projection='offset',
                                vmin_factor=-3.0, vmax_factor=1.0,
                                add_contours=False,
                                figsize=(8, 8), ax=ax1, fig=fig,
                                crop=True, box_size=int(4 * image_statistics[prefix]['C95radii']),
                                save_name=image_to_analyse.replace('.fits', '_map'),
                                plot_title=f'Residual Map', plot_colorbar=False,
                                # CM='magma',
                                add_beam=True)
            if self.config.show_figures:
                plt.show()
                plt.close(fig)
            else:
                plt.close(fig)
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
                                                            last_level=last_level,
                                                            #    mask=None,
                                                            mask=mask_dilated,
                                                            verbose=1,
                                                            save_csv=True,
                                                            do_fit_ellipse=False,
                                                        #    show_figure=self.config.show_figures,
                                                            show_figure = False)[-1]

                flux_density, flux_density_err = img_props['total_flux_mask'], img_props[
                    'total_flux_error']

                print('Flux density = ', flux_density)
                _FLUXES.append(flux_density)
                _FLUXES_err.append(flux_density_err)
                plt.close('all')
            FLUXES = mlibs.np.asarray(_FLUXES)
            FLUXES_err = mlibs.np.asarray(_FLUXES_err)
            freqlist = mlibs.getfreqs(sub_band_images)

            if self.config.show_figures:
                _verbose = 2
            else:
                _verbose = 0
            mini, result_1, param_dict = (
                mlibs.do_fit_spec_RC_linear(freqlist,
                                            FLUXES * 1000,
                                            FLUXES_err * 1000,
                                            basename_save=image_to_analyse,
                                            title_text=r'WSClean Sub-band Images',
                                            plot_errors_shade=True, do_mcmc_fit=True,
                                            verbose=_verbose))
            if self.config.show_figures:
                plt.show()
            # plt.clf()
            plt.close()
            plt.close('all')
            
            #  Append spectral fit results into image_statistics[prefix] 
            # Prefer MCMC param_dict (asymmetric bounds) when available;
            # fall back to lmfit result_1.params (symmetric stderr).
            spidx_props = {}
            for par_name in result_1.params:
                p = result_1.params[par_name]
                spidx_props[f'spidx_{par_name}_lmfit']  = p.value
                spidx_props[f'spidx_{par_name}_stderr']  = p.stderr if p.stderr is not None else np.nan
            if param_dict is not None:
                for par_name, vals in param_dict.items():
                    spidx_props[f'spidx_{par_name}_best']  = vals.get('best',        np.nan)
                    spidx_props[f'spidx_{par_name}_lower'] = vals.get('lower',       np.nan)
                    spidx_props[f'spidx_{par_name}_upper'] = vals.get('upper',       np.nan)
                    spidx_props[f'spidx_{par_name}_lo_bound'] = vals.get('lower_bound', np.nan)
                    spidx_props[f'spidx_{par_name}_hi_bound'] = vals.get('upper_bound', np.nan)
            spidx_props['spidx_freqs']  = list(freqlist)
            spidx_props['spidx_fluxes'] = list(FLUXES * 1000)
            spidx_props['spidx_fluxes_err'] = list(FLUXES_err * 1000)
            image_statistics[prefix].update(spidx_props)
            
        except:
            print('--==>> Some error found in the sub-bands images.')
            pass
        return (image_statistics, image_list)

    # def create_mask(self, imagename, rms_mask, sigma_mask, mask_grow_iterations,
    #                 PLOT=False, use_residual=True, sigma_residual=1.5,
    #                 use_snr_grow=True, sigma_grow=5.0):
    #     """
    #     Create a WSClean-compatible binary mask FITS file for a given image.

    #     Parameters
    #     ----------
    #     imagename : str
    #         Path to the cleaned image FITS file.
    #     rms_mask : float
    #         RMS noise level used for thresholding (typically mad_std of residual).
    #     sigma_mask : float
    #         Sigma threshold for seeds (use_snr_grow=True) or flat threshold
    #         (use_snr_grow=False). Decremented by 2 each iteration if no mask found.
    #     mask_grow_iterations : int
    #         Number of binary dilation iterations (beam-sized structuring element).
    #     PLOT : bool
    #         If True, pass to mlibs for diagnostic plots.
    #     use_residual : bool
    #         If True (default), locate and use the residual image for gating.
    #         Automatically disabled if the residual file does not exist.
    #     sigma_residual : float
    #         Reject pixels where |residual| > sigma_residual * rms_res (default 3.0).
    #         Raise to 4.0-5.0 for extended sources or early self-cal steps.
    #     use_snr_grow : bool
    #         If True (default), use mask_dilation_snr (seed-and-grow). This
    #         eliminates isolated artefact blobs by requiring spatial connectivity
    #         to high-confidence seeds. Falls back to mask_dilation if residual
    #         is unavailable.
    #     sigma_grow : float
    #         Lower grow threshold for mask_dilation_snr (default 5.0). Pixels
    #         above this level that are connected to seeds and pass the residual
    #         gate are included in the mask.
    #     """
    #     # -- Locate residual image (WSClean MFS naming, then plain naming) ----
    #     residual_image = imagename.replace('-MFS-image.fits', '-MFS-residual.fits')
    #     if not os.path.exists(residual_image):
    #         residual_image = imagename.replace('-image.fits', '-residual.fits')
    #     if not os.path.exists(residual_image):
    #         use_residual = False
    #         use_snr_grow = False
        
    #     if use_residual:
    #         print(" ++==>> Using residual image for masking emission.")
    #         print(f" ++==>> Residual image: {os.path.basename(residual_image)}")
    #     if use_snr_grow:
    #         print(" ++==>> Using snr_grow based masking method.")

    #     valid_sigma_mask = sigma_mask
    #     while True:
    #         if use_snr_grow and use_residual:
    #             mask_valid = mlibs.mask_dilation_snr(
    #                 imagename,
    #                 residual=residual_image,
    #                 rms=rms_mask,
    #                 sigma=valid_sigma_mask * 0.5,
    #                 sigma_seed=valid_sigma_mask,
    #                 sigma_residual=sigma_residual,
    #                 iterations=mask_grow_iterations,
    #                 PLOT=PLOT,
    #                 verbose = 2 #debug testing
    #             )[1]
    #         else:
    #             mask_valid = mlibs.mask_dilation(
    #                 imagename,
    #                 PLOT=PLOT,
    #                 rms=rms_mask,
    #                 dilation_size=None,
    #                 sigma=valid_sigma_mask,
    #                 iterations=mask_grow_iterations,
    #                 residual=residual_image if use_residual else None,
    #                 sigma_residual=sigma_residual,
    #                 min_size_beams=0.0,
    #             )[1]
    #         if mask_valid.sum() > 0:
    #             break
    #         print(' ++>> No mask found with sigma_mask:', valid_sigma_mask)
    #         print(' ++>> Reducing sigma_mask by 2 until valid mask is found...')
    #         valid_sigma_mask = valid_sigma_mask - 2.0
    #         if valid_sigma_mask <= 6:
    #             print("Reached minimum sigma threshold without finding a valid mask.")
    #             break

    #     mask_wsclean = mask_valid * 1.0
    #     mask_name = imagename.replace('.fits', '') + '_mask.fits'
    #     mlibs.pf.writeto(mask_name, mask_wsclean, overwrite=True)
    #     return mask_name
    
    # def create_mask(self, imagename, rms_mask, sigma_mask, mask_grow_iterations,
    #                     PLOT=False, use_residual=True, sigma_residual=3.0):

    #         # Locate the residual image (WSClean naming convention)
    #         residual_image = imagename.replace('-MFS-image.fits', '-MFS-residual.fits')
    #         # print("     >>>>>>>>>> USING NEW TESTING CREATE_MASK")
    #         print(f"     >>>>>>>>>> residual_image {residual_image}")
    #         if not os.path.exists(residual_image):
    #             residual_image = imagename.replace('-image.fits', '-residual.fits')
    #         if not os.path.exists(residual_image):
    #             use_residual = False
    #         # if use_residual:
    #             # print("     >>>>>>>>>> ALL OKAY")
    #         valid_sigma_mask = sigma_mask
    #         while True:
    #             mask_valid = mlibs.mask_dilation(
    #                 imagename,
    #                 PLOT=PLOT,
    #                 rms=rms_mask,
    #                 dilation_size=None,
    #                 sigma=valid_sigma_mask,
    #                 iterations=mask_grow_iterations,
    #                 residual=residual_image if use_residual else None,
    #                 sigma_residual=sigma_residual,
    #                 min_size_beams=0.0,
    #             )[1]
    #             if mask_valid.sum() > 0:
    #                 break
    #             print(' ++>> No mask found with sigma_mask:', valid_sigma_mask)
    #             print(' ++>> Reducing sigma_mask by 2 until valid mask is found...')
    #             valid_sigma_mask = valid_sigma_mask - 2.0
    #             if valid_sigma_mask <= 6:
    #                 print("Reached minimum sigma threshold without finding a valid mask.")
    #                 break

    #         # Final mask at sigma - 1 to slightly relax the threshold after
    #         # finding the valid level - gives a marginally larger mask.
    #         mask_valid = mlibs.mask_dilation(
    #             imagename,
    #             PLOT=PLOT,
    #             rms=rms_mask,
    #             dilation_size=None,
    #             sigma=valid_sigma_mask - 1,
    #             iterations=mask_grow_iterations,
    #             residual=residual_image if use_residual else None,
    #             sigma_residual=sigma_residual,
    #             min_size_beams=0.0,
    #         )[1]

    #         mask = mask_valid
    #         mask_wsclean = mask * 1.0
    #         mask_name = imagename.replace('.fits', '') + '_mask.fits'
    #         mlibs.pf.writeto(mask_name, mask_wsclean, overwrite=True)
    #         return mask_name

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

    # def plot_visibilities(self, g_vis, name, with_DATA=True, with_MODEL=False,
    #                       with_CORRECTED=False, with_RESIDUAL=False):

    #     if with_DATA == True:
    #         plotfile = os.path.dirname(g_vis) + '/selfcal/plots/' + name + '_uvwave_amp_data.jpg'
    #         if not os.path.isfile(plotfile):
    #             plotms(vis=g_vis, xaxis='UVwave', yaxis='amp', avgantenna=False, avgscan=False,
    #                    # antenna=self.config.antennas, spw=self.config.spws,
    #                    coloraxis='baseline',
    #                    ydatacolumn='data', avgchannel='9999', avgtime='9999',
    #                    correlation='RR', plotrange=[0, 0, 0, 0],
    #                    width=2000, height=800, showgui=False, overwrite=True, dpi=1200,
    #                    highres=True,
    #                    customsymbol=True, symbolsize=4, symbolshape='diamond',
    #                    # plotrange=[-1,-1,-1,0.3],
    #                    plotfile=plotfile)
    #         else:
    #             pass

    #         plotfile = os.path.dirname(g_vis) + '/selfcal/plots/' + name + '_freq_amp_data.jpg'
    #         if not os.path.isfile(plotfile):
    #             plotms(vis=g_vis, xaxis='freq', yaxis='amp', avgantenna=False, avgscan=False,
    #                    # antenna=self.config.antennas, spw=self.config.spws,
    #                    coloraxis='baseline',
    #                    ydatacolumn='data',
    #                    #    avgchannel='8',
    #                    avgtime='9999',
    #                    correlation='RR', plotrange=[0, 0, 0, 0],
    #                    width=2000, height=800, showgui=False, overwrite=True, dpi=1200,
    #                    highres=True,
    #                    customsymbol=True, symbolsize=4, symbolshape='diamond',
    #                    plotfile=plotfile)
    #         else:
    #             pass

    #     if with_CORRECTED == True:
    #         plotfile = os.path.dirname(
    #             g_vis) + '/selfcal/plots/' + name + '_uvwave_amp_corrected_div_model.jpg'
    #         if not os.path.isfile(plotfile):
    #             plotms(vis=g_vis, xaxis='UVwave', yaxis='amp',
    #                    # antenna=self.config.antennas, spw=self.config.spws,
    #                    coloraxis='baseline', avgantenna=False, avgscan=False,
    #                    ydatacolumn='corrected/model', avgchannel='9999', avgtime='9999',
    #                    correlation='RR',
    #                    width=2000, height=800, showgui=False, overwrite=True, dpi=1200,
    #                    highres=True,
    #                    customsymbol=True, symbolsize=4, symbolshape='diamond',
    #                    plotrange=[0, 0, 0, 10],
    #                    plotfile=plotfile)
    #         else:
    #             pass

    #         plotfile = os.path.dirname(
    #             g_vis) + '/selfcal/plots/' + name + '_uvwave_amp_corrected.jpg'
    #         if not os.path.isfile(plotfile):
    #             plotms(vis=g_vis, xaxis='UVwave', yaxis='amp', avgantenna=False, avgscan=False,
    #                    # antenna=self.config.antennas, spw=self.config.spws,
    #                    coloraxis='baseline',
    #                    # plotrange=[-1,-1,0,0.3],
    #                    ydatacolumn='corrected', avgchannel='9999', avgtime='9999',
    #                    correlation='RR', plotrange=[0, 0, 0, 0],
    #                    width=2000, height=800, showgui=False, overwrite=True, dpi=1200,
    #                    highres=True,
    #                    customsymbol=True, symbolsize=4, symbolshape='diamond',
    #                    plotfile=plotfile)
    #         else:
    #             pass

    #         plotfile = os.path.dirname(g_vis) + '/selfcal/plots/' + name + '_freq_amp_corrected.jpg'
    #         if not os.path.isfile(plotfile):
    #             plotms(vis=g_vis, xaxis='freq', yaxis='amp', avgantenna=False, avgscan=False,
    #                    # antenna=self.config.antennas, spw=self.config.spws,
    #                    coloraxis='baseline',
    #                    ydatacolumn='corrected',
    #                    #    avgchannel='8',
    #                    avgtime='9999',
    #                    correlation='RR', plotrange=[0, 0, 0, 0],
    #                    width=2000, height=800, showgui=False, overwrite=True, dpi=1200,
    #                    highres=True,
    #                    customsymbol=True, symbolsize=4, symbolshape='diamond',
    #                    plotfile=plotfile)
    #         else:
    #             pass

    #     if with_MODEL == True:
    #         plotfile = os.path.dirname(g_vis) + '/selfcal/plots/' + name + '_uvwave_amp_model.jpg'
    #         if not os.path.isfile(plotfile):
    #             plotms(vis=g_vis, xaxis='UVwave', yaxis='amp', avgantenna=False, avgscan=False,
    #                    # antenna=self.config.antennas, spw=self.config.spws,
    #                    coloraxis='baseline',
    #                    ydatacolumn='model', avgchannel='9999', avgtime='9999',
    #                    correlation='RR', plotrange=[0, 0, 0, 0],
    #                    width=2000, height=800, showgui=False, overwrite=True, dpi=1200,
    #                    highres=True,
    #                    customsymbol=True, symbolsize=4, symbolshape='diamond',
    #                    plotfile=plotfile)
    #         else:
    #             pass

    #         plotfile = os.path.dirname(g_vis) + '/selfcal/plots/' + name + '_freq_amp_model.jpg'
    #         if not os.path.isfile(plotfile):
    #             plotms(vis=g_vis, xaxis='freq', yaxis='amp', avgantenna=False, avgscan=False,
    #                    # antenna=self.config.antennas, spw=self.config.spws,
    #                    coloraxis='baseline',
    #                    ydatacolumn='model',
    #                    #    avgchannel='',
    #                    avgtime='9999',
    #                    correlation='RR', plotrange=[0, 0, 0, 0],
    #                    width=2000, height=800, showgui=False, overwrite=True, dpi=1200,
    #                    highres=True,
    #                    customsymbol=True, symbolsize=4, symbolshape='diamond',
    #                    plotfile=plotfile)
    #         else:
    #             pass

    #     if with_RESIDUAL == True:
    #         plotfile = os.path.dirname(
    #             g_vis) + '/selfcal/plots/' + name + '_uvwave_amp_corrected-model.jpg'
    #         if not os.path.isfile(plotfile):
    #             plotms(vis=g_vis, xaxis='UVwave', yaxis='amp',
    #                    # antenna=self.config.antennas, spw=self.config.spws,
    #                    coloraxis='baseline', avgantenna=False, avgscan=False,
    #                    ydatacolumn='corrected-model', avgchannel='9999', avgtime='9999',
    #                    correlation='RR', plotrange=[0, 0, 0, 0],
    #                    width=2000, height=800, showgui=False, overwrite=True, dpi=1200,
    #                    highres=True,
    #                    customsymbol=True, symbolsize=4, symbolshape='diamond',
    #                    plotfile=plotfile)
    #         else:
    #             pass

    #     pass

    # def plot_uvwave(self, g_vis, name):

    #     plotfile = os.path.dirname(g_vis) + '/selfcal/plots/' + name + '_uvwave.jpg'
    #     if not os.path.isfile(plotfile):
    #         plotms(vis=g_vis, xaxis='uwave', yaxis='vwave', avgantenna=False,
    #                # antenna=self.config.antennas, spw=self.config.spws,
    #                coloraxis='spw',
    #                ydatacolumn='data', avgchannel='16', avgtime='12',
    #                correlation='RR,LL',
    #                width=800, height=800, showgui=False, overwrite=True, dpi=1200, highres=True,
    #                # customsymbol=True,symbolsize=1,symbolshape='diamond',
    #                plotfile=plotfile)
    #     plotfile = os.path.dirname(g_vis) + '/selfcal/plots/' + name + '_uv.jpg'
    #     if not os.path.isfile(plotfile):
    #         plotms(vis=g_vis, xaxis='u', yaxis='v', avgantenna=False,
    #                # antenna=self.config.antennas, spw=self.config.spws,
    #                coloraxis='spw',
    #                ydatacolumn='data', avgchannel='9999', avgtime='2',
    #                correlation='RR,LL',
    #                width=800, height=800, showgui=False, overwrite=True, dpi=1200, highres=True,
    #                # customsymbol=True,symbolsize=1,symbolshape='diamond',
    #                plotfile=plotfile)
    #     pass

    # def plot_visibilities(self, g_vis, name, with_DATA=True, with_MODEL=False,
    #                       with_CORRECTED=False, with_RESIDUAL=False):
    #     """
    #     Temp function.
    #     """
    #     pass

    # def plot_uvwave(self, g_vis, name):
    #     """
    #     Temp function.
    #     """
    #     pass


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

    # def calibration_table_plot(self, table, stage='selfcal',
    #                            table_type='gain_phase', kind='',
    #                            xaxis='time', yaxis='phase',
    #                            fields=[''], showgui=True):
    #     if not os.path.exists(os.path.dirname(table) + '/plots/' + stage):
    #         os.makedirs(os.path.dirname(table) + '/plots/' + stage)

    #     if yaxis == 'phase':
    #         plotrange = [-1, -1, -180, 180]
    #     else:
    #         plotrange = [-1, -1, -1, -1]

    #     if fields == '':

    #         plotms(vis=table, xaxis=xaxis, yaxis=yaxis, field='',
    #                gridcols=1, gridrows=1, coloraxis='spw', antenna='', plotrange=plotrange,
    #                width=1000, height=400, dpi=600, overwrite=True, showgui=showgui,
    #                # correlation='LL,RR',
    #                plotfile=os.path.dirname(
    #                    table) + '/plots/' + stage + '/' + table_type + '_' + xaxis + '_' + yaxis + '_field_' + str(
    #                    'all') + '.jpg')

    #         plotms(vis=table, xaxis=xaxis, yaxis=yaxis, field='', avgbaseline=True,
    #                gridcols=1, gridrows=1, coloraxis='spw', antenna='', plotrange=plotrange,
    #                width=1000, height=400, dpi=600, overwrite=True, showgui=showgui,
    #                # correlation='LL,RR',
    #                plotfile=os.path.dirname(
    #                    table) + '/plots/' + stage + '/' + table_type + '_' + xaxis + '_' +
    #                         yaxis + '_avgbaseline_field_' + str(
    #                    'all') + '.jpg')

    #     else:

    #         for FIELD in fields:
    #             plotms(vis=table, xaxis=xaxis, yaxis=yaxis, field=FIELD,
    #                    # gridcols=4,gridrows=4,coloraxis='spw',antenna='',iteraxis='antenna',
    #                    # width=2048,height=1280,dpi=256,overwrite=True,showgui=False,
    #                    gridcols=1, gridrows=1, coloraxis='spw', antenna='',
    #                    plotrange=plotrange,
    #                    width=1000, height=400, dpi=600, overwrite=True, showgui=showgui,
    #                    # correlation='LL,RR',
    #                    plotfile=os.path.dirname(
    #                        table) + '/plots/' + stage + '/' + table_type + '_' + xaxis + '_' + yaxis + '_field_' + str(
    #                        FIELD) + '.jpg')

    #     pass

    def calibration_table_plot(self, table, stage='selfcal',
                               table_type='gain_phase', kind='',
                               xaxis='time', yaxis='phase',
                               fields=[''], showgui=True):
        """
        Plot a gain calibration table.

        Attempts to use the pure-Python gain_plots module first (no GUI,
        produces _all.png and _avgbaseline.png for each coloraxis variant).
        Falls back to plotms if gain_plots is unavailable or raises an error.

        Parameters
        ----------
        table      : str   path to the CASA .tb caltable
        stage      : str   subdirectory under plots/ (default 'selfcal')
        table_type : str   used as filename prefix for the output plots
        yaxis      : str   'phase' | 'amp' | 'amplitude'
        showgui    : bool  only relevant for the plotms fallback
        """
        out_dir = os.path.join(os.path.dirname(table), 'plots', stage)
        os.makedirs(out_dir, exist_ok=True)

        # -- Try pure-Python gain_plots first ----------------------------------
        # Import failure is reported once but does NOT trigger the fallback -
        # only a runtime failure during the actual plotting does.
        try:
            import gain_plots
        except ImportError:
            print('  [calibration_table_plot] WARNING: gain_plots module not '
                  'found. Add gain_plots.py to the same directory as ph4ser.py.')
            gain_plots = None

        if gain_plots is not None:
            try:
                gain_plots.plot_caltable_all_axes(
                    table_path       = table,
                    yaxis            = yaxis,
                    gap_pad_fraction = 0.1,
                    savefig          = True,
                    out_dir          = out_dir,
                    prefix           = table_type,
                )
                return   # success - skip plotms fallback
            except Exception as e:
                print(f'  [calibration_table_plot] gain_plots plotting failed '
                      f'({type(e).__name__}: {e}).')
                print('  [calibration_table_plot] Falling back to plotms.')

        # -- plotms fallback --------------------------------------------------─
        if yaxis == 'phase':
            plotrange = [-1, -1, -180, 180]
            _yaxis_pm = 'phase'
        else:
            plotrange  = [-1, -1, -1, -1]
            _yaxis_pm  = 'amp'   # plotms uses 'amp', not 'amplitude'

        _base = os.path.join(out_dir, f'{table_type}_{xaxis}_{yaxis}')

        if fields == '' or fields == ['']:
            plotms(vis=table, xaxis=xaxis, yaxis=_yaxis_pm, field='',
                   gridcols=1, gridrows=1, coloraxis='spw', antenna='',
                   plotrange=plotrange,
                   width=1000, height=400, dpi=600,
                   overwrite=True, showgui=showgui,
                   plotfile=f'{_base}_field_all.jpg')

            plotms(vis=table, xaxis=xaxis, yaxis=_yaxis_pm, field='',
                   avgbaseline=True,
                   gridcols=1, gridrows=1, coloraxis='spw', antenna='',
                   plotrange=plotrange,
                   width=1000, height=400, dpi=600,
                   overwrite=True, showgui=showgui,
                   plotfile=f'{_base}_avgbaseline_field_all.jpg')
        else:
            for FIELD in fields:
                plotms(vis=table, xaxis=xaxis, yaxis=_yaxis_pm, field=FIELD,
                       gridcols=1, gridrows=1, coloraxis='spw', antenna='',
                       plotrange=plotrange,
                       width=1000, height=400, dpi=600,
                       overwrite=True, showgui=showgui,
                       plotfile=f'{_base}_field_{FIELD}.jpg')



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

        if return_solution_stats is True:
            SNRs, percentiles_SNRs = make_plot_check(cut_off=cut_off,return_solution_stats=return_solution_stats)
        else:
            make_plot_check(cut_off=cut_off)
        compare_phase_variation()
        if calmode == 'ap':
            compare_amp_variation()

        if return_solution_stats is True:
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
                    with_multiscale=False,
                    scales="None", maxmscales='6',
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

        imaging_script_path = os.path.dirname(os.path.abspath(__file__))
        # print(' >> Imaging script path:', imaging_script_path)

        # os.system("export OPENBLAS_NUM_THREADS=1 && python imaging_with_wsclean.py --f " +
        # os.system("export OPENBLAS_NUM_THREADS=1 && python " + imaging_script_path + "/imaging_with_wsclean34.py --f " +
        # os.system("export OPENBLAS_NUM_THREADS=1 && python " + imaging_script_path + "/imaging_with_wsclean36.py --f " +
        os.system("export OPENBLAS_NUM_THREADS=1 && python " + imaging_script_path + "/imaging_with_wsclean.py --f " +
                  g_name + " --sx "
                  + str(imsize) + " --sy " + str(imsizey) + " --niter "
                  + str(niter) + " --data " + datacolumn + " --cellsize " + cell
                  + ' --nsigma_automask ' + nsigma_automask + ' --mask ' + str(mask)
                  + ' --nsigma_autothreshold ' + nsigma_autothreshold
                  # +' --opt_args '+ opt_args
                  + ' --quiet ' + str(quiet) + ' --with_multiscale ' + str(with_multiscale)
                  + ' --scales ' + scales + ' --maxmscales ' + str(maxmscales)
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
                    with_multiscale=False, scales="'0,5,20,40'", maxmscales='6',
                    uvtaper=[], PLOT=False, with_DATA=True, with_CORRECTED=True, with_MODEL=True):

        g_vis = g_name + '.ms'
        if imsizey is None:
            imsizey = imsize
        if base_name is None:
            base_name = str(n_interaction) + '_update_model_image_'
        else:
            base_name = base_name

        imaging_script_path = os.path.dirname(os.path.abspath(__file__))
        print(' >> Imaging script path:', imaging_script_path)

        # os.system("export OPENBLAS_NUM_THREADS=1 && python imaging_with_wsclean.py --f " +
        # os.system("export OPENBLAS_NUM_THREADS=1 && python " + imaging_script_path + "/imaging_with_wsclean34.py --f " +
        # os.system("export OPENBLAS_NUM_THREADS=1 && python " + imaging_script_path + "/imaging_with_wsclean36.py --f " +
        os.system("export OPENBLAS_NUM_THREADS=1 && python " + imaging_script_path + "/imaging_with_wsclean.py --f " +
                  g_name + " --sx "
                  + str(imsize) + " --sy " + str(imsizey) + " --niter "
                  + str(niter) + " --data " + datacolumn + " --cellsize " + cell
                  + ' --nsigma_automask ' + nsigma_automask + ' --mask ' + str(mask)
                  + ' --nsigma_autothreshold ' + nsigma_autothreshold
                  + ' --nc ' + str(nc) + ' --negative_arg ' + negative_arg
                  # +' --opt_args '+ opt_args
                  + ' --quiet ' + str(quiet) + ' --with_multiscale ' + str(with_multiscale)
                  + ' --scales ' + scales + ' --maxmscales ' + str(maxmscales)
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
                      action='apply', append_cal_table=False,
                      flagbackup=True, calwt=False,
                      PLOT=False, with_CORRECTED=True, with_MODEL=True, with_DATA=True,
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
        else:
            if append_cal_table:
                gain_tables.append(caltable)
            else:
                pass
    
        return (gain_tables, spwmap)
    
    @staticmethod
    def _pair_spwmaps(tables, spwmaps):
        """
        Return exactly one spwmap entry per gain table.

        Missing entries are padded with [] (i.e. no SPW mapping), so that the
        table <-> spwmap pairing stays consistent automatically, no matter how
        many tables are being carried over from a previous step.
        """
        paired = [list(sm) for sm in spwmaps[:len(tables)]]
        paired += [[]] * (len(tables) - len(paired))
        return paired

    def _get_initial_tables_and_spwmap(self, iteration, step_spwmap,
                                       keep_tables_from=None):
        """
        Return (gain_tables_start, spwmap_start) for a selfcal solve step.

        If 'delay_K' is in config.steps, computes a fresh delay calibration
        (gaintype='K') for this iteration (applied calonly as a diagnostic)
        and prepends it as the starting tables, with an explicit empty spwmap
        entry so that spwmaps_applied always carries exactly one entry per
        table in gain_tables_applied.

        Delays are re-computed independently for each iteration so that
        improved model/data quality after each phase step is exploited.

        If `keep_tables_from` names a previous step (e.g. 'p0'), that step's
        gain tables are also prepended, together with their spwmaps, so that
        the current step solves incrementally on top of them. This is what
        general_settings['keep_p0'] uses to turn the default chain
        `p0 | p1 > p2 > ap1` into `p0 > p1 > p2 > ap1`; p2/ap1 then inherit
        p0 automatically, since they copy the full cumulative lists of the
        previous step.

        iteration        : the current step's iteration string ('0', '1', ...)
        step_spwmap      : the spwmap for the current step only
                           ([] when no SPW combining, [[0,0,...]] when combining).
        keep_tables_from : step name whose tables/spwmaps to carry over, or None.

        The returned spwmap always holds one entry per returned gain table plus
        the current step's own entry last, as gaincal/applycal expect.
        """
        step_entry = step_spwmap[0] if step_spwmap else []
        tables_start, spwmaps_start = [], []

        if 'delay_K' in self.config.steps:
            delay_tables, _ = self.self_gain_cal(
                self.g_name,
                n_interaction=iteration,
                minsnr=2.0,
                solint='inf',
                flagbackup=False,
                gaintype='K',
                combine='scan',
                refant=self.refant,
                refantmode=self.config.refantmode,
                minblperant=self.config.minblperant,
                calmode='p',
                spwmap=[],
                solnorm=False,
                calwt=False,
                applymode='calflag',
                action='',
                append_cal_table=True,
                PLOT=False,
                gain_tables=[]
            )
            tables_start += delay_tables
            spwmaps_start += self._pair_spwmaps(delay_tables, [])

        if (keep_tables_from is not None
                and keep_tables_from in self.gain_tables_applied):
            # .copy() is required: self_gain_cal appends to this list in place.
            kept_tables = self.gain_tables_applied[keep_tables_from].copy()
            print(f' ++==> Keeping {keep_tables_from} gain table(s) as the root '
                  f'of the chain: {kept_tables}')
            tables_start += kept_tables
            spwmaps_start += self._pair_spwmaps(
                kept_tables, self.spwmaps_applied.get(keep_tables_from, []))

        if not tables_start:
            return [], step_spwmap

        return tables_start, spwmaps_start + [step_entry]

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
                     clipminmax=[0, 50.0],
                     # channelavg=True, chanbin=1, timeavg=True, timebin='24s',
                     # timedevscale=timedevscale, freqdevscale=freqdevscale,
                     action=action, flagbackup=False, savepars=False)
            flagdata(vis=g_vis, mode='extend', field='', spw='',
                    action='apply', datacolumn=datacolumn,
                    combinescans=False, flagbackup=False,
                    growtime=75.0, growfreq=75.0, extendpols=True)
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
        This function comes from the e-MERLIN CASA Pipeline (Javier Mold\'on).
        https://github.com/e-merlin/eMERLIN_CASA_pipeline/blob/master/functions/eMCP_functions.py#L1501
        """
        # Find phase solutions per scan:
        # tablename = calib_dir +
        # if not os.path.exists(tablename):
        if os.path.exists(tablename):
            os.system(f'rm -rf {tablename}') #must run on every find_refant call. 
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
            print('You may want to use refantmode = "flex".')
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
        """Initialize storage dictionaries for pipeline data"""
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
        """Set base name attributes."""
        self.g_name_ = f"{self.config.path}{self.config.vis_name}"
        self.g_vis_  = f"{self.g_name_}.ms"
        # g_name / g_vis start as the base; may be updated by _prepare_visibility or _phase_shift
        self.g_name  = self.g_name_
        self.g_vis   = self.g_vis_

    def _reconcile_chanbin(self, vis, channel_width, spw_selection):
        """
        Make a `channel_width` map consistent with the SPW selection that is
        about to be passed to mstransform.

        mstransform requires `chanbin` to be either a single value (applied to
        every selected SPW) or a list with exactly one entry per *selected*
        SPW. A map built over all rows of the SPECTRAL_WINDOW table -- which is
        what happens when the MS still describes leftover SPWs that carry no
        data -- would otherwise be silently misaligned or rejected. In that
        case the map is sliced down to the selected SPW ids here.

        Args:
            vis: Input visibility file (.ms).
            channel_width: Scalar or sequence of per-SPW channel bin widths.
            spw_selection: Comma-separated SPW selection string handed to
                mstransform (as returned by get_spwids(return_string=True)).

        Returns:
            list: Channel bin widths (>= 1) aligned with `spw_selection`.
        """
        selected = [int(spwid) for spwid in str(spw_selection).split(',')
                    if spwid.strip() != '']

        if np.isscalar(channel_width):
            chanbin = [max(1, int(channel_width))]
        else:
            chanbin = [max(1, int(width)) for width in np.atleast_1d(channel_width)]

        if len(chanbin) <= 1 or not selected or len(chanbin) == len(selected):
            return chanbin

        nspw_total = get_nspw(vis)
        if len(chanbin) == nspw_total and max(selected) < nspw_total:
            # map covers every described SPW, including ones with no data:
            # keep only the entries matching the selection actually used.
            aligned = [chanbin[spwid] for spwid in selected]
            print(f'     !!==> channel_width has {len(chanbin)} entries but only '
                  f'{len(selected)} SPWs carry data ({spw_selection}); the MS '
                  'metadata describes leftover SPWs.')
            print(f'     ==> Realigned channel_width to the selected SPWs: {aligned}')
            return aligned

        raise ValueError(
            f'channel_width has {len(chanbin)} entries, which matches neither the '
            f'{len(selected)} selected SPWs ({spw_selection}) nor the {nspw_total} '
            f'SPWs described by {vis}. Rebuild it with '
            'pipeline.get_chan_avg_map(vis, chan_out_avg=...), which follows the '
            'same SPW selection used here.')

    def _prepare_visibility(self):
        """
        Prepare the input MS before calibration.

        Always performs an mstransform to ensure correct SPW metadata and
        correlation selection, and to initialise/propagate WEIGHT_SPECTRUM.
        Optionally averages in time and/or frequency according to
        general_settings flags. The original MS is replaced in-place so that
        self.g_name / self.g_vis remain unchanged.

        Note: usewtspectrum=True propagates an existing WEIGHT_SPECTRUM column
        through the transform unchanged (or correctly averaged); if none exists
        it is initialised from WEIGHT.
        """
        gs      = self.config.general_settings
        vis     = self.g_vis
        vis_tmp = vis.replace('.ms', '.temp.ms')

        print('++==> Preparing visibility (SPW check / optional averaging).')

        # --- build mstransform kwargs incrementally ---
        mst_kwargs = dict(
            vis            = vis,
            outputvis      = vis_tmp,
            spw            = self.get_spwids(vis, return_string=True),
            datacolumn     = 'data',
            correlation    = gs['correlations'],
            usewtspectrum  = True,   # propagate/initialise WEIGHT_SPECTRUM
        )

        # optional time averaging
        if gs['do_average_time'] and gs['timebin'] is not None:
            print(f"     ++==> Time averaging requested: timebin = {gs['timebin']}")
            mst_kwargs['timeaverage'] = True
            mst_kwargs['timebin']     = gs['timebin']
        else:
            print('     --==> No time averaging requested (do_average_time=False or timebin=None).')

        # optional frequency averaging
        # NOTE: mstransform uses chanaverage/chanbin, NOT split's 'width' parameter
        if gs['do_average_freq'] and gs['channel_width'] is not None:
            print(f"     ++==> Frequency averaging requested: channel_width = {gs['channel_width']}")
            mst_kwargs['chanaverage'] = True
            mst_kwargs['chanbin']     = self._reconcile_chanbin(vis,
                                                                gs['channel_width'],
                                                                mst_kwargs['spw'])
        else:
            print('     --==> No frequency averaging requested (do_average_freq=False or channel_width=None).')

        print(f"     ==> Running mstransform with parameters: {mst_kwargs}")
        mstransform(**mst_kwargs)

        # replace original with the prepared MS (rename temp back to original name)
        print('     ==> Replacing original MS with prepared MS.')
        os.system(f'rm -rf {vis}')
        os.system(f'mv {vis_tmp} {vis}')
        print(f'     ++==> Visibility prepared: {vis}')
        
    def _phase_shift(self):
        """
        Optionally shift the phase centre of the MS.

        If general_settings['new_phasecentre'] is set, runs phaseshift and
        updates self.g_name / self.g_vis to point to the phase-shifted MS.
        """
        new_phasecentre = self.config.general_settings['new_phasecentre']

        if new_phasecentre is None:
            print('     --==> No phase shift requested.')
            return

        print(' ++==>> A shift for the phase centre was provided.')
        phs_name = self.g_name_ + '_phs'
        phs_vis  = phs_name + '.ms'

        _, self.rame, self.or_phc = self.get_phase_centre(self.g_vis_)

        if not os.path.exists(phs_vis):
            print(f' ++==>> Running phaseshift on input MS.')
            print(f'          From: {self.or_phc}')
            print(f'          To:   {new_phasecentre}')
            phaseshift(vis       = self.g_vis_,
                    outputvis = phs_vis,
                    phasecenter = f"{self.rame} {new_phasecentre}")
        else:
            print('--==>> Phase-shifted MS already exists. Skipping phaseshift.')

        # update working names to point at the phase-shifted MS
        self.g_name = phs_name
        self.g_vis  = phs_vis

    def _create_directories(self):
        # print(f"++==> Preparing to selfcalibrate {self.g_vis_}.")
        print('++==> Creating basic directory structure.')
        if not os.path.exists(self.config.path + 'selfcal/'):
            os.makedirs(self.config.path + 'selfcal/')
        if not os.path.exists(self.config.path + 'selfcal/plots'):
            os.makedirs(self.config.path + 'selfcal/plots')

    def _generate_listobs(self):
        print('++==> Generating listobs file.')
        listobs(vis=self.g_vis, listfile=self.g_name + '.listobs.txt', overwrite=True)
        if self.g_name != self.g_name_:
            listobs(vis=self.g_vis_, listfile=self.g_name_ + '.listobs.txt', overwrite=True)

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
        """Execute startup step."""
        print(f"++==> Preparing to selfcalibrate {self.config.path}{self.config.vis_name}")
        self._init_names()
        self._create_directories()
        self._prepare_visibility()   # fix SPW metadata + optional averaging (always runs)
        self._phase_shift()          # update phase centre (only if requested)
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
        initweights(vis=self.g_vis,dowtsp=False,wtmode='nyq')

    def _run_autoflag_init(self):
        """
        Run automatic rflag on the data before selfcalibration.
        """
        if self.config.plotting_verbosity > 1:
            self.plot_visibilities(g_vis=self.g_vis, name='before_initial_flag',
                                with_DATA=True,
                                with_MODEL=False, with_CORRECTED=False)

        # run_autoflag(g_vis, display='both', action='apply', mode = 'tfcrop',
        #              timecutoff=3.5, freqcutoff=3.5,
        #              winsize=5, datacolumn='data')
        self.run_autoflag(self.g_vis, display='report', action='apply', mode='rflag',
                          timedevscale=4.0, freqdevscale=4.0,
                          winsize=7, datacolumn='data')
        self.run_autoflag(self.g_vis, action='apply', 
                          mode = 'clip',
                          datacolumn='data')

        summary = flagdata(vis=self.g_vis, field='', mode='summary')
        self.flag_data_steps['autoflag_init'] = summary.copy()

        if self.config.plotting_verbosity > 1:
            self.plot_visibilities(g_vis=self.g_vis, name='after_initial_flag',
                                with_DATA=True,
                                with_MODEL=False, with_CORRECTED=False)

    def _run_select_refant(self):
        if self.config.refant == '':
            # if 'select_refant' in steps and 'select_refant' not in steps_performed:
            print(' ++==> Estimating order of best reference antennas...')
            tablename_refant = os.path.dirname(self.g_name) + '/selfcal/find_refant.phase'
            refant = self.find_refant(msfile=self.g_vis, field='',
                                      #   combine='spw',
                                      tablename=tablename_refant)
            print(' ++==> Preferential reference antenna order = ', refant)
            self.refant = refant
            # print(f"The reference antenna is {refant}")
        else:
            print(' ++==> Using provided reference antenna order.')
            print(f"     => {self.config.refant}")
            self.refant = self.config.refant

    def _run_fov_image(self):
        """
        Create a FOV dirty image.
        """
        # niter = 50#knowing the dirty image is enough.
        # robust = 0.5  # or 0.5 if lots of extended emission.
        self.run_wsclean(self.g_name_,
                         imsize=self.config.init_parameters['fov_image']['imsize'],
                        #  imsizey=self.config.init_parameters['fov_image']['imsizey'],
                         cell=self.config.init_parameters['fov_image']['cell'],
                         robust=self.config.init_parameters['fov_image']['robust'],
                         base_name=self.config.init_parameters['fov_image']['basename'],
                         nsigma_automask='6.0', nsigma_autothreshold='2.0',
                         n_interaction='0', savemodel=False, quiet=self.config.quiet,
                         datacolumn='DATA',
                         nc=self.config.nc,negative_arg=self.config.negative_arg,
                         with_multiscale=False,
                         shift = self.config.init_parameters['fov_image']['FIELD_SHIFT'],
                         niter=self.config.init_parameters['fov_image']['niter'],
                         PLOT=False)

        self.image_statistics, self.image_list = \
            self.compute_image_stats(path=self.config.path,
                                     image_list=self.image_list,
                                     image_statistics=self.image_statistics,
                                     prefix=self.config.init_parameters['fov_image']['basename'],
                                     sigma=10.0,
                                     selfcal_step='p0')
        plt.clf()
        plt.close('all')

        # mask_grow_iterations = global_parameters['mask_grow_iterations']

        # mask_name = create_mask(image_list['test_image_0'],
        #                         rms_mask=rms_mask,
        #                         sigma_mask=p0_params['sigma_mask'],
        #                         mask_grow_iterations=mask_grow_iterations)        
        
        file_list = glob.glob(f"{self.config.path}*MFS-image.fits")
        file_list.sort(key=os.path.getmtime, reverse=False)

        try:
            self.image_list['FOV_image'] = file_list[-1]
        except:
            self.image_list['FOV_image'] = file_list
        self.image_list['FOV_residual'] = self.image_list['FOV_image'].replace(
            'MFS-image.fits','MFS-residual.fits')
        self.image_list['FOV_model'] = self.image_list['FOV_image'].replace(
            'MFS-image.fits','MFS-model.fits')

    def _run_test_image(self):
        niter_test = self.config.init_parameters['test_image']['niter']
        robust = self.config.init_parameters['test_image']['robust']
        prefix = self.config.init_parameters['test_image']['prefix']
        if self.config.plotting_verbosity > 1:
            PLOT = True
        else:
            PLOT = False
        self.run_wsclean(self.g_name,
                         imsize=self.config.imsize,
                         imsizey=self.config.imsizey,
                         cell=self.config.cell_size,
                         robust=robust, base_name=prefix,
                         nsigma_automask='5.0', nsigma_autothreshold='2.5',
                         n_interaction='0', savemodel=False, quiet=self.config.quiet,
                         datacolumn='DATA', shift=self.config.global_parameters['FIELD_SHIFT'],
                         with_multiscale=False, scales='None', maxmscales='3',
                         nc=self.config.nc, negative_arg=self.config.negative_arg,
                         uvtaper=self.config.init_parameters['test_image']['uvtaper'],
                         niter=niter_test,
                         PLOT=PLOT, with_DATA=True, with_CORRECTED=False, with_MODEL=False)

        self.image_statistics, self.image_list = self.compute_image_stats(path=self.config.path,
                                                                          image_list=self.image_list,
                                                                          image_statistics=self.image_statistics,
                                                                          prefix=prefix,
                                                                          selfcal_step='test_image')
        plt.clf()
        plt.close('all')

        current_total_flux = abs(self.image_statistics['test_image']['total_flux_mask']) * 1000
        if current_total_flux > 10.0:
            self.image_statistics, self.image_list = self.compute_image_stats(path=self.config.path,
                                                                              image_list=self.image_list,
                                                                              image_statistics=self.image_statistics,
                                                                              prefix=prefix,
                                                                              sigma=15,
                                                                              selfcal_step='test_image')
            plt.clf()
            plt.close('all')
            
        current_total_flux = abs(self.image_statistics['test_image']['total_flux_mask']) * 1000
        if current_total_flux > 50.0:
            self.image_statistics, self.image_list = self.compute_image_stats(path=self.config.path,
                                                                              image_list=self.image_list,
                                                                              image_statistics=self.image_statistics,
                                                                              prefix=prefix,
                                                                              sigma=30,
                                                                              selfcal_step='test_image')
            plt.clf()
            plt.close('all')

        
        self.modified_robust = None
        # if self.config.multi_config == False:
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
                                with_multiscale=True, scales='None', maxmscales='3',
                                nc=self.config.nc, negative_arg=self.config.negative_arg,
                                uvtaper=self.config.init_parameters['test_image']['uvtaper'],
                                niter=self.config.init_parameters['test_image']['niter'],
                                PLOT=False)
            else:
                # self.modified_robust = robust + 0.5
                if self.config.multi_config == True:
                    self.modified_robust = robust + 0.25
                else:
                    if robust <= 0.75:
                    # if robust < 0.0: #TO-BE-REMOVED
                        # self.modified_robust = 0.75
                        self.modified_robust = robust + 0.25
                    else:
                        self.modified_robust = robust
                self.run_wsclean(self.g_name,
                                imsize=self.config.imsize,
                                imsizey=self.config.imsizey,
                                cell=self.config.cell_size,
                                robust=self.modified_robust, base_name=prefix,
                                nsigma_automask=self.config.global_parameters['nsigma_automask'],
                                nsigma_autothreshold=self.config.global_parameters['nsigma_autothreshold'],
                                n_interaction='0', savemodel=False, quiet=self.config.quiet,
                                datacolumn='DATA',
                                shift=self.config.global_parameters['FIELD_SHIFT'],
                                with_multiscale=True, scales='None', maxmscales='3',
                                nc=self.config.nc, negative_arg=self.config.negative_arg,
                                uvtaper=[''],
                                niter=self.config.init_parameters['test_image']['niter'],
                                PLOT=False)

            self.image_statistics, self.image_list = self.compute_image_stats(path=self.config.path,
                                                                            image_list=self.image_list,
                                                                            image_statistics=self.image_statistics,
                                                                            prefix=prefix,
                                                                            selfcal_step='test_image')
            plt.clf()
            plt.close('all')
            
        if self.config.plotting_verbosity >= 1:
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
                    abs(self.image_statistics['test_image']['total_flux_mask']) * 1000)
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
            if self.config.multi_config == True:
                self.selfcal_params = self.select_parameters(
                    abs(self.image_statistics['test_image']['total_flux_mask']) * 1000)
                self.parameter_selection['p0_pos'] = self.selfcal_params
                print(' ++++>> Template of Parameters to be used from now on:',
                      self.parameter_selection['p0_pos']['name'])
                if self.parameter_selection['p0_pos']['p0']['combine'] == 'spw':
                    self.parameter_selection['p0_pos']['p0']['spwmap'] = self.get_spwmap(self.g_vis)
            else:         
                try:
                    self.selfcal_params = self.select_parameters(
                        abs(self.image_statistics['selfcal_test_0']['total_flux_mask']) * 1000)
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

        # if self.config.multi_config == False and self.modified_robust is not None:
        if self.modified_robust is not None and self.p0_params['uvtaper'] == ['']:
            if self.p0_params['robust'] <= self.modified_robust:
                self.p0_params['robust'] = self.modified_robust
            else:
                pass
        
        # if self.config.params_trial_2 is None and self.modified_robust is not None:
        #     self.p0_params['robust'] = self.modified_robust

        # minblperant = 3
        # combine='spw'
        # if (self.p0_params['combine'] == 'spw' or
        #         self.p0_params['combine'] == 'scan,spw' or
        #         self.p0_params['combine'] == 'spw,scan'):
        #     self.p0_params['spwmap'] = self.get_spwmap(self.g_vis)
        if 'spw' in self.p0_params['combine']:
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
                            #  robust=0.5, #test2 Arp299 eM-C
                             base_name=prefix,
                             nsigma_automask=self.p0_params['nsigma_automask'],
                             nsigma_autothreshold=self.p0_params['nsigma_autothreshold'],
                             n_interaction=iteration, savemodel=False, quiet=self.config.quiet,
                             datacolumn='DATA',
                             shift=self.config.global_parameters['FIELD_SHIFT'],
                             uvtaper=self.p0_params['uvtaper'],
                             scales=self.p0_params['scales'],
                             with_multiscale=True, maxmscales='3',
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
            plt.clf()
            plt.close('all')

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
            
            if self.config.global_parameters['custom_mask'] is not None:
                mask_name = self.config.global_parameters['custom_mask']
                print(f" >> Using custom mask: {mask_name}")
            elif self.config.global_parameters['use_mask'] is True:
                print(f" >> Using auto-generated mask: {mask_name}")
            else:
                mask_name = None
                print(" >> No mask will be used.")

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
                             with_multiscale=self.p0_params['with_multiscale'],
                             scales=self.p0_params['scales'], 
                             maxmscales=self.p0_params['maxmscales'],
                             uvtaper=self.p0_params['uvtaper'],
                             nc=self.config.nc, negative_arg=self.config.negative_arg,
                             savemodel=True, mask=mask_name, PLOT=False,
                             robust=self.p0_params['robust'], datacolumn='DATA')

            self.image_statistics, self.image_list = self.compute_image_stats(path=self.config.path,
                                                                              image_list=self.image_list,
                                                                              image_statistics=self.image_statistics,
                                                                              prefix='start_image',
                                                                              sigma=self.p0_params[
                                                                                  'sigma_mask'],
                                                                              selfcal_step='test_image')
            plt.clf()
            plt.close('all')

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
            # print(' ++==> Estimating order of best reference antennas...')
            # self.tablename_refant = os.path.dirname(self.g_name) + '/selfcal/find_refant.phase'
            # self.refant = find_refant(msfile=self.g_vis, field='',
            #                         combine=self.p0_params['combine'],
            #                         tablename=tablename_refant)
            # print(' ++==> Preferential reference antenna order = ', refant)
            # # steps_performed.append('select_refant')
            # # print(f"The reference antenna is {refant}")
            # #     else:
            # #         refant = refant
            self._run_select_refant()

        if 'p0' not in self.steps_performed:
            # p0_params['solint'] = 'inf'
            # p0_params['combine'] = 'spw'
            if self.config.plotting_verbosity > 1:
                PLOT = True
            else:
                PLOT = False
            _gain_tables_start, _spwmap_start = self._get_initial_tables_and_spwmap(
                iteration, self.p0_params['spwmap'])
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
                                   spwmap=_spwmap_start,
                                   solnorm=self.config.solnorm,
                                   calwt=self.config.general_settings['calwt'],
                                   applymode=self.config.general_settings['applymode_p'],
                                   #  interp = 'cubicPD,'
                                   #           'cubicPD',
                                   #  interp='cubic,cubic',
                                   # interp='linearPD,'
                                   #        'linearflagrel',
                                   action='apply',
                                   PLOT=PLOT,
                                   gain_tables=_gain_tables_start
                                   )
            )

            summary = flagdata(vis=self.g_vis, field='', mode='summary')
            self.flag_data_steps['selfcal_p0'] = summary.copy()

            self.run_wsclean(self.g_name, 
                             robust=self.p0_params['robust'] if self.p0_params['robust']>0.0 else self.p0_params['robust']+0.25,
                            #  robust=0.5, #test2 Arp299 eM-C
                             imsize=self.config.imsize,
                             imsizey=self.config.imsizey,
                             cell=self.config.cell_size,
                             base_name='selfcal_test_0',
                             nsigma_automask=self.p0_params['nsigma_automask'],
                             nsigma_autothreshold=self.p0_params['nsigma_autothreshold'],
                             n_interaction='', savemodel=False, quiet=self.config.quiet,
                             with_multiscale=self.p0_params['with_multiscale'],
                             maxmscales=self.p0_params['maxmscales'],
                             datacolumn='CORRECTED_DATA',
                             uvtaper=self.p0_params['uvtaper'],
                             scales=self.p0_params['scales'],
                             nc=self.config.nc,
                             # negative_arg=negative_arg,
                            #  niter=self.config.global_parameters['niter'],
                             niter=50000,
                             shift=self.config.global_parameters['FIELD_SHIFT'],
                             PLOT=False)
            self.image_statistics, self.image_list = self.compute_image_stats(path=self.config.path,
                                                                              image_list=self.image_list,
                                                                              image_statistics=self.image_statistics,
                                                                              prefix='selfcal_test_0',
                                                                              sigma=10.0,
                                                                              selfcal_step='p0')
            plt.clf()
            plt.close('all')

            if self.config.params_trial_2 is None and self.config.multi_config == False:
                if abs(self.image_statistics['selfcal_test_0']['total_flux_mask']) * 1000 < 10.0:
                    """
                    After the first pass of selfcal [phase], the data quality may have
                    improved. If originaly, the total flux density was below 5 mJy,
                    now, after improved corrected phases, the flux may be above 5 mJy.

                    We can check if using a taper or a higher robust will increase the
                    flux density above 5 mJy. If Yes, the source will not considered as `very faint`,
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
                                     robust=0.5 if self.p0_params['uvtaper'] != [''] else 1.0, #using high robust with tapper is too aggresive
                                     base_name='selfcal_test_0',
                                     nsigma_automask=self.p0_params['nsigma_automask'],
                                     nsigma_autothreshold=self.p0_params['nsigma_autothreshold'],
                                     n_interaction='0', savemodel=False, quiet=self.config.quiet,
                                     datacolumn='CORRECTED_DATA',
                                     with_multiscale=self.p0_params['with_multiscale'],
                                     scales=self.p0_params['scales'],
                                     maxmscales=self.p0_params['maxmscales'],
                                     uvtaper=self.p0_params['uvtaper'],
                                     nc=self.config.nc, negative_arg=self.config.negative_arg,
                                    #  niter=self.config.global_parameters['niter'],
                                     niter=50000,
                                     shift=self.config.global_parameters['FIELD_SHIFT'],
                                     PLOT=False)
                    self.image_statistics, self.image_list = self.compute_image_stats(
                        path=self.config.path,
                        image_list=self.image_list,
                        image_statistics=self.image_statistics,
                        prefix='selfcal_test_0')
                    plt.clf()
                    plt.close('all')

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

            if self.config.global_parameters['custom_mask'] is not None:
                mask_name = self.config.global_parameters['custom_mask']
                print(f" >> Using custom mask: {mask_name}")
            elif self.config.global_parameters['use_mask'] is True:
                print(f" >> Using auto-generated mask: {mask_name}")
            else:
                mask_name = None
                print(" >> No mask will be used.")

            self.run_wsclean(self.g_name, robust=self.p1_params['robust'],
                             imsize=self.config.imsize,
                             imsizey=self.config.imsizey,
                             cell=self.config.cell_size,
                             nsigma_automask=self.p1_params['nsigma_automask'],
                             nsigma_autothreshold=self.p1_params['nsigma_autothreshold'],
                             n_interaction=iteration, savemodel=True, quiet=self.config.quiet,
                             with_multiscale=self.p1_params['with_multiscale'],
                             maxmscales=self.p1_params['maxmscales'],
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
            plt.clf()
            plt.close('all')

            self.steps_performed.append('update_model_1')

        self.phase_tables = []
        self.spwmaps = []

        if self.config.params_trial_2 is not None:
            # if (self.p1_params['combine'] == 'spw' or
            #         self.p1_params['combine'] == 'scan,spw' or
            #         self.p1_params['combine'] == 'spw,scan'):
            if 'spw' in self.p1_params['combine']:
                self.p1_spwmap = self.get_spwmap(self.g_vis)[0]
                self.spwmaps.append(self.p1_spwmap)
                self.p1_params['spwmap'] = self.spwmaps
            else:
                self.p1_spwmap = []
                self.spwmaps.append(self.p1_spwmap)
                self.p1_params['spwmap'] = self.spwmaps
        else:
            # if self.p1_params['combine'] == 'spw':
            if 'spw' in self.p1_params['combine']:
                self.p1_params['spwmap'] = self.get_spwmap(self.g_vis)

        if 'p1' not in self.steps_performed:
            if self.config.plotting_verbosity > 1:
                PLOT = True
            else:
                PLOT = False
            _gain_tables_start, _spwmap_start = self._get_initial_tables_and_spwmap(
                iteration, self.p1_params['spwmap'],
                keep_tables_from='p0' if self.config.general_settings.get(
                    'keep_p0', False) else None)
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
                                   spwmap=_spwmap_start,
                                   calmode=self.p1_params['calmode'],
                                   # interp='cubic,cubic',
                                   action='apply',
                                   PLOT=PLOT,
                                   gain_tables=_gain_tables_start,
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

        # if self.p2_params['combine'] == 'spw':
        if 'spw' in self.p2_params['combine']:
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
            self.run_wsclean(self.g_name, 
                            #  robust=self.p2_params['robust'],
                            robust=self.p2_params['robust'] if self.p2_params['robust']>0.0 else self.p2_params['robust']+0.25,
                            #  robust=0.75, #test2 Arp299 eM-C
                             imsize=self.config.imsize,
                             imsizey=self.config.imsizey,
                             cell=self.config.cell_size,
                             base_name='selfcal_test_1',
                             nsigma_automask=self.p2_params['nsigma_automask'],
                             nsigma_autothreshold=self.p2_params['nsigma_autothreshold'],
                             n_interaction='', savemodel=False, quiet=self.config.quiet,
                             with_multiscale=self.p2_params['with_multiscale'],
                             scales=self.p2_params['scales'],
                             maxmscales=self.p2_params['maxmscales'],
                             datacolumn='CORRECTED_DATA',
                             uvtaper=self.p2_params['uvtaper'],
                             nc=self.config.nc,
                             negative_arg=self.config.negative_arg,
                             shift=self.config.global_parameters['FIELD_SHIFT'],
                            #  niter=self.config.global_parameters['niter'],
                             niter=50000,
                             PLOT=False)

            self.image_statistics, self.image_list = self.compute_image_stats(path=self.config.path,
                                                                              image_list=self.image_list,
                                                                              image_statistics=self.image_statistics,
                                                                              sigma=self.p2_params[
                                                                                  'sigma_mask'],
                                                                              prefix='selfcal_test_1')
            plt.clf()
            plt.close('all')

            mask_name = self.create_mask(self.image_list['selfcal_test_1'],
                                         rms_mask=None,
                                         sigma_mask=self.p2_params['sigma_mask'],
                                         mask_grow_iterations=self.p2_params[
                                             'mask_grow_iterations'])
            
            if self.config.global_parameters['custom_mask'] is not None:
                mask_name = self.config.global_parameters['custom_mask']
                print(f" >> Using custom mask: {mask_name}")
            elif self.config.global_parameters['use_mask'] is True:
                print(f" >> Using auto-generated mask: {mask_name}")
            else:
                mask_name = None
                print(" >> No mask will be used.")

            if self.config.plotting_verbosity > 1:
                PLOT = True
            else:
                PLOT = False
            self.run_wsclean(self.g_name, robust=self.p2_params['robust'],
                             imsize=self.config.imsize,
                             imsizey=self.config.imsizey,
                             cell=self.config.cell_size,
                             nsigma_automask=self.p2_params['nsigma_automask'],
                             nsigma_autothreshold=self.p2_params['nsigma_autothreshold'],
                             n_interaction=iteration, savemodel=True, quiet=self.config.quiet,
                             with_multiscale=self.p2_params['with_multiscale'],
                             scales=self.p2_params['scales'],
                             maxmscales=self.p2_params['maxmscales'],
                             datacolumn='CORRECTED_DATA', mask=mask_name,
                             shift=self.config.global_parameters['FIELD_SHIFT'],
                             uvtaper=self.p2_params['uvtaper'],
                             nc=self.config.nc, negative_arg=self.config.negative_arg,
                             niter=self.config.global_parameters['niter'],
                             PLOT=PLOT, with_DATA=False, with_CORRECTED=False, with_MODEL=True)

            self.image_statistics, self.image_list = self.compute_image_stats(path=self.config.path,
                                                                              image_list=self.image_list,
                                                                              image_statistics=self.image_statistics,
                                                                              sigma=self.p2_params[
                                                                                  'sigma_mask'],
                                                                              prefix='2_update_model_image')
            plt.clf()
            plt.close('all')

            self.steps_performed.append('update_model_2')

        if 'p2' not in self.steps_performed:
            if self.config.plotting_verbosity > 1:
                PLOT = True
            else:
                PLOT = False
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
                                   PLOT=PLOT,
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
        self.vis_split_name_p_nostatwt = f"{os.path.dirname(self.g_name)}/{self.config.field}{self.config.savename}_p_nostatwt.ms"

        if not os.path.exists(self.vis_split_name_p_statwt):
            print(f' ++==> Splitting phase-only self-calibrated visibility.')
            print(f'       Filename is: {self.vis_split_name_p_statwt}')
            if self.config.general_settings['new_phasecentre'] is None:
                split(vis=self.g_name + '.ms',
                      outputvis=self.vis_split_name_p_statwt,
                      datacolumn='corrected', keepflags=True)
                # split(vis=self.g_name + '.ms',
                #       outputvis=self.vis_split_name_p_nostatwt,
                #       datacolumn='corrected', keepflags=True)
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

                # split(vis=self.g_name + '.ms',
                #       outputvis=self.vis_split_name_p_nostatwt.replace('.ms', '_temp.ms'),
                #       datacolumn='corrected', keepflags=True)
                # print(' ++==> Phase-shiftting to original phasecentre...')
                # phaseshift(vis=self.vis_split_name_p_nostatwt.replace('.ms', '_temp.ms'),
                #            outputvis=self.vis_split_name_p_nostatwt,
                #            phasecenter=f"J2000 {self.or_phc}"
                #            )
                # os.system(f"rm -r {self.vis_split_name_p_nostatwt.replace('.ms', '_temp.ms')}")



        if self.config.instrument == 'eM':
            minblperant = 3
        else:
            minblperant = self.config.minblperant

        # if self.config.params_trial_2 is not None or self.multi_config == True:
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
        # if (self.ap1_params['combine'] == 'spw' or
        #         self.ap1_params['combine'] == 'scan,spw' or
        #         self.ap1_params['combine'] == 'spw,scan'):
        if 'spw' in self.ap1_params['combine']:
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
            self.run_wsclean(self.g_name, 
                            #  robust=self.ap1_params['robust'],
                             robust=self.ap1_params['robust'] if self.ap1_params['robust']>0.0 else self.ap1_params['robust']+0.25,
                            #  robust=0.5, #test2 Arp299 eM-C
                             imsize=self.config.imsize,
                             imsizey=self.config.imsizey,
                             cell=self.config.cell_size,
                             base_name='selfcal_test_2',
                             nsigma_automask=self.ap1_params['nsigma_automask'],
                             nsigma_autothreshold=self.ap1_params['nsigma_autothreshold'],
                             n_interaction='', savemodel=False, quiet=self.config.quiet,
                             with_multiscale=self.ap1_params['with_multiscale'],
                             scales=self.ap1_params['scales'],
                             maxmscales=self.ap1_params['maxmscales'],
                             datacolumn='CORRECTED_DATA',
                             uvtaper=self.ap1_params['uvtaper'],
                             shift=self.config.global_parameters['FIELD_SHIFT'],
                             nc=self.config.nc,
                             # negative_arg=self.config.negative_arg,
                            #  niter=self.config.global_parameters['niter'],
                             niter=50000,
                             PLOT=False)

            self.image_statistics, self.image_list = self.compute_image_stats(path=self.config.path,
                                                                              image_list=self.image_list,
                                                                              image_statistics=self.image_statistics,
                                                                              sigma=self.ap1_params[
                                                                                  'sigma_mask'],
                                                                              prefix='selfcal_test_2')
            plt.clf()
            plt.close('all')

            # mask_grow_iterations = self.ap1_params['mask_grow_iterations']

            mask_name = self.create_mask(self.image_list['selfcal_test_2'],
                                         rms_mask=None,
                                         sigma_mask=self.ap1_params['sigma_mask'],
                                         mask_grow_iterations=self.ap1_params[
                                             'mask_grow_iterations'])

            if self.config.global_parameters['custom_mask'] is not None:
                mask_name = self.config.global_parameters['custom_mask']
                print(f" >> Using custom mask: {mask_name}")
            elif self.config.global_parameters['use_mask'] is True:
                print(f" >> Using auto-generated mask: {mask_name}")
            else:
                mask_name = None
                print(" >> No mask will be used.")

            self.run_wsclean(self.g_name, robust=self.ap1_params['robust'],
                             imsize=self.config.imsize,
                             imsizey=self.config.imsizey,
                             cell=self.config.cell_size,
                             nsigma_automask=self.ap1_params['nsigma_automask'],
                             nsigma_autothreshold=self.ap1_params['nsigma_autothreshold'],
                             n_interaction=iteration, savemodel=True, quiet=self.config.quiet,
                             with_multiscale=self.ap1_params['with_multiscale'],
                             scales=self.ap1_params['scales'],
                             maxmscales=self.ap1_params['maxmscales'],
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
            plt.clf()
            plt.close('all')

            self.steps_performed.append('update_model_3')

        if 'ap1' not in self.steps_performed:
            if self.config.plotting_verbosity >= 1:
                PLOT = True
            else:
                PLOT = False
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
                                   PLOT=PLOT,
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
        self.vis_split_name_1_nostatwt = f"{os.path.dirname(self.g_name)}/{self.config.field}{self.config.savename}_nostatwt.ms"

        if not os.path.exists(self.vis_split_name_1_statwt):
            print(f' ++==> Splitting final self-calibrated (p+ap) visibility.')
            print(f'       Filename is: {self.vis_split_name_1_statwt}')
            if self.config.general_settings['new_phasecentre'] is None:
                split(vis=self.g_name + '.ms',
                      outputvis=self.vis_split_name_1_statwt,
                      datacolumn='corrected', keepflags=True)
                split(vis=self.g_name + '.ms',
                      outputvis=self.vis_split_name_1_nostatwt,
                      datacolumn='corrected', keepflags=True)
                print(' ++==> Running statw on split data...')
                statwt(vis=self.vis_split_name_1_statwt,
                       statalg=self.config.general_settings['statwt_statalg'],
                       timebin=self.config.general_settings['timebin_statw'],
                       datacolumn='data')
            else:
                # Phase-shift happens AFTER all imaging in _run_phaseshift_final_ms(),
                # so that selfcal_image is consistent with all other test images
                # (same shifted phase centre) and PB correction is reliable.
                split(vis=self.g_name + '.ms',
                      outputvis=self.vis_split_name_1_statwt,
                      datacolumn='corrected', keepflags=True)
                split(vis=self.g_name + '.ms',
                      outputvis=self.vis_split_name_1_nostatwt,
                      datacolumn='corrected', keepflags=True)
                print(' ++==> Running statwt on split data...')
                statwt(vis=self.vis_split_name_1_statwt,
                       statalg=self.config.general_settings['statwt_statalg'],
                       timebin=self.config.general_settings['timebin_statw'],
                       datacolumn='data')

        # if os.path.exists(self.vis_split_name_p_statwt):
        #     print(' ++==> Running statw on split phase-only self-calibrated visibility...')
        #     statwt(vis=self.vis_split_name_p_statwt,
        #            statalg=self.config.general_settings['statwt_statalg'],
        #            timebin=self.config.general_settings['timebin_statw'],
        #            datacolumn='data')

        print(' ++==> Imaging visibilities after self-calibration...')

        niter = self.config.global_parameters['niter']

        if (self.config.params_trial_2 is not None) or (self.config.multi_config):
            ROBUSTS = [-0.25, 0.25]
        else:
            # ROBUSTS = [0.5]
            ROBUSTS = [self.ap1_params['robust'] - 0.25 if self.ap1_params['robust']>0.0 else self.ap1_params['robust']+0.25]
            # ROBUSTS = [-1.0]

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
                                     with_multiscale=self.config.global_parameters[
                                         'with_multiscale'],
                                     n_interaction='', savemodel=False, quiet=self.config.quiet,
                                     scales=self.config.global_parameters['scales'],
                                     maxmscales=self.config.global_parameters['maxmscales'],
                                     datacolumn='DATA',
                                     shift=self.config.global_parameters['FIELD_SHIFT'],
                                    #  uvtaper=self.config.global_parameters['uvtaper'],
                                     uvtaper=self.ap1_params['uvtaper'],
                                     nc=self.config.nc,
                                     negative_arg='negative',
                                     niter=niter,
                                     PLOT=False)

                    self.image_statistics, self.image_list = self.compute_image_stats(
                        path=self.config.path,
                        image_list=self.image_list,
                        image_statistics=self.image_statistics,
                        prefix='selfcal_image')
                    plt.clf()
                    plt.close('all')
                    
                    # rms_final = self.image_statistics['selfcal_image']
                    rms_final = mlibs.mad_std(mlibs.load_fits_data(self.image_list['selfcal_image_residual']))
                    mask_name = self.create_mask(self.image_list['selfcal_image'],
                                                rms_mask=rms_final,
                                                sigma_mask=8.0,
                                                mask_grow_iterations=1
                                                # mask_grow_iterations=self.config.global_parameters['mask_grow_iterations']
                                                )
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

        if self.config.plotting_verbosity > 1:
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
                # Phase-shift deferred to _run_phaseshift_final_ms() after all imaging.
                split(vis=self.g_vis,
                      outputvis=self.vis_split_name_flag,
                      datacolumn='corrected', keepflags=True)

            statwt(vis=self.vis_split_name_flag,
                   statalg=self.config.general_settings['statwt_statalg'],
                   timebin=self.config.general_settings['timebin_statw'],
                   datacolumn='data')
            
        if self.config.plotting_verbosity > 1:
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
                         nsigma_automask=self.config.global_parameters['nsigma_automask'],
                         nsigma_autothreshold=self.config.global_parameters['nsigma_autothreshold'],
                         with_multiscale=self.config.global_parameters['with_multiscale'],
                         n_interaction='', savemodel=False, quiet=True,
                        #  with_multiscale=True,
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
        plt.clf()
        plt.close('all')

        # steps_performed.append('run_autoflag_final')

    def _run_phaseshift_final_ms(self):
        """
        Phase-shift all final split MSes back to the original phase centre.
        Called automatically after all imaging steps when new_phasecentre was set,
        so that delivered MSes are at the original centre while all imaging was
        done consistently at the shifted centre.
        Each MS is shifted to a temp name, the original is deleted, and the temp
        is renamed back so the filenames remain unchanged downstream.
        """
        if self.config.general_settings['new_phasecentre'] is None:
            return

        candidates = [
            ('vis_split_name_1_statwt',  getattr(self, 'vis_split_name_1_statwt',  None)),
            ('vis_split_name_1_nostatwt', getattr(self, 'vis_split_name_1_nostatwt', None)),
            ('vis_split_name_flag',       getattr(self, 'vis_split_name_flag',       None)),
        ]

        for attr_name, ms_path in candidates:
            if ms_path is None or not os.path.exists(ms_path):
                continue
            temp_path = ms_path.replace('.ms', '_phs_temp.ms')
            print(f' ++==> Phase-shifting {os.path.basename(ms_path)} to original phase centre...')
            phaseshift(vis=ms_path,
                       outputvis=temp_path,
                       phasecenter=f"J2000 {self.or_phc}")
            os.system(f"rm -r {ms_path}")
            os.system(f"mv {temp_path} {ms_path}")
            print(f'       Done: {os.path.basename(ms_path)}')

    # def _run_organise_products(self):
    #     """
    #     Move all non-MS files produced during the self-calibration run into a
    #     subdirectory called selfcal_products/ inside the working directory.

    #     Files moved:
    #       - All FITS images, residuals, models, masks, pb-corrected cubes
    #       - All CSVs, JPGs, PDFs produced by the pipeline and mlibs
    #       - CASA log files (casa-*.log)
    #       - listobs files
    #       - Gain table plot files (selfcal/plots/ is already a subdirectory;
    #         its contents are left in place)

    #     Files NOT moved:
    #       - *.ms directories (measurement sets)
    #       - *.ms.flagversions directories
    #       - The selfcal/ subdirectory and its contents (plots, gain tables)
    #       - selfcal_products/ itself

    #     Symlinks to the key final products (selfcal_image MFS fits and its
    #     residual) are created in the working directory for quick access.
    #     """
    #     work_dir = os.path.dirname(os.path.abspath(self.g_name))
    #     products_dir = os.path.join(work_dir, 'selfcal_products')
    #     os.makedirs(products_dir, exist_ok=True)

    #     # Extensions to move
    #     move_extensions = {
    #         '.fits', '.jpg', '.jpeg', '.png', '.pdf',
    #         '.csv', '.log', '.listobs',
    #     }

    #     # Names / patterns to always leave in place
    #     skip_names = {'selfcal_products', 'selfcal'}

    #     moved, skipped, failed = [], [], []

    #     for entry in sorted(os.listdir(work_dir)):
    #         full_path = os.path.join(work_dir, entry)

    #         # Never touch directories (covers .ms, .ms.flagversions, selfcal/, etc.)
    #         if os.path.isdir(full_path):
    #             continue

    #         # Skip if already inside products_dir (shouldn't happen, but safe)
    #         if entry in skip_names:
    #             continue

    #         # Only move files with recognised extensions
    #         _, ext = os.path.splitext(entry)
    #         if ext.lower() not in move_extensions:
    #             skipped.append(entry)
    #             continue

    #         dest = os.path.join(products_dir, entry)
    #         try:
    #             shutil.move(full_path, dest)
    #             moved.append(entry)
    #         except Exception as exc:
    #             print(f'[organise_products] Could not move {entry}: {exc}')
    #             failed.append(entry)

    #     print(f'[organise_products] Moved {len(moved)} files -> selfcal_products/')
    #     if failed:
    #         print(f'[organise_products] Failed to move {len(failed)} files: '
    #               + ', '.join(failed))

    #     # -- Convenience symlinks in work_dir for the final selfcal image ------
    #     for key in ('selfcal_image', 'selfcal_image_residual'):
    #         src_fits = self.image_list.get(key)
    #         if src_fits is None:
    #             continue
    #         fname    = os.path.basename(src_fits)
    #         src_dest = os.path.join(products_dir, fname)
    #         link     = os.path.join(work_dir, fname)
    #         if os.path.exists(src_dest) and not os.path.exists(link):
    #             try:
    #                 os.symlink(src_dest, link)
    #                 print(f'[organise_products] Symlink -> {fname}')
    #             except Exception as exc:
    #                 print(f'[organise_products] Could not create symlink for '
    #                       f'{fname}: {exc}')

    def _run_organise_products(self):
        """
        Move all non-MS files into selfcal_products/, then copy a curated set
        of key files back to the working directory for quick access.

        Files moved to selfcal_products/:
          - All FITS, JPG, PNG, PDF, CSV, LOG, LISTOBS files

        Files NOT moved:
          - *.ms and *.ms.flagversions directories
          - selfcal/ subdirectory (gain tables, plots)
          - selfcal_products/ itself

        Files copied back to the working directory after the move:
          - selfcal_image MFS-image-pb.fits   (primary pb-corrected image)
          - selfcal_image MFS-residual-pb.fits (pb-corrected residual)
          - *_selfcal_results.pdf
          - *_selfcal_masks.pdf
          - *MFS-image-pb_map.jpg
          - *MFS-image-pb_Lgrow_levels.jpg
          - *MFS-image-pb__RC_alpha_fit_linear.jpg
        Each copy is guarded by os.path.exists() checks.
        """
        work_dir     = os.path.dirname(os.path.abspath(self.g_name))
        products_dir = os.path.join(work_dir, 'selfcal_products')
        os.makedirs(products_dir, exist_ok=True)

        move_extensions = {
            '.fits', '.jpg', '.jpeg', '.png', '.pdf',
            '.csv', '.log', '.listobs',
        }
        skip_names = {'selfcal_products', 'selfcal'}

        moved, failed = [], []

        for entry in sorted(os.listdir(work_dir)):
            full_path = os.path.join(work_dir, entry)
            if os.path.isdir(full_path) or entry in skip_names:
                continue
            _, ext = os.path.splitext(entry)
            if ext.lower() not in move_extensions:
                continue
            dest = os.path.join(products_dir, entry)
            try:
                shutil.move(full_path, dest)
                moved.append(entry)
            except Exception as exc:
                print(f'[organise_products] Could not move {entry}: {exc}')
                failed.append(entry)

        print(f'[organise_products] Moved {len(moved)} files -> selfcal_products/')
        if failed:
            print(f'[organise_products] Failed: ' + ', '.join(failed))

        # -- Copy key files back to the working directory ----------------------
        # All paths are derived exactly from self.image_list and self.g_name
        # so only the final selfcal_image products are copied, not intermediate
        # images from other self-cal steps.
        copy_back = []

        # 1. pb-corrected final image and residual (derived from selfcal_image)
        selfcal_fits = self.image_list.get('selfcal_image', '')
        mfs_base = selfcal_fits.replace('-MFS-image.fits', '')
        for suffix in ['-MFS-image-pb.fits', '-MFS-residual-pb.fits']:
            fname = os.path.basename(mfs_base + suffix)
            src   = os.path.join(products_dir, fname)
            if os.path.exists(src):
                copy_back.append(src)
            else:
                print(f'[organise_products] Not found, skipping: {fname}')

        # 2. Diagnostic JPGs tied to the final selfcal_image specifically
        for suffix in ['MFS-image-pb_map.jpg',
                       'MFS-image-pb_Lgrow_levels.jpg',
                       'MFS-image-pb__RC_alpha_fit_linear.jpg']:
            fname = os.path.basename(mfs_base + '-' + suffix)
            src   = os.path.join(products_dir, fname)
            if os.path.exists(src):
                copy_back.append(src)
            else:
                print(f'[organise_products] Not found, skipping: {fname}')

        # 3. Summary PDFs derived from self.g_name
        g_base = os.path.basename(self.g_name)
        for suffix in ['_selfcal_results.pdf', '_selfcal_masks.pdf']:
            fname = g_base + suffix
            src   = os.path.join(products_dir, fname)
            if os.path.exists(src):
                copy_back.append(src)
            else:
                print(f'[organise_products] Not found, skipping: {fname}')

        copied = []
        for src in copy_back:
            dest = os.path.join(work_dir, os.path.basename(src))
            if os.path.exists(dest):
                continue   # already there if run twice
            try:
                shutil.copy2(src, dest)
                copied.append(os.path.basename(src))
            except Exception as exc:
                print(f'[organise_products] Could not copy back '
                      f'{os.path.basename(src)}: {exc}')

        if copied:
            print(f'[organise_products] Copied back {len(copied)} key files '
                  'to working directory:')
            for f in copied:
                print(f'  {f}')

    def _run_report_results(self):
        """
        Save self-calibration statistics to CSV and produce a summary figure
        with (a) before/after radio maps and (b) per-step metric evolution.
        The mask progression figure is delegated to _run_report_masks().
        Flagging statistics use the existing plot_category_data machinery.
        """
        # -- 1. Save statistics tables ----------------------------------------
        self.df = pd.DataFrame.from_dict(self.image_statistics, orient='index')
        self.df.to_csv(self.g_name + '_selfcal_statistics.csv',
                       header=True, index=False)
        self.df_gt = pd.DataFrame.from_dict(self.gain_tables_applied, orient='index')
        self.df_gt.to_csv(self.g_name + '_tables_applied.csv',
                          header=True, index=False)
 
        # -- 2. Build ordered step sequence from what is actually present ----─
        # test_image_0, selfcal_test_0, selfcal_test_2, selfcal_image always
        # exist; selfcal_test_1 is only produced during the p2 stage.
        _candidates = [
            ('test_image',   'Test\nImage'),
            ('selfcal_test_0', 'After\np0'),
            ('selfcal_test_1', 'After\np1'),   # optional - only if p1/p2 ran
            ('selfcal_test_2', 'Before\nap1'),
            ('selfcal_image',  'Final\n(after ap1)'),
        ]
        step_sequence = [(k, lbl) for k, lbl in _candidates
                         if k in self.image_statistics]
        step_keys   = [k   for k, _ in step_sequence]
        step_labels = [lbl for _, lbl in step_sequence]
 
        # -- 3. Extract per-step metrics with error estimates ----------------─
        def _get(k, field, scale=1.0):
            return self.image_statistics[k].get(field, np.nan) * scale
 
        peak_vals     = [_get(k, 'peak_of_flux')            for k in step_keys]
        peak_err_vals = [_get(k, 'peak_error')              for k in step_keys]
        rms_raw_vals  = [_get(k, 'mad_std_residual')        for k in step_keys]
        rms_vals      = [v * 1e6 for v in rms_raw_vals]
 
        # SNR = peak / rms_residual;  σ_SNR ≈ SNR × (σ_peak / peak)
        snr_vals = []
        snr_err_vals = []
        for pk, pk_err, rms_r in zip(peak_vals, peak_err_vals, rms_raw_vals):
            snr = pk / rms_r if rms_r > 0 else np.nan
            snr_vals.append(snr)
            snr_err_vals.append(snr * (pk_err / pk) if (pk > 0 and np.isfinite(snr)) else np.nan)
 
        flux_vals     = [abs(_get(k, 'total_flux_mask')) * 1e3  for k in step_keys]
        flux_err_vals = [abs(_get(k, 'total_flux_error')) * 1e3 for k in step_keys]

        # Residual peak convergence: max and |min| residual per step
        res_max_vals = [_get(k, 'max_residual') * 1e6 for k in step_keys]
        res_min_vals = [_get(k, 'min_residual') * 1e6 for k in step_keys]
        
        # Spectral index - prefer MCMC best+asymmetric bounds, fall back to lmfit symmetric stderr.
        # Steps without sub-band fit carry np.nan and are shown as gaps in the plot.
        spidx_vals    = [_get(k, 'spidx_alpha_best')  for k in step_keys]
        spidx_lo_vals = [_get(k, 'spidx_alpha_lower') for k in step_keys]
        spidx_hi_vals = [_get(k, 'spidx_alpha_upper') for k in step_keys]
        # for i, k in enumerate(step_keys):
        #     if np.isnan(spidx_vals[i]):
        #         spidx_vals[i]    = _get(k, 'spidx_alpha_lmfit')
        #         spidx_lo_vals[i] = _get(k, 'spidx_alpha_stderr')
        #         spidx_hi_vals[i] = _get(k, 'spidx_alpha_stderr')
        for i, k in enumerate(step_keys):
            spidx_vals[i]    = _get(k, 'spidx_alpha_lmfit')
            spidx_lo_vals[i] = _get(k, 'spidx_alpha_stderr')
            spidx_hi_vals[i] = _get(k, 'spidx_alpha_stderr')
        has_spidx = any(np.isfinite(v) for v in spidx_vals)
 
        # -- 4. Shared spatial parameters for the image panels ----------------
        self.rms = mlibs.mad_std(mlibs.load_fits_data(self.image_list['selfcal_image_residual']))
        self.centre = mlibs.nd.maximum_position(
            np.nan_to_num(mlibs.load_fits_data(self.image_list['selfcal_image']), nan=0))[::-1]
        c95 = self.image_statistics['selfcal_image'].get('C95radii', 200)
        box_size = max(int(3.5 * c95), 200)
 
        # -- 5. Build figure --------------------------------------------------
        #  Layout - gridspec(4 rows × 8 cols, height_ratios 3:3:2:2):
        #    rows 0-1 : test_image_0 (cols 0-3)  |  selfcal_image (cols 4-7)
        #    rows 2-3 : SNR (0-1) | RMS (2-3) | Flux (4-5) | reserved \alpha (6-7)
        fig = plt.figure(figsize=(20, 11))
        gs  = gridspec.GridSpec(4, 10, figure=fig,
                                height_ratios=[3, 3, 2, 2],
                                hspace=0.50, wspace=0.65)
 
        # -- Radio maps ------------------------------------------------------─
        rms_test = mlibs.mad_std(mlibs.load_fits_data(self.image_list['test_image_residual']))
 
        ax_before = fig.add_subplot(gs[0:2, 0:5])
        ax_before = mlibs.eimshow(
            imagename=self.image_list['test_image'],
            center=self.centre,
            projection='offset',
            rms=rms_test,
            crop=True, box_size=box_size,
            ax=ax_before, fig=fig,
            plot_colorbar=True,
            vmin_factor=1.0, vmax_factor=1.0,
            plot_rms=True,
            add_beam=True,
            add_contours=True)
        ax_before.set_title('Test Image  (before self-cal)', fontsize=13, pad=6)
 
        ax_after = fig.add_subplot(gs[0:2, 5:10])
        ax_after = mlibs.eimshow(
            imagename=self.image_list['selfcal_image'],
            center=self.centre,
            projection='offset',
            rms=self.rms,
            crop=True, box_size=box_size,
            ax=ax_after, fig=fig,
            plot_colorbar=True,
            vmin_factor=1.0, vmax_factor=1.0,
            plot_rms=True,
            add_beam=True,
            add_frequency=True,
            add_contours=True)
        ax_after.set_title('Self-calibrated Image  (final)', fontsize=13, pad=6)
 
        # -- Helper for metric step-plots ------------------------------------─
        def _metric_ax(ax, yvals, ylabel, colour='#1f4e79', yerr=None):
            xs = list(range(len(step_labels)))
            ax.errorbar(xs, yvals,
                        yerr=yerr,
                        marker='o', linestyle='-.', linewidth=1.2, markersize=7,
                        capsize=3, capthick=1.0,
                        markeredgewidth=1.5, markerfacecolor='none',
                        color=colour)
            ax.set_xticks(xs)
            ax.set_xticklabels(step_labels, fontsize=10)
            ax.set_ylabel(ylabel, fontsize=12, labelpad=8)
            ax.set_ylim(np.nanmax([0,np.nanmin(yvals)*0.7]),
                        np.nanmax(yvals)*1.1)
            ax.tick_params(axis='y', labelsize=12)
            ax.grid(axis='y', linestyle=':', linewidth=0.6, alpha=0.7)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
 
        ax_snr = fig.add_subplot(gs[2:4, 0:2])
        _metric_ax(ax_snr, snr_vals,
                   r'SNR $= S_\mathrm{p} / \sigma_\mathrm{rms}$',
                   colour='#1f4e79',
                   yerr=snr_err_vals)
 
        ax_rms = fig.add_subplot(gs[2:4, 2:4])
        _metric_ax(ax_rms, rms_vals,
                   r'$\sigma_\mathrm{rms}\ [\mu\mathrm{Jy\,beam}^{-1}]$',
                   colour='#7b3f00')   # no error bars on rms
 
        ax_flux = fig.add_subplot(gs[2:4, 6:8])
        _metric_ax(ax_flux, flux_vals,
                   r'$S_\nu^\mathrm{mask}\ [\mathrm{mJy}]$',
                   colour='#145a32',
                   yerr=flux_err_vals)
 
        # In-band spectral index \alpha - asymmetric error bars from MCMC when available
        ax_spidx = fig.add_subplot(gs[2:4, 8:10])
        if has_spidx:
            xs = list(range(len(step_labels)))
            # errorbar requires (2, N) array for asymmetric errors
            yerr_spidx = np.array([
                [lo if np.isfinite(lo) else 0.0 for lo in spidx_lo_vals],
                [hi if np.isfinite(hi) else 0.0 for hi in spidx_hi_vals],
            ])
            ax_spidx.errorbar(
                xs, spidx_vals,
                yerr=yerr_spidx,
                marker='o', linestyle='-.', linewidth=1.2, markersize=7,
                capsize=3, capthick=1.0,
                markeredgewidth=1.5, markerfacecolor='none',
                color='#4a0e8f')
            ax_spidx.axhline(0, color='grey', linewidth=0.7, linestyle=':')
            ax_spidx.set_xticks(xs)
            ax_spidx.set_xticklabels(step_labels, fontsize=10)
            ax_spidx.set_ylabel(r'$\alpha$ (in-band)', fontsize=12)
            ax_spidx.tick_params(axis='y', labelsize=12)
            ax_spidx.grid(axis='y', linestyle=':', linewidth=0.6, alpha=0.7)
        else:
            ax_spidx.set_xticks([])
            ax_spidx.set_yticks([])
            ax_spidx.text(0.5, 0.5, 'no sub-band\nspectral fit',
                          ha='center', va='center', fontsize=9, color='grey',
                          transform=ax_spidx.transAxes)
        ax_spidx.spines['top'].set_visible(False)
        ax_spidx.spines['right'].set_visible(False)
 
        # Residual peak convergence: max(R), |min(R)| and rms per step
        ax_res = fig.add_subplot(gs[2:4, 4:6])
        xs = list(range(len(step_labels)))
        has_res = any(np.isfinite(v) for v in res_max_vals)
        if has_res:
            # raw min_residual is negative — plot as-is so asymmetry is visible
            res_min_raw = [_get(k, 'min_residual') * 1e6 for k in step_keys]
            ax_res.plot(xs, res_max_vals,
                        marker='o', linestyle='-.', linewidth=1.2, markersize=7,
                        markeredgewidth=1.5, markerfacecolor='none',
                        color='#c0392b', label=r'$R_\mathrm{max}$')
            ax_res.plot(xs, res_min_raw,
                        marker='s', linestyle='-.', linewidth=1.2, markersize=7,
                        markeredgewidth=1.5, markerfacecolor='none',
                        color='#2980b9', label=r'$R_\mathrm{min}$')
            # +/-rms band as reference: well-calibrated residuals should sit inside
            ax_res.plot(xs,  rms_vals,
                        marker='', linestyle='--', linewidth=2.0,
                        color='black', alpha=0.9, label=r'$\pm\sigma_\mathrm{rms}$')
            ax_res.plot(xs, [-v for v in rms_vals],
                        marker='', linestyle='--', linewidth=2.0,
                        color='black', alpha=0.9)
            ax_res.axhline(0, color='black', linewidth=0.5, linestyle=':')
            ax_res.set_xticks(xs)
            ax_res.set_xticklabels(step_labels, fontsize=10)
            ax_res.set_ylabel(r'$[\mu\mathrm{Jy\,beam}^{-1}]$', fontsize=12)
            ax_res.tick_params(axis='y', labelsize=12)
            ax_res.legend(fontsize=12, framealpha=0.6, loc='best')
            ax_res.grid(axis='y', linestyle=':', linewidth=0.6, alpha=0.7)
        else:
            ax_res.set_xticks([])
            ax_res.set_yticks([])
            ax_res.text(0.5, 0.5, 'no residual\ndata',
                        ha='center', va='center', fontsize=9, color='grey',
                        transform=ax_res.transAxes)
        ax_res.spines['top'].set_visible(False)
        ax_res.spines['right'].set_visible(False)
 
        fig.suptitle(f'{os.path.basename(self.g_name)}  -  self-calibration summary',
                     fontsize=13, y=1.01)
 
        plt.savefig(self.g_name + '_selfcal_results.pdf',
                    dpi=300, bbox_inches='tight')
        plt.clf()
        plt.close()
 
        # -- 6. Mask progression figure --------------------------------------─
        self._run_report_masks()
 
        # -- 7. Flagging statistics ------------------------------------------─
        percentages_over_steps = {}
        for label, step in self.flag_data_steps.items():
            percentages_over_steps[label] = self.calculate_percentages(step)
        categories = ['field', 'observation', 'spw', 'correlation', 'antenna']
        for category in categories:
            self.plot_category_data(percentages_over_steps, category)
            
    # def _run_report_masks(self):
    #     """
    #     Produce <g_name>_selfcal_masks.pdf: one panel per self-cal step showing
    #     the WSClean cleaning mask. All panels share the same centre and zoom,
    #     derived from the bounding box of the final (selfcal_image) mask so that
    #     the zoom-in automatically reflects the true extent of the recovered
    #     emission.
    #     """
    #     # -- Collect mask files for steps that were actually produced --------─
    #     _mask_candidates = [
    #         ('test_image_0',   'Test\nImage'),
    #         ('selfcal_test_0', 'p0'),
    #         ('selfcal_test_1', 'p1/p2'),   # optional
    #         ('selfcal_test_2', 'ap1'),
    #         ('selfcal_image',  'Final'),
    #     ]
    #     mask_panels = []
    #     for key, label in _mask_candidates:
    #         if key not in self.image_list:
    #             continue
    #         mask_path = self.image_list[key].replace('-image.fits', '-image_mask.fits')
    #         if os.path.exists(mask_path):
    #             mask_panels.append((label, mask_path))
 
    #     if not mask_panels:
    #         print('[_run_report_masks] No mask files found - skipping mask figure.')
    #         return
 
    #     # -- Determine zoom from the bounding box of the final mask ----------─
    #     final_mask_path = self.image_list['selfcal_image'].replace(
    #         '-image.fits', '-image_mask.fits')
    #     try:
    #         final_mask_2d = np.squeeze(mlibs.load_fits_data(final_mask_path))
    #         nz = np.argwhere(final_mask_2d > 0)
    #         if len(nz) == 0:
    #             raise ValueError('final mask is empty')
    #         y_min, x_min = nz.min(axis=0)
    #         y_max, x_max = nz.max(axis=0)
    #         extent   = max(y_max - y_min, x_max - x_min)
    #         box_size = max(int(extent * 1.6), 200)   # 60 % padding around mask extent
    #         # Centre on the centroid of the final mask bounding box (x, y order)
    #         centre = (int((x_min + x_max) / 2), int((y_min + y_max) / 2))
    #     except Exception as exc:
    #         print(f'[_run_report_masks] Could not determine mask extent ({exc}); '
    #               'falling back to selfcal_image peak position.')
    #         box_size = 400
    #         centre = mlibs.nd.maximum_position(
    #             np.nan_to_num(mlibs.load_fits_data(self.image_list['selfcal_image']), nan=0))[::-1]
 
    #     # -- Build figure ----------------------------------------------------─
    #     n     = len(mask_panels)
    #     ncols = min(n, 4)                    # at most 4 panels per row
    #     nrows = int(np.ceil(n / ncols))
    #     fig, axes = plt.subplots(nrows, ncols,
    #                              figsize=(4.2 * ncols, 4.5 * nrows),
    #                              squeeze=False)
 
    #     for idx, (label, mask_path) in enumerate(mask_panels):
    #         row, col = divmod(idx, ncols)
    #         ax = axes[row][col]
    #         ax = mlibs.eimshow(
    #             imagename=mask_path,
    #             center=centre,
    #             crop=True, box_size=box_size,
    #             ax=ax, fig=fig,
    #             vmin_factor=0.0, vmax_factor=1.0,
    #             add_contours=False,
    #             plot_colorbar=False,
    #             show_axis='off',
    #             CM='magma')
    #         ax.set_title(label, fontsize=10)
 
    #     # Hide any unused axes (when n is not a multiple of ncols)
    #     for idx in range(n, nrows * ncols):
    #         row, col = divmod(idx, ncols)
    #         axes[row][col].set_visible(False)
 
    #     fig.suptitle(f'{os.path.basename(self.g_name)}  -  cleaning masks per self-cal step',
    #                  fontsize=12, y=1.01)
    #     plt.tight_layout()
    #     plt.savefig(self.g_name + '_selfcal_masks.pdf',
    #                 dpi=300, bbox_inches='tight')
    #     plt.clf()
    #     plt.close()
    
    def _run_report_masks(self):
            """
            Produce <g_name>_selfcal_masks.pdf: one panel per self-cal step showing
            the WSClean cleaning mask, cropped to the largest connected masked region
            of the final (selfcal_image) mask.  All panels share the same centre and
            zoom so the evolution across steps is directly comparable.

            Manual numpy crop is used because WSClean mask FITS files lack
            WCS/pixel-scale headers, which prevents mlibs.eimshow from performing
            its internal crop correctly.
            """
            # -- Collect mask files for steps that were actually produced --------─
            _mask_candidates = [
                ('test_image_0',   'Mask for p0'),
                ('selfcal_test_0', 'Mask for p1'),
                ('selfcal_test_1', 'Mask for p2'),   # optional - only if p2 ran
                ('selfcal_test_2', 'Mask for ap1'),
                ('selfcal_image',  'Final Mask (end of self-cal)'),
            ]
            mask_panels = []
            for key, label in _mask_candidates:
                if key not in self.image_list:
                    continue
                mask_path = self.image_list[key].replace('-image.fits', '-image_mask.fits')
                if os.path.exists(mask_path):
                    mask_panels.append((label, mask_path))

            if not mask_panels:
                print('[_run_report_masks] No mask files found - skipping mask figure.')
                return

            # -- Determine crop from the LARGEST CONNECTED COMPONENT of final mask ─
            final_mask_path = self.image_list['selfcal_image'].replace(
                '-image.fits', '-image_mask.fits')
            try:
                final_mask_2d = np.squeeze(mlibs.load_fits_data(final_mask_path))
                ny, nx        = final_mask_2d.shape
                binary        = (final_mask_2d > 0).astype(np.int32)

                labeled, n_features = mlibs.nd.label(binary)
                if n_features == 0:
                    raise ValueError('final mask is entirely empty')

                # Find the largest component (ignore background label 0)
                comp_sizes    = np.bincount(labeled.ravel())
                comp_sizes[0] = 0
                largest_label = int(comp_sizes.argmax())
                largest_comp  = (labeled == largest_label)

                nz            = np.argwhere(largest_comp)
                y_min, x_min  = nz.min(axis=0)
                y_max, x_max  = nz.max(axis=0)
                extent        = max(y_max - y_min, x_max - x_min)
                half          = max(int(extent * 1.2), 30)   # 40% padding each side
                cy            = int((y_min + y_max) / 2)
                cx            = int((x_min + x_max) / 2)
                print(f'[_run_report_masks] Largest component: {nz.shape[0]} px, '
                    f'extent={extent} px, crop half-width={half} px, '
                    f'centre=({cy}, {cx}),  {n_features} components total.')
            except Exception as exc:
                print(f'[_run_report_masks] Could not determine mask extent ({exc}); '
                    'using centre of image with fallback size.')
                final_mask_2d = np.squeeze(mlibs.load_fits_data(final_mask_path))
                ny, nx        = final_mask_2d.shape
                cy, cx        = ny // 2, nx // 2
                half          = min(ny, nx) // 4

            def _crop(mask_2d):
                """Crop mask_2d to [cy+/-half, cx+/-half], clamped to array bounds."""
                y0 = max(0,  cy - half)
                y1 = min(ny, cy + half)
                x0 = max(0,  cx - half)
                x1 = min(nx, cx + half)
                return mask_2d[y0:y1, x0:x1]

            # -- Build figure ----------------------------------------------------─
            n     = len(mask_panels)
            ncols = min(n, 5)
            nrows = int(np.ceil(n / ncols))
            fig, axes = plt.subplots(nrows, ncols,
                                    figsize=(4.2 * ncols, 4.5 * nrows),
                                    squeeze=False)

            for idx, (label, mask_path) in enumerate(mask_panels):
                row, col = divmod(idx, ncols)
                ax       = axes[row][col]
                try:
                    mask_2d = np.squeeze(mlibs.load_fits_data(mask_path))
                    cropped = _crop(mask_2d)
                    ax.imshow(cropped, origin='lower', cmap='magma',
                            vmin=0, vmax=1, interpolation='nearest')
                except Exception as exc:
                    print(f'[_run_report_masks] Could not load {mask_path}: {exc}')
                    ax.text(0.5, 0.5, 'N/A', ha='center', va='center',
                            transform=ax.transAxes, fontsize=10, color='grey')
                ax.set_title(label, fontsize=10)
                ax.axis('off')

            # Hide any unused axes (when n is not a multiple of ncols)
            for idx in range(n, nrows * ncols):
                row, col = divmod(idx, ncols)
                axes[row][col].set_visible(False)

            fig.suptitle(f'{os.path.basename(self.g_name)}  -  Main/Larger (zoom-in) cleaning mask per self-cal step',
                        fontsize=12, y=1.01)
            plt.tight_layout()
            plt.savefig(self.g_name + '_selfcal_masks.pdf',
                        dpi=300, bbox_inches='tight')
            plt.clf()
            plt.close()


if __name__ == '__main__':
    config = Configuration()
    
    # config.general_settings['do_average_time'] = True
    # config.general_settings['timebin'] = '8s'

    # config.general_settings['do_average_freq'] = True
    # # config.general_settings['channel_width'] = [2]
    # # NOTE: channel_width (the per-SPW mstransform chanbin map) is computed from
    # # the MS itself, right after the Pipeline instance is created (see below).

    # # config.general_settings['new_phasecentre'] = None
    # config.general_settings['new_phasecentre'] = '09:55:50.684 +69.40.43.763'
    # # config.global_parameters['custom_mask'] = "/media/sagauga/void/astronomical-data/M82_v2/eM_C/sc_v12/CY2204/standard/M82_A_K_9216x5120_0.008asec_mask_v2.fits"
    
    # # config.cell_size = None
    # config.cell_size = '0.008arcsec'
    # # config.receiver = None
    # config.imsize = int(1024*9)
    # config.imsizey = int(1024*5)
    
    config.general_settings['do_average_time'] = True
    config.general_settings['timebin'] = '8s'

    config.general_settings['do_average_freq'] = True
    # config.general_settings['channel_width'] = [2]
    #or
    chan_out_avg = 32
    # NOTE: channel_width (the per-SPW mstransform chanbin map) is computed from
    # the MS itself, right after the Pipeline instance is created (see below).

    config.general_settings['new_phasecentre'] = None
    # config.general_settings['new_phasecentre'] = '09:55:50.684 +69.40.43.763'
    # config.global_parameters['custom_mask'] = "/media/sagauga/void/astronomical-data/M82_v2/eM_C/sc_v12/CY2204/standard/M82_A_K_9216x5120_0.008asec_mask_v2.fits"
    
    config.cell_size = None
    # config.cell_size = '0.008arcsec'
    # config.receiver = None
    config.imsize = int(1024*2)
    config.imsizey = int(1024*2)
    
    pipeline = Pipeline(config)

    # Per-SPW channel binning map for the frequency averaging done in
    # _prepare_visibility() (startup step). The Pipeline holds `config` by
    # reference, so setting it here is picked up by the startup step below.
    if config.general_settings['do_average_freq']:
        config.general_settings['channel_width'] = pipeline.get_chan_avg_map(
            config.path + config.vis_name + '.ms', chan_out_avg=chan_out_avg)

    # if 'startup' in config.steps:
    
    pipeline.run_step('startup')
    
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

    if config.general_settings['new_phasecentre'] is not None:
        pipeline._run_phaseshift_final_ms()

    if 'report_results' in config.steps:
        pipeline.run_step('report_results')

    pipeline._run_organise_products()

    # exit()