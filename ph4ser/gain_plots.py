"""
gain_plots.py
=============
Pure-Python gain calibration table plotter.
Replaces plotms-based calibration_table_plot calls in ph4ser.py.

All functions here are self-contained and can be used independently
from a Jupyter notebook or imported by ph4ser.

This was made with the help of Claude. This was adapted from previous testing code. 
All functions were checked for any bug/error. 

Public API
----------
plot_caltable(table_path, yaxis, coloraxis, ...)
    Plot one coloraxis variant. Produces two files:
      <out_dir>/<prefix>_<yaxis>_colorby_<coloraxis>_all[_figNN].png
      <out_dir>/<prefix>_<yaxis>_colorby_<coloraxis>_avgbaseline[_figNN].png

plot_caltable_all_axes(table_path, yaxis, ...)
    Convenience wrapper: calls plot_caltable for spw, antenna, correlation.
"""

import os
import warnings
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
from datetime import datetime, timezone
import datetime as dt

# ── CASA table tool (optional - raises ImportError if not available) ───────────
try:
    from casatools import table as _tb_tool
    _CASATOOLS_OK = True
except ImportError:
    _CASATOOLS_OK = False

# ── Layout constants ───────────────────────────────────────────────────────────
MAX_PANELS  = 6        # max observations per figure
FIG_WIDTH   = 13.0    # fixed total figure width (inches)

_FIG_HEIGHT = {1: 4.8, 2: 8.0, 3: 11.0}   # total height per row count

_FONTS = {
    # nrows: (axis_label, tick, panel_title, suptitle)
    1: (10, 9, 10, 10),
    2: ( 9, 8,  9,  9),
    3: ( 8, 7,  8,  9),
}

# Fraction of figure height reserved above panels for suptitle + legend
_TOP_RESERVED = {1: 0.22, 2: 0.16, 3: 0.12}

_CORR_COLORS = ['#2196F3', '#F44336']   # blue = pol 0 (RR/XX), red = pol 1 (LL/YY)

MJD_EPOCH = datetime(1858, 11, 17, tzinfo=timezone.utc)


# ── Helpers ────────────────────────────────────────────────────────────────────

def _normalise_yaxis(yaxis):
    """Accept 'amp' or 'amplitude' -> returns canonical 'amplitude'."""
    if yaxis in ('amp', 'amplitude'):
        return 'amplitude'
    if yaxis == 'phase':
        return 'phase'
    raise ValueError(f"yaxis must be 'phase' or 'amplitude' (got {yaxis!r})")


def _mjd_sec_to_datetime(mjd_sec):
    return MJD_EPOCH + dt.timedelta(days=mjd_sec / 86400.0)


def _circ_mean_deg(angles_deg):
    """Circular mean of an array of angles in degrees."""
    r = np.deg2rad(angles_deg)
    return np.rad2deg(np.arctan2(np.nanmean(np.sin(r)), np.nanmean(np.cos(r))))


# ── Table I/O ──────────────────────────────────────────────────────────────────

def read_caltable(table_path, yaxis='phase'):
    """
    Read a CASA gain calibration table via casatools.table.

    Parameters
    ----------
    table_path : str   path to the .tb directory
    yaxis      : 'phase' | 'amplitude' | 'amp'

    Returns
    -------
    dict with keys:
        time, obs_id, scan, antenna, spw  - (N,) arrays
        values  - (N, n_pol)  phase [deg] or amplitude
        flag    - (N, n_pol)  bool, True = flagged
        ant_names - list[str]
        n_pol     - int
        yaxis     - str  (normalised)

    Raises
    ------
    RuntimeError  if casatools is not importable
    """
    if not _CASATOOLS_OK:
        raise RuntimeError(
            'casatools is not available - cannot read caltable directly.')

    yaxis = _normalise_yaxis(yaxis)
    tb = _tb_tool()
    tb.open(table_path)
    time    = tb.getcol('TIME')
    obs_id  = tb.getcol('OBSERVATION_ID')
    scan    = tb.getcol('SCAN_NUMBER')
    antenna = tb.getcol('ANTENNA1')
    spw     = tb.getcol('SPECTRAL_WINDOW_ID')
    cparam  = tb.getcol('CPARAM')   # (n_pol, n_chan, N)
    flag    = tb.getcol('FLAG')     # (n_pol, n_chan, N)
    tb.close()

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        flag_3d       = flag.astype(bool)
        cparam_masked = np.where(flag_3d, np.nan + 1j * np.nan, cparam)
        cparam_avg    = np.nanmean(cparam_masked, axis=1)   # (n_pol, N)
        flag_avg      = np.all(flag_3d, axis=1)             # (n_pol, N)

    values = (np.angle(cparam_avg, deg=True) if yaxis == 'phase'
              else np.abs(cparam_avg)).T       # (N, n_pol)
    flag_T = flag_avg.T                        # (N, n_pol)
    n_pol  = values.shape[1]

    tb.open(table_path + '/ANTENNA')
    ant_names = list(tb.getcol('NAME'))
    tb.close()

    return dict(time=time, obs_id=obs_id, scan=scan, antenna=antenna,
                spw=spw, values=values, flag=flag_T,
                ant_names=ant_names, n_pol=n_pol, yaxis=yaxis)


# ── Time-axis compression ──────────────────────────────────────────────────────

def compress_time_axis(time, scan, gap_pad_fraction=0.04):
    """
    Replace inter-scan gaps with a fixed small padding.

    Returns
    -------
    t_comp      : (N,) compressed seconds starting at 0
    scan_bounds : list of (x_start, x_end) per scan (in compressed coords)
    scan_ids    : list of int scan numbers (same order as scan_bounds)
    """
    unique_scans = np.unique(scan)
    total_on_src = sum(
        time[scan == s].max() - time[scan == s].min() for s in unique_scans
    )
    gap_pad = max(gap_pad_fraction * total_on_src, 30.0)

    t_comp = np.zeros_like(time)
    cursor, scan_bounds, scan_ids = 0.0, [], []

    for s in unique_scans:
        mask  = scan == s
        t_rel = time[mask] - time[mask].min()
        t_comp[mask] = cursor + t_rel
        dur   = t_rel.max() if t_rel.size > 1 else 0.0
        scan_bounds.append((cursor, cursor + dur))
        scan_ids.append(int(s))
        cursor += dur + gap_pad

    return t_comp, scan_bounds, scan_ids


def _make_tick_formatter(t_comp, time_raw):
    """FuncFormatter: compressed-time value -> UTC HH:MM string."""
    sort_idx   = np.argsort(t_comp)
    t_c_sorted = t_comp[sort_idx]
    t_r_sorted = time_raw[sort_idx]

    def formatter(x, pos):
        idx = np.clip(np.searchsorted(t_c_sorted, x), 0, len(t_c_sorted) - 1)
        return _mjd_sec_to_datetime(t_r_sorted[idx]).strftime('%H:%M')
    return ticker.FuncFormatter(formatter)


# ── Colour maps ────────────────────────────────────────────────────────────────

def build_color_map(data, coloraxis):
    """
    Build a colour map from the full dataset (all obs) so that colours are
    consistent across panels within a figure.

    Parameters
    ----------
    data      : dict returned by read_caltable (full, not obs-sliced)
    coloraxis : 'spw' | 'antenna' | 'correlation'

    Returns
    -------
    color_map    : dict {group_id: matplotlib color}
    group_labels : dict {group_id: str}
    """
    if coloraxis == 'spw':
        ids    = np.sort(np.unique(data['spw']))
        cmap   = matplotlib.colormaps.get_cmap('tab20').resampled(len(ids))
        cmap_d = {s: cmap(i) for i, s in enumerate(ids)}
        labs   = {s: f'SPW {s}' for s in ids}

    elif coloraxis == 'antenna':
        ids       = np.sort(np.unique(data['antenna']))
        ant_names = data.get('ant_names', [])
        cmap      = matplotlib.colormaps.get_cmap('tab20').resampled(len(ids))
        cmap_d    = {a: cmap(i) for i, a in enumerate(ids)}
        labs      = {a: (ant_names[a] if a < len(ant_names) else f'Ant {a}')
                     for a in ids}

    elif coloraxis == 'correlation':
        n_pol  = data['n_pol']
        ids    = list(range(n_pol))
        pol_labels = ['RR', 'LL'] if n_pol == 2 else [f'P{p}' for p in ids]
        cmap_d = {p: _CORR_COLORS[p % len(_CORR_COLORS)] for p in ids}
        labs   = {p: pol_labels[p] for p in ids}

    else:
        raise ValueError(f"coloraxis must be 'spw', 'antenna', or 'correlation' "
                         f"(got {coloraxis!r})")

    return cmap_d, labs


# ── Axis decoration ────────────────────────────────────────────────────────────

def _decorate_ax(ax, t_comp, time_raw, yaxis, scan_bounds, x_lim,
                 obs_date, fs_label, fs_tick):
    ax.set_xlim(*x_lim)
    ax.xaxis.set_major_formatter(_make_tick_formatter(t_comp, time_raw))
    ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=5, prune='both'))
    plt.setp(ax.get_xticklabels(), rotation=25, ha='right', fontsize=fs_tick)
    ax.tick_params(axis='y', labelsize=fs_tick)
    ax.set_xlabel(f'Time UTC  [{obs_date}]', fontsize=fs_label)
    if yaxis == 'phase':
        ax.set_ylim(-185, 185)
        ax.set_ylabel('Phase (deg)', fontsize=fs_label)
        ax.axhline(0, color='k', lw=0.4, ls='--', alpha=0.4)
    else:
        ax.set_ylabel('Amplitude', fontsize=fs_label)
    ax.grid(True, alpha=0.2, lw=0.4)
    for xs, xe in scan_bounds:
        ax.axvspan(xs, xe, alpha=0.04, color='steelblue', lw=0)


# ── Data-to-axes renderers ─────────────────────────────────────────────────────

def _mean_per_bin(t_comp, mask_group, values, flag, yaxis, n_pol):
    """Mean (circular for phase) over a group mask at each time bin."""
    tc_unique = np.unique(t_comp[mask_group])
    tc_b, val_b = [], []
    for tv in tc_unique:
        at_t = mask_group & (t_comp == tv)
        vals_here = []
        for pol in range(n_pol):
            good = at_t & ~flag[:, pol]
            if good.any():
                vals_here.extend(values[good, pol].tolist())
        if not vals_here:
            continue
        tc_b.append(tv)
        arr = np.array(vals_here)
        val_b.append(_circ_mean_deg(arr) if yaxis == 'phase'
                     else float(np.nanmean(arr)))
    return tc_b, val_b


def _plot_all_antennas(ax, t_comp, data_obs, coloraxis, color_map):
    """Scatter all unflagged solutions, coloured by coloraxis group."""
    spw     = data_obs['spw']
    antenna = data_obs['antenna']
    values  = data_obs['values']
    flag    = data_obs['flag']
    n_pol   = data_obs['n_pol']

    if coloraxis in ('spw', 'antenna'):
        group_arr = spw if coloraxis == 'spw' else antenna
        for gid in np.sort(np.unique(group_arr)):
            c = color_map[gid]
            m = group_arr == gid
            for pol in range(n_pol):
                good = m & ~flag[:, pol]
                if good.any():
                    ax.scatter(t_comp[good], values[good, pol],
                               s=4, color=c, alpha=0.5,
                               rasterized=True, linewidths=0)

    elif coloraxis == 'correlation':
        for pol in range(n_pol):
            good = ~flag[:, pol]
            if good.any():
                ax.scatter(t_comp[good], values[good, pol],
                           s=4, color=color_map[pol], alpha=0.5,
                           rasterized=True, linewidths=0)


def _plot_avg_baseline(ax, t_comp, data_obs, yaxis, coloraxis, color_map):
    """
    Mean per (group, time_bin).

    spw         -> average over antennas + pols, scatter + line
    antenna     -> average over SPWs + pols, lines only (no scatter)
    correlation -> average over antennas + SPWs, scatter + line
    """
    spw     = data_obs['spw']
    antenna = data_obs['antenna']
    values  = data_obs['values']
    flag    = data_obs['flag']
    n_pol   = data_obs['n_pol']

    if coloraxis == 'spw':
        for gid in np.sort(np.unique(spw)):
            c = color_map[gid]
            tc_b, val_b = _mean_per_bin(t_comp, spw == gid,
                                        values, flag, yaxis, n_pol)
            if tc_b:
                ax.scatter(tc_b, val_b, s=8, color=c, alpha=0.8,
                           rasterized=True, linewidths=0)
                ax.plot(tc_b, val_b, color=c, alpha=0.35, lw=0.7)

    elif coloraxis == 'antenna':
        # Lines only - one trace per antenna (mean over SPWs + pols)
        for gid in np.sort(np.unique(antenna)):
            c = color_map[gid]
            tc_b, val_b = _mean_per_bin(t_comp, antenna == gid,
                                        values, flag, yaxis, n_pol)
            if tc_b:
                ax.plot(tc_b, val_b, color=c, alpha=0.75, lw=1.2,
                        rasterized=True)

    elif coloraxis == 'correlation':
        for pol in range(n_pol):
            c = color_map[pol]
            tc_unique = np.unique(t_comp)
            tc_b, val_b = [], []
            for tv in tc_unique:
                at_t = (t_comp == tv) & ~flag[:, pol]
                if not at_t.any():
                    continue
                tc_b.append(tv)
                arr = values[at_t, pol]
                val_b.append(_circ_mean_deg(arr) if yaxis == 'phase'
                             else float(np.nanmean(arr)))
            if tc_b:
                ax.scatter(tc_b, val_b, s=8, color=c, alpha=0.8,
                           rasterized=True, linewidths=0)
                ax.plot(tc_b, val_b, color=c, alpha=0.35, lw=0.7)


# ── Figure construction ────────────────────────────────────────────────────────

def _grid_shape(n):
    """
    Return (nrows, ncols) for n panels (always <= 2 cols).

    n=1->(1,1)  n=2->(1,2)  n=3,4->(2,2)  n=5,6->(3,2)
    """
    if n == 1:      return (1, 1)
    elif n == 2:    return (1, 2)
    elif n <= 4:    return (2, 2)
    else:           return (3, 2)


def _make_figure_with_legend(nrows, ncols, supertitle,
                              color_map, group_labels,
                              coloraxis, averaged):
    """
    Build a figure with the correct size, suptitle, and shared legend
    placed centred below the title without overlapping panels.

    Returns
    -------
    fig        : matplotlib Figure
    axes_flat  : list of Axes (length nrows*ncols, may include hidden extras)
    fonts      : (fs_label, fs_tick, fs_title)
    """
    fig_h   = _FIG_HEIGHT[nrows]
    top_res = _TOP_RESERVED[nrows]
    fs_label, fs_tick, fs_title, fs_sup = _FONTS[nrows]

    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(FIG_WIDTH, fig_h),
                             squeeze=False)

    fig.subplots_adjust(
        top    = 1.0 - top_res,
        bottom = 0.09,
        left   = 0.07,
        right  = 0.97,
        hspace = 0.55,
        wspace = 0.28,
    )

    # Suptitle at top of reserved band
    fig.suptitle(supertitle, fontsize=fs_sup,
                 y=1.0 - top_res * 0.05,
                 va='top')

    # Legend centred in middle of reserved band
    use_line = (coloraxis == 'antenna' and averaged)
    if use_line:
        handles = [Line2D([0], [0], color=color_map[g], lw=2,
                          label=group_labels[g])
                   for g in sorted(color_map)]
    else:
        handles = [Line2D([0], [0], marker='o', color='w',
                          markerfacecolor=color_map[g], markersize=6,
                          label=group_labels[g])
                   for g in sorted(color_map)]

    ncol     = max(1, int(np.ceil(len(handles) / 2)))
    legend_y = 1.0 - top_res * 0.55

    # fig.legend(
    #     handles=handles,
    #     title=coloraxis.capitalize(),
    #     loc='upper center',
    #     bbox_to_anchor=(0.5, legend_y),
    #     ncol=ncol,
    #     fontsize=fs_tick,
    #     title_fontsize=fs_tick,
    #     framealpha=0.8,
    #     handlelength=1.5,
    #     handletextpad=0.4,
    #     columnspacing=0.8,
    # )

    return fig, axes.flatten().tolist(), (fs_label, fs_tick, fs_title)


# ── Public API ─────────────────────────────────────────────────────────────────

def plot_caltable(table_path, yaxis='phase',
                  coloraxis='spw',
                  gap_pad_fraction=0.04,
                  savefig=True, out_dir=None, prefix=''):
    """
    Plot a CASA gain calibration table.

    Observations are grouped into figures of up to 6 panels (3 rows × 2 cols).
    Produces two files per figure batch:
      <out_dir>/<prefix>_<yaxis>_colorby_<coloraxis>_all[_figNN].png
      <out_dir>/<prefix>_<yaxis>_colorby_<coloraxis>_avgbaseline[_figNN].png

    Parameters
    ----------
    table_path      : str   path to the CASA .tb directory
    yaxis           : 'phase' | 'amplitude' | 'amp'
    coloraxis       : 'spw' | 'antenna' | 'correlation'
    gap_pad_fraction: float  inter-scan visual gap (fraction of on-source span)
    savefig         : bool
    out_dir         : str or None  (defaults to table's parent directory)
    prefix          : str  prepended to output filenames

    Returns
    -------
    list of matplotlib Figure objects
    """
    yaxis = _normalise_yaxis(yaxis)
    data  = read_caltable(table_path, yaxis=yaxis)

    unique_obs               = np.unique(data['obs_id'])
    color_map, group_labels  = build_color_map(data, coloraxis)

    if out_dir is None:
        out_dir = os.path.dirname(os.path.abspath(table_path))
    os.makedirs(out_dir, exist_ok=True)

    table_name = os.path.basename(table_path.rstrip('/'))
    sep        = '_' if prefix and not prefix.endswith('_') else ''
    base_fname = f'{prefix}{sep}{yaxis}_colorby_{coloraxis}'

    obs_batches = [
        unique_obs[i:i + MAX_PANELS]
        for i in range(0, len(unique_obs), MAX_PANELS)
    ]
    n_figs      = len(obs_batches)
    all_figures = []

    for fig_idx, batch in enumerate(obs_batches):
        n            = len(batch)
        nrows, ncols = _grid_shape(n)

        fig_suffix = (f'  -  figure {fig_idx+1}/{n_figs}' if n_figs > 1 else '')
        base_sup   = (f'{table_name}\n'
                      f'{yaxis.capitalize()} vs Time  '
                      f'[colour: {coloraxis}]{fig_suffix}')

        fig_all, axes_all, fonts = _make_figure_with_legend(
            nrows, ncols,
            supertitle   = base_sup + '  |  All antennas',
            color_map    = color_map,
            group_labels = group_labels,
            coloraxis    = coloraxis,
            averaged     = False,
        )
        fig_avg, axes_avg, _ = _make_figure_with_legend(
            nrows, ncols,
            supertitle   = base_sup + '  |  Baseline-averaged',
            color_map    = color_map,
            group_labels = group_labels,
            coloraxis    = coloraxis,
            averaged     = True,
        )
        fs_label, fs_tick, fs_title = fonts

        for panel_idx, obs in enumerate(batch):
            ax_all = axes_all[panel_idx]
            ax_avg = axes_avg[panel_idx]

            mask = data['obs_id'] == obs
            data_obs = {k: (v[mask] if isinstance(v, np.ndarray)
                            and v.shape[0] == mask.shape[0] else v)
                        for k, v in data.items()}

            time = data_obs['time']
            scan = data_obs['scan']

            t_comp, scan_bounds, _ = compress_time_axis(
                time, scan, gap_pad_fraction)

            x_min = scan_bounds[0][0]
            x_max = scan_bounds[-1][1]
            x_pad = max((x_max - x_min) * 0.02, 10.0)
            x_lim = (x_min - x_pad, x_max + x_pad)

            obs_date    = _mjd_sec_to_datetime(time.min()).strftime('%Y-%m-%d')
            panel_title = f'OBS {obs}  [{obs_date}]'

            _plot_all_antennas(ax_all, t_comp, data_obs, coloraxis, color_map)
            _decorate_ax(ax_all, t_comp, time, yaxis,
                         scan_bounds, x_lim, obs_date, fs_label, fs_tick)
            ax_all.set_title(panel_title, fontsize=fs_title, pad=3)

            _plot_avg_baseline(ax_avg, t_comp, data_obs, yaxis, coloraxis, color_map)
            _decorate_ax(ax_avg, t_comp, time, yaxis,
                         scan_bounds, x_lim, obs_date, fs_label, fs_tick)
            ax_avg.set_title(panel_title, fontsize=fs_title, pad=3)

        # Hide unused panels when batch < grid capacity
        for ax in axes_all[n:] + axes_avg[n:]:
            ax.set_visible(False)

        if savefig:
            suf = f'_fig{fig_idx+1:02d}' if n_figs > 1 else ''
            for fig, variant in [(fig_all, 'all'), (fig_avg, 'avgbaseline')]:
                fpath = os.path.join(out_dir, f'{base_fname}_{variant}{suf}.png')
                fig.savefig(fpath, dpi=150, bbox_inches='tight')
                print(f'  [gain_plots] Saved -> {fpath}')

        plt.close(fig_all)
        plt.close(fig_avg)
        all_figures.extend([fig_all, fig_avg])

    return all_figures


def plot_caltable_all_axes(table_path, yaxis='phase',
                            gap_pad_fraction=0.04,
                            savefig=True, out_dir=None, prefix=''):
    """
    Convenience wrapper: run plot_caltable for all three coloraxis options
    (spw, antenna, correlation) in sequence.

    Parameters are identical to plot_caltable except coloraxis is omitted.

    Returns
    -------
    dict {coloraxis: list_of_figures}
    """
    results = {}
    # for cax in ('spw', 'antenna', 'correlation'):
    for cax in ('spw', 'correlation'):
        print(f'  [gain_plots] Plotting coloraxis={cax!r} ...')
        try:
            figs = plot_caltable(
                table_path,
                yaxis=yaxis,
                coloraxis=cax,
                gap_pad_fraction=gap_pad_fraction,
                savefig=savefig,
                out_dir=out_dir,
                prefix=prefix,
            )
            results[cax] = figs
        except Exception as e:
            print(f'  [gain_plots] WARNING: coloraxis={cax!r} failed: {e}')
            results[cax] = []
    return results
