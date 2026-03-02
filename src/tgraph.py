#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import numpy as np

# -----------------------
#  0) Python version check
# -----------------------
if sys.version_info[0] < 3:
    print("Please use Python 3 or above.")
    sys.exit(1)

# -----------------------
#  1) set matplotlib backend to GTK3Agg
# -----------------------
import matplotlib
matplotlib.use("GTK3Agg")

import matplotlib.cm as cm
from matplotlib.figure import Figure
from matplotlib.backends.backend_gtk3agg import (
    FigureCanvasGTK3Agg as FigureCanvas
)
from matplotlib.backends.backend_gtk3 import NavigationToolbar2GTK3 as NavigationToolbar


# -----------------------
#  2) import PyGObject
# -----------------------
import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, GLib

# -----------------------
#  3) import tdata.py (ensure script directory is on import path)
# -----------------------
# when installed the script may live in a bindir; we want to load the
# sibling tdata.py that will be placed next to the script.  ``realpath``
# resolves symlinks in case the user installed via a symlink.
_script_dir = os.path.dirname(os.path.realpath(__file__))
if _script_dir not in sys.path:
    sys.path.insert(0, _script_dir)

import tdata

# -----------------------------------------------------------------------------------
#          following section contains data handling, command-line parsing,
#          global variables, etc. from the original tgraph.py
# -----------------------------------------------------------------------------------

tgraph_version = "GTK-ported-1.0"

print("GTK tgraph version:", tgraph_version)

# default columns
xcol = 0
ycol = 1
zcol = 2
vcol = 1

graph_stride = 1


# parse command line and construct filelist; argparse simplifies argument handling
import argparse

def _parse_cols(s: str):
    parts = s.split(":")
    if len(parts) not in (2, 3):
        raise argparse.ArgumentTypeError("-c/--cols must be x:v or x:y:v")
    cols = [int(p) - 1 for p in parts]
    if len(cols) == 2:
        return cols[0], None, cols[1]
    else:
        return cols[0], cols[1], cols[2]


def _parse_range(s: str):
    parts = s.split(":")
    if len(parts) != 2:
        raise argparse.ArgumentTypeError("range must be MIN:MAX")
    return float(parts[0]), float(parts[1])


def parse_cmdline(argv=None):
    """Parse command line arguments.

    Parameters
    ----------
    argv : list or None
        List of arguments (excluding program name). If ``None`` the function
        will use ``sys.argv[1:]``.  Providing ``argv`` is useful for tests.
    """
    parser = argparse.ArgumentParser(
        description="GTK tgraph plotting tool (GTK front-end for tdata)")
    parser.add_argument("-c", "--cols", metavar="X[:Y]:V",
                        help="select columns, e.g. 1:2 or 1:2:3",
                        type=_parse_cols)
    parser.add_argument("-x", "--xrange", metavar="MIN:MAX",
                        help="x axis limits", type=_parse_range)
    parser.add_argument("-y", "--yrange", metavar="MIN:MAX",
                        help="y axis limits", type=_parse_range)
    parser.add_argument("-v", "--vrange", metavar="MIN:MAX",
                        help="value (v) limits", type=_parse_range)
    parser.add_argument("-s", "--stride", metavar="N", type=int,
                        help="data stride", default=1)
    parser.add_argument("-t", "--time-label", metavar="STR",
                        help="time column label", default="time")
    parser.add_argument("-m", "--marker", action="store_true",
                        help="plot markers")
    parser.add_argument("files", nargs="*",
                        help="input files; supports { } and [ ] syntax")
    return parser.parse_args(argv)


# ---------------------------------------------------------------------------
#  glue logic that used to execute on import; factor it into a callable
# ---------------------------------------------------------------------------
def initialize(argv=None):
    """Perform command-line parsing and build global state.

    This function used to run at import time; tests rely on importing the
    module after adjusting ``sys.argv``.  ``argv`` may also be passed
    explicitly for convenience.
    """
    global args, filelist, xcol, ycol, zcol, vcol
    global graph_stride, timelabel_str
    global got_xrange, got_yrange, got_vrange
    global graph_xmin, graph_xmax, graph_ymin, graph_ymax
    global graph_vmin, graph_vmax
    global graph_time, graph_timelist, graph_timeindex

    args = parse_cmdline(argv)

    filelist = tdata.tFileList()
    print("Trying to open files:")

    # apply simple options
    if args.cols is not None:
        xcol, ycol_tmp, vcol_tmp = args.cols
        if ycol_tmp is not None:
            ycol = ycol_tmp
        vcol = vcol_tmp

    if args.marker:
        matplotlib.rcParams['lines.marker'] = 'o'

    graph_stride = args.stride

    timelabel_str = args.time_label.lower()

    # report the columns that will be used (defaults if not overridden)
    print("cols:", xcol+1, ycol+1, ":", vcol+1)

    got_xrange = got_yrange = got_vrange = 0
    if args.xrange is not None:
        graph_xmin, graph_xmax = args.xrange
        got_xrange = 1
    if args.yrange is not None:
        graph_ymin, graph_ymax = args.yrange
        got_yrange = 1
    if args.vrange is not None:
        graph_vmin, graph_vmax = args.vrange
        got_vrange = 1

    # bracket handling state
    openCBrack = inCBrack = openSBrack = inSBrack = endSBrack = 0

    for token in args.files:
        if token == "{":
            openCBrack = 1
            inCBrack = 0
            continue
        elif token == "}":
            openCBrack = 0
            inCBrack = 0
            continue
        elif token == "[":
            openSBrack = 1
            inSBrack = 0
            endSBrack = 0
            continue
        elif token == "]":
            openSBrack = 0
            inSBrack = 0
            endSBrack = 1
            continue

        # treat as filename
        filelist.add(token, timelabel_str)
        print(filelist.file[-1].filename)
        filelist.file[-1].data.set_cols(xcol=xcol, ycol=ycol, zcol=2, vcol=vcol)
        if inSBrack:
            filelist.append_file_i2_to_i1(-2, -1)
        if openSBrack:
            inSBrack = 1
            openSBrack = 0
        if inCBrack:
            filelist.merge_file_i2_into_i1(-2, -1)
        if openCBrack:
            inCBrack = 1
            openCBrack = 0

    # no files error
    if len(filelist.file) == 0:
        print("No files given on command line.\nUse -h for help.")
        sys.exit(1)

    # set columns and initialize for all files
    for i in range(len(filelist.file)):
        filelist.file[i].data.set_cols(xcol=xcol, ycol=ycol, zcol=2, vcol=vcol)

    graph_time = filelist.mintime()
    graph_timelist = filelist.get_timelist()
    graph_timeindex = tdata.geti_from_t(graph_timelist, graph_time)

    # if command-line did not specify x/y/v range, derive from data
    if not hasattr(sys.modules[__name__], 'graph_xmin'):
        graph_xmin = tdata.inf_to_1e300(filelist.minx())
        graph_xmax = tdata.inf_to_1e300(filelist.maxx())
    if not hasattr(sys.modules[__name__], 'graph_ymin'):
        graph_ymin = tdata.inf_to_1e300(filelist.miny())
        graph_ymax = tdata.inf_to_1e300(filelist.maxy())
    if not hasattr(sys.modules[__name__], 'graph_vmin'):
        graph_vmin = tdata.inf_to_1e300(filelist.minv())
        graph_vmax = tdata.inf_to_1e300(filelist.maxv())


# execute initialization whenever the module is imported (or
# for backwards compatibility when it is reloaded during tests).
initialize()

# -----------------------
#   some global state
# -----------------------
graph_3dOn = 0
graph_axis_on = 1
graph_plot_surface = 0
graph_plot_scatter = 0
graph_clear_on_replot = 1
graph_plot_closest_t = 1
graph_plot_grid = 1
# logarithmic scales
graph_xscale = 'linear'
graph_yscale = 'linear'
# colormap
graph_colormap_str = "coolwarm"  # e.g. "jet" or "coolwarm"
graph_colormap = getattr(cm, graph_colormap_str)

# -----------------------
#   dictionary stores related settings
# -----------------------
graph_labelsOn = 0
graph_labels = {
    'title': '',
    'x-axis': '',
    'y-axis': '',
    'v-axis': '',
    'fontsize': matplotlib.rcParams['font.size'],
    'timeformat': '%g'
}

graph_legendOn = 0
graph_legend = {
    'fontsize': matplotlib.rcParams['font.size'],
    'loc': 'upper right',
    'ncol': 1,
    'fancybox': matplotlib.rcParams['legend.fancybox'],
    'shadow': matplotlib.rcParams['legend.shadow']
}
if hasattr(matplotlib.rcParams, 'legend.frameon'):
    graph_legend['frameon'] = matplotlib.rcParams['legend.frameon']
if hasattr(matplotlib.rcParams, 'legend.framealpha'):
    graph_legend['framealpha'] = matplotlib.rcParams['legend.framealpha']
graph_legend['handlelength'] = matplotlib.rcParams['legend.handlelength']

# global pointer to main application window (used as transient parent)
main_window = None

def _get_parent(menuitem):
    """Return a suitable parent window for dialogs.

    If we have stored the main window in :data:`main_window` we always use
    that.  Otherwise fall back to ``menuitem.get_toplevel()`` (the behaviour
    used earlier).  ``menuitem`` may be ``None`` when callbacks aren't
    triggered from a menu.
    """
    if main_window is not None:
        return main_window
    if menuitem:
        try:
            return menuitem.get_toplevel()
        except Exception:
            pass
    return None

graph_settings = {
    'colormap': graph_colormap_str,
    'linewidth': matplotlib.rcParams['lines.linewidth'],
    'antialiased': 0,
    'shade': 1,
    'edgecolor': 'none',
    'stride': graph_stride
}

graph_limits = {
    'xmin': graph_xmin,
    'xmax': graph_xmax,
    'ymin': graph_ymin,
    'ymax': graph_ymax,
    'vmin': graph_vmin,
    'vmax': graph_vmax
}

graph_linecolors = {}
graph_linestyles = {}
graph_linemarkers = {}
graph_linemarkersizes = {}
graph_linewidths = {}
graph_coltrafos = {}

# assign a default color to each file at initialization
def set_graph_globals_for_file_i(filelist, i):
    # similar to tgraph.py's use of color_cycle
    if matplotlib.__version__ < '1.5.1':
        # older matplotlib versions used axes.color_cycle
        color_cycle = matplotlib.rcParams.get('axes.color_cycle', [
            '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
            '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
            '#bcbd22', '#17becf'
        ])
        c = color_cycle[i % len(color_cycle)]
    else:
        # axes.prop_cycle
        cycle_list = list(matplotlib.rcParams['axes.prop_cycle'])
        c = cycle_list[i % len(cycle_list)]['color']

    graph_linecolors[f"#{i}"] = c
    graph_linestyles[f"#{i}"] = "-"
    graph_linemarkers[f"#{i}"] = matplotlib.rcParams['lines.marker']
    graph_linemarkersizes[f"#{i}"] = matplotlib.rcParams['lines.markersize']
    graph_linewidths[f"#{i}"] = ""
    graph_coltrafos[f"#{i}"] = ""
    # store the file's name in graph_legend
    graph_legend[f"#{i}"] = filelist.file[i].name

for i in range(len(filelist.file)):
    set_graph_globals_for_file_i(filelist, i)

print("(tmin, tmax) =", filelist.mintime(), filelist.maxtime())
print("(xmin, xmax) =", graph_xmin, graph_xmax)
print("(ymin, ymax) =", graph_ymin, graph_ymax)
print("(vmin, vmax) =", graph_vmin, graph_vmax)
print("stride =", graph_stride)

# -----------------------------------------------------------------------------------
#   2D/3D plotting functions
# -----------------------------------------------------------------------------------
def axplot2d_at_time(filelist, ax, t):
    if graph_clear_on_replot:
        ax.clear()
    # restore axis scales (clear() resets them)
    ax.set_xscale(graph_xscale)
    ax.set_yscale(graph_yscale)
    for i in range(len(filelist.file)):
        f = filelist.file[i]
        xs = f.data.getx(t, graph_plot_closest_t)
        ys = f.data.getv(t, graph_plot_closest_t)  # ycol -> vcol
        color = graph_linecolors.get(f"#{i}", "blue")
        style = graph_linestyles.get(f"#{i}", "-")
        marker = graph_linemarkers.get(f"#{i}", "")
        markersz = graph_linemarkersizes.get(f"#{i}", 6)
        lw = graph_linewidths.get(f"#{i}", "")
        if lw != "":
            lw = float(lw)

        if graph_plot_scatter:
            if not marker:
                marker = "o"
            ax.scatter(xs, ys, label=f.name, color=color, marker=marker)
        else:
            ax.plot(xs, ys, label=f.name, color=color,
                    linestyle=style, marker=marker, markersize=markersz,
                    linewidth=lw if lw else None)

    ax.set_xlim(graph_limits['xmin'], graph_limits['xmax'])
    ax.set_ylim(graph_limits['vmin'], graph_limits['vmax'])

    if graph_labelsOn:
        ax.set_xlabel(graph_labels['x-axis'], fontsize=graph_labels['fontsize'])
        ax.set_ylabel(graph_labels['v-axis'], fontsize=graph_labels['fontsize'])

        # title and time
        title = graph_labels['title']
        tf = graph_labels['timeformat']
        if tf:
            tstr = tf % t
            # append time to title, right-aligned
            title += "    " + tstr
        ax.set_title(title)
    else:
        # hide labels/title by clearing text
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_title("")

    if graph_plot_grid:
        ax.grid(True)
    if graph_legendOn:
        ax.legend(fontsize=graph_legend['fontsize'],
                  loc=graph_legend['loc'],
                  ncol=graph_legend['ncol'])

    if not graph_axis_on:
        ax.set_axis_off()

def axplot3d_at_time(filelist, ax, t):
    """
    Like the original tgraph.py: x and y form surface coordinates while v
    serves as the z value.
    """
    from mpl_toolkits.mplot3d import axes3d
    if graph_clear_on_replot:
        ax.clear()

    for i in range(len(filelist.file)):
        f = filelist.file[i]
        # read blocks to determine data line breaks
        blocks = f.data.getblocks(t)
        # reshape
        x = np.reshape(f.data.getx(t, graph_plot_closest_t), (blocks, -1))
        y = np.reshape(f.data.gety(t, graph_plot_closest_t), (blocks, -1))
        z = np.reshape(f.data.getv(t, graph_plot_closest_t), (blocks, -1))

        color = graph_linecolors.get(f"#{i}", "blue")
        marker = graph_linemarkers.get(f"#{i}", "")
        lw = graph_linewidths.get(f"#{i}", "")
        if lw != "":
            lw = float(lw)

        if graph_plot_surface:
            if graph_settings['colormap']:
                colormap_obj = getattr(cm, graph_settings['colormap'])
                ax.plot_surface(x, y, z,
                                rstride=graph_stride, cstride=graph_stride,
                                cmap=colormap_obj,
                                linewidth=float(graph_settings['linewidth']),
                                antialiased=int(graph_settings['antialiased']),
                                shade=int(graph_settings['shade']),
                                edgecolor=graph_settings['edgecolor'])
            else:
                ax.plot_surface(x, y, z,
                                rstride=graph_stride, cstride=graph_stride,
                                color=color,
                                linewidth=float(graph_settings['linewidth']),
                                antialiased=int(graph_settings['antialiased']),
                                shade=int(graph_settings['shade']),
                                edgecolor=graph_settings['edgecolor'])
        else:
            if graph_plot_scatter:
                if not marker:
                    marker = "o"
                ax.scatter(x, y, z, color=color, marker=marker, label=f.name)
            else:
                ax.plot_wireframe(x, y, z,
                                  rstride=graph_stride, cstride=graph_stride,
                                  color=color, linewidth=lw if lw else None,
                                  label=f.name)

    ax.set_xlim(graph_limits['xmin'], graph_limits['xmax'])
    ax.set_ylim(graph_limits['ymin'], graph_limits['ymax'])
    ax.set_zlim(graph_limits['vmin'], graph_limits['vmax'])

    if graph_labelsOn:
        ax.set_xlabel(graph_labels['x-axis'], fontsize=graph_labels['fontsize'])
        ax.set_ylabel(graph_labels['y-axis'], fontsize=graph_labels['fontsize'])
        ax.set_zlabel(graph_labels['v-axis'], fontsize=graph_labels['fontsize'])

        title = graph_labels['title']
        tf = graph_labels['timeformat']
        if tf:
            tstr = tf % t
            title += "    " + tstr
        ax.set_title(title)
    else:
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_zlabel("")
        ax.set_title("")

    # note: 3D surfaces do not support a legend
    if graph_legendOn and not graph_plot_surface:
        ax.legend()

    if not graph_axis_on:
        ax.set_axis_off()

graph_delay = 1  # added: controls playback delay (ms)

# added: functions to increment/decrement time, jump to first/last frame, playback controls
def min_graph_time():
    global graph_time, graph_timeindex
    graph_timeindex = 0
    graph_time = graph_timelist[graph_timeindex]
    replot()

def max_graph_time():
    global graph_time, graph_timeindex
    graph_timeindex = len(graph_timelist) - 1
    graph_time = graph_timelist[graph_timeindex]
    replot()

def inc_graph_time():
    global graph_time, graph_timeindex
    if graph_timeindex < len(graph_timelist) - 1:
        graph_timeindex += 1
    graph_time = graph_timelist[graph_timeindex]
    replot()

def dec_graph_time():
    global graph_time, graph_timeindex
    if graph_timeindex > 0:
        graph_timeindex -= 1
    graph_time = graph_timelist[graph_timeindex]
    replot()

# use a global variable to store timer ID for playback control
playback_id = None

def play_graph_time():
    """
    Each call increments the time frame and replots.  If not at the last
    frame, schedule another timeout; otherwise stop playback.
    """
    global graph_timeindex, playback_id, graph_delay

    inc_graph_time()  # inc_graph_time already replots
    if graph_timeindex < len(graph_timelist) - 1:
        # continue playback
        playback_id = GLib.timeout_add(int(graph_delay), play_graph_time)
    else:
        playback_id = None  # playback finished

def start_play_graph_time():
    """
    Start playback from current time frame.  If already playing, stop first.
    """
    global playback_id
    if playback_id:
        # if already playing -> stop
        GLib.source_remove(playback_id)
        playback_id = None
    else:
        # if not playing -> start
        if graph_timeindex >= len(graph_timelist) - 1:
            # if at last frame, go back to start
            min_graph_time()
        playback_id = GLib.timeout_add(int(graph_delay), play_graph_time)

# set graph_time when user presses Enter or leaves the Entry
def set_graph_time_from_entry(entry):
    global graph_time, graph_timeindex
    try:
        val = float(entry.get_text().strip())
        # find the time index closest to val
        graph_timeindex = tdata.geti_from_t(graph_timelist, val)
        graph_time = graph_timelist[graph_timeindex]
        replot()
    except ValueError:
        pass

def set_delay_from_entry(entry):
    global graph_delay
    try:
        val = float(entry.get_text().strip())
        graph_delay = val
    except ValueError:
        pass

# update display of the time Entry
def update_time_entry(entry):
    global graph_time
    entry.set_text(f"{graph_time}")

# -----------------------------------------------------------------------------------
#   some utility functions (replot, toggles, etc. from original tgraph)
# -----------------------------------------------------------------------------------
def setup_axes(fig, is3d, old_ax=None):
    """
    Recreate Axes on a figure (2D or 3D).  Remove the old axes if provided.
    """
    if old_ax is not None:
        fig.delaxes(old_ax)
    if not is3d:
        return fig.add_subplot(111)
    else:
        from mpl_toolkits.mplot3d import Axes3D
        return fig.add_subplot(111, projection='3d')

# global references for menu callbacks
global_fig = None
global_ax = None
global_canvas = None

def replot():
    global global_fig, global_ax, global_canvas
    global graph_time, graph_3dOn

    # in headless/test scenarios we may not have canvas/axes set yet
    if global_ax is None or global_canvas is None:
        return

    if graph_3dOn == 0:
        axplot2d_at_time(filelist, global_ax, graph_time)
    else:
        axplot3d_at_time(filelist, global_ax, graph_time)

    global_canvas.draw()

# menu callbacks: toggles
def toggle_timeframe_update(_menuitem):
    global graph_clear_on_replot
    graph_clear_on_replot = 0 if graph_clear_on_replot else 1
    replot()

def toggle_timeframe_closest(_menuitem):
    global graph_plot_closest_t
    graph_plot_closest_t = 0 if graph_plot_closest_t else 1
    replot()

def toggle_axis_on(_menuitem):
    global graph_axis_on
    graph_axis_on = 0 if graph_axis_on else 1
    replot()

def toggle_log_xscale(_menuitem):
    global global_ax, graph_xscale
    # flip stored state
    graph_xscale = 'log' if graph_xscale == 'linear' else 'linear'
    if global_ax is not None:
        global_ax.set_xscale(graph_xscale)
    replot()

def toggle_log_yscale(_menuitem):
    global global_ax, graph_yscale
    graph_yscale = 'log' if graph_yscale == 'linear' else 'linear'
    if global_ax is not None:
        global_ax.set_yscale(graph_yscale)
    replot()

def toggle_grid(_menuitem):
    global graph_plot_grid
    graph_plot_grid = 0 if graph_plot_grid else 1
    replot()

def toggle_line_scatter(_menuitem):
    global graph_plot_surface, graph_plot_scatter
    if graph_plot_surface:
        # if it was surface mode before, turn it off first
        graph_plot_surface = 0
    graph_plot_scatter = 0 if graph_plot_scatter else 1
    replot()

def toggle_2d_3d(_menuitem):
    global graph_3dOn, global_ax, global_fig
    graph_3dOn = 0 if graph_3dOn else 1
    global_ax = setup_axes(global_fig, graph_3dOn, global_ax)
    replot()

def toggle_surface(_menuitem):
    global graph_plot_surface, graph_plot_scatter
    if graph_plot_surface == 1:
        graph_plot_surface = 0
    else:
        graph_plot_surface = 1
        graph_plot_scatter = 0
    replot()

def toggle_labels(_menuitem):
    global graph_labelsOn
    graph_labelsOn = 0 if graph_labelsOn else 1
    # if labels are off we could clear them; instead simplify to hide/show only
    replot()

def toggle_legend(_menuitem):
    global graph_legendOn
    graph_legendOn = 0 if graph_legendOn else 1
    replot()

# -----------------------------------------------------------------------------------
#   A 'dialog' used to edit several key->value pairs
# -----------------------------------------------------------------------------------
class GTKDialogKeyValue(Gtk.Dialog):
    """
    Similar to the original WTdialog; used to edit string/numeric values in a
    dictionary.
    """
    def __init__(self, title, keyvalues: dict, parent=None):
        super().__init__(title=title, transient_for=parent, modal=True)
        # force dialog type hint so Wayland compositors treat it as a floating dialog
        try:
            self.set_type_hint(Gtk.WindowTypeHint.DIALOG)
        except Exception:
            pass
        self.set_default_size(400, 300)

        self.keyvalues = keyvalues
        self.entries = {}

        vbox = self.get_content_area()
        grid = Gtk.Grid(column_spacing=8, row_spacing=6, margin=10)
        vbox.add(grid)

        row = 0
        # create a row (Label + Entry) for each key in the dict
        for key in sorted(keyvalues.keys()):
            label = Gtk.Label(label=key, xalign=1.0)
            entry = Gtk.Entry()
            entry.set_text(str(keyvalues[key]))
            self.entries[key] = entry
            grid.attach(label, 0, row, 1, 1)
            grid.attach(entry, 1, row, 1, 1)
            row += 1

        # bottom button
        self.ok_button = Gtk.Button(label="Apply")
        self.ok_button.connect("clicked", self.on_apply)
        grid.attach(self.ok_button, 0, row, 2, 1)

        self.show_all()

    def on_apply(self, button):
        # fetch values from the dialog entries
        for k, e in self.entries.items():
            text = e.get_text().strip()
            self.keyvalues[k] = text
        self.destroy()

# -----------------------------------------------------------------------------------
#   menu callbacks: various input dialogs -> update globals -> replot
# -----------------------------------------------------------------------------------
def input_graph_limits(_menuitem):
    global graph_limits
    parent = _get_parent(_menuitem)
    dlg = GTKDialogKeyValue("Edit Limits", graph_limits, parent=parent)
    dlg.run()  # block until dialog closes
    dlg.destroy()

    # write edited results back to globals
    # convert to float
    graph_limits['xmin'] = float(graph_limits['xmin'])
    graph_limits['xmax'] = float(graph_limits['xmax'])
    graph_limits['ymin'] = float(graph_limits['ymin'])
    graph_limits['ymax'] = float(graph_limits['ymax'])
    graph_limits['vmin'] = float(graph_limits['vmin'])
    graph_limits['vmax'] = float(graph_limits['vmax'])

    replot()

def input_graph_labels(_menuitem):
    global graph_labels, graph_labelsOn
    parent = _get_parent(_menuitem)
    dlg = GTKDialogKeyValue("Edit Labels", graph_labels, parent=parent)
    dlg.run()
    dlg.destroy()

    # convert some values
    matplotlib.rcParams['font.size'] = float(graph_labels['fontsize'])
    # enable display
    graph_labelsOn = 1
    replot()

def input_graph_legend(_menuitem):
    global graph_legend, graph_legendOn, filelist
    parent = _get_parent(_menuitem)
    dlg = GTKDialogKeyValue("Edit Legend", graph_legend, parent=parent)
    dlg.run()
    dlg.destroy()

    # handle #0, #1, #2... corresponding file names
    for i in range(len(filelist.file)):
        k = f"#{i}"
        if k in graph_legend:
            filelist.file[i].name = graph_legend[k]
    # parse loc
    s = str(graph_legend['loc']).strip()
    if s.isdigit():
        graph_legend['loc'] = int(s)
    else:
        # try parsing as "x, y"
        if "," in s:
            parts = s.replace("(", "").replace(")", "").split(",")
            graph_legend['loc'] = (float(parts[0]), float(parts[1]))
    graph_legend['ncol'] = int(graph_legend['ncol'])

    # enable display
    graph_legendOn = 1
    replot()

def input_graph_settings(_menuitem):
    global graph_settings, graph_stride, graph_colormap, graph_colormap_str
    parent = _get_parent(_menuitem)
    dlg = GTKDialogKeyValue("Graph Settings", graph_settings, parent=parent)
    dlg.run()
    dlg.destroy()

    # apply changes
    matplotlib.rcParams['lines.linewidth'] = float(graph_settings['linewidth'])
    graph_stride = int(graph_settings['stride'])
    cmap_name = graph_settings['colormap']
    graph_colormap_str = cmap_name
    if cmap_name:
        try:
            graph_colormap = getattr(cm, cmap_name)
        except:
            pass

    replot()

def input_graph_xcolumns(_menuitem):
    """
    Used to modify the x-column for each file
    """
    parent = _get_parent(_menuitem)
    # prepare a dict
    d = {}
    for i in range(len(filelist.file)):
        d[f"#{i}"] = str(filelist.file[i].data.get_xcol0() + 1)
    dlg = GTKDialogKeyValue("Select x-Columns", d, parent=parent)
    dlg.run()
    dlg.destroy()
    # apply
    for i in range(len(filelist.file)):
        newcol = int(d[f"#{i}"]) - 1
        filelist.file[i].data.set_xcols(newcol)

    # recompute x min/max
    global graph_xmin, graph_xmax
    graph_xmin = tdata.inf_to_1e300(filelist.minx())
    graph_xmax = tdata.inf_to_1e300(filelist.maxx())
    graph_limits['xmin'] = graph_xmin
    graph_limits['xmax'] = graph_xmax
    replot()

def input_graph_ycolumns(_menuitem):
    parent = _get_parent(_menuitem)
    d = {}
    for i in range(len(filelist.file)):
        d[f"#{i}"] = str(filelist.file[i].data.get_ycol0() + 1)
    dlg = GTKDialogKeyValue("Select y-Columns", d, parent=parent)
    dlg.run()
    dlg.destroy()
    for i in range(len(filelist.file)):
        newcol = int(d[f"#{i}"]) - 1
        filelist.file[i].data.set_ycols(newcol)

    global graph_ymin, graph_ymax
    graph_ymin = tdata.inf_to_1e300(filelist.miny())
    graph_ymax = tdata.inf_to_1e300(filelist.maxy())
    graph_limits['ymin'] = graph_ymin
    graph_limits['ymax'] = graph_ymax
    replot()

def input_graph_vcolumns(_menuitem):
    parent = _get_parent(_menuitem)
    d = {}
    for i in range(len(filelist.file)):
        d[f"#{i}"] = str(filelist.file[i].data.get_vcol0() + 1)
    dlg = GTKDialogKeyValue("Select v-Columns", d, parent=parent)
    dlg.run()
    dlg.destroy()
    for i in range(len(filelist.file)):
        newcol = int(d[f"#{i}"]) - 1
        filelist.file[i].data.set_vcols(newcol)

    global graph_vmin, graph_vmax
    graph_vmin = tdata.inf_to_1e300(filelist.minv())
    graph_vmax = tdata.inf_to_1e300(filelist.maxv())
    graph_limits['vmin'] = graph_vmin
    graph_limits['vmax'] = graph_vmax
    replot()

def input_graph_linecolors(_menuitem):
    parent = _get_parent(_menuitem)
    d = {}
    for i in range(len(filelist.file)):
        d[f"#{i}"] = graph_linecolors[f"#{i}"]
    dlg = GTKDialogKeyValue("Edit Line Colors", d, parent=parent)
    dlg.run()
    dlg.destroy()
    for i in range(len(filelist.file)):
        graph_linecolors[f"#{i}"] = d[f"#{i}"]
    replot()

def input_graph_linestyles(_menuitem):
    parent = _get_parent(_menuitem)
    d = {}
    for i in range(len(filelist.file)):
        d[f"#{i}"] = graph_linestyles[f"#{i}"]
    dlg = GTKDialogKeyValue("Edit Line Styles", d, parent=parent)
    dlg.run()
    dlg.destroy()
    for i in range(len(filelist.file)):
        graph_linestyles[f"#{i}"] = d[f"#{i}"]
    replot()

def input_graph_linemarkers(_menuitem):
    parent = _get_parent(_menuitem)
    d = {}
    for i in range(len(filelist.file)):
        d[f"#{i}"] = graph_linemarkers[f"#{i}"]
    dlg = GTKDialogKeyValue("Edit Line Markers", d, parent=parent)
    dlg.run()
    dlg.destroy()
    for i in range(len(filelist.file)):
        graph_linemarkers[f"#{i}"] = d[f"#{i}"]
    replot()

def input_graph_linemarkersizes(_menuitem):
    parent = _get_parent(_menuitem)
    d = {}
    for i in range(len(filelist.file)):
        d[f"#{i}"] = str(graph_linemarkersizes[f"#{i}"])
    dlg = GTKDialogKeyValue("Edit Line Markersizes", d, parent=parent)
    dlg.run()
    dlg.destroy()
    for i in range(len(filelist.file)):
        graph_linemarkersizes[f"#{i}"] = float(d[f"#{i}"])
    replot()

def input_graph_linewidths(_menuitem):
    parent = _get_parent(_menuitem)
    d = {}
    for i in range(len(filelist.file)):
        d[f"#{i}"] = str(graph_linewidths[f"#{i}"])
    dlg = GTKDialogKeyValue("Edit Line Widths", d, parent=parent)
    dlg.run()
    dlg.destroy()
    for i in range(len(filelist.file)):
        val = d[f"#{i}"].strip()
        graph_linewidths[f"#{i}"] = val
    replot()

def input_graph_coltrafos(_menuitem):
    """
    Like the original script, allows expression-based column transforms.
    This example keeps only the UI; the actual transform_col behavior is not
    demonstrated here.
    """
    parent = _get_parent(_menuitem)
    d = {}
    for i in range(len(filelist.file)):
        d[f"#{i}"] = graph_coltrafos[f"#{i}"]
    dlg = GTKDialogKeyValue("Transform Columns", d, parent=parent)
    dlg.run()
    dlg.destroy()
    for i in range(len(filelist.file)):
        expr = d[f"#{i}"].strip()
        graph_coltrafos[f"#{i}"] = expr
        if expr != "":
            print(f"Transform file #{i}: {expr}")
            filelist.file[i].data.transform_col(expr, c_index_shift=1)
    replot()

def about_dialog(_menuitem):
    parent = _get_parent(_menuitem)
    md = Gtk.MessageDialog(
        transient_for=parent,
        modal=True,
        text=f"tgraph {tgraph_version}\n\n"
                       "Copyright (C) 2015 Wolfgang Tichy.\n"
                       "Ported to GTK by Yingjie Wang.\n",
        buttons=Gtk.ButtonsType.OK
    )
    # hint as a dialog for Wayland
    try:
        md.set_type_hint(Gtk.WindowTypeHint.DIALOG)
    except Exception:
        pass
    md.run()
    md.destroy()

def _find_help_file():
    """Return the path to ``tgraph.txt`` if it exists.

    We prefer a copy next to ``tdata.py`` (useful during development); if
    that file is missing look under ``$(prefix)/share/doc/tgraph`` where
    ``make install`` now places it.  ``prefix`` is assumed to be the parent
    directory of the directory containing this script (i.e. ``bin``).
    """
    # development / source tree location
    candidate = os.path.join(os.path.dirname(tdata.__file__), "tgraph.txt")
    if os.path.exists(candidate):
        return candidate

    # installed location: ../share/doc/tgraph/tgraph.txt
    script_dir = os.path.dirname(os.path.realpath(__file__))
    prefix = os.path.dirname(script_dir)
    candidate = os.path.join(prefix, "share", "doc", "tgraph", "tgraph.txt")
    if os.path.exists(candidate):
        return candidate

    return None


def help_dialog(_menuitem):
    parent = _get_parent(_menuitem)
    helptext = "tgraph.txt not found!"
    path = _find_help_file()
    if path is not None:
        with open(path, "r") as f:
            helptext = f.read()

    # display using TextView
    dialog = Gtk.Dialog(title="tgraph help", transient_for=parent, modal=True)
    try:
        dialog.set_type_hint(Gtk.WindowTypeHint.DIALOG)
    except Exception:
        pass
    dialog.set_default_size(600, 400)
    box = dialog.get_content_area()
    scrolled_win = Gtk.ScrolledWindow()
    scrolled_win.set_hexpand(True)
    scrolled_win.set_vexpand(True)
    textview = Gtk.TextView()
    textview.set_editable(False)
    textbuffer = textview.get_buffer()
    textbuffer.set_text(helptext)
    scrolled_win.add(textview)
    box.pack_start(scrolled_win, True, True, 0)

    button_ok = Gtk.Button(label="Close")
    button_ok.connect("clicked", lambda w: dialog.destroy())
    box.pack_end(button_ok, False, False, 5)

    dialog.show_all()
    dialog.run()
    dialog.destroy()

# -----------------------------
# Added: Open File dialog function
# -----------------------------
def open_file_dialog(_menuitem):
    """
    Allow the user to choose multiple new files via a FileChooserDialog and
    add them to the filelist.
    """
    parent = _get_parent(_menuitem)
    dialog = Gtk.FileChooserDialog(
        title="Open File(s)",
        parent=parent,
        action=Gtk.FileChooserAction.OPEN
    )
    try:
        dialog.set_type_hint(Gtk.WindowTypeHint.DIALOG)
    except Exception:
        pass
    dialog.add_buttons(
        Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
        Gtk.STOCK_OPEN,   Gtk.ResponseType.OK
    )
    dialog.set_select_multiple(True)  # allow multiple selection

    # filters can be added here if desired
    # file_filter = Gtk.FileFilter()
    # file_filter.set_name("All files")
    # file_filter.add_pattern("*")
    # dialog.add_filter(file_filter)

    response = dialog.run()
    if response == Gtk.ResponseType.OK:
        # retrieve all selected filenames
        filenames = dialog.get_filenames()
        print("User selected:", filenames)

        # for each file, hand it over to filelist
        from math import isinf
        for fname in filenames:
            if not os.path.isfile(fname):
                continue
            # assume default timelabel_str
            timelabel_str = "time"
            filelist.add(fname, timelabel_str)
            print("Added:", filelist.file[-1].filename)
            # set default columns
            # xcol, ycol, zcol, vcol etc. may be taken from globals
            filelist.file[-1].data.set_cols(xcol=0, ycol=1, zcol=2, vcol=1)

        # update graph_timelist, graph_time, etc.
        global graph_timelist, graph_time, graph_timeindex
        graph_timelist = filelist.get_timelist()
        # if initially empty, we can now compute a new time range
        # graph_time = filelist.mintime()
        # graph_timeindex = tdata.geti_from_t(graph_timelist, graph_time)

        # update x/y/v ranges
        global graph_xmin, graph_xmax, graph_ymin, graph_ymax, graph_vmin, graph_vmax
        from tdata import inf_to_1e300
        graph_xmin = inf_to_1e300(filelist.minx())
        graph_xmax = inf_to_1e300(filelist.maxx())
        graph_ymin = inf_to_1e300(filelist.miny())
        graph_ymax = inf_to_1e300(filelist.maxy())
        graph_vmin = inf_to_1e300(filelist.minv())
        graph_vmax = inf_to_1e300(filelist.maxv())

        # sync into graph_limits dictionary
        graph_limits['xmin'] = graph_xmin
        graph_limits['xmax'] = graph_xmax
        graph_limits['ymin'] = graph_ymin
        graph_limits['ymax'] = graph_ymax
        graph_limits['vmin'] = graph_vmin
        graph_limits['vmax'] = graph_vmax

        # assign line colors/names for the new files
        for i in range(len(filelist.file)):
            # if not initialized previously, call set_graph_globals_for_file_i
            # simplest approach: apply set_graph_globals_for_file_i to the most
            # recently added files
            pass

        # redraw
        replot()

    dialog.destroy()

# -----------------------------
# Added: save sequence of frames function
# -----------------------------
def save_movieframes_dialog(_menuitem):
    """
    Pop up a "Save As" dialog to pick an output filename (e.g. "frame.png").
    Then save every time frame as "basename_{index}.ext".
    """
    parent = _get_parent(_menuitem)
    dialog = Gtk.FileChooserDialog(
        title="Save Movie Frames",
        parent=parent,
        action=Gtk.FileChooserAction.SAVE
    )
    try:
        dialog.set_type_hint(Gtk.WindowTypeHint.DIALOG)
    except Exception:
        pass
    dialog.add_buttons(
        Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
        Gtk.STOCK_SAVE,   Gtk.ResponseType.OK
    )
    dialog.set_do_overwrite_confirmation(True)

    response = dialog.run()
    if response == Gtk.ResponseType.OK:
        outname = dialog.get_filename()
        # outname might be /home/user/frame.png
        print("Saving frames to base name:", outname)

        # handle extension and base name
        base, ext = os.path.splitext(outname)
        if not ext:
            # if user didn't supply an extension, default to .png
            ext = ".png"

        # all time frames
        global graph_timeindex, graph_time, graph_timelist
        i1 = 0
        i2 = len(graph_timelist)
        # compute number of digits needed for zero padding
        digits = len(str(i2-1))
        fmt_str = f"%0{digits}d"

        # save current index
        old_index = graph_timeindex
        old_time = graph_time

        # loop and save each frame
        for idx in range(i2):
            graph_timeindex = idx
            graph_time = graph_timelist[idx]
            replot()  # allow the display to update for the current time

            # construct filename
            frame_idx = fmt_str % idx
            filename = f"{base}_{frame_idx}{ext}"
            global_canvas.print_figure(filename)
            if idx == 0:
                first_file = filename
            if idx == i2-1:
                last_file = filename

        # restore previous time
        graph_timeindex = old_index
        graph_time = old_time
        replot()

        # inform the user
        info_dialog = Gtk.MessageDialog(
            parent=None,
            flags=0,
            message_format=f"Frames saved:\n  {first_file} ... {last_file}",
            buttons=Gtk.ButtonsType.OK
        )
        info_dialog.run()
        info_dialog.destroy()

    dialog.destroy()

# -----------------------------------------------------------------------------------
#   GTK main window
# -----------------------------------------------------------------------------------
class GTKGraphWindow(Gtk.Window):
    def __init__(self):
        super().__init__(title="GTK tgraph")
        global main_window
        main_window = self
        self.set_default_size(900, 700)

        # create a vertical layout
        vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=0)
        self.add(vbox)

        # 1) menubar
        menubar = self.build_menubar()
        vbox.pack_start(menubar, False, False, 0)

        # 2) add a row of time playback controls at the top
        time_controls_box = self.build_time_controls()
        vbox.pack_start(time_controls_box, False, False, 4)

        # 3) matplotlib Figure + Canvas
        global global_fig, global_ax, global_canvas
        global_fig = Figure(figsize=(6,5), dpi=100)
        # global_ax = setup_axes(global_fig, graph_3dOn, None)
        global_ax = global_fig.add_subplot(111)
        global_canvas = FigureCanvas(global_fig)
        global_canvas.set_hexpand(True)
        global_canvas.set_vexpand(True)

        vbox.pack_start(global_canvas, True, True, 0)

        toolbar = NavigationToolbar(global_canvas)
        vbox.pack_start(toolbar, False, False, 0)

        # initial drawing
        replot()
        self.show_all()

    def build_menubar(self):
        menubar = Gtk.MenuBar()

        # File
        file_menu = Gtk.Menu()
        file_item = Gtk.MenuItem(label="File")
        file_item.set_submenu(file_menu)
        menubar.append(file_item)

        mi_open = Gtk.MenuItem(label="Open File")
        mi_open.connect("activate", open_file_dialog)
        file_menu.append(mi_open)

        mi_save = Gtk.MenuItem(label="Save Movie Frames")
        mi_save.connect("activate", save_movieframes_dialog)
        file_menu.append(mi_save)

        mi_exit = Gtk.MenuItem(label="Exit")
        mi_exit.connect("activate", lambda w: self.destroy())
        file_menu.append(mi_exit)

        # Options
        options_menu = Gtk.Menu()
        options_item = Gtk.MenuItem(label="Options")
        options_item.set_submenu(options_menu)
        menubar.append(options_item)

        add_opt = lambda lbl, cb: (
            lambda i: (i.connect("activate", cb), options_menu.append(i))
        )(Gtk.MenuItem(label=lbl))

        add_opt("Toggle Timeframe update/add", toggle_timeframe_update)
        add_opt("Toggle Timeframe closest/exact", toggle_timeframe_closest)
        add_opt("Toggle Axis on/off", toggle_axis_on)
        add_opt("Toggle log/lin x", toggle_log_xscale)
        add_opt("Toggle log/lin y", toggle_log_yscale)
        add_opt("Toggle Grid on/off", toggle_grid)
        add_opt("Toggle Line/Scatter", toggle_line_scatter)
        add_opt("Toggle 2D/3D", toggle_2d_3d)
        add_opt("Toggle 3D-Surface", toggle_surface)
        add_opt("Toggle Labels", toggle_labels)
        add_opt("Toggle Legend", toggle_legend)

        # Settings
        settings_menu = Gtk.Menu()
        settings_item = Gtk.MenuItem(label="Settings")
        settings_item.set_submenu(settings_menu)
        menubar.append(settings_item)

        def add_set(lbl, cb):
            mi = Gtk.MenuItem(label=lbl)
            mi.connect("activate", cb)
            settings_menu.append(mi)

        add_set("Select x-Columns", input_graph_xcolumns)
        add_set("Select y-Columns", input_graph_ycolumns)
        add_set("Select v-Columns", input_graph_vcolumns)
        add_set("Edit Limits", input_graph_limits)
        add_set("Edit Labels", input_graph_labels)
        add_set("Edit Legend", input_graph_legend)
        add_set("Graph Settings", input_graph_settings)

        # Lines
        lines_menu = Gtk.Menu()
        lines_item = Gtk.MenuItem(label="Lines")
        lines_item.set_submenu(lines_menu)
        menubar.append(lines_item)

        def add_line_set(lbl, cb):
            mi = Gtk.MenuItem(label=lbl)
            mi.connect("activate", cb)
            lines_menu.append(mi)

        add_line_set("Edit Line Colors", input_graph_linecolors)
        add_line_set("Edit Line Styles", input_graph_linestyles)
        add_line_set("Edit Line Markers", input_graph_linemarkers)
        add_line_set("Edit Line Markersizes", input_graph_linemarkersizes)
        add_line_set("Edit Line Widths", input_graph_linewidths)

        # Transformations
        transform_menu = Gtk.Menu()
        transform_item = Gtk.MenuItem(label="Transformations")
        transform_item.set_submenu(transform_menu)
        menubar.append(transform_item)

        mi_transform = Gtk.MenuItem(label="Transform Columns")
        mi_transform.connect("activate", input_graph_coltrafos)
        transform_menu.append(mi_transform)

        # Help
        help_menu = Gtk.Menu()
        help_item = Gtk.MenuItem(label="Help")
        help_item.set_submenu(help_menu)
        menubar.append(help_item)

        mi_about = Gtk.MenuItem(label="About")
        mi_about.connect("activate", about_dialog)
        help_menu.append(mi_about)

        mi_help = Gtk.MenuItem(label="Read help in tgraph.txt")
        mi_help.connect("activate", help_dialog)
        help_menu.append(mi_help)

        return menubar

    def build_time_controls(self):
            """
            Construct a row: Time [Entry] << < Play > >> Delay [Entry]
            """
            global graph_delay
            hbox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=6)

            # 'Time' Label
            lbl_time = Gtk.Label(label="Time:")
            hbox.pack_start(lbl_time, False, False, 0)

            # Time entry
            self.time_entry = Gtk.Entry()
            self.time_entry.set_width_chars(8)
            # initial value
            self.time_entry.set_text(f"{graph_time}")
            # parse and set time when user presses Enter or focus leaves
            self.time_entry.connect("activate", lambda w: set_graph_time_from_entry(w))
            self.time_entry.connect("focus-out-event", lambda w, e: (set_graph_time_from_entry(w), False))
            hbox.pack_start(self.time_entry, False, False, 0)

            # buttons: <<, <
            btn_min = Gtk.Button(label="<<")
            btn_min.connect("clicked", lambda w: (min_graph_time(), update_time_entry(self.time_entry)))
            hbox.pack_start(btn_min, False, False, 0)

            btn_dec = Gtk.Button(label="<")
            btn_dec.connect("clicked", lambda w: (dec_graph_time(), update_time_entry(self.time_entry)))
            hbox.pack_start(btn_dec, False, False, 0)

            # Play
            btn_play = Gtk.Button(label="Play")
            btn_play.connect("clicked", lambda w: start_play_graph_time())
            hbox.pack_start(btn_play, False, False, 0)

            # button: >
            btn_inc = Gtk.Button(label=">")
            btn_inc.connect("clicked", lambda w: (inc_graph_time(), update_time_entry(self.time_entry)))
            hbox.pack_start(btn_inc, False, False, 0)

            # button: >>
            btn_max = Gtk.Button(label=">>")
            btn_max.connect("clicked", lambda w: (max_graph_time(), update_time_entry(self.time_entry)))
            hbox.pack_start(btn_max, False, False, 0)

            # delay
            lbl_delay = Gtk.Label(label="  Delay:")
            hbox.pack_start(lbl_delay, False, False, 0)

            self.delay_entry = Gtk.Entry()
            self.delay_entry.set_width_chars(4)
            self.delay_entry.set_text(str(graph_delay))
            self.delay_entry.connect("activate", lambda w: set_delay_from_entry(w))
            self.delay_entry.connect("focus-out-event", lambda w, e: (set_delay_from_entry(w), False))
            hbox.pack_start(self.delay_entry, False, False, 0)

            return hbox

# -----------------------------------------------------------------------------------
#   main entry point
# -----------------------------------------------------------------------------------
def main():
    win = GTKGraphWindow()
    win.connect("destroy", Gtk.main_quit)
    Gtk.main()

if __name__ == "__main__":
    main()
