import sys
import os
import importlib
import matplotlib
matplotlib.use("Agg")  # headless backend for tests

import pytest

# we will import tgraph after adjusting argv to avoid GTK initialization errors

def reload_with_argv(argv_list):
    # ensure there's at least one existing filename so module doesn't exit
    args = list(argv_list)
    if not any(os.path.exists(a) for a in args if isinstance(a, str)):
        args.append("example-data/psi.X0")
    sys.argv = ["tgraph.py"] + args
    if "tgraph" in sys.modules:
        importlib.reload(sys.modules["tgraph"])
        return sys.modules["tgraph"]
    else:
        import tgraph
        return tgraph


def test_parse_cols_and_ranges():
    tg = reload_with_argv(["-c", "1:2:3", "-x", "0:10", "-y", "5:15", "-v", "1:2"])
    # ensure arguments parsed correctly
    assert tg.xcol == 0
    assert tg.ycol == 1
    assert tg.vcol == 2
    assert tg.graph_xmin == 0.0
    assert tg.graph_xmax == 10.0
    assert tg.graph_ymin == 5.0
    assert tg.graph_ymax == 15.0
    assert tg.graph_vmin == 1.0
    assert tg.graph_vmax == 2.0


def test_parse_marker_and_stride():
    tg = reload_with_argv(["-m", "-s", "5"])
    assert tg.matplotlib.rcParams['lines.marker'] == 'o'
    assert tg.graph_stride == 5


def test_bad_col_option():
    with pytest.raises(SystemExit):
        # argparse should exit on invalid format
        reload_with_argv(["-c", "notvalid"])


class DummyAx:
    """Minimal stand-in for a matplotlib Axes.  """
    def __init__(self):
        self._xscale = 'linear'
        self._yscale = 'linear'
        self._xlabel = ''
        self._ylabel = ''
        self._title = ''

    def set_xscale(self, v):
        self._xscale = v

    def set_yscale(self, v):
        self._yscale = v

    def get_xscale(self):
        return self._xscale

    def get_yscale(self):
        return self._yscale

    def clear(self):
        # clear resets scales to default
        self._xscale = 'linear'
        self._yscale = 'linear'

    def set_xlabel(self, txt, fontsize=None):
        self._xlabel = txt

    def set_ylabel(self, txt, fontsize=None):
        self._ylabel = txt

    def set_title(self, txt):
        self._title = txt

    # methods called by axplot2d_at_time
    def plot(self, *args, **kwargs):
        pass

    def scatter(self, *args, **kwargs):
        pass


@pytest.fixture
def tg_module():
    # import with at least one file so filelist is non-empty
    tg = reload_with_argv(["example-data/psi.X0"])
    return tg


def test_toggle_labels(tg_module):
    tg = tg_module
    # give the graph some nonempty labels so toggle has visible effect
    tg.graph_labels['x-axis'] = 'X'
    tg.graph_labels['v-axis'] = 'V'
    fig = matplotlib.figure.Figure()
    ax = fig.add_subplot(111)
    tg.graph_labelsOn = 1
    tg.axplot2d_at_time(tg.filelist, ax, tg.graph_time)
    # labels should now appear
    assert ax.get_xlabel() == 'X'
    assert ax.get_ylabel() == 'V'
    # turn off
    tg.graph_labelsOn = 0
    ax.clear()
    tg.axplot2d_at_time(tg.filelist, ax, tg.graph_time)
    assert ax.get_xlabel() == '' and ax.get_ylabel() == '' and ax.get_title() == ''


def test_toggle_scales(tg_module):
    tg = tg_module
    fig = matplotlib.figure.Figure()
    ax = fig.add_subplot(111)
    # initial linear
    tg.graph_xscale = 'linear'
    tg.graph_yscale = 'linear'
    # flip x
    tg.toggle_log_xscale(None)
    assert tg.graph_xscale == 'log'
    tg.axplot2d_at_time(tg.filelist, ax, tg.graph_time)
    assert ax.get_xscale() == 'log'
    # flip y
    tg.toggle_log_yscale(None)
    assert tg.graph_yscale == 'log'
    ax.clear()
    tg.axplot2d_at_time(tg.filelist, ax, tg.graph_time)
    assert ax.get_yscale() == 'log'
