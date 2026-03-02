import sys
import os
import importlib
import matplotlib
matplotlib.use("Agg")  # headless backend for tests

# ensure the top-level project directory is on the Python path so we can
# import the `tgraph` module (pytest may change cwd to tests/ during collection)
_here = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
if _here not in sys.path:
    sys.path.insert(0, _here)

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


def test_dialog_transient_and_hint(tg_module):
    # the GTKDialogKeyValue should honour the parent and set a dialog hint
    from gi.repository import Gtk
    Gtk.init_check()
    try:
        parent = Gtk.Window()
    except RuntimeError:
        pytest.skip("GTK not usable in test environment")
    dlg = tg_module.GTKDialogKeyValue("t", {"k": "v"}, parent=parent)
    assert dlg.get_transient_for() is parent
    assert dlg.get_type_hint() == Gtk.WindowTypeHint.DIALOG
    dlg.destroy()
    parent.destroy()


def test_main_window_global(tg_module):
    from gi.repository import Gtk
    Gtk.init_check()
    try:
        win = tg_module.GTKGraphWindow()
    except RuntimeError:
        pytest.skip("GTK not usable in test environment")
    assert tg_module.main_window is win
    win.destroy()


def test_transform_parent_fallback(tg_module):
    """Ensure the transform dialog uses the stored main window when called.

    This catches the earlier NameError bug and also exercises the global
    parent lookup logic.
    """
    from gi.repository import Gtk
    Gtk.init_check()
    try:
        win = tg_module.GTKGraphWindow()
    except RuntimeError:
        pytest.skip("GTK not usable in test environment")
    # monkeypatch GTKDialogKeyValue to record what parent was passed
    called = {}
    class Recorder(tg_module.GTKDialogKeyValue):
        def __init__(self, title, keyvalues, parent=None):
            called['parent'] = parent
            super().__init__(title, keyvalues, parent=parent)
    monkeypatch = pytest.MonkeyPatch()
    monkeypatch.setattr(tg_module, 'GTKDialogKeyValue', Recorder)
    tg_module.input_graph_coltrafos(None)
    assert called.get('parent') is tg_module.main_window
    monkeypatch.undo()
    win.destroy()


def test_get_parent_logic(tg_module):
    # _get_parent should prefer global main_window when set, and handle None
    tg_module.main_window = None
    class Dummy:
        def __init__(self, w):
            self._w = w
        def get_toplevel(self):
            return self._w
    a = object()
    assert tg_module._get_parent(Dummy(a)) is a
    assert tg_module._get_parent(None) is None
    w = object()
    tg_module.main_window = w
    assert tg_module._get_parent(Dummy(object())) is w
    assert tg_module._get_parent(None) is w
