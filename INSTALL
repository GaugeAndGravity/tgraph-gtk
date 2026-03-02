DEPENDENCIES

The GTK frontend requires a recent Python 3 installation along with
these packages:

  * numpy
  * matplotlib
  * pygobject (PyGObject / GTK 3)

Most distributions provide the above; for example on Debian/Ubuntu:

    sudo apt install python3-matplotlib python3-gi

Older versions of tgraph also used Tk and tkinter, but the GTK port does
not require them.

INSTALLATION

A simple ``make install`` is provided to mimic the behaviour of typical
GNU/Linux packages.  By default it installs into ``~/.local`` for a
regular user and ``/usr/local`` when run as root.  The resulting layout
is:

    $(prefix)/bin/tgraph          # executable script
    $(prefix)/bin/tdata.py        # helper module
    $(prefix)/bin/tgraph.txt      # descriptive text file

The installer merely copies the two Python files into ``$(bindir)``
and places the help text under ``$(prefix)/share/doc/tgraph``.  ``tgraph``
is made executable and will locate ``tdata.py`` next to itself at
runtime, so nothing needs to be added to ``$PYTHONPATH``.

Usage after installation is as simple as running ``tgraph`` from your
PATH.  Uninstalling is just as easy (``make uninstall``) – the help
text is removed along with the binary and module.
If you prefer to run from the source tree, there is no need to install
anything; make sure ``src`` is on ``PYTHONPATH`` or call the script
directly (``python3 src/tgraph.py``).

For development the tests exercise the module directly and already
adjust ``sys.path`` so that ``import tgraph`` works from the repository
root.
