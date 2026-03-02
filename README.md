# tgraph.py — plot 1D and 2D data files

This is a GTK port of [tgraph](https://github.com/wofti/tgraph). The original script is using TK, which will not work with Wayland, and looks bad on HiDPi monitors. This working in progress GTK3 port will solve that, because GTK widgets:
- automatically works natively with either Xorg or Wayland,
- automatically follows your global HiDPi setting,
- and fits your global desktop GTK theme.

## Installation
A simple `make install` target is provided; see the accompanying
`INSTALL` file for details about prefixes and where the script will be
copied.  You can also run directly from the source tree (`python3 src/tgraph.py`) if you prefer.


## Help
More help and documentation is in the file tgraph.txt.


## Plotting examples
Example to plot columns 1 and 2 of a text data file in ygraph format:
```
tgraph.py -c 1:2 example-data/psi.X0
```

Example to plot columns 1, 2, and 4 of VTK data files:
```
tgraph.py -c 1:2:4 [ example-data/rho_*.vtk ]
```


The data file examples live in `example-data/`:

* **1D data in ygraph format:** `psi.X0`

* **2D data for gnuplot's `splot`:** `psi.XY0`

* **2D VTK (STRUCTURED_POINTS):** `rho_0862.vtk`, `rho_0942.vtk`
