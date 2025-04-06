#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import numpy as np

# -----------------------
#  0) Python 版本检查
# -----------------------
if sys.version_info[0] < 3:
    print("Please use Python 3 or above.")
    sys.exit(1)

# -----------------------
#  1) 设置 matplotlib 后端为 GTK3Agg
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
#  2) 导入 PyGObject
# -----------------------
import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, GLib

# -----------------------
#  3) 导入 tdata.py
# -----------------------
import tdata

# -----------------------------------------------------------------------------------
#          以下部分为原 tgraph.py 中的数据处理、命令行解析、全局变量等
# -----------------------------------------------------------------------------------

tgraph_version = "GTK-ported-1.0"

print("GTK tgraph version:", tgraph_version)

# 默认列
xcol = 0
ycol = 1
zcol = 2
vcol = 1
print("Default cols:", xcol+1, ycol+1, ":", vcol+1)

graph_stride = 1

# 解析命令行，构造 filelist
filelist = tdata.tFileList()
argvs = sys.argv[1:]
print("Trying to open files:")

gotC = 0
gotx = 0
goty = 0
gotv = 0
gots = 0
gott = 0
timelabel_str = "time"
got_xrange = 0
got_yrange = 0
got_vrange = 0
openCBrack = 0
inCBrack = 0
openSBrack = 0
inSBrack = 0
endSBrack = 0

for argv in argvs:
    if argv.startswith("-c"):
        gotC = 1
        continue
    elif argv.startswith("-x"):
        gotx = 1
        continue
    elif argv.startswith("-y"):
        goty = 1
        continue
    elif argv.startswith("-v"):
        gotv = 1
        continue
    elif argv.startswith("-s"):
        gots = 1
        continue
    elif argv.startswith("-t"):
        gott = 1
        continue
    elif argv.startswith("-m"):
        matplotlib.rcParams['lines.marker'] = 'o'
        continue

    # 处理 -c 后跟的 “xcol:vcol” 或 “xcol:ycol:vcol”
    if gotC == 1:
        cols = argv.split(":")
        xcol = int(cols[0]) - 1
        if len(cols) == 2:
            vcol = int(cols[1]) - 1
        elif len(cols) == 3:
            ycol = int(cols[1]) - 1
            vcol = int(cols[2]) - 1
        print("cols:", xcol+1, ycol+1, ":", vcol+1)
        gotC = 0
        continue

    # 处理 -x/-y/-v 后跟的 “min:max”
    if gotx == 1 or goty == 1 or gotv == 1:
        parts = argv.split(":")
        if len(parts) == 2:
            mn = float(parts[0])
            mx = float(parts[1])
            if gotx:
                graph_xmin, graph_xmax = mn, mx
                got_xrange = 1
            elif goty:
                graph_ymin, graph_ymax = mn, mx
                got_yrange = 1
            else:
                graph_vmin, graph_vmax = mn, mx
                got_vrange = 1
        gotx = goty = gotv = 0
        continue

    # 处理 -s
    if gots == 1:
        graph_stride = int(argv)
        gots = 0
        continue

    # 处理 -t
    if gott == 1:
        timelabel_str = str(argv).lower()
        gott = 0
        continue

    # 处理 { } [ ] 等文件合并语法
    if argv == "{":
        openCBrack = 1
        inCBrack = 0
        continue
    elif argv == "}":
        openCBrack = 0
        inCBrack = 0
        continue
    elif argv == "[":
        openSBrack = 1
        inSBrack = 0
        endSBrack = 0
        continue
    elif argv == "]":
        openSBrack = 0
        inSBrack = 0
        endSBrack = 1
    else:
        # 视作文件名
        filelist.add(argv, timelabel_str)
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

# 若未传入任何文件，打印帮助退出
if len(filelist.file) == 0:
    print("No files given on command line.\nExample usage:")
    print("  ./gtk_tgraph_full.py file1 file2 -c 1:2")
    sys.exit(1)

# 为所有文件设置列并初始化
for i in range(len(filelist.file)):
    filelist.file[i].data.set_cols(xcol=xcol, ycol=ycol, zcol=2, vcol=vcol)

graph_time = filelist.mintime()
graph_timelist = filelist.get_timelist()
graph_timeindex = tdata.geti_from_t(graph_timelist, graph_time)

# 如果命令行未指定 x/y/v 范围，则从数据中自动获取
if not hasattr(sys.modules[__name__], 'graph_xmin'):
    graph_xmin = tdata.inf_to_1e300(filelist.minx())
    graph_xmax = tdata.inf_to_1e300(filelist.maxx())
if not hasattr(sys.modules[__name__], 'graph_ymin'):
    graph_ymin = tdata.inf_to_1e300(filelist.miny())
    graph_ymax = tdata.inf_to_1e300(filelist.maxy())
if not hasattr(sys.modules[__name__], 'graph_vmin'):
    graph_vmin = tdata.inf_to_1e300(filelist.minv())
    graph_vmax = tdata.inf_to_1e300(filelist.maxv())

# -----------------------
#   一些全局状态
# -----------------------
graph_3dOn = 0
graph_axis_on = 1
graph_plot_surface = 0
graph_plot_scatter = 0
graph_clear_on_replot = 1
graph_plot_closest_t = 1
graph_plot_grid = 1
# colormap
graph_colormap_str = "coolwarm"  # 比如 "jet" 或 "coolwarm"
graph_colormap = getattr(cm, graph_colormap_str)

# -----------------------
#   字典存储相关设置
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

# 初始化时，为每个文件分配一个默认颜色
def set_graph_globals_for_file_i(filelist, i):
    # 类似 tgraph.py 对 color_cycle 的使用
    if matplotlib.__version__ < '1.5.1':
        # 老版本 color_cycle
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
    # graph_legend 中存储该文件的名字
    graph_legend[f"#{i}"] = filelist.file[i].name

for i in range(len(filelist.file)):
    set_graph_globals_for_file_i(filelist, i)

print("(tmin, tmax) =", filelist.mintime(), filelist.maxtime())
print("(xmin, xmax) =", graph_xmin, graph_xmax)
print("(ymin, ymax) =", graph_ymin, graph_ymax)
print("(vmin, vmax) =", graph_vmin, graph_vmax)
print("stride =", graph_stride)

# -----------------------------------------------------------------------------------
#   2D/3D 绘图函数
# -----------------------------------------------------------------------------------
def axplot2d_at_time(filelist, ax, t):
    if graph_clear_on_replot:
        ax.clear()
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
    ax.set_xlabel(graph_labels['x-axis'], fontsize=graph_labels['fontsize'])
    ax.set_ylabel(graph_labels['v-axis'], fontsize=graph_labels['fontsize'])

    # 标题与时间
    title = graph_labels['title']
    tf = graph_labels['timeformat']
    if tf:
        tstr = tf % t
        # 标题右对齐显示时间
        title += "    " + tstr
    ax.set_title(title)

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
    与原 tgraph.py 类似，以 x, y 为面上坐标，v 为 z。
    """
    from mpl_toolkits.mplot3d import axes3d
    if graph_clear_on_replot:
        ax.clear()

    for i in range(len(filelist.file)):
        f = filelist.file[i]
        # 读取 blocks 以判断数据换行
        blocks = f.data.getblocks(t)
        # 变形
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

    ax.set_xlabel(graph_labels['x-axis'], fontsize=graph_labels['fontsize'])
    ax.set_ylabel(graph_labels['y-axis'], fontsize=graph_labels['fontsize'])
    ax.set_zlabel(graph_labels['v-axis'], fontsize=graph_labels['fontsize'])

    title = graph_labels['title']
    tf = graph_labels['timeformat']
    if tf:
        tstr = tf % t
        title += "    " + tstr
    ax.set_title(title)

    # 注意 3D surface 不支持 legend
    if graph_legendOn and not graph_plot_surface:
        ax.legend()

    if not graph_axis_on:
        ax.set_axis_off()

graph_delay = 1  # 新增，用于控制播放的延时(ms)

# 新增：对时间值进行增/减、首/末帧、播放相关函数
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

# 播放时要有个全局变量存放计时器ID，以便停止
playback_id = None

def play_graph_time():
    """
    每次调用时，把时间帧 +1 并 replot。如果还没到最后，则继续下一次定时。
    如果到最后一帧，则停止。
    """
    global graph_timeindex, playback_id, graph_delay

    inc_graph_time()  # 已经会 replot
    if graph_timeindex < len(graph_timelist) - 1:
        # 继续播放
        playback_id = GLib.timeout_add(int(graph_delay), play_graph_time)
    else:
        playback_id = None  # 播放结束

def start_play_graph_time():
    """
    从当前时间帧开始播放。如果已经在播，则先停止。
    """
    global playback_id
    if playback_id:
        # 正在播放 -> 停止
        GLib.source_remove(playback_id)
        playback_id = None
    else:
        # 没在播放 -> 启动
        if graph_timeindex >= len(graph_timelist) - 1:
            # 若在末帧，则从头开始
            min_graph_time()
        playback_id = GLib.timeout_add(int(graph_delay), play_graph_time)

# 当用户在 Entry 中敲回车或离开时，设置 graph_time
def set_graph_time_from_entry(entry):
    global graph_time, graph_timeindex
    try:
        val = float(entry.get_text().strip())
        # 找到与 val 最接近的时间索引
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

# 更新时间 Entry 的显示
def update_time_entry(entry):
    global graph_time
    entry.set_text(f"{graph_time}")

# -----------------------------------------------------------------------------------
#   一些工具函数（原 tgraph 中的 replot, toggle 等）
# -----------------------------------------------------------------------------------
def setup_axes(fig, is3d, old_ax=None):
    """
    在 fig 上重新创建 Axes（2D 或 3D）。删除旧轴。
    """
    if old_ax is not None:
        fig.delaxes(old_ax)
    if not is3d:
        return fig.add_subplot(111)
    else:
        from mpl_toolkits.mplot3d import Axes3D
        return fig.add_subplot(111, projection='3d')

# 全局引用，方便菜单回调
global_fig = None
global_ax = None
global_canvas = None

def replot():
    global global_fig, global_ax, global_canvas
    global graph_time, graph_3dOn

    if graph_3dOn == 0:
        axplot2d_at_time(filelist, global_ax, graph_time)
    else:
        axplot3d_at_time(filelist, global_ax, graph_time)

    global_canvas.draw()

# 菜单回调：toggles
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
    global global_ax
    if global_ax.get_xscale() == 'linear':
        global_ax.set_xscale('log')
    else:
        global_ax.set_xscale('linear')
    replot()

def toggle_log_yscale(_menuitem):
    global global_ax
    if global_ax.get_yscale() == 'linear':
        global_ax.set_yscale('log')
    else:
        global_ax.set_yscale('linear')
    replot()

def toggle_grid(_menuitem):
    global graph_plot_grid
    graph_plot_grid = 0 if graph_plot_grid else 1
    replot()

def toggle_line_scatter(_menuitem):
    global graph_plot_surface, graph_plot_scatter
    if graph_plot_surface:
        # 若原先是 surface，先关掉
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
    # 若 labels 关闭，可把 label 全清空；也可以仅不显示
    # 这里简化仅视作“显示/隐藏”
    replot()

def toggle_legend(_menuitem):
    global graph_legendOn
    graph_legendOn = 0 if graph_legendOn else 1
    replot()

# -----------------------------------------------------------------------------------
#   “对话框” 用于编辑若干 key->value 的界面
# -----------------------------------------------------------------------------------
class GTKDialogKeyValue(Gtk.Dialog):
    """
    类似原先 WTdialog，用来修改一个 dict 中的字符串 / 数值。
    """
    def __init__(self, title, keyvalues: dict, parent=None):
        super().__init__(title=title, transient_for=parent, modal=True)
        self.set_default_size(400, 300)

        self.keyvalues = keyvalues
        self.entries = {}

        vbox = self.get_content_area()
        grid = Gtk.Grid(column_spacing=8, row_spacing=6, margin=10)
        vbox.add(grid)

        row = 0
        # 为 dict 中每个键创建一行 (Label + Entry)
        for key in sorted(keyvalues.keys()):
            label = Gtk.Label(label=key, xalign=1.0)
            entry = Gtk.Entry()
            entry.set_text(str(keyvalues[key]))
            self.entries[key] = entry
            grid.attach(label, 0, row, 1, 1)
            grid.attach(entry, 1, row, 1, 1)
            row += 1

        # 底部按钮
        self.ok_button = Gtk.Button(label="Apply")
        self.ok_button.connect("clicked", self.on_apply)
        grid.attach(self.ok_button, 0, row, 2, 1)

        self.show_all()

    def on_apply(self, button):
        # 从对话框的 entries 中获取值
        for k, e in self.entries.items():
            text = e.get_text().strip()
            self.keyvalues[k] = text
        self.destroy()

# -----------------------------------------------------------------------------------
#   菜单回调：各种输入对话框 -> 改变全局设置 -> replot
# -----------------------------------------------------------------------------------
def input_graph_limits(_menuitem):
    global graph_limits
    dlg = GTKDialogKeyValue("Edit Limits", graph_limits)
    dlg.run()  # 阻塞等待
    dlg.destroy()

    # 将编辑结果写回全局
    # 转换为 float
    graph_limits['xmin'] = float(graph_limits['xmin'])
    graph_limits['xmax'] = float(graph_limits['xmax'])
    graph_limits['ymin'] = float(graph_limits['ymin'])
    graph_limits['ymax'] = float(graph_limits['ymax'])
    graph_limits['vmin'] = float(graph_limits['vmin'])
    graph_limits['vmax'] = float(graph_limits['vmax'])

    replot()

def input_graph_labels(_menuitem):
    global graph_labels, graph_labelsOn
    dlg = GTKDialogKeyValue("Edit Labels", graph_labels)
    dlg.run()
    dlg.destroy()

    # 转换一些值
    matplotlib.rcParams['font.size'] = float(graph_labels['fontsize'])
    # 打开显示
    graph_labelsOn = 1
    replot()

def input_graph_legend(_menuitem):
    global graph_legend, graph_legendOn, filelist
    dlg = GTKDialogKeyValue("Edit Legend", graph_legend)
    dlg.run()
    dlg.destroy()

    # 处理 #0, #1, #2... 对应 file 名
    for i in range(len(filelist.file)):
        k = f"#{i}"
        if k in graph_legend:
            filelist.file[i].name = graph_legend[k]
    # 解析 loc
    s = str(graph_legend['loc']).strip()
    if s.isdigit():
        graph_legend['loc'] = int(s)
    else:
        # 尝试解析成 "x, y"
        if "," in s:
            parts = s.replace("(", "").replace(")", "").split(",")
            graph_legend['loc'] = (float(parts[0]), float(parts[1]))
    graph_legend['ncol'] = int(graph_legend['ncol'])

    # 打开显示
    graph_legendOn = 1
    replot()

def input_graph_settings(_menuitem):
    global graph_settings, graph_stride, graph_colormap, graph_colormap_str
    dlg = GTKDialogKeyValue("Graph Settings", graph_settings)
    dlg.run()
    dlg.destroy()

    # 应用更改
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
    用来修改每个 file 的 xcol
    """
    # 准备一个 dict
    d = {}
    for i in range(len(filelist.file)):
        d[f"#{i}"] = str(filelist.file[i].data.get_xcol0() + 1)
    dlg = GTKDialogKeyValue("Select x-Columns", d)
    dlg.run()
    dlg.destroy()
    # 应用
    for i in range(len(filelist.file)):
        newcol = int(d[f"#{i}"]) - 1
        filelist.file[i].data.set_xcols(newcol)

    # 重新计算 x min max
    global graph_xmin, graph_xmax
    graph_xmin = tdata.inf_to_1e300(filelist.minx())
    graph_xmax = tdata.inf_to_1e300(filelist.maxx())
    graph_limits['xmin'] = graph_xmin
    graph_limits['xmax'] = graph_xmax
    replot()

def input_graph_ycolumns(_menuitem):
    d = {}
    for i in range(len(filelist.file)):
        d[f"#{i}"] = str(filelist.file[i].data.get_ycol0() + 1)
    dlg = GTKDialogKeyValue("Select y-Columns", d)
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
    d = {}
    for i in range(len(filelist.file)):
        d[f"#{i}"] = str(filelist.file[i].data.get_vcol0() + 1)
    dlg = GTKDialogKeyValue("Select v-Columns", d)
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
    d = {}
    for i in range(len(filelist.file)):
        d[f"#{i}"] = graph_linecolors[f"#{i}"]
    dlg = GTKDialogKeyValue("Edit Line Colors", d)
    dlg.run()
    dlg.destroy()
    for i in range(len(filelist.file)):
        graph_linecolors[f"#{i}"] = d[f"#{i}"]
    replot()

def input_graph_linestyles(_menuitem):
    d = {}
    for i in range(len(filelist.file)):
        d[f"#{i}"] = graph_linestyles[f"#{i}"]
    dlg = GTKDialogKeyValue("Edit Line Styles", d)
    dlg.run()
    dlg.destroy()
    for i in range(len(filelist.file)):
        graph_linestyles[f"#{i}"] = d[f"#{i}"]
    replot()

def input_graph_linemarkers(_menuitem):
    d = {}
    for i in range(len(filelist.file)):
        d[f"#{i}"] = graph_linemarkers[f"#{i}"]
    dlg = GTKDialogKeyValue("Edit Line Markers", d)
    dlg.run()
    dlg.destroy()
    for i in range(len(filelist.file)):
        graph_linemarkers[f"#{i}"] = d[f"#{i}"]
    replot()

def input_graph_linemarkersizes(_menuitem):
    d = {}
    for i in range(len(filelist.file)):
        d[f"#{i}"] = str(graph_linemarkersizes[f"#{i}"])
    dlg = GTKDialogKeyValue("Edit Line Markersizes", d)
    dlg.run()
    dlg.destroy()
    for i in range(len(filelist.file)):
        graph_linemarkersizes[f"#{i}"] = float(d[f"#{i}"])
    replot()

def input_graph_linewidths(_menuitem):
    d = {}
    for i in range(len(filelist.file)):
        d[f"#{i}"] = str(graph_linewidths[f"#{i}"])
    dlg = GTKDialogKeyValue("Edit Line Widths", d)
    dlg.run()
    dlg.destroy()
    for i in range(len(filelist.file)):
        val = d[f"#{i}"].strip()
        graph_linewidths[f"#{i}"] = val
    replot()

def input_graph_coltrafos(_menuitem):
    """
    与原脚本类似，可对列进行表达式变换。此处示例仅保留 UI，不详细演示 transform_col。
    """
    d = {}
    for i in range(len(filelist.file)):
        d[f"#{i}"] = graph_coltrafos[f"#{i}"]
    dlg = GTKDialogKeyValue("Transform Columns", d)
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
    md = Gtk.MessageDialog(
        text=f"tgraph {tgraph_version}\n\n"
                       "Copyright (C) 2015 Wolfgang Tichy.\n"
                       "Ported to GTK by Yingjie Wang.\n",
        buttons=Gtk.ButtonsType.OK
    )
    md.run()
    md.destroy()

def help_dialog(_menuitem):
    # 简单读一下 tgraph.txt
    tgraph_txt_file = os.path.join(os.path.dirname(tdata.__file__), "tgraph.txt")
    helptext = "tgraph.txt not found!"
    if os.path.exists(tgraph_txt_file):
        with open(tgraph_txt_file, "r") as f:
            helptext = f.read()

    # 用 TextView 显示
    dialog = Gtk.Dialog(title="tgraph help", modal=True)
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
# 新增：Open File 对话框函数
# -----------------------------
def open_file_dialog(_menuitem):
    """
    允许用户通过 FileChooserDialog 选择多个新文件，并添加到 filelist。
    """
    dialog = Gtk.FileChooserDialog(
        title="Open File(s)",
        parent=None,
        action=Gtk.FileChooserAction.OPEN
    )
    dialog.add_buttons(
        Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
        Gtk.STOCK_OPEN,   Gtk.ResponseType.OK
    )
    dialog.set_select_multiple(True)  # 允许多选

    # 可以根据需要设置过滤器
    # file_filter = Gtk.FileFilter()
    # file_filter.set_name("All files")
    # file_filter.add_pattern("*")
    # dialog.add_filter(file_filter)

    response = dialog.run()
    if response == Gtk.ResponseType.OK:
        # 获取所有选定文件
        filenames = dialog.get_filenames()
        print("User selected:", filenames)

        # 对每个文件，交给 filelist 处理
        from math import isinf
        for fname in filenames:
            if not os.path.isfile(fname):
                continue
            # 假设使用默认 timelabel_str
            timelabel_str = "time"
            filelist.add(fname, timelabel_str)
            print("Added:", filelist.file[-1].filename)
            # 设置默认列
            # xcol, ycol, zcol, vcol 等可从全局拿
            filelist.file[-1].data.set_cols(xcol=0, ycol=1, zcol=2, vcol=1)

        # 更新 graph_timelist, graph_time 等
        global graph_timelist, graph_time, graph_timeindex
        graph_timelist = filelist.get_timelist()
        # 若只有初始空的话，现在可以拿到新的 time 范围
        # graph_time = filelist.mintime()
        # graph_timeindex = tdata.geti_from_t(graph_timelist, graph_time)

        # 更新 x/y/v 范围
        global graph_xmin, graph_xmax, graph_ymin, graph_ymax, graph_vmin, graph_vmax
        from tdata import inf_to_1e300
        graph_xmin = inf_to_1e300(filelist.minx())
        graph_xmax = inf_to_1e300(filelist.maxx())
        graph_ymin = inf_to_1e300(filelist.miny())
        graph_ymax = inf_to_1e300(filelist.maxy())
        graph_vmin = inf_to_1e300(filelist.minv())
        graph_vmax = inf_to_1e300(filelist.maxv())

        # 同步到 graph_limits 字典
        graph_limits['xmin'] = graph_xmin
        graph_limits['xmax'] = graph_xmax
        graph_limits['ymin'] = graph_ymin
        graph_limits['ymax'] = graph_ymax
        graph_limits['vmin'] = graph_vmin
        graph_limits['vmax'] = graph_vmax

        # 给新文件分配线条颜色/名字
        for i in range(len(filelist.file)):
            # 如果之前没初始化过，就 set_graph_globals_for_file_i
            # 可以简单点，对最后几个新的 file 做 set_graph_globals_for_file_i
            pass

        # 重绘
        replot()

    dialog.destroy()

# -----------------------------
# 新增：保存序列帧函数
# -----------------------------
def save_movieframes_dialog(_menuitem):
    """
    弹出一个 “Save As” 对话框，指定输出文件名。例如 “frame.png”。
    然后把所有时间帧以 "basename_{index}.ext" 的形式保存下来。
    """
    dialog = Gtk.FileChooserDialog(
        title="Save Movie Frames",
        action=Gtk.FileChooserAction.SAVE
    )
    dialog.add_buttons(
        Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
        Gtk.STOCK_SAVE,   Gtk.ResponseType.OK
    )
    dialog.set_do_overwrite_confirmation(True)

    response = dialog.run()
    if response == Gtk.ResponseType.OK:
        outname = dialog.get_filename()
        # outname 可能是 /home/user/frame.png
        print("Saving frames to base name:", outname)

        # 处理扩展名, 基础名
        base, ext = os.path.splitext(outname)
        if not ext:
            # 若用户没写扩展名，给个默认
            ext = ".png"

        # 全部时间帧
        global graph_timeindex, graph_time, graph_timelist
        i1 = 0
        i2 = len(graph_timelist)
        # 打印需要的零填充位数
        digits = len(str(i2-1))
        fmt_str = f"%0{digits}d"

        # 存储当前 index
        old_index = graph_timeindex
        old_time = graph_time

        # 循环保存
        for idx in range(i2):
            graph_timeindex = idx
            graph_time = graph_timelist[idx]
            replot()  # 让画面更新到对应时间

            # 构造文件名
            frame_idx = fmt_str % idx
            filename = f"{base}_{frame_idx}{ext}"
            global_canvas.print_figure(filename)
            if idx == 0:
                first_file = filename
            if idx == i2-1:
                last_file = filename

        # 恢复之前的时间
        graph_timeindex = old_index
        graph_time = old_time
        replot()

        # 给用户提示
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
#   GTK 主窗口
# -----------------------------------------------------------------------------------
class GTKGraphWindow(Gtk.Window):
    def __init__(self):
        super().__init__(title="GTK tgraph")
        self.set_default_size(900, 700)

        # 创建一个垂直布局
        vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=0)
        self.add(vbox)

        # 1) 菜单栏
        menubar = self.build_menubar()
        vbox.pack_start(menubar, False, False, 0)

        # 2) 最上方增加一行“时间播放控件”
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

        # 初次绘制
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
            构造一行：Time [Entry] << < Play > >> Delay [Entry]
            """
            global graph_delay
            hbox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=6)

            # “Time” Label
            lbl_time = Gtk.Label(label="Time:")
            hbox.pack_start(lbl_time, False, False, 0)

            # 时间 Entry
            self.time_entry = Gtk.Entry()
            self.time_entry.set_width_chars(8)
            # 初始显示
            self.time_entry.set_text(f"{graph_time}")
            # 当用户敲回车或失焦时，解析并设置 time
            self.time_entry.connect("activate", lambda w: set_graph_time_from_entry(w))
            self.time_entry.connect("focus-out-event", lambda w, e: (set_graph_time_from_entry(w), False))
            hbox.pack_start(self.time_entry, False, False, 0)

            # 按钮: <<, <
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

            # 按钮: >
            btn_inc = Gtk.Button(label=">")
            btn_inc.connect("clicked", lambda w: (inc_graph_time(), update_time_entry(self.time_entry)))
            hbox.pack_start(btn_inc, False, False, 0)

            # 按钮: >>
            btn_max = Gtk.Button(label=">>")
            btn_max.connect("clicked", lambda w: (max_graph_time(), update_time_entry(self.time_entry)))
            hbox.pack_start(btn_max, False, False, 0)

            # 间隔
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
#   主入口
# -----------------------------------------------------------------------------------
def main():
    win = GTKGraphWindow()
    win.connect("destroy", Gtk.main_quit)
    Gtk.main()

if __name__ == "__main__":
    main()
