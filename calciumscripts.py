import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors
import matplotlib.patches as patches
import plotly.graph_objects as go
import pandas as pd
from scipy.ndimage import gaussian_filter

def calciumplot(
    stim=False,
    size=(8, 6),
    figdpi=800,
    plotmin=0,
    plotmax=3000,
    cellrange=[],
    offsetnbr=250,
    figname="",
    F=None,
    Fneu=None,
    iscell=None,
):
    """ Plot selected (cellrange) calcium traces, in a specified range (plotmin, plotmax). """

    fig, ax = plt.subplots(figsize=size)
    cm = plt.get_cmap("viridis")

    # Extract plotmin/max range, if iscell == 1.
    Flist = []
    Fneulist = []
    traceindex = []  # keep track of trace index.
    for trace in cellrange:
        if iscell[trace][0] == 1:
            Flist.append(F[trace][plotmin:plotmax])
            Fneulist.append(Fneu[trace][plotmin:plotmax])
            traceindex.append(trace)

    # plot with Y offset.
    offset = 0  # y offset per plot.

    # normalize colormap.
    colors = [
        cm(x / (len(Flist) - 0.6)) for x in list(range(len(Flist)))
    ]
    for index in range(len(Flist)):
        # Subtract neuropil*0.7, normalize to lower fifteenth of full trace.
        y = (Flist[index] - (0.7 * Fneulist[index])) / lowerperc(Flist[index])
        y *= 100
        x = np.arange(0, len(Flist[index]))
        x = x * (1 / 30.9)  # frames to seconds.
        plt.plot(
            x,
            y + offset,
            label="Trace #" + str(traceindex[index]),
            color=colors[index],
            lw=1.35,
            alpha=0.85,
        )
        offset += offsetnbr
        # plt.plot(z, label="Normalized to btm 10%")

    # Add rectangle stimulus marker. TODO: add stimulus marker from file.
    if stim:
        for i in x[67::309]:  # frames(index)
            stim = patches.Rectangle(
                (i, -75), 2, offset + 125, fill=True, alpha=0.15, color="r", linewidth=0
            )
            # rest = patches.Rectangle((i+3, -75), 8, 800, fill=True, alpha=0.3, color="b")
            ax.add_patch(stim)
            # ax.add_patch(rest)

    # Add L-shaped scalebar.
    scalebar(fontSize=12, scaleXsize=5, scaleYsize=100, lineWidth=3)

    plt.savefig(figname + ".png", dpi=figdpi)
    plt.show()


def calciumplotaxon(
    stim=False,
    size=(8, 6),
    figdpi=800,
    plotmin=0,
    plotmax=3000,
    cell_axon_pair=[],
    offsetnbr=250,
    figname="",
    F=None,
    F_axons=None,
    Fneu_axons=None,
    Fneu=None,
    iscell=None,
):
    """ Plot selected (cellrange) calcium traces from soma and axon pairs, in a specified range (plotmin, plotmax). """

    fig, ax = plt.subplots(figsize=size)
    cm = plt.get_cmap("viridis")

    # Extract plotmin/max range, if iscell == 1.
    Flist = []
    F_axonslist = []
    F_axonsigma = []
    Fneulist = []
    Fneu_axonslist = []
    traceindex = []  # keep track of trace index.

    for pair in cell_axon_pair:
        Flist.append(F[pair[0]][plotmin:plotmax])
        Fneulist.append(Fneu[pair[0]][plotmin:plotmax])
        F_axonslist.append(F_axons[pair[1]][plotmin:plotmax])
        Fneu_axonslist.append(Fneu_axons[pair[1]][plotmin:plotmax])
        F_axonsigma.append(pair[2])
        traceindex.append(pair[:2])

    # normalize colormap.
    colors = [
        cm(x / (len(Flist) - 0.6)) for x in list(range(len(Flist)))
    ]

    # Normalize and plot with Y offset.
    offset = 0  # y offset per plot.
    for index, F in enumerate(Flist):
        # Subtract neuropil*0.7, normalize to lower fifteenth of full trace.
        y = (F - (0.7 * Fneulist[index])) / lowerperc(F)
        y *= 100
        z = (F_axonslist[index] - (0.7 * Fneu_axonslist[index])) / lowerperc(F_axonslist[index])
        z *= 100
        x = np.arange(0, len(F_axonslist[index]))
        x = x * (1 / 30.9)  # frames to seconds.
        plt.plot(
            x,
            y + offset,
            label="Trace #" + str(traceindex[index]),
            color="black",
            lw=1.4,
            alpha=0.45,
        )
        plt.plot(
            x,
            gaussian_filter(z, sigma=F_axonsigma[index]) + offset,
            label="Normalized to btm 10%",
            color=colors[index],
            alpha=0.65,
            lw=1.3,
        )
        offset += offsetnbr

    # Add rectangle stimulus marker.
    if stim:
        for i in x[67::309]:  # frames(index)
            stim = patches.Rectangle(
                (i, -75), 2, offset + 125, fill=True, alpha=0.15, color="r", linewidth=0
            )
            # rest = patches.Rectangle((i+3, -75), 8, 800, fill=True, alpha=0.3, color="b")
            ax.add_patch(stim)
            # ax.add_patch(rest)

    # Add L-shaped scalebar.
    scalebar(fontSize=12, scaleXsize=5, scaleYsize=100, lineWidth=3)

    plt.savefig(figname + ".png", dpi=figdpi)
    plt.show()

def roistotext(stats=None, ops=None, iscell=None, outputdir=None):
    """ Convert ROIs from stat to text image readable by ImageJ. """
    neumask = np.zeros((ops['Ly'], ops['Lx']))
    f, axarr = plt.subplots(1,2)

    #for i in range(len(stat)):
    #    neuy, neux = np.unravel_index(stat[i]['neuropil_mask'], (ops['Ly'], ops['Lx']))
    #    neumask[neuy, neux] = 1

    im = np.zeros((ops['Ly'], ops['Lx']), dtype=int)
    for i, stat in enumerate(stats):
        if iscell[i][0] == 1:
            ypix = stat['ypix'][~stat['overlap']]
            xpix = stat['xpix'][~stat['overlap']]
            im[ypix,xpix] = i+1

    with open(outputdir + "roi_textimage.txt", "w") as outfile:
        np.savetxt(outfile, im, fmt="%i")

    axarr[0].imshow(neumask)
    axarr[1].imshow(im)

    #plt.savefig("foo.png")
    plt.show()


def calciumplotlive(
    cellrange=[], plotmin=0, plotmax=28000, offsetnbr=250, F=[], Fneu=[], iscell=[]
):
    """Plot selected (cellrange) calcium traces in a specified time-range (plotmin, plotmax)
    in live, browsable plotly.
    """

    # Extract plotmin/max range, if iscell == 1.
    Flist = []
    Fneulist = []
    traceindex = []  # keep track of trace index.
    for trace in cellrange:
        if iscell[trace][0] == 1:
            Flist.append(F[trace][plotmin:plotmax])
            Fneulist.append(Fneu[trace][plotmin:plotmax])
            traceindex.append(trace)

    offset = 0  # y offset per plot.
    fig = go.Figure()
    for index in range(len(Flist)):
        # Subtract neuropil*0.7, normalize to lower tenth of full trace.
        y = (Flist[index] - (0.7 * Fneulist[index])) / lowerperc(Flist[index])
        y = (y - y.mean()) * 100
        x = np.arange(0, len(Flist[index]))
        x = x * (1 / 30.9)  # frames to seconds.
        fig.add_trace(
            go.Scatter(
                x=x, y=y + offset, mode="lines", name="Trace " + str(traceindex[index])
            )
        )
        offset += offsetnbr

    fig.update_layout(autosize=False, width=1000, height=1000)
    fig.show()


def somaoverneuropil(
    plotmin=100,
    plotmax=28000,
    cellrange=False,
    F=None,
    Fneu=None,
    iscell=None,
    figname="foverfneu",
):
    """ Calculate DF/F0soma / DF/F0neuropil. """

    Flist = []
    Fneulist = []
    traceindex = []  # keep track of trace index.

    if not cellrange:
        for index, trace in enumerate(F):
            if iscell[index][0] == 1:
                Flist.append(trace[plotmin:plotmax])
                Fneulist.append(Fneu[index][plotmin:plotmax])
                traceindex.append(index)
    else:
        for trace in cellrange:
            if iscell[trace][0] == 1:
                Flist.append(F[trace][plotmin:plotmax])
                Fneulist.append(Fneu[trace][plotmin:plotmax])
                traceindex.append(trace)

    stons = []
    somadfs = []
    neudfs = []
    somaovernpil_raw = []
    if (len(Flist)) >= 80:
        limit = 80
    else:
        limit = len(Flist)

    for index in range(limit):
        somadf = Flist[index] / lowerperc(Flist[index])
        neudf = Fneulist[index] / lowerperc(Fneulist[index])
        stons.append(upperperc(somadf) / upperperc(neudf))
        somadfs.append(upperperc(somadf))
        neudfs.append(upperperc(neudf))
        somaovernpil_raw.append(Flist[index].mean() / Fneulist[index].mean())
        #print("Max df/f0 in trace #" + str(traceindex[index]) + " " + str(max(somadf)))

    print("Mean max df/f0", sum(somadfs) / len(somadfs))
    print("Mean max neuropil df/f0", sum(neudfs) / len(neudfs))
    print("DF/F0soma / DF/F0neuropil:", sum(stons) / len(stons))

    soverneuframe = pd.DataFrame(
        {
            "Trace": traceindex[:limit],
            "Soma DF/F0": somadfs,
            "Npil DF/F0": neudfs,
            "Soma/Npil DF/F0": stons,
            "Soma/Npil raw": somaovernpil_raw,
        }
    )
    soverneuframe.to_csv(figname + "somaoverneuropil.txt")


def lowerperc(data=None):
    """
    Return mean of lower fifteenth of data.
    """
    percent = 15 # of data to include in mean
    data = np.array(sorted(data))
    limit = int(percent * data.size / 100.0)
    return data[0:limit].mean()


def upperperc(data=None):
    """
    Return mean of upper .1% of data.
    """
    percent = 99.9  # of data to include in mean
    data = np.array(sorted(data))
    limit = int(percent * data.size / 100.0)
    return data[limit:].mean()


def scalebar(
    hideTicks=True,
    hideFrame=True,
    fontSize=7,
    scaleXsize=None,
    scaleYsize=None,
    lineWidth=1,
):
    """
    Add an L-shaped scalebar to the current figure.
    This removes current axis labels, ticks, and the figure frame.
    *** Adapted from Scott W. Harden's pyABF package. ***
    """

    scaleXunits = "Seconds"
    scaleYunits = "%\u0394F/F\u2080"

    # calculate the current data area
    x1, x2, y1, y2 = plt.axis()  # bounds
    xc, yc = (x1 + x2) / 2, (y1 + y2) / 2  # center point
    xs, ys = abs(x2 - x1), abs(y2 - y1)  # span

    # determine how big we want the scalebar to be
    if not scaleXsize:
        scaleXsize = abs(plt.xticks()[0][1] - plt.xticks()[0][0]) / 2.5
    if not scaleYsize:
        scaleYsize = abs(plt.yticks()[0][1] - plt.yticks()[0][0]) / 2

    # create the scale bar labels
    lblX = str(scaleXsize)
    lblY = str(scaleYsize)

    # prevent units unecessarially ending in ".0"
    if lblX.endswith(".0"):
        lblX = lblX[:-2]
    if lblY.endswith(".0"):
        lblY = lblY[:-2]

    # add units to the labels
    lblX = lblX + " " + scaleXunits
    lblY = lblY + "" + scaleYunits
    lblX = lblX.strip()
    lblY = lblY.strip()

    # determine the dimensions of the scalebar
    scaleBarPadX = 0
    scaleBarPadY = 0
    scaleBarX = x2 - scaleBarPadX * xs
    scaleBarX2 = scaleBarX - scaleXsize
    scaleBarY = y1 + scaleBarPadY * ys
    scaleBarY2 = scaleBarY + scaleYsize

    # determine the center of the scalebar (where text will go)
    scaleBarXc = (scaleBarX + scaleBarX2) / 2
    scaleBarYc = (scaleBarY + scaleBarY2) / 2

    # create a scalebar point array suitable for plotting as a line
    scaleBarXs = [scaleBarX2, scaleBarX, scaleBarX]
    scaleBarYs = [scaleBarY, scaleBarY, scaleBarY2]

    # the text shouldn't touch the scalebar, so calculate how much to pad it
    lblPadMult = 0.005
    lblPadMult += 0.002 * lineWidth
    lblPadX = xs * lblPadMult
    lblPadY = ys * lblPadMult

    # hide the old tick marks
    if hideTicks:
        plt.gca().get_yaxis().set_visible(False)
        plt.gca().get_xaxis().set_visible(False)

    # hide the square around the image
    if hideFrame:
        plt.gca().spines["top"].set_visible(False)
        plt.gca().spines["right"].set_visible(False)
        plt.gca().spines["bottom"].set_visible(False)
        plt.gca().spines["left"].set_visible(False)

    # now do the plotting
    plt.plot(scaleBarXs, scaleBarYs, "k-", lw=lineWidth)
    plt.text(
        scaleBarXc, scaleBarY - lblPadY, lblX, ha="center", va="top", fontsize=fontSize
    )
    plt.text(
        scaleBarX + lblPadX, scaleBarYc, lblY, ha="left", va="center", fontsize=fontSize
    )
