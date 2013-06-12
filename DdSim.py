#!/usr/bin/env python
"""Directional D simulation modules"""
__author__ = "Martin Zackrisson"
__copyright__ = "Swedish copyright laws apply"
__credits__ = ["Martin Zackrisson"]
__license__ = "GPL v3.0"
__version__ = "1.0"
__maintainer__ = "Martin Zackrisson"
__email__ = "martin.zackrisson@gu.se"
__status__ = "Development"

#
#       DEPENDENCIES
#

#SIMULATION
try:
    import simuPOP
except:
    print "WARNING: No simulation can be run until simuPOP is installed"

#PLOTTING
#import matplotlib
#matplotlib.use('Agg')
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt

#DATA
import numpy as np
from scipy.stats import pearsonr
from scipy.interpolate import spline

#UTILS
import time
import string

#
#       INTERACTIVE VIEW PLOTTING
#

plt.ion()

#
#       SIMULATION OUTPUT FILES
#
#
# If simulation is called without fpath argument, then
# FH is asumed to be a filepointer.
# A Fake_File instance can also be used instead of the
# string in the fpath to replace the filepointer.

FH = None


class Fake_File(object):

    def __init__(self):
        """An instance of Fake_File can be used instead of true file-pointers
        if user wishes to store the interim data in memory, however due to
        compatibility with downstream analysis this solution is actually slower
        than writing to a file,

        The class supports the following commands:

        write(str)
        read()
        close()

        The class can also be iterated:

            fh = Fake_File()
            for row in fh:
                pass

        """

        self._data = []

    def write(self, line):
        """Emulates the writing of a string.
        The string will be considered a line for the purpose of iteration.
        If this is not true, consider using read() followed by splitting on
        newline.

        :param line:
            A line to be written
        """
        self._data.append(line)

    def read(self):
        """Emulates reading the entire file contents.
        This disregards the emulated rows and returns a single string.
        The string can be split on "\n" to get lines back, given that
        writing was done with explicit line-breaks.

        :return:
            The entire contents of the data-stor as a string
        """
        return "\n".join(map(str, self._data))

    def __iter__(self):
        """To simulate the iteration described in the init-docstring.
        """
        for row in self._data:
            yield str(row)

    def close(self):
        """Placeholder function. Note that nothing changes to the
        Fake_File on calling this function. Further data can be
        added and all data is still in memory.
        The purpose if to have the same interface as a true file-pointer
        """

        pass

#
#       WRITING SIMULATION OUTPUTS IN VARIOUS FORMATS
#
# Note
#
# Currently writeMigrate only works for simulations with one locus
# due to unknown format in multi-loci migrate files
#


def writeCSV(fpath, Y, generation=-1):
    """Writes the simulation results (Y) for one generation as a
    comma separated vector file.

    :param fpath:
        Path to output-file. If file exists at path it
        will be overwritten without promptin question

    :param Y:
        A numpy array containing the allel frequencies,
        the second value in the return tuple of makeNParray

    :param generation:
        (Default value -1), specifice which generation's
        results to output. Negative numbers refers to
        positions for last, -1 being the last generation.
    """

    A = Y[generation, ...]
    nPop, nLoci, nAlleles = A.shape

    fh = open(fpath, 'w')

    #This makes the header row using numbers for pops starting at 1
    header = "\t".join(['pop'] + map(str, range(1, nPop + 1)))
    fh.write("{0}\n".format(header))

    #For each row(allele) write the data
    for l in range(nLoci):
        for a in range(nAlleles):
            data = "\t".join(['allele {0}:{1}'.format(l + 1, a + 1)] +
                             map(str, A[:, l, a].tolist()))
            fh.write("{0}\n".format(data))

    fh.close()


def writeMigrate(fpath, pop):
    """Writes the simulated results (Y) for one generation in
    migrate-n compatible format.

    NOTE: This only currently works for one locus simulations.

    :param fpath:
        Path to output-file. If file exists at path it
        will be overwritten without promptin question

    :param pop:
        The population object as returned from the simulation
    """
    nPop = len(pop.subPopSizes())
    nLoci = len(pop.individual(0).genotype(0))
    title = "simuPOP simulation {0}".format(time.ctime(time.time()))

    fh = open(fpath, 'w')

    fh.write("{0}\t{1}\t.  {2}\n".format(nPop, nLoci, title))

    popheader = "{0}\tPopulation{1}\n"
    a_entry = "{0}{1}\t{2}\n"
    locus_entry = "{0}.{1}"

    #HERE SHOULD ITERATE OVER LOCI....

    for popI in range(nPop):

        #GET POP NAME AS A, B, C etc
        popName = string.uppercase[popI]
        #GET POP SIZE
        popSize = pop.subPopSize(popI)

        #WRITE POPULATION HEADER
        fh.write(popheader.format(popSize, popName))

        for i, ind in enumerate(pop.individuals(popI)):
            fh.write(a_entry.format(
                popName, i,
                "\t".join([locus_entry.format(a, b) for a, b in
                           zip(ind.genotype(0), ind.genotype(1))])))

    fh.close()


def print_M(M):
    """Outputs M in simple to read format to stdout

    :param M:
        A 2D migration array
    """
    for row in M:
        for val in row:
            print "{0:.2f}\t".format(val),
        print

#
#       LOADING SIMULATION RESULTS
#
# Helper function for loading the stored output of the simulation
# into numpy arrays so it can be used with the write, plot and
# statistics functions in this module.
#


def makeNParray(fpath, alleles):
    """makeNParray takes two arguments and returns three arrays

    :param fpath:
        A path to the file containing the simulation
        information, or an instance of Fake_File contiannign
        a simulation.

    :param alleles:
        The number of alleles possible per loci as specified
        when executing the simulation.

    :return X, Y, Z:

        X - A numpy 1D array which contains the generations (1, 2, ... N)
        Y - A 4D numpy array with allele frequencies. Axis specifying:
            (generation, population, locus, allel)
        Z - A 2D numpy array with population sizes, axis being:
            (generation, population)
    """
    if isinstance(fpath, Fake_File):
        fh = fpath
    else:
        fh = open(fpath, 'r')

    X = []
    Y = []
    Z = []
    nLoci = None

    for line in fh:
        line = line.strip()
        if line != "":
            d = eval(line)
            X.append(d[0])
            Z.append(d[1])
            sp_y = []
            for sp in d[2:]:
                if nLoci is None:
                    nLoci = len(sp.keys())
                sp_loc_y = [None] * nLoci
                for loc in sp.keys():
                    y = [0] * alleles
                    for a, v in sp[loc].items():
                        y[a] = v
                    sp_loc_y[loc] = y
                sp_y.append(sp_loc_y)
            Y.append(sp_y)

    fh.close()
    return np.array(X), np.array(Y), np.array(Z)

#
#       MAKING PLOTS
#
# Three plots are implemented to evaluate single simulations
#


def plotTop10(X, Y):
    """Plots the on average most common alleles per loci.
    In total 10 subplots are generated.
    If more than one locus in simulation, then top alleles per loci is drawn.

    :param X:
        Numpy array for X-axis
    :param Y:
        Numpy array for the different Y-datas
    :return:
        Matplotlib.Figure
    """
    loci = Y.shape[-2]
    graph_per_loci = int(10 / loci)

    fig = plt.figure()

    loc = -1
    a = 0
    for g in range(10):

        if g % graph_per_loci == 0 and loc < loci - 1:

            loc += 1
            top10 = Y[..., loc, :].mean(axis=0).mean(
                axis=0).argsort()[-10:].tolist()
            top10.reverse()
            a = 0

        ax = fig.add_subplot(5, 2, g + 1)
        ax.plot(X, Y[..., loc, top10[a]])
        ax.set_ylim(0, 1)
        ax.set_title("Locus {0} Allele {1}".format(loc, top10[a]))
        a += 1

    fig.tight_layout()
    return fig


def plotGeneralPopInfo(X, Y, Z, step=10,
                       shift_square_threshold=0.0001):
    """Plots general information about the populations, such as
    how many alleles are present in each and how many of those
    still have significant movement over time at each generation.

    Finally it plots the number of members post migration at each
    generation. This may fluctuate somewhat even when populations
    are kept stable as migration happens after mating.

    :param X:
        Numpy array as returned from makeNParray
    :param Y:
        Numpy array as returned from makeNParray
    :param Z:
        Numpy array as returned from makeNParray
    :param step:
        The delta generations for investigating moving allele frequenciesself.
    :param shift_square_threshold:
        The threshold for reporting an allele as moving.
    :return:
        Matplotlib.Figure
    """

    delta = (Y[step:] - Y[:-step]) ** 2
    fig = plt.figure()

    #Alleles in pop
    ax = fig.add_subplot(3, 1, 1)
    ax.set_title("Number of alleles in population")
    ax.plot(X, (Y > 0).sum(axis=-1).sum(axis=-1))

    #Moving alleles
    ax = fig.add_subplot(3, 1, 2)
    ax.set_title("Alleles with shifted freq after {0} generations".format(
        step))
    ax.plot(X[:-step], (delta > shift_square_threshold).sum(axis=-1).sum(
        axis=-1))

    #Moving alleles
    ax = fig.add_subplot(3, 1, 3)
    ax.set_title("Population sizes")
    ax.plot(X, Z)

    fig.tight_layout()
    return fig


def plotDd(X, Y):
    """Plots the Directional D over time.

    The same data is plotted twice, but half is made transparent
    in each subplot.

    :param X:
        The X as returned from makeNParray
    :param Y:
        The Y as returned from makeNParray
    :return:
        Matplotlib.Figure
    """

    data = get_Dd(Y)

    fig = plt.figure()
    axes = [fig.add_subplot(2, 1, 2), fig.add_subplot(2, 1, 1)]

    nPop = Y.shape[1]
    legend = []
    label_pattern = "{0} -> {1}"
    markers = "o<+*sphx"
    markers = ['-', '--', '-.', ':']
    colors = "bgrmky"

    for source_Pop in range(nPop):
        for sink_Pop in range(nPop):
            if source_Pop != sink_Pop:
                label = label_pattern.format(source_Pop + 1, sink_Pop + 1)
                legend.append(label)
                tmpY = data[source_Pop][sink_Pop]
                tmpY[np.isnan(tmpY)] = 1
                axes[source_Pop < nPop / 2].plot(
                    X, tmpY,
                    markers[sink_Pop],
                    color=colors[source_Pop],
                    label=label)
                axes[not(source_Pop < nPop / 2)].plot(
                    X, tmpY,
                    markers[sink_Pop],
                    color=colors[source_Pop],
                    label=label,
                    alpha=0.15)

    for ax in axes:
        ax.set_ylim((0, 1.1))
        ax.legend(legend, loc='best', ncol=nPop)

    fig.tight_layout()
    return fig


def plotHistograms(data, startPos=0, bins=10, label="{0} -> {1}",
                   smooth=False, subplots=True):
    """
    Plots a composite graph of superimposed histograms.

    :param data:
        A 3D array with the first dimension describing the source
        population, the second the target population and the third
        the experiment iteration.
    :param startPos:
        The index of the thrid dimension for where to start presenting
        the data.
    :param bins:
        A bins argument to Numpy.histogram, either an int for the number
        of bins for each series, a tuple for the bounds or a list of
        specific bin boundries. See Numpy.histogram for more info.
    :param label:
        A string pattern to put as legend label. \{0\} will be replaced
        with a source integer starting at 1, and \{1\} target int.
    :param smooth:
        If histogram line should be smoothened (BETA/Works poorly)
    :param subplots:
        If figure should have subplots (one per source in dataset)
    :return:
        Matplotlib.Figure
    """

    fig = plt.figure()
    XX = []
    YY = []
    Labels = []
    Colors = [
        (142, 56, 142),
        (113, 113, 198),
        (125, 158, 192),
        (56, 142, 142),
        (113, 198, 113),
        (142, 142, 56),
        (197, 193, 170),
        (198, 113, 113),
        (85, 85, 85),
        (170, 170, 170)]

    Axes = []
    pltOnAx = []

    if subplots:
        for source in range(data.shape[0]):
            Axes.append(fig.add_subplot(
                1, data.shape[0], source + 1,
                sharey=(len(Axes) > 0 and Axes[0] or None)))

    else:
        Axes.append(fig.gca())

    for source in range(data.shape[0]):

        if subplots:
            ax = Axes[source]
        else:
            ax = Axes[0]

        for target in range(data.shape[1]):

            if source != target:
                tmpD = data[source, target, :]
                Y, X = np.histogram(tmpD[np.isfinite(tmpD)])
                XX.append(X)
                YY.append(Y)
                Labels.append(label.format(source + 1, target + 1))
                pltOnAx.append(ax)

    minX = (min([x.min() for x in XX]), )
    maxX = (max([x.max() for x in XX]), )
    maxY = max([y.max() for y in YY])

    zeroY = (0, )

    for i in range(len(XX)):

        X = XX[i]
        Y = YY[i]

        X = (X[1:] + X[:-1]) / 2.0

        X = np.r_[minX, X, maxX]
        Y = np.r_[zeroY, Y, zeroY]

        if smooth:

            vX = X
            X = np.linspace(X.min(), X.max(), bins * 100)
            Y = spline(vX, Y, X, order=10)

        c = [v / 256.0 for v in Colors[i % len(Colors)]]

        pltOnAx[i].plot(X, Y, '-', color=c, label=Labels[i])
        pltOnAx[i].fill_between(X, Y, 0, color=c, alpha=0.4)

    for ax in Axes:
        ax.legend(ncol=(subplots and 1 or 2), prop={'size': 'x-small'})
        ax.set_xlim(minX[0], maxX[0])
        ax.set_ylim(0, maxY)
        cl = plt.getp(ax, 'xmajorticklabels')
        plt.setp(cl, fontsize='x-small')

    fig.tight_layout()

    return fig


def plotBoxes(data, startPos=0):
    """
    Plots a composite graph of boxplots and heatmap for data.

    :param data:
        A 3D array with the first dimension describing the source
        population, the second the target population and the third
        the experiment iteration.
    :param startPos:
        The index of the thrid dimension for where to start presenting
        the data.
    :return:
        Matplotlib.Figure
    """
    subSize = 0.8
    xOff = 0.165
    yOff = 0.135
    W = 1 - xOff - 0.1  # 0.11
    H = 1 - yOff - 0.1  # 0.105

    #SETTING UP AXES
    fig = plt.figure()
    nPop = data.shape[1]
    if nPop != 4:
        raise Exception("Not implemented for other than 4 pops")

    majorAx = plt.gca()
    majorAx.set_xlim(0, nPop)
    majorAx.set_ylim(0, nPop)
    #majorAx = AA.Axes(fig, [0, 0, nPop + 1, nPop + 1])
    majorAx.set_ylabel("Source Population")
    majorAx.set_yticks([.5, 1.5, 2.45, 3.4])
    majorAx.set_yticklabels(['4', '3', '2', '1'])
    majorAx.set_xlabel("Target Population")
    majorAx.set_xticks([.6, 1.55, 2.5, 3.45])
    majorAx.set_xticklabels(['1', '2', '3', '4'])
    fig.add_axes(majorAx)
    axes = []
    divVal = float(nPop)
    for target in range(nPop)[::-1]:
        for source in range(nPop):
            ax = plt.axes(
                [xOff + W * (source / divVal),
                 yOff + H * (target / divVal),
                 W * subSize / divVal, H * subSize / divVal])
            fig.add_axes(ax)
            axes.append(ax)

    #axes = [fig.add_subplot(nPop, nPop, i + 1) for i in range(nPop ** 2)]
    #data = get_Dd(Y)

    heatMapVals = []
    plotI = 0
    for sourcePop in range(nPop):
        heatRow = []
        for targetPop in range(nPop):

            if sourcePop != targetPop:

                tmpD = data[sourcePop][targetPop][startPos:]
                tmpCleanD = tmpD[np.isfinite(tmpD)]
                if tmpCleanD.size == 0:
                    r = plt.Rectangle((0, 0), 1, 1)
                    r.set_facecolor('white')
                    axes[plotI].add_patch(r)
                    axes[plotI].set_axis_off()
                    tmpMean = np.nan
                else:
                    axes[plotI].boxplot(tmpCleanD)
                    axes[plotI].set_ylim(0, 1.0)
                    axes[plotI].set_xticklabels([])
                    tl = plt.getp(axes[plotI], 'ymajorticklabels')
                    plt.setp(tl, fontsize='x-small')
                    """
                    axes[plotI].set_xticklabels(["{0} -> {1}".format(
                        sourcePop + 1, targetPop + 1)])
                    """
                    tmpMean = tmpCleanD.mean()

            else:
                tmpMean = np.nan
                if plotI > 0:
                    r = plt.Rectangle((0, 0), 1, 1)
                    r.set_facecolor('grey')
                    axes[plotI].add_patch(r)
                    axes[plotI].set_axis_off()

            print "{0} -> {1}: {2}".format(
                sourcePop + 1, targetPop + 1, tmpMean)
            heatRow.append(tmpMean)
            plotI += 1

        heatMapVals.append(heatRow)

    heatVals = np.array(heatMapVals)
    imax = axes[0].imshow(heatVals, interpolation='nearest', cmap=plt.cm.RdBu,
                          vmin=0, vmax=1)
    il = plt.getp(imax.axes, 'ymajorticklabels')
    plt.setp(il, fontsize='x-small')
    il = plt.getp(imax.axes, 'xmajorticklabels')
    plt.setp(il, fontsize='x-small')
    divider = make_axes_locatable(axes[0])
    cax = divider.append_axes("right", "5%", pad="3%")
    cb = plt.colorbar(imax, cax=cax)
    cb.set_ticks([0, 0.5, 1])
    cb.set_ticklabels(['0', '.5', '1'])
    cl = plt.getp(cb.ax, 'ymajorticklabels')
    plt.setp(cl, fontsize='x-small')

    axes[0].set_axis_on()
    axes[0].set_yticks(np.arange(nPop))
    axes[0].set_xticks(np.arange(nPop))
    axes[0].set_xticklabels(np.arange(nPop) + 1)
    axes[0].set_yticklabels(np.arange(nPop) + 1)
    #axes[0].set_ylabel("Source")
    #axes[0].set_xlabel("Target")

    #fig.tight_layout()

    return fig
#
#       POPULATION GENETICS FUNCTIONS
#
# Population genetics functions used to calculate the
# directional D and migration.


def _get_vps(v):
    """Helper function to emulate vector algebra of equation 2.

    Note: Not the full equation 2 and as implemented takes serveral
    vectors in parallel when dealing with several loci.

    :param v:
        A numpy-array with allele frequencies
    :return:
        A scalar per locus in a numpy array.
    """

    return (v ** 2).sum(axis=-1)


def _hmean(v):
    """Calculates harmonic mean

    :param v:
        A numpy array with loci as last dimension.
    :return:
        The harmonic mean of several loci's D.
    """
    return v.shape[-1] / (1 / v).sum(axis=-1)


def get_D(Y, pop1, pop2, generation=None):
    """Equation 4c, calculates Jost's D on vector-like allele arrays.
    Note that this implimentation works for several loci in parallell.

    :param Y:
        Numpy.array Y as returned from makeNParray
    :param pop1:
        An int specifying the first population
    :param pop2:
        Either a number specifying the second population or
        A numpy array such that it has the information
        on the second population.
    :param generation:
        (Default None), if specified, D will be returned for the
        specified generation(s), else it will be returned for all
        generations.
    :return:
        A Numpy.array holding D per locus
    """

    if generation is None:
        generation = np.s_[:]

    vPop1 = Y[generation, pop1, ...]
    if isinstance(pop2, int):
        vPop2 = Y[generation, pop2, ...]
    else:
        vPop2 = pop2

    D = _get_vps(vPop1 - vPop2) / (_get_vps(vPop1) + _get_vps(vPop2))
    return D


def VMP_f(pop1, pop2):
    """f for calculating the Virtual Migrant Pool as specified in
    equation 7.

    This function should only be used as an argmuent to
    get_Virual_Migrant_Pool

    :param pop1:
        A Numpy.array of allele frequencies
    :param pop2:
        A Numpy.array of allele frequencies
    :return:
        The geometric mean of the composing populations normalized.
    """
    pool = np.sqrt(pop1 * pop2)
    gamma = pool.sum(axis=-1)
    gamma = gamma.reshape(gamma.shape + (1,))
    return pool / gamma


def get_Virual_Migrant_Pool(Y, pop1, pop2, generation=None,
                            f=VMP_f):
    """Common interface for obtaining the Virtual Migrant Pool (VMP)
    independent on which f-function used.

    :param Y:
        A numpy array as returned by makeNParray
    :param pop1:
        An integer specifying population 1
    :param pop2:
        An integer specifying population 2
    :param generation:
        (Default None). If specified, VMP will be constructed
        for the specific generation(s), else it will be
        constructed for all generations.
    :param f:
        (Default VMP_f) A VMP construction fucntion. It should share interface
        with the method VMP_f
    :return:
        A numpy array representing the VMP
    """
    if generation is None:
        generation = np.s_[:]

    P1 = Y[generation, pop1, ...]
    P2 = Y[generation, pop2, ...]

    return f(P1, P2)


def get_Dd(Y, generation=None, D=get_D, **kwargs):
    """Generates a Directional D matrix as default.

    :param Y:
        A numpy array as returned by makeNParray
    :param generation:
        (Default None). If specified, Dd will be constructed
        for the specific generation(s), else it will be
        constructed for all generations.
    :param D:
        (Default get_D). A method for calculating distances.
        It should share the same interface as get_D
    :param **kwargs:
        Key-word arguments are passed along to get_Virual_Migrant_Pool
    :return:
        Dd numpy array
    """

    if not hasattr(Y, 'shape'):
        Y = np.array(Y)

    nPop = Y.shape[1]
    #comparisons = (nPop - 1) * nPop / 2

    data = [[0] * nPop for i in range(nPop)]

    for pop1 in range(nPop):

        for pop2 in range(nPop):

            if pop1 > pop2:

                pool = get_Virual_Migrant_Pool(
                    Y, pop1, pop2, generation=generation, **kwargs)

                data[pop1][pop2] = _hmean(D(Y, pop1, pool,
                                          generation=generation))
                data[pop2][pop1] = _hmean(D(Y, pop2, pool,
                                          generation=generation))

    return data


def get_normed_M(Y, generation=None, D=get_D, **kwargs):
    """Constructs a normalized M array.

    :param Y:
        A numpy array as returned by makeNParray
    :param generation:
        (Default None). If specified, M will be constructed
        for the specific generation(s), else it will be
        constructed for all generations.
    :param D:
        (Default get_D). A method for calculating distances
        for get_Dd to use.
        It should share the same interface as get_D.
    :param **kwargs:
        Further arguments are passed along to get_Dd
    :return:
        A normalized M array
    """
    D = np.array(get_Dd(Y, generation=generation,
                        D=D, **kwargs))

    M = (1 - D) / D
    finM = M[np.isfinite(M)]
    if finM.size > 0:
        M /= M[np.isfinite(M)].max()

    return M

#
#       SIMULATION FUNCTIONS
#
# Functions directly relating to running the simulation and
# gathering data from the simulation


def _getSimulationStep(generations, simulationStep=None):
    """Returning a suggested simulation step if non was
    specified

    :param generations:
        The number of generations being run
    :param simulationStep:
        If not None, this will be the step size, overriding the step
        size calibration
    :return:
        Step size
    """

    if simulationStep is None:
        if generations <= 200:
            s = 10
        else:
            s = 20
    else:
        s = simulationStep

    return s


def _myOUT(pop, param):
    """Method responsible at recording the current state of the
    population at each simulation step

    :param pop:
        The simuPOP simulation
    :param param:
        The number of loci
    :return:
        True
    """
    nLoci = param
    #nLoci = len(pop.lociPos())
    d = [pop.dvars().gen, pop.subPopSizes()]

    for sp in range(len(d[-1])):
        subPop = {}
        for locus in range(nLoci):
            simuPOP.stat(pop, alleleFreq=locus, vars=['alleleFreq_sp'])
            sd = pop.dvars(sp).alleleFreq
            k = sd.keys()[0]
            subPop[k] = sd[k]

        d.append(subPop)

    FH.write("{0}\n".format(d))

    return True


def simuMigration(subPopSizes, migrationMatrix, generations, nAlleles,
                  mutationRates, simulationStep=None,
                  pre_migration_generations=0, loci=[1],
                  stable_pop_sizes=True, microsats=True, fpath=None):
    """Simulate the change of allele frequencies among subpopulations
    as a result of migration.

    The following is specified in simuMigration but not possible to
    change with arguments:

    The simulation is run as a diploid, sexually randoom mating  species.

    :param subPopSizes:
        If a scalar, all subpopulations will have equal starting
        size. The number of subpopulations will be determined from
        the shape of the migrationMatrix.
        If a list or tuple, the subpopulations and sized will be
        set from the values in the tuple.
    :param migrationMatrix:
        A square numpy array of migration frequencies.
        Postions i,j signifies migration from i to j.
    :param generations:
        Number of generations with migration
    :param nAlleles:
        The number of potential alleles at each loci.
    :param mutationRates:
        A scalar determining mutation rates.
    :param simulationStep:
        (Default None). If supplied, this will be the step
        length for the simulation. Else step length will
        be determined by _getSimulationStep and depend on
        the number of generations.
    :param pre_migration_generations:
        (Default 0). If set to positive value, an initial
        phase of evolution without migration will be invoked
        for the number of generations specified. This does
        not affect the number of generations with migration.
        Also note that the pre migration generations will be
        part of the output data.
    :param loci
        (Default [1]) The number of loci per chromosome as a list.
        Example one chromosome, 8 loci:
            loci=[8]
        Example two chromosomes with 2 and 3 loci respectively:
            loci=[2, 3]
    :param stable_pop_sizes:
        (Default True). If true, mating is done such that the
        next generation has the same number as the initial
        generation.
        If false, each parent pair produces two offsprings
        and sizes of the populations will drift due to migration.
    :param microsats:
        (Default True). If true, the StepwiseMutator is invoked
        which emulates microsatelites. That is, allele #9 may
        only mutate to #8 or #10.
        If false KAlleleMutator is invoked and an allele may
        mutate into any other allele.
    :param fpath:
        (default None). If not supplied, then the modules
        FH-object is assumed to be a file pointer and used
        to output simulation data.
        If a Fake_File instance, this will be used to output
        simulation data.
        If a string, a new file pointer will be created
        overwriting any data at the fpath without warnning.
    :return:
        simuPOP Population as it is after the simulation.
    """
    sTime = time.time()
    if fpath:
        global FH
        if isinstance(fpath, Fake_File):
            FH = fpath
        else:
            FH = open(fpath, 'w')

    #DEAL WITH IMPUT (ALLOW BOTH LIST OF POP-SIZES AS WELL AS
    #SCALAR FOR ALL SUBPOP SIZES
    numOfSubPops = migrationMatrix.shape[0]
    if type(subPopSizes) != list:
        subPopSizes = [subPopSizes] * numOfSubPops

    if stable_pop_sizes:
        nextGenSize = subPopSizes
    else:
        nextGenSize = []

    if microsats:
        preOps = [simuPOP.StepwiseMutator(rates=mutationRates,
                                          maxAllele=nAlleles - 1)]
    else:
        preOps = [simuPOP.KAlleleMutator(k=nAlleles, rates=mutationRates)]

    #if type(mutationRates) != list:
    #    mutationRates = [mutationRates] * nAlleles

    nLoci = sum(loci)

    #CREATE POPULATION
    pop = simuPOP.Population(
        size=subPopSizes,
        subPopNames=map(str, range(numOfSubPops)),
        ploidy=2,
        loci=loci,  # Loci on each chromosome
        #chromTypes=[AUTOSOME],
        infoFields=['migrate_to'])

    #MAKE A MIGRATION MATRIX
    a = migrationMatrix.tolist()

    #DETERMINE SIMULATION STEP
    s = _getSimulationStep(generations, simulationStep)

    if pre_migration_generations > 0:

        pop.evolve(
            initOps=[
                simuPOP.InitSex(),
                simuPOP.InitGenotype(freq=[1.0 / nAlleles
                                           for i in range(nAlleles)])
            ],
            preOps=preOps,
            matingScheme=simuPOP.RandomMating(subPopSize=nextGenSize),
            postOps=[
                simuPOP.PyOperator(func=_myOUT, param=nLoci, step=s),
            ],
            gen=pre_migration_generations,
        )

        print "DONE PRE MIGRATION PART", (time.time() - sTime)

        pop.evolve(
            preOps=preOps,
            matingScheme=simuPOP.RandomMating(subPopSize=nextGenSize),
            postOps=[
                simuPOP.Migrator(rate=a),
                simuPOP.PyOperator(func=_myOUT, param=nLoci, step=s),
            ],
            gen=generations,
        )

    else:

        pop.evolve(
            initOps=[
                simuPOP.InitSex(),
                simuPOP.InitGenotype(freq=[1.0 / nAlleles
                                           for i in range(nAlleles)])
            ],
            preOps=preOps,
            matingScheme=simuPOP.RandomMating(subPopSize=nextGenSize),
            postOps=[
                simuPOP.Migrator(rate=a),
                simuPOP.PyOperator(func=_myOUT, param=nLoci, step=s),
            ],
            gen=generations,
        )

    if fpath:
        FH.close()

    print "SIM TOOK", (time.time() - sTime)
    return pop


def runMultiRunStoreFinalMig(nSimulations, evalF=get_Dd, **kwargs):
    """Helper function to repeat the same simulation and extract
    final migration pattern at last generation.

    :param nSimulations:
        Integer determining the number of repeats
    :param evalF:
        The funcion to evaluate the allel frequencies of the last
        generation.
    :param **kwargs:
        Further arguments and keyword arguments are passed along
        to simuMigration, some of these such as 'fpath' and
        'nAlleles' are used by evalSimuSettings as well.
    :return:
        The result of the evalF-method such that the experiment
        iterations are on the last axis.
    """
    if not 'fpath' in kwargs:
        print "Don't know where you are saving, cant do anything"
        return
    if not 'nAlleles' in kwargs:
        print "Don't know the number of alleles, cant do a thing"
        return

    fpath = kwargs['fpath']
    nAlleles = kwargs['nAlleles']
    data = []

    for i in range(nSimulations):
        print "STARTING", i
        simuMigration(**kwargs)
        X, Y, Z = makeNParray(fpath, nAlleles)
        M = evalF(Y, generation=-1)
        data.append(M)
        print "DONE WITH", i

    D = np.array(data)
    D = D.ravel().reshape(D.shape[1:] + (D.shape[0], ), order='F')

    return D


def evalSimuSettings(nSimulations, *args, **kwargs):
    """Helper function to repeat the same simulation and extract
    correlation between the observed and expected migration.

    :param nSimulations:
        Integer determining the number of repeats
    :param *args / **kwargs:
        Further arguments and keyword arguments are passed along
        to simuMigration, some of these such as 'fpath' and
        'nAlleles' are used by evalSimuSettings as well.
    :return:
        A list of pearson correlation r's for each iteration
    """
    if not 'fpath' in kwargs:
        print "Don't know where you are saving, cant do anything"
        return
    if not 'nAlleles' in kwargs:
        print "Don't know the number of alleles, cant do a thing"
        return

    Mig = None
    pearson = []

    fpath = kwargs['fpath']
    nAlleles = kwargs['nAlleles']
    migR = migrationMatrix.ravel()

    for i in range(nSimulations):
        print "STARTING", i
        simuMigration(*args, **kwargs)
        X, Y, Z = makeNParray(fpath, nAlleles)
        M = get_normed_M(Y, generation=-1).ravel()
        M[np.isnan(M)] = 0
        finM = np.isfinite(M)
        #if Mig is None:
        Mig = migR[finM]
        #if migrationMatrix.size - Mig.size != np.sqrt(migrationMatrix.size):
        #    print "You were unlucky, run again"
        #    return

        print "DONE WITH", i

        pearson.append(pearsonr(M[finM], Mig)[0])

    return pearson


def metaEvalSimu(iterations=1000, fpath='pop.data'):
    """Helper function to evaluate a set of different parameters
    with evalSimuSettings.

    A numpy-array of the pearson values for each simulation is
    saved out as well as a histogram for each settings combination.

    :param iterations:
        (Default 1000). The number of iterations for each
        simulation.

    :param fpath:
        (Default 'pop.data'). Where data will be stored throughout
        the meta-evaluation.
    """

    global migrationMatrix
    global migrationMatrix2

    """
    constant_kwargs = {'mutationRates': 0.0001, 'simulationStep': 1,
                       'pre_migration_generations': 500,
                       'fpath': fpath}
    """

    def makeHist(pstats, i, title, fpath, iterations):
        """Method outputs the data as npy file and histogram svg"""
        P = np.array(pstats)
        np.save("{0}.{1}.sim{2}.npy".format(fpath, iterations, i), P)
        P = P[np.isfinite(P)]
        plt.clf()
        plt.hist(P, bins=100)
        plt.gca().set_title("Sim {0} (N={1}): {2}".format(i, iterations, title))
        plt.tight_layout()
        plt.savefig('{0}.{1}.sim{2}.svg'.format(fpath, iterations, i))

    """
    pstats = evalSimuSettings(iterations,
                              subPopSizes=1000,
                              migrationMatrix=migrationMatrix,
                              generations=1000,
                              nAlleles=256,
                              loci=[4],
                              **constant_kwargs)

    makeHist(pstats, 0, '4 loci on one chromosome', fpath, iterations)

    pstats = evalSimuSettings(iterations,
                              subPopSizes=1000,
                              migrationMatrix=migrationMatrix,
                              generations=1000,
                              nAlleles=256,
                              loci=[2, 2],
                              **constant_kwargs)

    makeHist(pstats, 1, '2 by 2 loci', fpath, iterations)

    pstats = evalSimuSettings(iterations,
                              subPopSizes=1000,
                              migrationMatrix=migrationMatrix,
                              generations=1000,
                              nAlleles=256,
                              loci=[8],
                              **constant_kwargs)

    makeHist(pstats, 2, '8 loci', fpath, iterations)

    pstats = evalSimuSettings(iterations,
                              subPopSizes=5000,
                              migrationMatrix=migrationMatrix,
                              generations=1000,
                              nAlleles=256,
                              loci=[4],
                              **constant_kwargs)

    makeHist(pstats, 3, '5000 individiuals per subpop', fpath, iterations)


    pstats = evalSimuSettings(iterations,
                              subPopSizes=1000,
                              migrationMatrix=migrationMatrix,
                              generations=5000,
                              nAlleles=256,
                              loci=[4],
                              **constant_kwargs)

    makeHist(pstats, 4, '5000 generations', fpath, iterations)

    pstats = evalSimuSettings(iterations,
                              subPopSizes=1000,
                              migrationMatrix=migrationMatrix2,
                              generations=1000,
                              nAlleles=256,
                              loci=[4],
                              **constant_kwargs)

    makeHist(pstats, 5, 'Stronger migration', fpath, iterations)

    pstats = evalSimuSettings(iterations,
                              subPopSizes=1000,
                              migrationMatrix=migrationMatrix3,
                              generations=1000,
                              nAlleles=256,
                              loci=[4],
                              **constant_kwargs)

    makeHist(pstats, 6, 'Second migration pattern', fpath, iterations)

    pstats = evalSimuSettings(iterations,
                              subPopSizes=1000,
                              migrationMatrix=migrationMatrix,
                              generations=1000,
                              nAlleles=128,
                              loci=[4],
                              **constant_kwargs)

    makeHist(pstats, 7, '128 Possible alleles per loci', fpath, iterations)

    pstats = evalSimuSettings(iterations,
                              subPopSizes=10000,
                              migrationMatrix=migrationMatrix2,
                              generations=10000,
                              nAlleles=32,
                              loci=[4],
                              **constant_kwargs)

    makeHist(pstats, 8, 'Stronger Migration, 32 Possible alleles per loci, 10000 generations, 10000 individuals / pop', fpath, iterations)

    pstats = evalSimuSettings(iterations,
                              subPopSizes=1000,
                              migrationMatrix=migrationMatrix2,
                              generations=10000,
                              nAlleles=32,
                              loci=[4],
                              **constant_kwargs)

    makeHist(pstats, 9, 'Stronger Migration, 32 alleles, 4 loci, 10000 gen, 1000 individuals', fpath, iterations)
    """


#
#       MIGRATION MATRICES
#
# Ready made migration matrices
#

migrationMatrix = np.array(
    [[0, .0001, .0002, 0],
    [0, 0, .0001, 0],
    [0, 0, 0, 0],
    [0, 0, 0, 0]])

migrationMatrix2 = migrationMatrix * 5

migrationMatrix3 = np.array(
    [[0, .0002, .0001, 0],
    [0, 0, .0002, 0],
    [0, 0, 0, 0],
    [0, 0, 0, 0]])

#
#       RUN-TIME BEHAVIOUR
#
# If file is run from shell it will run the metaEvalSimu method
#

if __name__ == "__main__":

    metaEvalSimu()
