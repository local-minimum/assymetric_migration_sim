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

import simuPOP
import matplotlib
matplotlib.use('Agg')

import numpy as np
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
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
        """
        self._data.append(line)

    def read(self):
        """Emulates reading the entire file contents.
        This disregards the emulated rows and returns a single string.
        The string can be split on "\n" to get lines back, given that
        writing was done with explicit line-breaks.
        """
        return "".join(map(str, self._data))

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

    ** fpath    Path to output-file. If file exists at path it
                will be overwritten without promptin question

    ** Y        A numpy array containing the allel frequencies,
                the second value in the return tuple of makeNParray

    ** generation
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


def writeMigrate(fpath, Y, generation=-1):
    """Writes the simulated results (Y) for one generation in
    migrate-n compatible format.

    NOTE: This only currently works for one locus simulations.

    ** fpath    Path to output-file. If file exists at path it
                will be overwritten without promptin question

    ** Y        A numpy array containing the allel frequencies,
                the second value in the return tuple of makeNParray

    ** generation
                (Default value -1), specifice which generation's
                results to output. Negative numbers refers to
                positions for last, -1 being the last generation.

    """
    generations, nPop, nAlleles = Y.shape
    print "Generation {0} ({1}), Populations {2}, Alleles {3}".format(
        generation, generations, nPop, nAlleles)

    loci = Y.shape[-2]
    title = "simuPOP simulation {0}".format(time.ctime(time.time()))

    fh = open(fpath, 'w')

    fh.write("{0} {1} . {2}\n".format(nPop, loci, title))

    popheader = "{0} Population{1}\n"
    a_entry = "{0}{1}\t{2}\n"
    zpos = len(str(nAlleles + 1))

    #HERE SHOULD ITERATE OVER LOCI....

    loc = 1
    for popI in range(nPop):

        #GET POP NAME AS A, B, C etc
        pop = string.uppercase[popI]

        #WRITE POPULATION HEADER
        fh.write(popheader.format(nAlleles, pop))

        #WRITE ALL ALLELES
        for alleleI in range(nAlleles):
            fh.write(a_entry.format(pop, str(alleleI + 1).zfill(zpos),
                                    Y[generation, popI, loc, alleleI]))

    fh.close()


def print_M(M):
    """Outputs M in simple to read format to stdout"""
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

    ** fpath    A path to the file containing the simulation
                information, or an instance of Fake_File contiannign
                a simulation.

    ** alleles  The number of alleles possible per loci as specified
                when executing the simulation.

    Outputs:

    ** X        A numpy 1D array which contains the generations
                (1, 2, ... N)

    ** Y        A 4D numpy array with allele frequencies.
                Axis specifying:
                (generation, population, locus, allel).

    ** Z        A 2D numpy array with population sizes, axis being:
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

    Arguments X and Y are numpy arrays as returned by makeNParray

    Returns a matplotlib figure.
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


def plotGeneralPopInfo(
        X, Y, Z, step=10, shift_square_threshold=0.0001):
    """Plots general information about the populations, such as
    how many alleles are present in each and how many of those
    still have significant movement over time at each generation.

    Finally it plots the number of members post migration at each
    generation. This may fluctuate somewhat even when populations
    are kept stable as migration happens after mating.

    Arguments X, Y and Z are numpy arrays as returned by makeNParray.

    ** step     The delta generations for investigating moving
                allele frequencies.

    ** shift_square_threshold
                The threshold for reporting an allele as moving.
    Returns a matplotlib figure.
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

    Arguments X and Y are numpy arrays as returned by makeNParray.

    Returns a matplotlib figure.
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
                axes[source_Pop < nPop / 2].plot(
                    X, data[source_Pop][sink_Pop],
                    markers[sink_Pop],
                    color=colors[source_Pop],
                    label=label)
                axes[not(source_Pop < nPop / 2)].plot(
                    X, data[source_Pop][sink_Pop],
                    markers[sink_Pop],
                    color=colors[source_Pop],
                    label=label,
                    alpha=0.15)

    for ax in axes:
        ax.set_ylim(top=1.001)
        ax.legend(legend, loc='best', ncol=nPop)

    fig.tight_layout()
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

    Arguments:

    ** v    A numpy-array with allele frequencies

    Returns a scalar per locus in a numpy array.
    """

    return (v ** 2).sum(axis=-1)


def _hmean(v):
    """Calculates harmonic mean

    Arguments:

    ** v    A numpy array with loci as last dimension.

    Returns the harmonic mean of several loci's D.
    """
    return v.shape[-1] / (1 / v).sum(axis=-1)


def get_D(Y, pop1, pop2, generation=None):
    """Equation 4c, calculates Jost's D on vector-like allele arrays.
    Note that this implimentation works for several loci in parallell.

    Arguments:

    ** Y    Numpy array Y as returned from makeNParray

    ** pop1 An int specifying the first population

    ** pop2 Either a number specifying the second population or
            A numpy array such that it has the information
            on the second population.

    ** generation
            (Default None), if specified, D will be returned for the
            specified generation(s), else it will be returned for all
            generations.

    Returns a numpy array holding D per locus
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

    Arguments:

    ** pop1 A numpy array of allele frequencies

    ** pop2 A numpy array of allele frequencies

    Returns the geometric mean of the composing populations
    normalized.
    """
    pool = np.sqrt(pop1 * pop2)
    gamma = pool.sum()
    return pool / gamma


def get_Virual_Migrant_Pool(Y, pop1, pop2, generation=None,
                            f=VMP_f):
    """Common interface for obtaining the Virtual Migrant Pool (VMP)
    independent on which f-function used.

    Arguments:

    ** Y    A numpy array as returned by makeNParray

    ** pop1 An integer specifying population 1

    ** pop2 An integer specifying population 2

    ** generation
            (Default None). If specified, VMP will be constructed
            for the specific generation(s), else it will be
            constructed for all generations.

    ** f    (Default VMP_f) A VMP construction fucntion. It should share interface
            with the method VMP_f

    Returns a numpy array representing the VMP
    """
    if generation is None:
        generation = np.s_[:]

    P1 = Y[generation, pop1, ...]
    P2 = Y[generation, pop2, ...]

    return f(P1, P2)


def get_Dd(Y, generation=None, D=get_D, **kwargs):
    """Generates a Directional D matrix as default.

    Arguments:

    ** Y    A numpy array as returned by makeNParray

    ** generation
            (Default None). If specified, Dd will be constructed
            for the specific generation(s), else it will be
            constructed for all generations.

    ** D    (Default get_D). A method for calculating distances.
            It should share the same interface as get_D

    further key-word arguments are passed along to get_Virual_Migrant_Pool

    Returns a Dd numpy array
    """
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

    Arguments:

    ** Y    A numpy array as returned by makeNParray

    ** generation
            (Default None). If specified, M will be constructed
            for the specific generation(s), else it will be
            constructed for all generations.

    ** D    (Default get_D). A method for calculating distances
            for get_Dd to use.
            It should share the same interface as get_D.

    Further arguments are passed along to get_Dd

    Returns a normalized M array
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
    specified"""

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

    Arguments:

    ** pop  The simuPOP simulation

    ** param
            The number of loci
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

    Arguments:

    ** subPopSizes
            If a scalar, all subpopulations will have equal starting
            size. The number of subpopulations will be determined from
            the shape of the migrationMatrix.
            If a list or tuple, the subpopulations and sized will be
            set from the values in the tuple.

    ** migrationMatrix
            A square numpy array of migration frequencies.
            Postions i,j signifies migration from i to j.

    ** generations
            Number of generations with migration

    ** nAlleles
            The number of potential alleles at each loci.

    ** mutationRates
            A scalar determining mutation rates.

    ** simulationStep
            (Default None). If supplied, this will be the step
            length for the simulation. Else step length will
            be determined by _getSimulationStep and depend on
            the number of generations.

    ** pre_migration_generations
            (Default 0). If set to positive value, an initial
            phase of evolution without migration will be invoked
            for the number of generations specified. This does
            not affect the number of generations with migration.
            Also note that the pre migration generations will be
            part of the output data.

    ** loci (Default [1]) The number of loci per chromosome as
            a list.
            Example one chromosome, 8 loci:
                loci=[8]
            Example two chromosomes with 2 and 3 loci respectively:
                loci=[2, 3]

    ** stable_pop_sizes
            (Default True). If true, mating is done such that the
            next generation has the same number as the initial
            generation.
            If false, each parent pair produces two offsprings
            and sizes of the populations will drift due to migration.

    ** microsats
            (Default True). If true, the StepwiseMutator is invoked
            which emulates microsatelites. That is, allele #9 may
            only mutate to #8 or #10.
            If false KAlleleMutator is invoked and an allele may
            mutate into any other allele.

    ** fpath
            (default None). If not supplied, then the modules
            FH-object is assumed to be a file pointer and used
            to output simulation data.
            If a Fake_File instance, this will be used to output
            simulation data.
            If a string, a new file pointer will be created
            overwriting any data at the fpath without warnning.

    Returns a simuPOP Population as it is after the simulation.
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


def evalSimuSettings(nSimulations, *args, **kwargs):
    """Helper function to repeat the same simulation and extract
    correlation between the observed and expected migration.

    Arguments:

    ** nSimulations
            Integer determining the number of repeats

    Further arguments and keyword arguments are passed along
    to simuMigration, some of these such as 'fpath' and
    'nAlleles' are used by evalSimuSettings as well.

    Returns a list of pearson correlation r's for each iteration
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

    Arguments:

    ** iterations
            (Default 1000). The number of iterations for each
            simulation.

    ** fpath
            (Default 'pop.data'). Where data will be stored throughout
            the meta-evaluation.
    """
    global migrationMatrix
    global migrationMatrix2
    constant_kwargs = {'mutationRates': 0.0001, 'simulationStep': 1,
                       'pre_migration_generations': 500,
                       'fpath': fpath}

    def makeHist(pstats, i, title, fpath, iterations):
        """Method outputs the data as npy file and histogram svg"""
        P = np.array(pstats)
        np.save("{0}.{1}.sim{2}.npy".format(fpath, iterations, i), P)
        plt.clf()
        plt.hist(P, bins=100)
        plt.gca().set_title("Sim {0} (N={1}): {2}".format(i, iterations, title))
        plt.tight_layout()
        plt.savefig('{0}.{1}.sim{2}.svg'.format(fpath, iterations, i))

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
                              migrationMatrix=migrationMatrix,
                              generations=1000,
                              nAlleles=128,
                              loci=[4],
                              **constant_kwargs)

    makeHist(pstats, 5, '128 Possible alleles per loci', fpath, iterations)

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

#
#       RUN-TIME BEHAVIOUR
#
# If file is run from shell it will run the metaEvalSimu method
#

if __name__ == "__main__":

    metaEvalSimu()
