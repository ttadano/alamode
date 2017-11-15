#!/usr/bin/env python
#
# plotband.py
#
# Simple script to visualize phonon dispersion relations
#
# Copyright (c) 2014 Terumasa Tadano
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#

import numpy as np
import optparse
import matplotlib as mpl
try:
    mpl.use("Qt5agg")
except:
    pass
import matplotlib.pyplot as plt

# parser options
usage = "usage: %prog [options] file1.bands file2.bands ... "
parser = optparse.OptionParser(usage=usage)

parser.add_option("--nokey", action="store_false", dest="print_key", default=True,
                  help="don't print the key in the figure")
parser.add_option("-u", "--unit", action="store", type="string", dest="unitname", default="kayser",
                  help="print the band dispersion in units of UNIT. Available options are kayser, meV, and THz", metavar="UNIT")
parser.add_option("--emin", action="store", type="float", dest="emin",
                  help="minimum value of the energy axis")
parser.add_option("--emax", action="store", type="float", dest="emax",
                  help="maximum value of the energy axis")
parser.add_option("--normalize", action="store_true", dest="normalize_xaxis", default=False,
                  help="normalize the x axis to unity.")


# font styles 
mpl.rc('font', **{'family': 'Times New Roman', 'sans-serif': ['Helvetica']})
mpl.rc('xtick', labelsize = 16)
mpl.rc('ytick', labelsize = 16)
mpl.rc('axes' , labelsize = 16)
mpl.rc('lines', linewidth = 1.5)
mpl.rc('legend', fontsize='small')
#mpl.rc('text', usetex=True)
# line colors and styles
color = ['b', 'g', 'r', 'm', 'k', 'c', 'y', 'r']
lsty = ['-', '-', '-', '-', '--', '--', '--', '--']


def get_kpath_and_kval(file_in):

    ftmp = open(file_in, 'r')
    kpath = ftmp.readline().rstrip('\n').split()
    kval = ftmp.readline().rstrip('\n').split()
    ftmp.close()

    if kpath[0] == '#' and kval[0] == '#':
        kval_float = [float(val) for val in kval[1:]]
        kpath_mod = []
        for i in range(len(kpath[1:])):
            if kpath[i + 1] == 'G':
                kpath_mod.append('$\Gamma$')
            else:
                kpath_mod.append('$\mathrm{' + kpath[i + 1] + '}$')
        return kpath_mod, kval_float
    else:
        return [], []


def change_scale(array, str_scale):

    str_tmp = str_scale.lower()

    if str_tmp == 'kayser':
        print("Band structure will be shown in units of cm^{-1}")
        return array

    elif str_tmp == 'mev':
        print("Band structure will be shown in units of meV")
        kayser_to_mev = 0.0299792458 * 1.0e+12 * \
            6.62606896e-34 / 1.602176565e-19 * 1000

        for i in range(len(array)):
            for j in range(len(array[i])):
                for k in range(1, len(array[i][j])):
                    array[i][j][k] *= kayser_to_mev

        return array

    elif str_tmp == 'thz':
        print("Band structure will be shown in units of THz")
        kayser_to_thz = 0.0299792458

        for i in range(len(array)):
            for j in range(len(array[i])):
                for k in range(1, len(array[i][j])):
                    array[i][j][k] *= kayser_to_thz

        return array

    else:
        print("Unrecognizable option for --unit %s" % str_scale)
        print("Band structure will be shown in units of cm^{-1}")
        return array


def normalize_to_unity(array, array_axis):

    for i in range(len(array)):
        max_val = array[i][-1][0]

        factor_normalize = 1.0 / max_val

        for j in range(len(array[i])):
            array[i][j][0] *= factor_normalize

    max_val = array_axis[-1]
    factor_normalize = 1.0 / max_val

    for i in range(len(array_axis)):
        array_axis[i] *= factor_normalize

    return array, array_axis


def get_xy_minmax(array):

    xmin, xmax, ymin, ymax = [0, 0, 0, 0]

    for i in range(len(array)):
        xtmp = array[i][-1][0]
        xmax = max(xmax, xtmp)

    for i in range(len(array)):
        for j in range(len(array[i])):
            for k in range(1, len(array[i][j])):
                ytmp = array[i][j][k]
                ymin = min(ymin, ytmp)
                ymax = max(ymax, ytmp)

    return xmin, xmax, ymin, ymax


if __name__ == '__main__':
    '''
    Simple script for visualizing phonon dispersion relations.
    Usage:
    $ python plot_band.py [options] file1.bands file2.bands ...

    For details of available options, please type
    $ python plot_band.py -h
    '''

    options, args = parser.parse_args()
    files = args[0:]
    nfiles = len(files)

    if nfiles == 0:
        print("Usage: plotband.py [options] file1.bands file2.bands ...")
        print("For details of available options, please type\n$ python plotband.py -h")
        exit(1)
    else:
        print("Number of files = %d" % nfiles)

    xtickslabels, xticksvars = get_kpath_and_kval(files[0])
    data_merged = []

    for file in files:
        data_tmp = np.loadtxt(file, dtype=float)
        data_merged.append(data_tmp)

    data_merged = change_scale(data_merged, options.unitname)

    if options.normalize_xaxis:
        data_merged, xticksvars = normalize_to_unity(data_merged, xticksvars)

    xmin, xmax, ymin, ymax = get_xy_minmax(data_merged)
    fig, ax = plt.subplots()
    
    for i in range(len(data_merged)):
        ax.plot(data_merged[i][0:, 0], data_merged[i][0:, 1],
                 linestyle=lsty[i], color=color[i], label=files[i])

        for j in range(2, len(data_merged[i][0][0:])):
            ax.plot(data_merged[i][0:, 0], data_merged[i][0:, j],
                     linestyle=lsty[i], color=color[i])

    if options.unitname.lower() == "mev":
        ax.set_ylabel("Frequency (meV)", labelpad=20)
    elif options.unitname.lower() == "thz":
        ax.set_ylabel("Frequency (THz)", labelpad=20)
    else:
        ax.set_ylabel("Frequency (cm${}^{-1}$)", labelpad=10)

    if options.emin == None and options.emax == None:
        factor = 1.05
        ymin *= factor
        ymax *= factor
    else:
        if options.emin != None:
            ymin = options.emin
        if options.emax != None:
            ymax = options.emax

        if ymin > ymax:
            print("Warning: emin > emax")

    plt.axis([xmin, xmax, ymin, ymax])

    ax.set_xticks(xticksvars[0:])
    ax.set_xticklabels(xtickslabels[0:])
    ax.xaxis.grid(True, linestyle='-')

    if options.print_key:
        plt.legend(loc='best', prop={'size': 10})

#	plt.savefig('band_tmp.png', dpi=300, transparent=True)
    plt.show()
