#!/usr/bin/env python

import sys
import os
import shutil
import subprocess
import numpy as np


def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


def gen_alminput_si(fname, norder=1, prefix='si222', dfset='DFSET'):

    pos = [[0.000, 0.000, 0.000],
           [0.000, 0.000, 0.500],
           [0.000, 0.250, 0.250],
           [0.000, 0.250, 0.750],
           [0.000, 0.500, 0.000],
           [0.000, 0.500, 0.500],
           [0.000, 0.750, 0.250],
           [0.000, 0.750, 0.750],
           [0.125, 0.125, 0.125],
           [0.125, 0.125, 0.625],
           [0.125, 0.375, 0.375],
           [0.125, 0.375, 0.875],
           [0.125, 0.625, 0.125],
           [0.125, 0.625, 0.625],
           [0.125, 0.875, 0.375],
           [0.125, 0.875, 0.875],
           [0.250, 0.000, 0.250],
           [0.250, 0.000, 0.750],
           [0.250, 0.250, 0.000],
           [0.250, 0.250, 0.500],
           [0.250, 0.500, 0.250],
           [0.250, 0.500, 0.750],
           [0.250, 0.750, 0.000],
           [0.250, 0.750, 0.500],
           [0.375, 0.125, 0.375],
           [0.375, 0.125, 0.875],
           [0.375, 0.375, 0.125],
           [0.375, 0.375, 0.625],
           [0.375, 0.625, 0.375],
           [0.375, 0.625, 0.875],
           [0.375, 0.875, 0.125],
           [0.375, 0.875, 0.625],
           [0.500, 0.000, 0.000],
           [0.500, 0.000, 0.500],
           [0.500, 0.250, 0.250],
           [0.500, 0.250, 0.750],
           [0.500, 0.500, 0.000],
           [0.500, 0.500, 0.500],
           [0.500, 0.750, 0.250],
           [0.500, 0.750, 0.750],
           [0.625, 0.125, 0.125],
           [0.625, 0.125, 0.625],
           [0.625, 0.375, 0.375],
           [0.625, 0.375, 0.875],
           [0.625, 0.625, 0.125],
           [0.625, 0.625, 0.625],
           [0.625, 0.875, 0.375],
           [0.625, 0.875, 0.875],
           [0.750, 0.000, 0.250],
           [0.750, 0.000, 0.750],
           [0.750, 0.250, 0.000],
           [0.750, 0.250, 0.500],
           [0.750, 0.500, 0.250],
           [0.750, 0.500, 0.750],
           [0.750, 0.750, 0.000],
           [0.750, 0.750, 0.500],
           [0.875, 0.125, 0.375],
           [0.875, 0.125, 0.875],
           [0.875, 0.375, 0.125],
           [0.875, 0.375, 0.625],
           [0.875, 0.625, 0.375],
           [0.875, 0.625, 0.875],
           [0.875, 0.875, 0.125],
           [0.875, 0.875, 0.625]]

    with open(fname, 'w') as f:
        f.write(
            "&general\n PREFIX = %s; MODE = opt; NAT = 64; NKD = 1; KD = Si\n/\n" % prefix)
        f.write("&optimize\nDFSET = %s\n/\n" % dfset)
        f.write("&interaction\nNORDER = %d\n/\n" % norder)
        f.write("&cell\n 20.406\n 1.0 0.0 0.0\n 0.0 1.0 0.0\n 0.0 0.0 1.0\n/\n")
        f.write("&cutoff\n Si-Si None 7.6\n/\n")
        f.write("&position\n")
        for x in pos:
            f.write("1 %f %f %f\n" % (x[0], x[1], x[2]))
        f.write("/\n")


def run_alm_si(almbin, project_root):

    shutil.copy('%s/example/Si/reference/DFSET_harmonic' %
                project_root, 'DFSET_harmonic')
    shutil.copy('%s/example/Si/reference/DFSET_cubic' %
                project_root, 'DFSET_cubic')
    with open('DFSET_merged', 'w') as f:
        subprocess.run(['cat', 'DFSET_harmonic', 'DFSET_cubic'], stdout=f)
    gen_alminput_si('ALM1.in', 1, dfset='DFSET_harmonic', prefix='si222')
    gen_alminput_si('ALM2.in', 2, dfset='DFSET_merged', prefix='si222_cubic')
    try:
        with open('ALM1.log', 'w') as f:
            subprocess.run([almbin, 'ALM1.in'], stdout=f)
    except:
        return 1
    try:
        with open('ALM2.log', 'w') as f:
            subprocess.run([almbin, 'ALM2.in'], stdout=f)
    except:
        return 1

    return 0


def check_consistency_alm(project_root, abs_tol=0.01, rel_tol=1.0e-9):

    fname_ref = '%s/example/Si/reference/si222.fcs' % project_root
    data_ref = np.loadtxt(fname_ref,
                          comments=['#', '*'],
                          skiprows=15,
                          max_rows=25,
                          usecols=[2, ])
    data_now = np.loadtxt('si222.fcs',
                          comments=['#', '*'],
                          skiprows=15,
                          max_rows=25,
                          usecols=[2, ])
    isclose_all = True
    for val1, val2 in zip(data_ref, data_now):
        isclose_all = isclose_all & isclose(val1, val2, abs_tol)

    if not isclose_all:
        return 1

    fname_ref = '%s/example/Si/reference/si222_cubic.fcs' % project_root
    data_ref = np.loadtxt(fname_ref,
                          comments=['#', '*'],
                          skiprows=15,
                          max_rows=64,
                          usecols=[2, ])
    data_now = np.loadtxt('si222_cubic.fcs',
                          comments=['#', '*'],
                          skiprows=15,
                          max_rows=64,
                          usecols=[2, ])
    isclose_all = True
    for val1, val2 in zip(data_ref, data_now):
        isclose_all = isclose_all & isclose(val1, val2, rel_tol, abs_tol)
#        print(isclose_all, isclose(val1, val2, rel_tol, abs_tol), val1, val2, val1-val2)

    if not isclose_all:
        return 1

    return 0


def gen_anphoninput_si(prefix, mode, fname):

    with open(fname, 'w') as f:
        f.write("&general\n PREFIX = %s; MODE = %s; FCSXML = si222_cubic.xml; NKD = 1; KD = Si\n/\n" % (prefix, mode))
        f.write("&cell\n  10.203\n  0.0 0.5 0.5\n  0.5 0.0 0.5\n  0.5 0.5 0.0\n/\n")

        if mode == "phonons":
            f.write("&kpoint\n  1\n  G 0.0 0.0 0.0 X 0.5 0.5 0.0 51\n  X 0.5 0.5 1.0 G 0.0 0.0 0.0 51\n  G 0.0 0.0 0.0 L 0.5 0.5 0.5 51\n/\n")
        elif mode == "RTA":
            f.write("&kpoint\n  2\n  10 10 10\n/\n")


def run_anphon_si(anphonbin, project_root):

    gen_anphoninput_si(prefix='si222', mode="phonons", fname='phband.in')
    gen_anphoninput_si(prefix='si222', mode='RTA', fname='RTA.in')
    try:
        with open('phband.log', 'w') as f:
            subprocess.run([anphonbin, 'phband.in'], stdout=f)
    except:
        return 1
    try:
        with open('RTA.log', 'w') as f:
            subprocess.run([anphonbin, 'RTA.in'], stdout=f)
    except:
        return 1

    return 0


def check_consistency_anphon(project_root, abs_tol=0.01, rel_tol=1.0e-9):

    fname_ref = '%s/example/Si/reference/si222.bands' % project_root
    data_ref = np.loadtxt(fname_ref)
    data_now = np.loadtxt('si222.bands')

    m1, n1 = np.shape(data_ref)
    m2, n2 = np.shape(data_now)

    if not (m1 == m2 and n1 == n2):
        return 1

    isclose_all = True
    for i in range(m1):
        for j in range(n1):
            isclose_all = isclose_all & isclose(
                data_ref[i, j], data_now[i, j], abs_tol)

    if not isclose_all:
        return 1

    fname_ref = '%s/example/Si/reference/si222.kl' % project_root
    data_ref = np.loadtxt(fname_ref)
    data_now = np.loadtxt('si222.kl')

    m1, n1 = np.shape(data_ref)
    m2, n2 = np.shape(data_now)

    if not (m1 == m2 and n1 == n2):
        return 1

    isclose_all = True
    for i in range(m1):
        for j in range(n1):
            isclose_all = isclose_all & isclose(
                data_ref[i, j], data_now[i, j], abs_tol)
            #print(data_ref[i,j], data_now[i,j], data_ref[i,j]- data_now[i,j])

    if not isclose_all:
        return 1

    return 0


def runtest_silicon(almbin, anphonbin, project_root):

    info = run_alm_si(almbin, project_root)

    if info > 0:
        print("ALM code failed to execute.\nPlease check if the alm binary exist at %s" % almbin)
        return 1

    info = check_consistency_alm(project_root, abs_tol=0.01, rel_tol=1.0e-6)

    if info == 0:
        print("Silicon ALM --> pass")
    else:
        print("Silicon ALM --> failed")

    if info > 0:
        return 1

    info = run_anphon_si(anphonbin, project_root)

    if info > 0:
        print("ANPHON code failed to execute.\nPlease check if the alm binary exist at %s" % anphonbin)
        return 1

    info = check_consistency_anphon(project_root, abs_tol=0.1, rel_tol=0.01)

    if info == 0:
        print("Silicon ANPHON --> pass")
    else:
        print("Silicon ANPHON --> failed")

    return info


if __name__ == '__main__':

    #project_root = sys.argv[1]

    build_dir = os.getcwd()
    project_root = os.path.dirname(build_dir)

    dirname = '%s/test/si' % project_root
    if not os.path.exists(dirname):
        os.mkdir(dirname)
    os.chdir(dirname)

    almbin = "%s/_build/alm/alm" % project_root
    anphonbin = "%s/_build/anphon/anphon" % project_root

    info = runtest_silicon(almbin, anphonbin, project_root)

    if info == 0:
        sys.exit(0)
    else:
        sys.exit(1)
