########### Particle Reading Routines ###############

from __future__ import print_function

import itertools
import datetime
import os

import numpy as np

import reflexible.conv2netcdf4
from .FortFlex import readnumpart, readparticles

def read_particles(H, **kwargs):
    """
    Accepts a header object as input, returns dictionary of particle values
    keyed by datestring from H['available_dates'].

    **DEPENDENCY**
        Requires FortFlex.so module compiled using f2py. See FortFlex.f for more details.

    Usage::

        > FLEXDATA = read_particles(H, **kwargs)

    Returns:

        A particle dictionary keyed by date strings.

        FLEXDATA[datestring]['itime']
        FLEXDATA[datestring]['numpart']
        FLEXDATA[datestring]['xmass']
        FLEXDATA[datestring]['xlon']
        FLEXDATA[datestring]['ylat']
        FLEXDATA[datestring]['z']
        FLEXDATA[datestring]['npoint']
        FLEXDATA[datestring]['itramem']
        FLEXDATA[datestring]['topo']
        FLEXDATA[datestring]['pv']
        FLEXDATA[datestring]['qv']
        FLEXDATA[datestring]['rho']
        FLEXDATA[datestring]['hmix']
        FLEXDATA[datestring]['tr']
        FLEXDATA[datestring]['tt']

    Arguments

      .. tabularcolumns::  |l|L|

      =============         ========================================
      keyword               Description [default]
      =============         ========================================
      date                  which yyyymmddhhmmss from available_dates
                            or use (time_ret)
      time_ret              index to time
      ipout                 1 or 2 for particle output frequency
      BinaryFile            Use BinaryFile vs. FortFlex [False]
      calcfoot              Will cause footprint to be calculated
                            [False]
      verbose               more output
      =============         ========================================


    .. note::
        most arguments are able to be extracted from the header "H"

    """
    # set up the return dictionary (FLEXDATA updates fd, fd is returned)
    FLEXDATA = {}
    fd = reflexible.conv2netcdf4.Structure()

    # OPS is the options Structure, sets defaults, then update w/ kwargs
    fd.options = OPS = reflexible.conv2netcdf4.Structure()
    OPS.unit = H.unit
    # allows to select an index of npsec when calling readgrid
    OPS.time_ret = 0
    OPS.date = None
    OPS.verbose = False
    OPS.BinaryFile = False
    # add keyword overrides and options to header
    OPS.update(kwargs)
    # H.update(OPS)

    # get times to return
    get_dates = None
    if OPS.time_ret is not None:
        get_dates = []
        time_ret = OPS.time_ret
        if isinstance(time_ret, (int, np.int64)):
            time_ret = [time_ret]

        if time_ret[0] < 0:
            time_ret = np.arange(len(H.available_dates))

        for t in time_ret:
            get_dates.append(H.available_dates[t])

    # define what dates to extract if user has explicitly defined a 'date'
    if OPS.date is not None:
        date = OPS.date
        if time_ret is not None:
            Warning("overwriting time_ret variable, date was requested")
        get_dates = []
        if not isinstance(date, list):
            date = date.strip().split(',')
        for d in date:
            try:
                get_dates.append(H.available_dates[H.available_dates.index(d)])
                time_ret = None
            except:
                print("Cannot find date: %s in H['available_dates']\n" % d)

    if get_dates is None:
        raise ValueError("Must provide either time_ret or date value.")
    else:
        # assign dates for indexing fd
        fd.dates = get_dates[:]

    if OPS.ipout is not None:
        ipout = OPS.ipout

    print('getting particles for: ', get_dates)
    # Some pre-definitions
    fail = 0
    # set filename prefix
    prefix = 'partposit_'

    # Determine what module to read, try to use FortFlex before pure Python
    # import the FortFlex / Fortran module
    try:
        from .FortFlex import readnumpart, readparticles
        useFortFlex = True
        print('Using readgrid from FortFlex')
    except:
        useFortFlex = False
        print('Cannot load FortFlex, reverting to BinaryFile.')
    if not useFortFlex:
        # TODO: write _readpartBF
        #readpart = _readpartBF
        OPS.BinaryFile = True
        print('Using BinaryFile')

    # --------------------------------------------------
    # Loop over all times, given in field H['available_dates']
    # --------------------------------------------------

    for datestring in get_dates:
        print(datestring)
        FLEXDATA[datestring] = fdc = dict()
        if ipout == 1:
            filename = os.path.join(H['pathname'], prefix + datestring)
        elif ipout == 2:
            if datestring != H.available_dates[-1]:
                continue
            else:
                filename = os.path.join(H['pathname'], prefix + 'end')
        if os.path.exists(filename):
            H.filename = filename
            print('reading: ' + filename)
            if OPS.verbose:
                print('with values:')
                inputvars = ['filename']
                for v in inputvars:
                    print(v, " ==> ", H[v])

            if OPS.BinaryFile:
                print("Reading {0} with BinaryFile".format(filename))
            else:
                itime, numpart = readnumpart(filename)
                itime, npoint, xlon, ylat, z, itramem, topo, pv, qv, rho, hmix, tr, tt, xmass = readparticles(filename, numpart, H.nspec)
                print("Read %d particles" % numpart)

                fdc['itime'] = itime

                fdc['timestamp'] = \
                    datetime.datetime.strptime(datestring, '%Y%m%d%H%M%S')
                fdc['numpart'] = numpart
                fdc['xlon'] = xlon
                fdc['ylat'] = ylat
                fdc['z'] = z
                fdc['xmass'] = xmass
                fdc['npoint'] = npoint
                fdc['itramem'] = itramem
                fdc['topo'] = topo
                fdc['pv'] = pv
                fdc['qv'] = qv
                fdc['rho'] = rho
                fdc['hmix'] = hmix
                fdc['tr'] = tr
                fdc['tt'] = tt
        else:
            print('***ERROR: file %s not found! \n' % filename)
            fail = 1

        fd.set_with_dict(FLEXDATA)

    return fd
