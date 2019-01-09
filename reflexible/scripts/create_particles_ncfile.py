#!/usr/bin/env python
from __future__ import print_function

"""Script to convert FLEXPART particle dumps into a NetCDF4 file.

:Author: Matthew Krupcale
:Contact:  mkrupcale@matthewkrupcale.com
:Created:  2017-01-15

This script demonstrates how to convert FLEXDATA partucle dump files
into a file with netCDF4 format.  Help to use this script can be get with::

  $ create_particles_ncfile.py -h

"""

import sys
import warnings
import platform
import getpass
import datetime
import os.path

import netCDF4 as nc
import numpy as np

from reflexible.scripts import (
    Header, read_command, read_releases, get_fpdirs, read_conffiles,
    read_species)

UNITS = ['conc', 'pptv', 'time', 'footprint', 'footprint_total']
"""Used in combination with H.nested to determine the value of the
   ncid.iout attribute when the COMMAND file is not available.
"""

# Global variables for hosting different command line arguments.
# The default values here are not relevant.
COMPLEVEL = 9
MIN_SIZE = False

def output_units(ncid):
    """The ncid.units attribute computation.

    This function computes the value of the ncid.units attribute required
    when concentration output, wet and dry deposition variables are created.

    Parameters
    ----------
    ncid : Python object
        the object associated to the netCDF4 file

    Return
    ------
    string
        the value of the ncid.units attributte
    """
    # TODO: Check the units here
    if ncid.ldirect == 1:
        # forward simulation
        if ncid.ind_source == 1:
            units = 'kg'
        else:
            units = 'kg m-3'
    else:
        # backward simulation
        if ncid.ind_receptor == 1:
            units = 'kg2 m-3'
        else:
            units = ''
    return units

def write_header(H, ncid, write_releases, species):
    """Create netCDF4 dimensions and variables.

    Create the netCDF4 variables (and the required dimensions) that will be
    stored in the netCDF4 file.

    Parameters
    ----------
    H : Python object
      The Header object.
    ncid : Python object
      The netCDF4 file object.
    write_releases : boolean
      True if releases have to be written into the netCDF4 file.
    species :
      Species data to be written into the netCDF4 file.

    """
    iout = ncid.iout

    # Create dimensions

    # time
    ncid.createDimension('time', None)
    adate, atime = str(H.ibdate), str(H.ibtime).zfill(6)
    timeunit = 'seconds since ' + adate[:4] + '-' + adate[4:6] + \
        '-' + adate[6:8] + ' ' + atime[:2] + ':' + atime[2:4]

    # observation
    ncid.createDimension('obs', None)
    # number of species
    ncid.createDimension('numspec', H.nspec)
    # number of release points
    ncid.createDimension('pointspec', H.numpointspec)
    # number of age classes
    ncid.createDimension('nageclass', H.nageclass)
    # number of actual release points
    ncid.createDimension('numpoint', H.numpoint)

    # Create variables

    # time
    tID = ncid.createVariable('time', 'i4', ('time',),
                              zlib=True, complevel=COMPLEVEL)
    tID.axis = 'T'
    tID.units = timeunit
    tID.calendar = 'proleptic_gregorian'

    # number of particles
    num_particles_ID = ncid.createVariable('num_particles', 'i4', ('time',),
                                           zlib=True, complevel=COMPLEVEL)
    num_particles_ID.long_name = 'number of particles'
    num_particles_ID.description = 'number of particles recorded at this time'
    num_particles_ID.sample_dimension = 'obs'

    # TODO: Fix the 'standard_names' here and elsewhere
    # lon
    lonID = ncid.createVariable('longitude', 'f4', ('obs',),
                                zlib=True, complevel=COMPLEVEL)
    lonID.long_name = 'longitude in degree east'
    lonID.axis = 'X'
    lonID.units = 'degrees_east'
    lonID.standard_name = 'grid_longitude'
    lonID.description = 'particle longitude'

    # lat
    latID = ncid.createVariable('latitude', 'f4', ('obs',),
                                zlib=True, complevel=COMPLEVEL)
    latID.long_name = 'latitude in degree north'
    latID.axis = 'Y'
    latID.units = 'degrees_north'
    latID.standard_name = 'grid_latitude'
    latID.description = 'particle latitude'

    # height
    levID = ncid.createVariable('height', 'f4', ('obs',),
                                zlib=True, complevel=COMPLEVEL)
    levID.axis = 'Z'
    levID.units = 'm'
    levID.positive = 'up'
    levID.standard_name = 'height'
    levID.long_name = 'particle height above ground'

    mass_units = output_units(ncid)

    # Particle mass outputs
    massID = ncid.createVariable('mass', 'f4', ('obs','numspec'),
                                 zlib=True, complevel=COMPLEVEL)
    massID.units = mass_units
    massID.long_name = 'mass of particle'

    if not MIN_SIZE:
        # particle release points
        npointID = ncid.createVariable('npoint', 'i4', ('obs',),
                                       zlib=True, complevel=COMPLEVEL)
        npointID.long_name = 'numpoint'
        npointID.description = 'release point of each particle'

        # particle release times
        itramemID = ncid.createVariable('itramem', 'i4', ('obs',),
                                       zlib=True, complevel=COMPLEVEL)
        itramemID.units = timeunit
        itramemID.calendar = 'proleptic_gregorian'
        itramemID.long_name = 'release time'
        itramemID.description = 'release time of each particle'

        # topography
        topoID = ncid.createVariable('topo', 'f4', ('obs',),
                                     zlib=True, complevel=COMPLEVEL)
        topoID.units = 'm'
        topoID.positive = 'up'
        topoID.long_name = 'topography'
        topoID.description = 'topography at each particle location'

        # potential vorticity
        pvID = ncid.createVariable('pv', 'f4', ('obs',),
                                   zlib=True, complevel=COMPLEVEL)
        pvID.units = 'K m2 kg-1 s-1'
        pvID.long_name = 'potential vorticity'
        pvID.description = 'potential vorticity at each particle location'

        # specific humidity
        qvID = ncid.createVariable('qv', 'f4', ('obs',),
                                    zlib=True, complevel=COMPLEVEL)
        qvID.units = 'kg kg-1'
        qvID.long_name = 'specific humidity'
        qvID.description = 'specific humidity at each particle location'

        # temperature
        ttID = ncid.createVariable('tt', 'f4', ('obs',),
                                   zlib=True, complevel=COMPLEVEL)
        ttID.units = 'K'
        ttID.long_name = 'temperature'
        ttID.description = 'temperature at each particle location'

        # air density
        rhoID = ncid.createVariable('rho', 'f4', ('obs',),
                                    zlib=True, complevel=COMPLEVEL)
        rhoID.units = 'kg m-3'
        rhoID.long_name = 'air density'
        rhoID.description = 'air density at each particle location'

        # planetary boundary layer height
        hmixID = ncid.createVariable('hmix', 'f4', ('obs',),
                                     zlib=True, complevel=COMPLEVEL)
        hmixID.units = 'm'
        hmixID.positive = 'up'
        hmixID.long_name = 'planetary boundary layer height'
        hmixID.description = 'planetary boundary layer height at each' + \
                             'particle location'

        # tropopause
        trID = ncid.createVariable('tr', 'f4', ('obs',),
                                   zlib=True, complevel=COMPLEVEL)
        trID.units = 'm'
        trID.positive = 'up'
        trID.long_name = 'altitude of thermal tropopause'
        trID.description = 'altitude of thermal tropopause at each particle location'

    if write_releases:
        # RELCOM
        relcomID = ncid.createVariable('RELCOM', 'S45', ('numpoint',),
                                       zlib=True, complevel=COMPLEVEL)
        # Fill RELCOM with default values ("NA" means Non-Available)
        relcomID[:] = np.array(["NA"] * H.numpoint, dtype="S45")
        relcomID.long_name = 'release point name'

        # RELLNG1
        rellng1ID = ncid.createVariable('RELLNG1', 'f4', ('numpoint',),
                                        zlib=True, complevel=COMPLEVEL)
        rellng1ID.units = 'degrees_east'
        rellng1ID.long_name = 'release longitude lower left corner'

        # RELLNG2
        rellng2ID = ncid.createVariable('RELLNG2', 'f4', ('numpoint',),
                                        zlib=True, complevel=COMPLEVEL)
        rellng2ID.units = 'degrees_east'
        rellng2ID.long_name = 'release longitude upper right corner'

        # RELLAT1
        rellat1ID = ncid.createVariable('RELLAT1', 'f4', ('numpoint',),
                                        zlib=True, complevel=COMPLEVEL)
        rellat1ID.units = 'degrees_north'
        rellat1ID.long_name = 'release latitude lower left corner'

        # RELLAT2
        rellat2ID = ncid.createVariable('RELLAT2', 'f4', ('numpoint',),
                                        zlib=True, complevel=COMPLEVEL)
        rellat2ID.units = 'degrees_north'
        rellat2ID.long_name = 'release latitude upper right corner'

        # RELZZ1
        relzz1ID = ncid.createVariable('RELZZ1', 'f4', ('numpoint',),
                                       zlib=True, complevel=COMPLEVEL)
        relzz1ID.units = 'm'
        relzz1ID.long_name = 'release height bottom'

        # RELZZ2
        relzz2ID = ncid.createVariable('RELZZ2', 'f4', ('numpoint',),
                                       zlib=True, complevel=COMPLEVEL)
        relzz2ID.units = 'm'
        relzz2ID.long_name = 'release height top'

        # RELKINDZ
        relkindzID = ncid.createVariable('RELKINDZ', 'i4', ('numpoint',),
                                         zlib=True, complevel=COMPLEVEL)
        relkindzID.long_name = 'release kind'

        # RELSTART
        relstartID = ncid.createVariable('RELSTART', 'i4', ('numpoint',),
                                         zlib=True, complevel=COMPLEVEL)
        relstartID.units = 's'
        relstartID.long_name = 'release start relative to simulation start'

        # RELEND
        relendID = ncid.createVariable('RELEND', 'i4', ('numpoint',),
                                       zlib=True, complevel=COMPLEVEL)
        relendID.units = 's'
        relendID.long_name = 'release end relative to simulation start'

        # RELPART
        relpartID = ncid.createVariable('RELPART', 'i4', ('numpoint',),
                                        zlib=True, complevel=COMPLEVEL)
        relpartID.long_name = 'number of release particles'

        # RELXMASS
        relxmassID = ncid.createVariable('RELXMASS', 'f4',
                                         ('numspec', 'numpoint'),
                                         zlib=True, complevel=COMPLEVEL)
        relxmassID.long_name = 'total release particles mass'

    # LAGE
    lageID = ncid.createVariable('LAGE', 'i4', ('nageclass',))
    lageID.units = 's'
    lageID.long_name = 'age class'

    return iout


def write_variables(H, ncid, write_releases, releases):
    """Fill netCDF4 variables with data.

    The netCDF4 variables created in the ``write_header`` function are filled
    with data.

    Parameters
    ----------
    H : Python object
      the header object
    ncid : Python object
      the netCDF4 file object
    write_releases : boolean
      True if releases have to be written into the netCDF4 file.
    """
    iout = ncid.iout

    # Fill variables with data.

    # time
    time_var = ncid.variables['time']
    if ncid.ipout == 1:
        time_var[:] = nc.date2num(H.available_dates_dt, units=time_var.units,
                                  calendar=time_var.calendar)
        dates = H.available_dates
    else:
        time_var[:] = nc.date2num(H.available_dates_dt[-1], units=time_var.units,
                                  calendar=time_var.calendar)
        dates = [H.available_dates[-1],]

    for idt, date in enumerate(dates):
        # read particle data
        try:
            H.read_particles(time_ret=idt)
            particles = H.particles[date]
        except IOError:
            # Oops, we have got an error while reading, so close the file
            ncid.close()
            # and re-raise the error
            raise

    if not MIN_SIZE:
        pass

    if write_releases:
        # release point information
        ncid.variables['RELSTART'][:] = H.ireleasestart
        ncid.variables['RELEND'][:] = H.ireleaseend
        ncid.variables['RELKINDZ'][:] = H.kindz
        ncid.variables['RELLNG1'][:] = H.xp1
        ncid.variables['RELLNG2'][:] = H.xp2
        ncid.variables['RELLAT1'][:] = H.yp1
        ncid.variables['RELLAT2'][:] = H.yp2
        ncid.variables['RELZZ1'][:] = H.zpoint1
        ncid.variables['RELZZ2'][:] = H.zpoint2
        ncid.variables['RELPART'][:] = H.npart
        ncid.variables['RELXMASS'][:, :] = H.xmass.T
        if "release_point_names" in releases:
            relnames = releases["release_point_names"]
            ncid.variables['RELCOM'][:len(relnames)] = relnames

    # Age classes
    ncid.variables['LAGE'][:] = H.lage

def create_ncfile(pathnames, command_path=None, releases_path=None,
                  write_releases=True, dirout=None, outfile=None):
    """Main function that create a netCDF4 file from a FLEXPART output.

    Parameters
    ----------
    pathnames : string
      the file where the FLEXDATA <options> and <output> are specified.
    command_path : string
      path for the associated COMMAND file.
    releases_path : string
      path for the associated RELEASES file.
    write_releases : string
      whether output of release point information.
    dirout : string
      the dir where the netCDF4 file will be created.
    outfile : string
      the complete path of the output file (overrides the ``dirout`` argument)

    Return
    ------
    tuple
      (the path of the netCDF4 file, the options dir, the output dir).
    """

    options_dir, output_dir = get_fpdirs(pathnames)
    H = Header(output_dir)

    if options_dir:
        command = read_conffiles("COMMAND", options_dir, command_path)
        if command['IPOUT'] == 0:
            raise ValueError("Command file indicates that particles are not output. Ensure IPOUT>0.")
        elif command['IPOUT'] == 1:
            fprefix = 'partposit_'
        else:
            fprefix = 'partposit_end_'
        releases = read_conffiles("RELEASES", options_dir, releases_path)
        species = read_species(options_dir, H.nspec)

    if outfile:
        # outfile has priority over previous flags
        ncfname = outfile
    else:
        if dirout is None:
            fprefix = os.path.join(output_dir, fprefix)
        else:
            fprefix = os.path.join(dirout, fprefix)
        ncfname = fprefix + "%s%s" % (H.ibdate, H.ibtime) + ".nc"

    print("About to create new netCDF4 file: '%s'" % ncfname)
    ncid = nc.Dataset(ncfname, 'w')
    write_metadata(H, command, ncid)
    write_header(H, ncid, write_releases, species)
    write_variables(H, ncid, write_releases, releases)
    ncid.close()
    return (ncfname, options_dir, output_dir)


def main():
    """Parses the passed command line arguments.

    The passed arguments will be used to create the netCDF4 file.
    """
    import argparse
    global MIN_SIZE, COMPLEVEL

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d", "--dirout",
        help=("The dir where the netCDF4 file will be created.  "
              "If not specified, then the <output> dir is used.")
        )
    parser.add_argument(
        "-o", "--outfile",
        help=("The complete path for the output file."
              "This overrides the --dirout flag."))
    parser.add_argument(
        "-C", "--command-path",
        help=("The path for the associated COMMAND file.  "
              "If not specified, then the <options>/COMMAND is used.")
        )
    parser.add_argument(
        "-R", "--releases-path",
        help=("The path for the associated RELEASES file.  "
              "If not specified, then the <options>/RELEASES is used.")
        )
    parser.add_argument(
        "-r", "--dont-write-releases", action="store_true",
        help=("Don't write release point information.")
        )
    parser.add_argument(
        "--min-size", dest="min_size", action="store_true",
        help=("Do not write redundant fields (orography) so as to reduce "
              "netCDF4 file size.")
        )
    parser.add_argument(
        "--complevel", type=int, default=9,
        help="Compression level for the netCDF4 file."
        )
    parser.add_argument(
        "pathnames", nargs="?",
        help="The Flexpart pathnames file stating where options and output are."
        )

    args = parser.parse_args()

    MIN_SIZE = args.min_size
    COMPLEVEL = args.complevel

    if args.pathnames is None:
        # The FLEXDATA pathnames file is mandatory
        parser.print_help()
        sys.exit(1)

    ncfname, options_dir, output_dir = create_ncfile(
        args.pathnames, args.command_path, args.releases_path,
        not args.dont_write_releases, args.dirout, args.outfile)

    print("New netCDF4 file is available in: '%s'" % ncfname)

if __name__ == '__main__':
    main()
