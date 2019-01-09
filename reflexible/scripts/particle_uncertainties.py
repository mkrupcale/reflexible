#!/usr/bin/env python
from __future__ import print_function

"""Script to use FLEXPART particle dumps to calculate particle mass
means and uncertainties in the FLEXPART grid and output to a NetCDF4 file.

:Author: Matthew Krupcale
:Contact:  mkrupcale@matthewkrupcale.com
:Created:  2017-03-04

This script demonstrates how to use FLEXPART partucle dump files to calculate
mass means and errors in each grid cell and output to a NetCDF4 file. Help to
use this script can be gotten with::

  $ particle_uncertainties.py -h

"""

import sys
import os.path
from datetime import datetime

import netCDF4 as nc
import numpy as np

from reflexible.conv2netcdf4 import (
    Header, read_command, read_releases, get_fpdirs)

from reflexible.scripts import (
    read_conffiles, read_species, write_metadata)

DATE_FORMAT_STR = "%Y-%m-%d"
TIME_FORMAT_STR = "%H:%M:%S"
DATETIME_FORMAT_STR = DATE_FORMAT_STR + " " + TIME_FORMAT_STR

# Global variables for hosting different command line arguments.
DEFAULT_COMPLEVEL = 9

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

def write_header(H, ncid, write_releases, species, nlon, nlat, nz,
                 min_size=False, complevel=DEFAULT_COMPLEVEL):
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

    # level
    ncid.createDimension('height', nz)

    # lat
    ncid.createDimension('latitude', nlat)

    # lon
    ncid.createDimension('longitude', nlon)

    # number of species
    ncid.createDimension('numspec', H.nspec)
    # number of age classes
    ncid.createDimension('nageclass', H.nageclass)
    # dimension for release point characters
    ncid.createDimension('nchar', 45)
    # number of actual release points
    ncid.createDimension('numpoint', H.numpoint)

    # Create variables

    # time
    tID = ncid.createVariable('time', 'i4', ('time',),
                              zlib=True, complevel=complevel)
    tID.axis = 'T'
    tID.units = timeunit
    tID.calendar = 'proleptic_gregorian'

    # height
    levID = ncid.createVariable('height', 'f4', ('height',),
                                zlib=True, complevel=complevel)
    levID.axis = 'Z'
    levID.units = 'm'
    levID.positive = 'up'
    levID.standard_name = 'height'
    levID.long_name = 'grid cell upper height above ground'

    # lat
    latID = ncid.createVariable('latitude', 'f4', ('latitude',),
                                zlib=True, complevel=complevel)
    latID.long_name = 'latitude in degree north'
    latID.axis = 'Y'
    latID.units = 'degrees_north'
    latID.standard_name = 'grid_latitude'
    latID.description = 'grid cell center latitudes'

    # lon
    lonID = ncid.createVariable('longitude', 'f4', ('longitude',),
                                zlib=True, complevel=complevel)
    lonID.long_name = 'longitude in degree east'
    lonID.axis = 'X'
    lonID.units = 'degrees_east'
    lonID.standard_name = 'grid_longitude'
    lonID.description = 'grid cell center longitudes'

    # Mass and mass uncertainty outputs for each cell, time, and species
    gIDs = ('time', 'height', 'latitude', 'longitude')
    chunksizes = (1, nz, nlat, nlon)

    numID = ncid.createVariable('num_particles', 'i4', gIDs,
                                chunksizes=chunksizes,
                                zlib=True, complevel=complevel,
                                fill_value=0)
    numID.long_name = 'number of particles'

    mass_units = output_units(ncid)

    for i in range(0, H.nspec):
        anspec = "%3.3d" % (i + 1)

        var_name = "spec" + anspec + "_mass"
        massID = ncid.createVariable(var_name, 'f4', gIDs,
                                     chunksizes=chunksizes,
                                     zlib=True, complevel=complevel,
                                     fill_value=0.)
        massID.units = mass_units
        massID.long_name = 'mean particle mass of species ' + str(H.species[i])

        var_name = "spec" + anspec + "_mass_error"
        mass_error_ID = ncid.createVariable(var_name, 'f4', gIDs,
                                            chunksizes=chunksizes,
                                            zlib=True, complevel=complevel,
                                            fill_value=0.)
        mass_error_ID.units = mass_units
        mass_error_ID.long_name = 'standard error in the mean particle mass of species ' + str(H.species[i])

    if not min_size:
        # ORO
        oroID = ncid.createVariable('ORO', 'i4', ('latitude', 'longitude'),
                                    chunksizes=(H.numygrid, H.numxgrid),
                                    zlib=True, complevel=complevel)
        oroID.standard_name = 'surface altitude'
        oroID.long_name = 'outgrid surface altitude'
        oroID.units = 'm'

    if write_releases:
        # RELCOM
        relcomID = ncid.createVariable('RELCOM', 'S45', ('numpoint',),
                                       zlib=True, complevel=complevel)
        # Fill RELCOM with default values ("NA" means Non-Available)
        relcomID[:] = np.array(["NA"] * H.numpoint, dtype="S45")
        relcomID.long_name = 'release point name'

        # RELLNG1
        rellng1ID = ncid.createVariable('RELLNG1', 'f4', ('numpoint',),
                                        zlib=True, complevel=complevel)
        rellng1ID.units = 'degrees_east'
        rellng1ID.long_name = 'release longitude lower left corner'

        # RELLNG2
        rellng2ID = ncid.createVariable('RELLNG2', 'f4', ('numpoint',),
                                        zlib=True, complevel=complevel)
        rellng2ID.units = 'degrees_east'
        rellng2ID.long_name = 'release longitude upper right corner'

        # RELLAT1
        rellat1ID = ncid.createVariable('RELLAT1', 'f4', ('numpoint',),
                                        zlib=True, complevel=complevel)
        rellat1ID.units = 'degrees_north'
        rellat1ID.long_name = 'release latitude lower left corner'

        # RELLAT2
        rellat2ID = ncid.createVariable('RELLAT2', 'f4', ('numpoint',),
                                        zlib=True, complevel=complevel)
        rellat2ID.units = 'degrees_north'
        rellat2ID.long_name = 'release latitude upper right corner'

        # RELZZ1
        relzz1ID = ncid.createVariable('RELZZ1', 'f4', ('numpoint',),
                                       zlib=True, complevel=complevel)
        relzz1ID.units = 'm'
        relzz1ID.long_name = 'release height bottom'

        # RELZZ2
        relzz2ID = ncid.createVariable('RELZZ2', 'f4', ('numpoint',),
                                       zlib=True, complevel=complevel)
        relzz2ID.units = 'm'
        relzz2ID.long_name = 'release height top'

        # RELKINDZ
        relkindzID = ncid.createVariable('RELKINDZ', 'i4', ('numpoint',),
                                         zlib=True, complevel=complevel)
        relkindzID.long_name = 'release kind'

        # RELSTART
        relstartID = ncid.createVariable('RELSTART', 'i4', ('numpoint',),
                                         zlib=True, complevel=complevel)
        relstartID.units = 's'
        relstartID.long_name = 'release start relative to simulation start'

        # RELEND
        relendID = ncid.createVariable('RELEND', 'i4', ('numpoint',),
                                       zlib=True, complevel=complevel)
        relendID.units = 's'
        relendID.long_name = 'release end relative to simulation start'

        # RELPART
        relpartID = ncid.createVariable('RELPART', 'i4', ('numpoint',),
                                        zlib=True, complevel=complevel)
        relpartID.long_name = 'number of release particles'

        # RELXMASS
        relxmassID = ncid.createVariable('RELXMASS', 'f4',
                                         ('numspec', 'numpoint'),
                                         zlib=True, complevel=complevel)
        relxmassID.long_name = 'total release particles mass'

    # LAGE
    lageID = ncid.createVariable('LAGE', 'i4', ('nageclass',))
    lageID.units = 's'
    lageID.long_name = 'age class'

    return iout

def write_variables(H, ncid, write_releases, releases,
                    llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat, nlon, nlat, nz,
                    start_datetime=None, end_datetime=None, min_size=False):
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
    H.available_dates_dt = np.asarray(H.available_dates_dt)
    time_mask = (H.available_dates_dt>=start_datetime) & (H.available_dates_dt<=end_datetime)
    if np.count_nonzero(time_mask) == 0:
        raise RuntimeError("No available datetimes selected for output")
    if ncid.ipout == 1:
        time_var[:] = nc.date2num(H.available_dates_dt[time_mask], units=time_var.units,
                                  calendar=time_var.calendar)
        dates = np.asarray(H.available_dates)[time_mask]
    else:
        if H.available_dates_dt[-1] not in H.available_dates_dt[time_mask]:
            raise RuntimeError("IPOUT=%d, but end datetime not selected for output" % ncid.ipout)
        time_var[:] = nc.date2num(H.available_dates_dt[-1], units=time_var.units,
                                  calendar=time_var.calendar)
        dates = [H.available_dates[-1],]

    # levels
    if nz is None:
        ncid.variables['height'][:] = H.outheight
    else:
        ncid.variables['height'][:] = np.linspace(
            H.outheight[-1]/nz,
            H.outheight[-1],
            nz)
    lower_heights = np.append([0,], ncid.variables['height'][:-1])

    # latitudes (grid cell centers)
    if not llcrnrlat:
        llcrnrlat = ncid.outlat0
    if not urcrnrlat:
        urcrnrlat = llcrnrlat + H.numygrid * ncid.dyout
    if not nlat:
        nlat = H.numygrid
        dy = ncid.dyout
    else:
        dy = (urcrnrlat - llcrnrlat) / nlat
    ncid.variables['latitude'][:] = np.linspace(
        llcrnrlat + 0.5 * dy,
        llcrnrlat + (nlat-0.5) * dy,
        nlat)

    # longitudes (grid cell centers)
    if not llcrnrlon:
        llcrnrlon = ncid.outlon0
    if not urcrnrlon:
        urcrnrlon = llcrnrlon + H.numxgrid * ncid.dxout
    if not nlon:
        nlon = H.numxgrid
        dx = ncid.dxout
    else:
        dx = (urcrnrlon - llcrnrlon) / nlon
    ncid.variables['longitude'][:] = np.linspace(
        llcrnrlon + 0.5 * dx,
        llcrnrlon + (nlon-0.5) * dx,
        nlon)

    for idt, date in enumerate(dates):
        # read particle data
        try:
            H.read_particles(time_ret=idt, ipout=ncid.ipout)
            particles = H.particles[date]
        except IOError:
            # Oops, we have got an error while reading, so close the file
            ncid.close()
            # and re-raise the error
            raise

        numpart = particles['numpart']
        xlon = particles['xlon'][:].T
        ylat = particles['ylat'][:].T
        z = particles['z'][:].T

        for k, upper_cell_height in enumerate(ncid.variables['height']):
            for j, cell_lat in enumerate(ncid.variables['latitude']):
                for i, cell_lon in enumerate(ncid.variables['longitude']):
                    cell_particles_mask = ((xlon-(cell_lon-0.5*dx))*((cell_lon+0.5*dx)-xlon) >= 0) & ((ylat-(cell_lat-0.5*dy))*((cell_lat+0.5*dy)-ylat) >= 0) & ((z-lower_heights[k])*(upper_cell_height-z) >= 0)

                    n = ncid.variables['num_particles']
                    n[idt, k, j, i] = np.count_nonzero(cell_particles_mask)

                    for ispec in range(H.nspec):
                        anspec = "%3.3d" % (ispec + 1)

                        # Calculate values in each cell
                        mass_name = "spec" + anspec + "_mass"
                        mass = ncid.variables[mass_name]
                        mass_error_name = "spec" + anspec + "_mass_error"
                        mass_error = ncid.variables[mass_error_name]

                        xmass = particles['xmass'][:, ispec].T

                        if n[idt, k, j, i] > 0:
                            mass[idt, k, j, i] = np.mean(xmass[cell_particles_mask])
                        if n[idt, k, j, i] > 1:
                            mass_error[idt, k, j, i] = np.std(xmass[cell_particles_mask], ddof=1)/np.sqrt(n[idt, k, j, i])

    if not min_size:
        # Orography
        ncid.variables['ORO'][:, :] = H.oro

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
                  write_releases=True, dirout=None, outfile=None,
                  llcrnrlon=None, llcrnrlat=None, urcrnrlon=None, urcrnrlat=None,
                  nlon=None, nlat=None, nz=None,
                  start_datetime=None, end_datetime=None,
                  min_size=False, complevel=DEFAULT_COMPLEVEL):
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
            fprefix = 'partposit_grid_mass_'
        else:
            fprefix = 'partposit_end_grid_mass_'
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

    if not llcrnrlon:
        llcrnrlon = H.outlon0
    if not urcrnrlon:
        urcrnrlon = llcrnrlon + H.numxgrid * H.dxout
    if not llcrnrlat:
        llcrnrlat = H.outlat0
    if not urcrnrlat:
        urcrnrlat = llcrnrlat + H.numygrid * H.dyout
    if not nlon:
        nlon = H.numxgrid
    if not nlat:
        nlat = H.numygrid
    if not nz:
        nz = H.numzgrid

    cache_size = 16 * nlon * nlat * nz

    print("About to create new netCDF4 file: '%s'" % ncfname)
    ncid = nc.Dataset(ncfname, 'w', chunk_cache=cache_size)
    write_metadata(H, command, ncid)
    write_header(H, ncid, write_releases, species, nlon, nlat, nz, min_size,
                 complevel)
    write_variables(H, ncid, write_releases, releases, llcrnrlon, llcrnrlat,
                    urcrnrlon, urcrnrlat, nlon, nlat, nz, start_datetime,
                    end_datetime, min_size)
    ncid.close()
    return (ncfname, options_dir, output_dir)


def main():
    """Parses the passed command line arguments.

    The passed arguments will be used to create the netCDF4 file.
    """
    import argparse

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
        "--complevel", type=int, default=DEFAULT_COMPLEVEL,
        help="Compression level for the netCDF4 file."
    )
    parser.add_argument(
        '-llcrnrlon', type=float,
        help="Lower-left longitude"
    )
    parser.add_argument(
        '-llcrnrlat', type=float,
        help="Lower-left latitude"
    )
    parser.add_argument(
        '-urcrnrlon', type=float,
        help="Upper-right longitude"
    )
    parser.add_argument(
        '-urcrnrlat', type=float,
        help="Upper-right latitude"
    )
    parser.add_argument(
        '-nlon', '--number-longitudes', dest='nlon', type=int,
        help="Number of cells in longitude  direction"
    )
    parser.add_argument(
        '-nlat', '--number-latitudes', dest='nlat', type=int,
        help="Number of cells in latitude  direction"
    )
    parser.add_argument(
        '-nz', '--number-heights', dest='nz', type=int,
        help="Number of cells in vertical  direction"
    )
    parser.add_argument(
        '--start-datetime', nargs='?', default='0001-1-1 00:00:00',
        help="Minimum datetime for which to process particle masses (Format 'Y-m-d [H:M:S]')",
    )
    parser.add_argument(
        '--end-datetime', nargs='?', default='9999-1-1 00:00:00',
        help="Maximum datetime for which to process particle masses (Format 'Y-m-d [H:M:S]')",
    )
    parser.add_argument(
        'pathnames',
        help="The Flexpart pathnames file stating where options and output are."
    )

    args = parser.parse_args()

    try:
        args.start_datetime = datetime.strptime(args.start_datetime, DATETIME_FORMAT_STR)
    except ValueError:
        args.start_datetime = datetime.strptime(args.start_datetime, DATE_FORMAT_STR)
    try:
        args.end_datetime = datetime.strptime(args.end_datetime, DATETIME_FORMAT_STR)
    except ValueError:
        args.end_datetime = datetime.strptime(args.end_datetime, DATE_FORMAT_STR)

    ncfname, options_dir, output_dir = create_ncfile(
        args.pathnames, args.command_path, args.releases_path,
        not args.dont_write_releases, args.dirout, args.outfile,
        args.llcrnrlon, args.llcrnrlat, args.urcrnrlon, args.urcrnrlat,
        args.nlon, args.nlat, args.nz, args.start_datetime, args.end_datetime,
        args.min_size, args.complevel)

    print("New netCDF4 file is available in: '%s'" % ncfname)

if __name__ == '__main__':
    main()
