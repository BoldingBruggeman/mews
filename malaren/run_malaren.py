#!/usr/bin/env python

import datetime
from pathlib import Path
from typing import Optional
import numpy as np

import cftime

import pygetm

setup = "malaren"
nz = 20
ddu = 0.75
ddl = 0.75
Dgamma = 10.0
timestep = 5.0
use_adaptive = False


def create_domain(
    runtype: int,
    rivers: bool,
    **kwargs,
):
    import netCDF4
    import glob
    import os

    with netCDF4.Dataset(args.bathymetry_file) as nc:
        nc.set_auto_mask(False)
        domain = pygetm.domain.create_cartesian(
            nc["x"][:],
            nc["y"][:],
            lon=nc["lon"],
            lat=nc["lat"],
            H=nc[args.bathymetry_name][:, :],
            mask=np.where(nc[args.bathymetry_name][...] == -9999.0, 0, 1),
            z0=0.01,
        )
    domain.mask_indices(261, 270 + 1, 43, 43 + 1)
    domain.mask_indices(334, 335 + 1, 49, 63 + 1)
    domain.mask_indices(335, 335 + 1, 49, 63 + 1)
    domain.mask_indices(304, 304 + 1, 38, 44 + 1)

    domain.limit_velocity_depth()
    domain.cfl_check()

    if rivers:
        river_list = []
        # Inflows
        for river in glob.glob("Rivers/inflow_q*.nc"):
            #name = river.replace("Rivers/inflow_q_", "")
            name = river.split('_')[-1][:-3]
            with netCDF4.Dataset(river) as r:
                lon = r["lon"][:]
                lat = r["lat"][:]
                river_list.append(
                    domain.rivers.add_by_location(
                        name,
                        float(lon),
                        float(lat),
                        coordinate_type=pygetm.CoordinateType.LONLAT,
                    )
                )
        
        # Outflows
        for river in glob.glob("Rivers/outflow*_q.nc"):
            #name = river.replace("Rivers/inflow_q_", "")
            name = river.split('\\')[1]
            name = name.replace("_q.nc", "")
            with netCDF4.Dataset(river) as r:
                lon = r["lon"][:]
                lat = r["lat"][:]
                river_list.append(
                    domain.rivers.add_by_location(
                        name, float(lon), float(lat), spherical=True
                    )
                )

    return domain


def create_simulation(
    domain: pygetm.domain.Domain,
    runtype: pygetm.RunType,
    **kwargs,
) -> pygetm.simulation.Simulation:
    global use_adaptive
    if False:
        internal_pressure_method = pygetm.internal_pressure.BlumbergMellor()
    else:
        internal_pressure = pygetm.internal_pressure.ShchepetkinMcwilliams()

    if True:
        vertical_coordinates = pygetm.vertical_coordinates.GVC(
            nz, ddl=ddl, ddu=ddu, Dgamma=Dgamma, gamma_surf=True
        )
    elif False:
        try:
            use_adaptive = True
            vertical_coordinates = pygetm.vertical_coordinates.Adaptive(
                nz,
                timestep,
                cnpar=1.0,
                ddu=ddu,
                ddl=ddl,
                gamma_surf=True,
                Dgamma=Dgamma,
                csigma=0.001,
                cgvc=-0.001,
                hpow=3,
                chsurf=-0.001,
                hsurf=1.5,
                chmidd=-0.1,
                hmidd=0.5,
                chbott=-0.001,
                hbott=1.5,
                cneigh=-0.1,
                rneigh=0.25,
                decay=2.0 / 3.0,
                # cNN=1.0,
                cNN=0.1,
                drho=0.2,
                cSS=-1.0,
                dvel=0.1,
                chmin=0.1,
                hmin=0.5,
                nvfilter=1,
                vfilter=0.1,
                nhfilter=1,
                hfilter=0.2,
                split=1,
                timescale=3.0 * 3600.0,
            )
        except:
            print("Error: can not initialize Adaptive-coordinates")
            quit()
    else:
        vertical_coordinates = pygetm.vertical_coordinates.Sigma(nz, ddl=ddl, ddu=ddu)

    final_kwargs = dict(
        advection_scheme=pygetm.AdvectionScheme.SUPERBEE,
        # gotm=os.path.join(setup_dir, "gotmturb.nml"),
        # airsea=airsea,
        internal_pressure=internal_pressure,
        vertical_coordinates=vertical_coordinates,
        delay_slow_ip=True,
    )
    final_kwargs.update(kwargs)
    sim = pygetm.Simulation(domain, runtype=runtype, fabm = "fabm-selma.yaml", **final_kwargs)

    if sim.runtype < pygetm.RunType.BAROCLINIC:
        sim.sst = sim.airsea.t2m
    if sim.runtype == pygetm.RunType.BAROCLINIC:
        sim.radiation.set_jerlov_type(pygetm.Jerlov.Type_II)

    if args.initial and sim.runtype == pygetm.RunType.BAROCLINIC:
        if True:
            print("egon")
            # river["salt"].set(0.1)
            # river["temp"].set(0.5)
        else:
            print("Read froom files")
            # sim.salt.set(
            # pygetm.input.from_nc(
            #    os.path.join(args.setup_dir, "Input/initialConditions.nc"), "salt"
            # ),
            # on_grid=True,
        # )
        sim.temp.set(2)
        sim.salt.set(0.1)
        sim.density.convert_ts(sim.salt, sim.temp)
        sim.temp[..., domain.T.mask == 0] = pygetm.constants.FILL_VALUE
        sim.salt[..., domain.T.mask == 0] = pygetm.constants.FILL_VALUE

    ERA_path = "ERA5/era5_????.nc"
    sim.airsea.u10.set(pygetm.input.from_nc(ERA_path, "u10"))
    sim.airsea.v10.set(pygetm.input.from_nc(ERA_path, "v10"))
    sim.airsea.t2m.set(pygetm.input.from_nc(ERA_path, "t2m") - 273.15)
    sim.airsea.d2m.set(pygetm.input.from_nc(ERA_path, "d2m") - 273.15)
    sim.airsea.sp.set(pygetm.input.from_nc(ERA_path, "sp"))
    sim.airsea.tcc.set(pygetm.input.from_nc(ERA_path, "tcc"))
    ERA_path = "ERA5/precip_????.nc"
    sim.airsea.tp.set(pygetm.input.from_nc(ERA_path, "tp") / 3600.0)
    
    
    for river in sim.domain.rivers.values():
        if "outflow" in river.name:
            # Outflow
            river.flow.set(pygetm.input.from_nc(f"Rivers/{river.name}_q.nc", river.name))
        else:
            # Inflow
            river.flow.set(pygetm.input.from_nc(f"Rivers/inflow_q_{river.name}.nc", river.name)) # .mean() ## Error here if I remove .mean()
            river["pfas_c"].set(1.0)
    
    # sim["age_age_of_water"].river_follow[:] = False # By default, precipitation also has age 0
    
    return sim


def create_output(
    output_dir: str,
    sim: pygetm.simulation.Simulation,
    **kwargs,
):
    sim.logger.info("Setting up output")

    path = Path(output_dir, "meteo.nc")
    output = sim.output_manager.add_netcdf_file(
        str(path),
        interval=datetime.timedelta(days = 7),
        sync_interval=None,
    )
    output.request(
        "u10",
        "v10",
        "sp",
        "t2m",
        "tcc",
        "tp",
    )

    path = Path(output_dir, setup + "_2d.nc")
    output = sim.output_manager.add_netcdf_file(
        str(path),
        interval=datetime.timedelta(days=7),
        sync_interval=None,
    )
    output.request("Ht", "zt", "u1", "v1", "tausxu", "tausyv")
    if args.debug_output:
        output.request("maskt", "masku", "maskv")
        output.request("U", "V")
        # output.request("Du", "Dv", "dpdx", "dpdy", "z0bu", "z0bv", "z0bt")
        # output.request("ru", "rru", "rv", "rrv")

    if sim.runtype > pygetm.RunType.BAROTROPIC_2D:
        path = Path(output_dir, setup + "_3d.nc")
        output = sim.output_manager.add_netcdf_file(
            str(path),
            interval=datetime.timedelta(days=7),
            sync_interval=None,
        )
    output.request("Ht", "uk", "vk", "ww", "SS", "num")
    if args.debug_output:
        output.request("fpk", "fqk", "advpk", "advqk")  # 'diffpk', 'diffqk')

    if sim.runtype == pygetm.RunType.BAROCLINIC:
        output.request("temp", "salt", "rho", "NN", "rad", "sst", "hnt", "nuh")
        if args.debug_output:
            output.request("idpdx", "idpdy")
        if use_adaptive:
            output.request("nug", "ga", "dga")

    if sim.fabm:
        output.request("pfas_c", "selma_po", "total_chlorophyll_calculator_result") # , "age_age_of_water")


def run(
    sim: pygetm.simulation.Simulation,
    start: cftime.datetime,
    stop: cftime.datetime,
    dryrun: bool = False,
    **kwargs,
):
    if dryrun:
        print(f"")
        print(f"Making a dryrun - skipping sim.advance()")
        print(f"")
    else:
        sim.start(
            simstart,
            timestep=timestep,
            split_factor=20,
            **kwargs,
        )
        while sim.time < simstop:
            sim.advance()
        sim.finish()


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("start", help="simulation start time - yyyy-mm-dd hh:mi:ss")
    parser.add_argument("stop", help="simulation stop time - yyyy-mm-dd hh:mi:ss")
    parser.add_argument(
        "--setup_dir",
        type=Path,
        help="Path to configuration files - not used yet",
        default=".",
    )

    parser.add_argument(
        "--bathymetry_file",
        type=str,
        help="Name of bathymetry file",
        default="Bathymetry/bathymetry.nc",
    )

    parser.add_argument(
        "--bathymetry_name",
        type=str,
        help="Name of bathymetry variable",
        default="bathymetry",
    )

    parser.add_argument(
        "--output_dir", type=str, help="Path to save output files", default="."
    )

    parser.add_argument(
        "--initial",
        action="store_true",
        help="Initial run salinity and temerature are specified",
    )
    parser.add_argument(
        "--runtype",
        type=int,
        choices=(pygetm.BAROTROPIC_2D, pygetm.BAROTROPIC_3D, pygetm.BAROCLINIC),
        help="Run type",
        default=pygetm.BAROCLINIC,
    )
    parser.add_argument(
        "--no_rivers", action="store_false", dest="rivers", help="No river input"
    )
    parser.add_argument(
        "--no_output",
        action="store_false",
        dest="output",
        help="Do not save any results to NetCDF",
    )
    parser.add_argument(
        "--debug_output",
        action="store_true",
        help="Save additional variables for debugging",
    )
    parser.add_argument("--save_restart", help="File to save restart to")
    parser.add_argument("--load_restart", help="File to load restart from")
    parser.add_argument("--profile", help="File to save profiling report to")
    parser.add_argument("--dryrun", action="store_true", help="Do a dry run")
    parser.add_argument(
        "--plot_domain", action="store_true", help="Plot the calculation domain"
    )
    
    # # Testing only ->
    # import sys
    # sys.argv.extend(["start 1995-02-01 00:00:00",
    #                   "stop 1995-03-01 00:00:00"])
    # # Empty sys.argv ->
    # sys.argv = [sys.argv[0]]
    
    args = parser.parse_args()
    
    # # Testing only
    # args.load_restart = "restart_malaren_19950301.nc"
    # args.save_restart = "restart_malaren_19950401.nc"
    # args.output_dir = "19950301"
    
    # Debugging: %debug line_to_debug
    
    if args.output_dir != ".":
        p = Path(args.output_dir)
        if not p.is_dir():
            print(f"Folder {args.output_dir} does not exist - create and run again")
            exit()

    domain = create_domain(args.runtype, args.rivers)

    sim = create_simulation(domain, args.runtype)

    # for plot options see:
    # https://github.com/BoldingBruggeman/getm-rewrite/blob/fea843cbc78bd7d166bdc5ec71c8d3e3ed080a35/python/pygetm/domain.py#L1943
    if args.plot_domain:
        f = domain.plot(show_mesh=False, show_subdomains=False)
        if f is not None:
            f.savefig("domain_mesh.png")
        f = domain.plot(show_mesh=False, show_mask=True)
        if f is not None:
            f.savefig("domain_mask.png")

    if args.output and not args.dryrun:
        create_output(args.output_dir, sim)

    if args.save_restart and not args.dryrun:
        sim.output_manager.add_restart(args.save_restart)

    if args.load_restart and not args.dryrun:
        simstart = sim.load_restart(args.load_restart)

    simstart = datetime.datetime.strptime(args.start, "%Y-%m-%d %H:%M:%S")
    simstop = datetime.datetime.strptime(args.stop, "%Y-%m-%d %H:%M:%S")
    profile = setup if args.profile is not None else None
    run(
        sim,
        simstart,
        simstop,
        dryrun=args.dryrun,
        report=datetime.timedelta(days=7),
        report_totals=datetime.timedelta(days=7),
        profile=profile,
    )
