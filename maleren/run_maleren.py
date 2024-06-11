#!/usr/bin/env python

import datetime
import os.path
import pathlib
from typing import Optional

import cftime

import pygetm


def create_domain(
    runtype: int,
    rivers: bool = False,
    **kwargs,
):
    import netCDF4
    import numpy as np
    import glob

    if runtype > pygetm.BAROTROPIC_2D:
        final_kwargs = dict(
            nz=20,
            vertical_coordinate_method=pygetm.VerticalCoordinates.GVC,
            Dgamma=10.0,
            ddu=0.75,
            ddl=0.5,
            Dmin=0.2,
            Dcrit=1.0,
        )
    else:
        final_kwargs = dict(
            nz=1,
        )

    final_kwargs.update(kwargs)

    with netCDF4.Dataset("Bathymetry/bathymetry.nc") as nc:
        nc.set_auto_mask(False)
        domain = pygetm.domain.create_cartesian(
            nc["x"][:],
            nc["y"][:],
            lon=nc["lon"],
            lat=nc["lat"],
            H=nc["bathymetry"][:, :],
            mask=np.where(nc["bathymetry"][...] == -9999.0, 0, 1),
            z0=0.01,
            **final_kwargs,
        )
    domain.mask_indices(261, 270 + 1, 43, 43 + 1)
    domain.mask_indices(334, 335 + 1, 49, 63 + 1)
    domain.mask_indices(335, 335 + 1, 49, 63 + 1)
    domain.mask_indices(304, 304 + 1, 38, 44 + 1)

    domain.limit_velocity_depth()
    domain.cfl_check()

    if rivers:
        river_list = []
        for river in glob.glob("Rivers/inflow_q*.nc"):
            name = river.replace("Rivers/inflow_q_", "")
            name = name.replace(".nc", "")
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
    runtype: int,
    **kwargs,
) -> pygetm.simulation.Simulation:
    final_kwargs = dict(
        advection_scheme=pygetm.AdvectionScheme.SUPERBEE,
        # gotm=os.path.join(setup_dir, "gotmturb.nml"),
        # airsea=airsea,
        internal_pressure_method=pygetm.InternalPressure.SHCHEPETKIN_MCWILLIAMS,
        delay_slow_ip=True,
    )
    final_kwargs.update(kwargs)
    sim = pygetm.Simulation(domain, runtype=runtype, **final_kwargs)

    if sim.runtype < pygetm.BAROCLINIC:
        sim.sst = sim.airsea.t2m
    if sim.runtype == pygetm.BAROCLINIC:
        sim.radiation.set_jerlov_type(pygetm.Jerlov.Type_II)

    if args.initial and sim.runtype == pygetm.BAROCLINIC:
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

    return sim


def create_output(
    sim: pygetm.simulation.Simulation,
    **kwargs,
):
    sim.logger.info("Setting up output")
    output = sim.output_manager.add_netcdf_file(
        "meteo.nc", interval=datetime.timedelta(hours=1), sync_interval=None
    )
    output.request(
        "u10",
        "v10",
        "sp",
        "t2m",
        "tcc",
        "tp",
    )

    output = sim.output_manager.add_netcdf_file(
        "maleren_2d.nc", interval=datetime.timedelta(hours=1), sync_interval=None
    )
    output.request("zt", "u1", "v1", "tausxu", "tausyv")
    if args.debug_output:
        output.request("maskt", "masku", "maskv")
        output.request("U", "V")
        # output.request("Du", "Dv", "dpdx", "dpdy", "z0bu", "z0bv", "z0bt")
        # output.request("ru", "rru", "rv", "rrv")

    if sim.runtype > pygetm.BAROTROPIC_2D:
        output = sim.output_manager.add_netcdf_file(
            "maleren_3d.nc", interval=datetime.timedelta(hours=6), sync_interval=None
        )
    output.request("uk", "vk", "ww", "SS", "num")
    if args.debug_output:
        output.request("fpk", "fqk", "advpk", "advqk")  # 'diffpk', 'diffqk')

    if sim.runtype == pygetm.BAROCLINIC:
        output.request("temp", "salt", "rho", "NN", "rad", "sst", "hnt", "nuh")
        if args.debug_output:
            output.request("idpdx", "idpdy")

    if sim.fabm:
        output.request("par", "med_ergom_o2", "med_ergom_OFL", "med_ergom_dd")


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
            timestep=5.0,
            split_factor=20,
            **kwargs,
        )
        while sim.time < simstop:
            sim.advance()
        sim.finish()


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("start", help="simulation start time - yyyy-mm-dd hh:mi:ss")
    parser.add_argument("stop", help="simulation stop time - yyyy-mm-dd hh:mi:ss")
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
    args = parser.parse_args()

    domain = create_domain(args.runtype)
    if args.plot_domain:
        f = domain.plot()
        if f is not None:
            f.savefig("domain_mesh.png")
        f = domain.plot(show_mesh=False, show_mask=True)
        if f is not None:
            f.savefig("domain_mask.png")

    sim = create_simulation(domain, args.runtype)

    if args.output and not args.dryrun:
        create_output(sim)

    if args.save_restart and not args.dryrun:
        sim.output_manager.add_restart(args.save_restart)

    if args.load_restart and not args.dryrun:
        simstart = sim.load_restart(args.load_restart)

    simstart = datetime.datetime.strptime(args.start, "%Y-%m-%d %H:%M:%S")
    simstop = datetime.datetime.strptime(args.stop, "%Y-%m-%d %H:%M:%S")
    profile = "maleren" if args.profile is not None else None
    run(
        sim,
        simstart,
        simstop,
        dryrun=args.dryrun,
        report=datetime.timedelta(hours=1),
        report_totals=datetime.timedelta(days=1),
        profile=profile,
    )
