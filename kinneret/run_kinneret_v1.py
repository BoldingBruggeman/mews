#!/usr/bin/env python
# run example: mpiexec -n 20 python run_kinneret_v0.py "2022-02-01 12:00:00" "2022-02-02 12:00:00" --initial
import datetime
from pathlib import Path
from typing import Optional
import numpy as np

import cftime

import pygetm
import yaml

setup = "kinneret"
nz = 30
ddu = 0.75
ddl = 0.5
Dgamma = 10.0
timestep = 5.0
use_adaptive = False
use_tracers = False

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
            nc["X"][:],
            nc["Y"][:],
            lon=nc["longitude"][:],
            lat=nc["latitude"][:],
            H=-1.0 * nc["positive_rx0"][:, :], 
            mask=np.where(nc["positive_rx0"][...] == 9.96921e+36, 0, 1),
            #H=-1.0 * nc[args.bathymetry_name][:, :],
            #mask=np.where(nc[args.bathymetry_name][...] == 1.0e37, 0, 1),
            z0=0.01,
            #            **final_kwargs,
        )

    domain.limit_velocity_depth()
    domain.cfl_check()
    domain.smooth(rx0=0.08)
    domain.mask_shallow(1.0)

    if rivers:
        river_list = []
        for river in glob.glob(os.path.join("Rivers/Rivers2010_2020", "inflow_q*.nc")):
            name = os.path.basename(river)
            name = name.replace("inflow_q_", "").replace(".nc", "")
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
    # Tracers
        if use_tracers and os.path.isdir("Tracers"):
            for tracer_file in glob.glob("Tracers/*.nc"):
                name = os.path.basename(tracer_file)
                name = name.replace("Tracer_file_", "nctracer_").replace(".nc", "")
                with netCDF4.Dataset(tracer_file) as r:
                    lon = r["lon"][:]
                    lat = r["lat"][:]
                    river_list.append(
                        domain.rivers.add_by_location(
                            name, float(lon), float(lat), coordinate_type=pygetm.CoordinateType.LONLAT
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
    elif True:
        try:
            use_adaptive = True
            vertical_coordinates = pygetm.vertical_coordinates.Adaptive(
                nz,
                timestep,
                cnpar=5.0,
                ddu=ddu,
                ddl=ddl,
                gamma_surf=True,
                Dgamma=Dgamma,
                csigma=0.001,
                cgvc=-0.001,
                hpow=3,
                chsurf=0.001,
                hsurf=1.0,
                chmidd=-0.1,
                hmidd=0.5,
                chbott=-0.001,
                hbott=1.5,
                cneigh=0.001,
                rneigh=0.25,
                decay=2.0 / 3.0,
                # cNN=1.0,
                cNN=0.05,
                drho=0.5,
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
    calculate_heat_flux=True,   
    airsea = pygetm.airsea.FluxesFromMeteo(
            shortwave_method=pygetm.DOWNWARD_FLUX,
            calculate_evaporation=True,
    )
        

    final_kwargs = dict(
        advection_scheme=pygetm.AdvectionScheme.SUPERBEE,
        # gotm=os.path.join(setup_dir, "gotmturb.nml"),
        gotm="gotm.yaml",
        airsea=airsea,
        internal_pressure=internal_pressure,
        vertical_coordinates=vertical_coordinates,
        delay_slow_ip=False,
    )
    final_kwargs.update(kwargs)
    sim = pygetm.Simulation(
        domain,
        runtype=runtype,
        fabm="fabm.yaml",
        **final_kwargs,
    )

    if sim.runtype < pygetm.RunType.BAROCLINIC:
        sim.sst = sim.airsea.t2m
    if sim.runtype == pygetm.RunType.BAROCLINIC:
        #sim.radiation.set_jerlov_type(pygetm.Jerlov.Type_II)
        sim.radiation.A.set(0.7474)
        sim.radiation.kc1.set(0.4659) #1/g1 in gotm from shjar
        sim.radiation.kc2.set(0.2039)
     
    
    if not args.load_restart and sim.runtype == pygetm.RunType.BAROCLINIC:
        if True:
            sim.temp.set(16.0)
            sim.salt.set(0.4)
        else:
            sim.salt.set(0.4)
            # sim.salt.set(
            #    pygetm.input.from_nc(
            #        os.path.join(args.setup_dir, "Input/initial_TS.nc"),
            #        "S",
            #    ),
            #    on_grid=True,
            # )
            # sim.salt[...] = np.flip(sim.salt[...], axis=0)
            sim.temp.set(
                pygetm.input.from_nc(
                    str(Path(args.setup_dir, "Input/initial_TS.nc")),
                    "Temp",
                ),
                on_grid=True,
            )
            sim.temp[...] = np.flip(sim.temp[...], axis=0)
        sim.temp[..., sim.T.mask == 0] = pygetm.constants.FILL_VALUE
        sim.salt[..., sim.T.mask == 0] = pygetm.constants.FILL_VALUE
        sim.density.convert_ts(sim.salt, sim.temp)

    # sim["diatoms_c"] set(pygetm.input.from_nc("some.nc","dia")) example of initial conditions of diatoms
    ERA_path = "ERA5/era5_????.nc"
    sim.airsea.u10.set(pygetm.input.from_nc(ERA_path, "u10") * 2 )
    sim.airsea.v10.set(pygetm.input.from_nc(ERA_path, "v10") * 2)
    sim.airsea.t2m.set(pygetm.input.from_nc(ERA_path, "t2m") - 273.15)
    sim.airsea.d2m.set(pygetm.input.from_nc(ERA_path, "d2m") - 273.15)
    sim.airsea.sp.set(pygetm.input.from_nc(ERA_path, "sp"))
    sim.airsea.tcc.set(pygetm.input.from_nc(ERA_path, "tcc"))
    sim.airsea.tp.set(pygetm.input.from_nc(ERA_path, "tp") / 3600.0)
    ERA_path_ssrd = "ERA5_ssrd/era5_ssrd_????.nc"
    sim.airsea.swr_downwards.set(pygetm.input.from_nc(ERA_path_ssrd, "ssrd") * (0.85/3600.0))
    
    for river in sim.rivers.values():
        river.flow.set(pygetm.input.from_nc(f"Rivers/Claude_Rivers2010_2020/inflow_q_{river.name}.nc", "Flow"))
        
        name_lower = river.name.lower()
        is_withdrawal = any(tag in name_lower for tag in ["withdrawal", "out", "pump"])

        if is_withdrawal:
           print(f"[INFO] '{river.name}' identified as a withdrawal — skipping Temp and Salt.")
           continue
        
        
        river["temp"].set(
            pygetm.input.from_nc(f"Rivers/Claude_Rivers2010_2020/inflow_q_{river.name}.nc", "Temp")
        )
        river["salt"].set(
            pygetm.input.from_nc(f"Rivers/Claude_Rivers2010_2020/inflow_q_{river.name}.nc", "Salt")
        )

    _flux_path = "Input/dust.nc"
    sim.fabm.get_dependency("ammonium_dust_flux/flux").set(
        pygetm.input.from_nc(str(_flux_path),"aa_dust_flux")) #, preprocess=_add_coord
    print(sim.fabm.has_dependency("ammonium_dust_flux/flux"))
    sim.fabm.get_dependency("nitrate_dust_flux/flux").set(
        pygetm.input.from_nc(str(_flux_path),"nn_dust_flux"))
    sim.fabm.get_dependency("phosphate_dust_flux/flux").set(
        pygetm.input.from_nc(str(_flux_path),"po_dust_flux"))
    sim.fabm.get_dependency("detritus_c_dust_flux/flux").set(
        pygetm.input.from_nc(str(_flux_path),"dd_n_dust_flux")) #must fix in dust.nc
    sim.fabm.get_dependency("detritus_p_dust_flux/flux").set(
        pygetm.input.from_nc(str(_flux_path),"dd_p_dust_flux"))
    sim.fabm.get_dependency("detritus_n_dust_flux/flux").set(
        pygetm.input.from_nc(str(_flux_path),"dd_n_dust_flux"))
    sim.fabm.get_dependency("PFe_dust_flux/flux").set(
        pygetm.input.from_nc(str(_flux_path),"PFe_dust_flux"))

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
        interval=datetime.timedelta(hours=24),
        sync_interval=None,
    )
    output.request(
        "u10",
        "v10",
        "sp",
        "t2m",
        "tcc",
        "tp",
        "swr",
        "pe",
    )

    path = Path(output_dir, setup + "_2d.nc")
    output = sim.output_manager.add_netcdf_file(
        str(path),
        interval=datetime.timedelta(days=1),
        sync_interval=None,
    )
    output.request("Ht", "zt", "Dt", "u1", "v1", "tausxu", "tausyv", "pe")
    if args.debug_output:
        output.request("maskt", "masku", "maskv")
        output.request("U", "V")
        # output.request("Du", "Dv", "dpdx", "dpdy", "z0bu", "z0bv", "z0bt")
        # output.request("ru", "rru", "rv", "rrv")

    if sim.runtype > pygetm.RunType.BAROTROPIC_2D:
        path = Path(output_dir, setup + "_3d.nc")
        output = sim.output_manager.add_netcdf_file(
            str(path),
            interval=datetime.timedelta(days=1),
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
        with open('output.yaml', 'r') as f:
            config = yaml.safe_load(f)
        fabm_outputs = config.get('fabm_outputs', [])
        output.request(*fabm_outputs)

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
        default="Bathymetry/bathymetry_400m_smoothed.nc",
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
    args = parser.parse_args()

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
        report=datetime.timedelta(hours=24),
        report_totals=datetime.timedelta(days=7),
        profile=profile,
    )