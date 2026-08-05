#!/usr/bin/env python



import datetime
from pathlib import Path
from typing import Optional
import numpy as np
import xarray as xr

import cftime

import pygetm



setup = "ohra"
nz = 20
ddu = 1.5
ddl = 1.5
Dgamma = 15.0
timestep = 0.8
use_adaptive = False
light_A = 0.55
light_kc1 = 2
light_kc2 = 0.15

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
    domain.limit_velocity_depth()
    domain.cfl_check()
    domain.mask_shallow(1.0)
    if rivers:
        river_list = []
        # Inflows
        for river in glob.glob("Rivers/Inflow_file*.nc"):
            name = os.path.basename(river)
            name = name.replace("Inflow_file_", "").replace(".nc", "")
            with netCDF4.Dataset(river) as r:
                lon = r["lon"][:]
                lat = r["lat"][:]
                river_list.append(
                    domain.rivers.add_by_location(
                        name,
                        float(lon),
                        float(lat),
                        coordinate_type=pygetm.CoordinateType.LONLAT
                    )
                )
        # Outflows
        for river in glob.glob("Rivers/Outflow_file*.nc"):
            name = os.path.basename(river)
            name = name.replace("Outflow_file_", "").replace(".nc", "")
            with netCDF4.Dataset(river) as r:
                lon = r["lon"][:]
                lat = r["lat"][:]
            river_list.append(
                domain.rivers.add_by_location(
                    name,
                    float(lon),
                    float(lat),
                    coordinate_type=pygetm.CoordinateType.LONLAT,
                    zu=40.17,   # upper depth limit
                    zl=41.17    # lower depth limit
                )
            )

    return domain

# Initialize airsea with DOWNWARD_FLUX for SSRD

def rx0_to_supergrid(domain, rx0_u, rx0_v, fill_value=np.nan):
    """
    Map rx0 arrays returned by domain.get_rx0() onto supergrid arrays so they can be plotted
    with domain.plot(field=...).
    """
    # Supergrid shape
    ny_sup, nx_sup = domain.H.shape  # (ny*2+1, nx*2+1)

    # Put rx0_u onto U points: [1::2, 2:-2:2]
    rx0_u_sup = np.full((ny_sup, nx_sup), fill_value, dtype=float)
    rx0_u_sup[1::2, 2:-2:2] = rx0_u

    # Put rx0_v onto V points: [2:-2:2, 1::2]
    rx0_v_sup = np.full((ny_sup, nx_sup), fill_value, dtype=float)
    rx0_v_sup[2:-2:2, 1::2] = rx0_v

    return rx0_u_sup, rx0_v_sup


def create_simulation(
    domain: pygetm.domain.Domain,
    runtype: pygetm.RunType,
    **kwargs,
) -> pygetm.simulation.Simulation:
    import numpy as np
    import pandas as pd
    import os

    global use_adaptive
    if True:
        internal_pressure = pygetm.internal_pressure.ShchepetkinMcwilliams()
    else:
        internal_pressure = pygetm.internal_pressure.BlumbergMellor()

    if True:
        airsea = pygetm.airsea.FluxesFromMeteo(shortwave_method=pygetm.DOWNWARD_FLUX,     calculate_evaporation=False )

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
        airsea=airsea,
        internal_pressure=internal_pressure,
        vertical_coordinates=vertical_coordinates,
    #    delay_slow_ip=True,
        Dmin = 0.2,
        Dcrit = 1,
    )
    final_kwargs.update(kwargs)
    sim = pygetm.Simulation(domain, runtype=runtype, gotm = "gotm.yaml", **final_kwargs)
  #  sim = pygetm.Simulation(domain, runtype=runtype, gotm = "gotm.yaml", fabm = "fabm.yaml", **final_kwargs)
    if sim.runtype < pygetm.RunType.BAROCLINIC:
        sim.sst = sim.airsea.t2m
    if sim.runtype == pygetm.RunType.BAROCLINIC:
     #   sim.radiation.set_jerlov_type(pygetm.Jerlov.Type_II)
        sim.radiation.A.fill(light_A)
        sim.radiation.kc1.fill(light_kc1)
        sim.radiation.kc2.fill(light_kc2)


    if not args.load_restart and sim.runtype == pygetm.RunType.BAROCLINIC:
        if True:
           sim.temp.set(4)
           sim.salt.set(0)
           sim["zt"].all_values[:] = -2.64
           sim["zint"].all_values[:] = -2.64
           sim["ziot"].all_values[:] = -2.64
           sim["zot"].all_values[:] = -2.64

          # sim.zt.all_values = 
           #sim["zt"].all_values = -7 you might have to do the same with zint, ziot, zot as well - to the same value - but do test that , in 2018 it is -2.64

        else:
            print("Read froom files")
            # sim.salt.set(
            # pygetm.input.from_nc(
            #    os.path.join(args.setup_dir, "Input/initialConditions.nc"), "salt"
            # ),
            # on_grid=True,
        # )
        sim.density.convert_ts(sim.salt, sim.temp)
        sim.temp[..., sim.T.mask == 0] = pygetm.constants.FILL_VALUE
        sim.salt[..., sim.T.mask == 0] = pygetm.constants.FILL_VALUE

    met_path = "ERA5/era5_????.nc"
    


    sim.airsea.u10.set(pygetm.input.from_nc(met_path, "u10")*1)
    sim.airsea.v10.set(pygetm.input.from_nc(met_path, "v10")*1)
    sim.airsea.t2m.set(pygetm.input.from_nc(met_path, "t2m") - 273.15)
    sim.airsea.d2m.set(pygetm.input.from_nc(met_path, "d2m") - 273.15)
    sim.airsea.sp.set(pygetm.input.from_nc(met_path, "sp"))
    sim.airsea.tcc.set(pygetm.input.from_nc(met_path, "tcc"))
    sim.airsea.tp.set(pygetm.input.from_nc(met_path, "tp") * 0)
    if sim.airsea.shortwave_method == pygetm.DOWNWARD_FLUX:
        sim.airsea.swr_downwards.set(
            pygetm.input.from_nc(met_path, "ssrd")* ( 0.75/3600.0 )
        )
    for river in sim.rivers.values():
        if "outflow" in river.name:
           ### Outflow
            river.flow.set(pygetm.input.from_nc(f"Rivers/Outflow_file_{river.name}.nc", "flow"))
        else:
          ### Inflow
           river.flow.set(pygetm.input.from_nc(f"Rivers/Inflow_file_{river.name}.nc", "flow"))         
           river["temp"].follow_target_cell = False #(True makes it use the value from the basin)
           river["temp"].set(pygetm.input.from_nc(f"Rivers/Inflow_file_{river.name}.nc", "temp"))
           # Nutrients
         #  river["selmaprotbas_po"].follow_target_cell = False #(True makes it use the value from the basin)
         #  river["selmaprotbas_po"].set(pygetm.input.from_nc(f"Rivers/Inflow_file_{river.name}.nc", "selmaprotbas_po"))  
         #  river["selmaprotbas_nn"].follow_target_cell = False
         #  river["selmaprotbas_nn"].set(pygetm.input.from_nc(f"Rivers/Inflow_file_{river.name}.nc", "selmaprotbas_nn"))
         #  river["selmaprotbas_dd_p"].follow_target_cell = False
         #  river["selmaprotbas_dd_p"].set(pygetm.input.from_nc(f"Rivers/Inflow_file_{river.name}.nc", "selmaprotbas_ddp"))
         #  river["selmaprotbas_dd_n"].follow_target_cell = False
         #  river["selmaprotbas_dd_n"].set(pygetm.input.from_nc(f"Rivers/Inflow_file_{river.name}.nc", "selmaprotbas_ddn"))
         #  river["selmaprotbas_dd_c"].follow_target_cell = False
         #  river["selmaprotbas_dd_c"].set(pygetm.input.from_nc(f"Rivers/Inflow_file_{river.name}.nc", "selmaprotbas_ddc"))
         #  river["selmaprotbas_pw"].follow_target_cell = False
         #  river["selmaprotbas_pw"].set(pygetm.input.from_nc(f"Rivers/Inflow_file_{river.name}.nc", "selmaprotbas_pw"))
         #  river["selmaprotbas_o2"].follow_target_cell = False
         #  river["selmaprotbas_o2"].set(pygetm.input.from_nc(f"Rivers/Inflow_file_{river.name}.nc", "selmaprotbas_o2"))
         #  river["selmaprotbas_pw"].follow_target_cell = False
         #  river["selmaprotbas_aa"].set(0.0)
         #  river["selmaprotbas_pw"].follow_target_cell = False
         #  river["selmaprotbas_si"].set(200)
         #  river["selmaprotbas_pw"].follow_target_cell = False
         #  river["selmaprotbas_dd_si"].set(200)


    
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
        interval=datetime.timedelta(hours=240),
        sync_interval=None,
    )
    output.request(
        "u10",
        "v10",
        "sp",
        "t2m",
        "tcc",
        "tp",
  #      "ssrd",
    )

    path = Path(output_dir, setup + "_2d.nc")
    output = sim.output_manager.add_netcdf_file(
        str(path),
        interval=datetime.timedelta(hours=6),
        sync_interval=None,
    )
    output.request("Ht", "zt","Dt", "u1", "v1", "tausxu", "tausyv")
    if args.debug_output:
        output.request("maskt", "masku", "maskv")
        output.request("U", "V")
        # output.request("Du", "Dv", "dpdx", "dpdy", "z0bu", "z0bv", "z0bt")
        # output.request("ru", "rru", "rv", "rrv")

    if sim.runtype > pygetm.RunType.BAROTROPIC_2D:
        path = Path(output_dir, setup + "_3d.nc")
        output = sim.output_manager.add_netcdf_file(
            str(path),
            interval=datetime.timedelta(hours=6),
            sync_interval=None,
        )
    output.request("Ht", "uk", "vk", "ww", "SS", "num")
    if args.debug_output:
        #output.request("fpk", "fqk", "advpk", "advqk")  # 'diffpk', 'diffqk')
        output.request("fpk", "fqk","advpk", "advqk")  # 'diffpk', 'diffqk')


    if sim.runtype == pygetm.RunType.BAROCLINIC:
        output.request("temp", "salt", "rho", "swr", "NN", "rad", "sst", "hnt", "nuh")
        if args.debug_output:
            output.request("idpdx", "idpdy")
        if use_adaptive:
            output.request("nug", "ga", "dga")

    if sim.fabm:
        output.request( "selmaprotbas_po", "total_chlorophyll_calculator_result", "selmaprotbas_o2")


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
        debug_time = datetime.datetime.strptime(
            "2018-10-01 00:00:00", "%Y-%m-%d %H:%M:%S"
        )
        while sim.time < simstop:
            if sim.time < debug_time:
                x= False
            else:
                x = True
            sim.advance(check_finite=x)
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
        default="Bathymetry/bathymetry_smoothed_fixed_rx_525_75_rx_025.nc",
     #   default="Bathymetry/bathymetry.nc",
    )

    parser.add_argument(
        "--bathymetry_name",
        type=str,
        help="Name of bathymetry variable",
        default="bathymetry_rx01_local_v",
     #   default="bathymetry",
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
#        domain = create_domain(args.runtype, args.rivers, H = -2)
    sim = create_simulation(domain, args.runtype)

    # for plot options see:
    # https://github.com/BoldingBruggeman/getm-rewrite/blob/fea843cbc78bd7d166bdc5ec71c8d3e3ed080a35/python/pygetm/domain.py#L1943
   # if args.plot_domain:
   #     f = domain.plot(show_mesh=False, show_subdomains=True, tiling=domain.create_tiling())
   #     if f is not None:
   #         f.savefig("domain_mesh.png")
   #     f = domain.plot(show_mesh=False, show_mask=True)
   #     if f is not None:
   #         f.savefig("domain_mask.png")
    def rx0_u_to_field(domain, rx0_u, fill=np.nan):
        f = np.full(domain.H.shape, fill, dtype=float)
        f[1::2, 2:-2:2] = rx0_u          # U points
        f[1::2, 1:-3:2] = rx0_u          # fill gaps next to U points (visual only)
        return f

    def rx0_v_to_field(domain, rx0_v, fill=np.nan):
        f = np.full(domain.H.shape, fill, dtype=float)
        f[2:-2:2, 1::2] = rx0_v          # V points
        f[1:-3:2, 1::2] = rx0_v          # fill gaps next to V points (visual only)
        return f

    if args.plot_domain:
        domain.plot(show_mesh=True, show_subdomains=True, tiling=domain.create_tiling(), show_rivers=False).savefig("domain_mesh.png")
        domain.plot(show_mesh=False, show_mask=True, show_subdomains=True, tiling=domain.create_tiling(), show_rivers=False).savefig("domain_mask.png")

        rx0_u, rx0_v = domain.get_rx0(zmin=0.0, Dmin=0.0)

        domain.plot(
            field=rx0_u_to_field(domain, rx0_u),
            show_bathymetry=False, show_rivers=False,
            label="rx0_u", cmap="viridis",
        ).savefig("rx0_u.png", dpi=200)

        domain.plot(
            field=rx0_v_to_field(domain, rx0_v),
            show_bathymetry=False, show_rivers=False,
            label="rx0_v", cmap="viridis",
        ).savefig("rx0_v.png", dpi=200)

        domain.plot(show_mesh=False, show_mask=True, show_subdomains=True, tiling=domain.create_tiling(), show_rivers=True).savefig("domain_mask.png")
        
        

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
        report=datetime.timedelta(hours=3),
        report_totals=datetime.timedelta(days=7),
        profile=profile,
    )
