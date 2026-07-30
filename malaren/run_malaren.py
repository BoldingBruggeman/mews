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
use_tracers = False
bathy_name = "Bathymetry/bathymetry.nc"
depth_name = "bathy_smooth_fixedpoints"

# Met scaling
scl_wind = 1.0 # Multiply this number by u10 and v10

# Light extinction
light_A = 0.58 # non-visible fraction of shortwave radiation (-)
light_kc1 = 2.5 # attenuation of non-visible fraction of shortwave radiation (m-1)
light_kc2 = 1.0 # attenuation of visible fraction of shortwave radiation (m-1)

# Outflow regulation (i.e. dynamic outflow based on water level)
outflow_reg = True

### Output
# Frequency
interval_output = 7 # days
interval_report = 7 # days
interval_meteo = 140 # days
interval_highfreq = 24 # hours

# Local high-frequency output sites
output_local = True
local_sites = {"prastfjarden": [637650, 6590250],
               "gorvaln": [652950, 6591150],
               "fiskarfjarden": [664650, 6578250],
               "rodstensfjarden": [661050, 6572550],
               "granfjarden": [601350, 6598950],
               "galten": [567450, 6590250],
               "blacken": [587550, 6593850],
               "vasterasfjarden": [589050, 6604650],
               "arnofjarden": [626250, 6593250],
               "ekoln": [646650, 6626550], 
               "skarven": [658650, 6604950],
               "ullfjarden": [642750, 6608250],
               "riddarfjarden": [671850, 6580050]}

# Format: name as first entry, then X and Y coordinates between []                                                                  
def create_domain(
    runtype: int,
    rivers: bool,
    **kwargs,
):
    import netCDF4
    import glob
    import os

    with netCDF4.Dataset(bathy_name) as nc:
        nc.set_auto_mask(False)
        domain = pygetm.domain.create_cartesian(
            nc["x"][:],
            nc["y"][:],
            lon=nc["lon"],
            lat=nc["lat"],
            H=nc[depth_name][:, :],
            mask=np.where(nc[depth_name][...] == -9999.0, 0, 1),
            z0=0.01,
        )
    # domain.mask_indices(261, 270 + 1, 43, 43 + 1)
    # domain.mask_indices(334, 335 + 1, 49, 63 + 1)
    # domain.mask_indices(335, 335 + 1, 49, 63 + 1)
    # domain.mask_indices(304, 304 + 1, 38, 44 + 1)

    domain.limit_velocity_depth()
    domain.cfl_check()

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
        if outflow_reg:
            river_list.append(domain.rivers.add_by_location("norrstrom", 674250, 6580050, coordinate_type=pygetm.CoordinateType.XY))
            river_list.append(domain.rivers.add_by_location("hammarby", 674550, 6577950, coordinate_type=pygetm.CoordinateType.XY))
            river_list.append(domain.rivers.add_by_location("sodertalje", 649650, 6565950, coordinate_type=pygetm.CoordinateType.XY))
        else:
            for river in glob.glob("Rivers/Outflow_file*.nc"):
                name = os.path.basename(river)
                name = name.replace("Outflow_file_", "").replace(".nc", "")
                with netCDF4.Dataset(river) as r:
                    lon = r["lon"][:]
                    lat = r["lat"][:]
                    river_list.append(
                        domain.rivers.add_by_location(
                            name, float(lon), float(lat), coordinate_type=pygetm.CoordinateType.LONLAT
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
    import netCDF4
    import numpy as np
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
        # delay_slow_ip=False,
    )
    final_kwargs.update(kwargs)
    sim = pygetm.Simulation(domain, runtype=runtype, gotm = "gotm.yaml", fabm = "fabm-selmaprotbas.yaml", **final_kwargs)
    
    # for name, field in sorted(sim.output_manager.fields.items()):
        # print(f"{name} - {getattr(field, 'long_name', '')}")
    # raise(Exception("!!!"))
    
    if sim.runtype < pygetm.RunType.BAROCLINIC:
        sim.sst = sim.airsea.t2m
    if sim.runtype == pygetm.RunType.BAROCLINIC:
        sim.radiation.set_jerlov_type(pygetm.Jerlov.Type_III) # sim.radiation.set_jerlov_type(pygetm.Jerlov.Type_II)
        # Overwrite Jerlov type by user settings
        sim.radiation.A.fill(light_A)
        sim.radiation.kc1.fill(light_kc1)
        sim.radiation.kc2.fill(light_kc2)

    if not args.load_restart and sim.runtype == pygetm.RunType.BAROCLINIC:
        if True:
            sim.temp.set(2)
            sim.salt.set(0.1)
        else:
            print("Read from files")
            # sim.salt.set(
            # pygetm.input.from_nc(
            #    os.path.join(args.setup_dir, "Input/initialConditions.nc"), "salt"
            # ),
            # on_grid=True,
        # )
        sim.density.convert_ts(sim.salt, sim.temp)
        sim.temp[..., sim.T.mask == 0] = pygetm.constants.FILL_VALUE
        sim.salt[..., sim.T.mask == 0] = pygetm.constants.FILL_VALUE

    ERA_path = "ERA5/era5_????.nc"
    sim.airsea.u10.set(pygetm.input.from_nc(ERA_path, "u10") * scl_wind)
    sim.airsea.v10.set(pygetm.input.from_nc(ERA_path, "v10") * scl_wind)
    sim.airsea.t2m.set(pygetm.input.from_nc(ERA_path, "t2m") - 273.15)
    sim.airsea.d2m.set(pygetm.input.from_nc(ERA_path, "d2m") - 273.15)
    sim.airsea.sp.set(pygetm.input.from_nc(ERA_path, "sp"))
    sim.airsea.tcc.set(pygetm.input.from_nc(ERA_path, "tcc"))
    ERA_path = "ERA5/precip_????.nc"
    sim.airsea.tp.set(pygetm.input.from_nc(ERA_path, "tp") / 3600.0)
    for river in sim.rivers.values():
        if "outflow" in river.name and not(outflow_reg):
            ### Prescribed outflow
            river.flow.set(pygetm.input.from_nc(f"Rivers/Outflow_file_{river.name}.nc", "q"))
        elif river.name in ["norrstrom", "hammarby", "sodertalje"] and outflow_reg:
            ### Regulated outflow
            river.flow.set(0.0) # To be changed later at runtime
        elif "nctracer_" in river.name:
            ### Tracer
            filename_tracer = river.name.replace("nctracer_", "")
            river.flow.set(pygetm.input.from_nc(f"Tracers/Tracer_file_{filename_tracer}.nc", "flow"))

            ncfile = netCDF4.Dataset(f"Tracers/Tracer_file_{filename_tracer}.nc")
            tracer_names = np.array(list(ncfile.variables.keys()))
            tracer_names = tracer_names[~np.isin(tracer_names, ["lon", "lat", "time", "flow"])]

            for trc in tracer_names:
                river[trc].set(pygetm.input.from_nc(f"Tracers/Tracer_file_{filename_tracer}.nc", trc))

            ncfile.close()
        else:
            ### Inflow
            river.flow.set(pygetm.input.from_nc(f"Rivers/Inflow_file_{river.name}.nc", "q"))
            
            # Temperature and oxygen
            river["temp"].follow_target_cell = False
            river["temp"].set(pygetm.input.from_nc(f"Rivers/Inflow_file_{river.name}.nc", "temp"))
            river["selmaprotbas_o2"].follow_target_cell = False
            river["selmaprotbas_o2"].set(pygetm.input.from_nc(f"Rivers/Inflow_file_{river.name}.nc", "o2"))
            
            # Nutrients
            river["selmaprotbas_po"].follow_target_cell = False #(True makes it use the value from the basin)
            river["selmaprotbas_po"].set(pygetm.input.from_nc(f"Rivers/Inflow_file_{river.name}.nc", "PO4"))
            river["selmaprotbas_aa"].follow_target_cell = False
            river["selmaprotbas_aa"].set(pygetm.input.from_nc(f"Rivers/Inflow_file_{river.name}.nc", "NH4"))
            river["selmaprotbas_nn"].follow_target_cell = False
            river["selmaprotbas_nn"].set(pygetm.input.from_nc(f"Rivers/Inflow_file_{river.name}.nc", "NO3"))
            river["selmaprotbas_si"].follow_target_cell = False
            river["selmaprotbas_si"].set(pygetm.input.from_nc(f"Rivers/Inflow_file_{river.name}.nc", "Si"))
            river["selmaprotbas_dd_c"].follow_target_cell = False
            river["selmaprotbas_dd_c"].set(pygetm.input.from_nc(f"Rivers/Inflow_file_{river.name}.nc", "dd_n") * 6.625) # N input adjusted with Redfield ratio
            river["selmaprotbas_dd_n"].follow_target_cell = False
            river["selmaprotbas_dd_n"].set(pygetm.input.from_nc(f"Rivers/Inflow_file_{river.name}.nc", "dd_n"))
            river["selmaprotbas_dd_p"].follow_target_cell = False
            river["selmaprotbas_dd_p"].set(pygetm.input.from_nc(f"Rivers/Inflow_file_{river.name}.nc", "dd_p"))
            river["dom_dom_a"].follow_target_cell = False
            river["dom_dom_a"].set(0.0)
            river["dom_dom_b"].follow_target_cell = False
            river["dom_dom_b"].set(pygetm.input.from_nc(f"Rivers/Inflow_file_{river.name}.nc", "DOM"))
            
            river["selmaprotbas_pw"].follow_target_cell = False
            river["selmaprotbas_pw"].set(0.0) # Unlikely this is high in the river water, and concentration must be set in order to stabilise the simulation
            river["selmaprotbas_dd_si"].follow_target_cell = False
            river["selmaprotbas_dd_si"].set(pygetm.input.from_nc(f"Rivers/Inflow_file_{river.name}.nc", "Si")) # Needed to stabilise the simulation. Assume same as Si (For NO3 -> dd_n, it is also in a similar order of magnitude)
    
    sim["age_age_of_water"].river_follow[:] = False # By default, precipitation also has age 0
    
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
        interval=datetime.timedelta(days = interval_meteo),
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
        interval=datetime.timedelta(days=interval_output),
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
            interval=datetime.timedelta(days=interval_output),
            sync_interval=None,
        )
    output.request("Ht", "uk", "vk", time_average = True) # , "ww", "SS", "num"
    if args.debug_output:
        output.request("fpk", "fqk", "advpk", "advqk")  # 'diffpk', 'diffqk')

    if sim.runtype == pygetm.RunType.BAROCLINIC:
        output.request("temp", "sst", "hnt", "nuh", "rad") # "salt", "rho", "NN", "rad", , "nuh", "taub" # "rad" should add "zft"
        if args.debug_output:
            output.request("idpdx", "idpdy")
        if use_adaptive:
            output.request("nug", "ga", "dga")

    if sim.fabm:
        output.request("selmaprotbas_po", "selmaprotbas_nn", "selmaprotbas_aa", "selmaprotbas_o2", "phytoplankton_chla", "dom_dom_a", "dom_dom_b", "age_age_of_water")

    # Generate output for high-frequency single locations
    if output_local:
        for sitename,loc in local_sites.items():
            path = Path(output_dir, setup + "_local_" + sitename + ".nc")
            output = sim.output_manager.add_netcdf_file(
                str(path),
                interval=datetime.timedelta(hours = interval_highfreq),
                sync_interval=None,
            )
            transform = pygetm.output.operators.IndexXY.parameterize(x = loc[0], y = loc[1], coordinate_type=pygetm.CoordinateType.XY)
            output.request("temp", "selmaprotbas_po", "selmaprotbas_nn", "selmaprotbas_aa", "selmaprotbas_o2", "phytoplankton_chla", "dom_dom_a", "dom_dom_b", transforms=[transform])

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
        if outflow_reg:
            import math
            # Based on measurements
            # Saturating exponential equation above a certain wlvl
            def outflow_norrstrom(wlvl):
                background_flow = 21.5 + 5.8 * wlvl
                high_flow = 300 - 300 * math.exp(-8 * wlvl)
                return(-1 * max(0.0, background_flow, high_flow))
            def outflow_hammarby(wlvl):
                background_flow = -0.94 + 20.5 * wlvl
                high_flow = 95 - 95 * math.exp(-20 * (wlvl - 0.05))
                return(-1 * max(0.0, background_flow, high_flow))
            def outflow_sodertalje(wlvl):
                background_flow = 16 + 4 * wlvl
                high_flow = 35 - 35 * math.exp(-50 * (wlvl - 0.2))
                return(-1 * max(0.0, background_flow, high_flow))
            
            local_river_names = [var.name for var in sim.rivers.values()]
            if "norrstrom" in local_river_names:
                norrstrom_i = sim.rivers["norrstrom"].i
                norrstrom_j = sim.rivers["norrstrom"].j
            if "hammarby" in local_river_names:
                hammarby_i = sim.rivers["hammarby"].i
                hammarby_j = sim.rivers["hammarby"].j
            if "sodertalje" in local_river_names:
                sodertalje_i = sim.rivers["sodertalje"].i
                sodertalje_j = sim.rivers["sodertalje"].j
        
        sim.start(
            simstart,
            timestep=timestep,
            split_factor=20,
            **kwargs,
        )
        
        while sim.time < simstop:
            if outflow_reg:
                if "norrstrom" in local_river_names:
                    sim.rivers["norrstrom"].flow[...] = (outflow_norrstrom(sim["zt"].all_values[norrstrom_j, norrstrom_i]))
                if "hammarby" in local_river_names:
                    sim.rivers["hammarby"].flow[...] = (outflow_hammarby(sim["zt"].all_values[hammarby_j, hammarby_i]))
                if "sodertalje" in local_river_names:
                    sim.rivers["sodertalje"].flow[...] = (outflow_sodertalje(sim["zt"].all_values[sodertalje_j, sodertalje_i]))
            sim.advance() # For debugging mode: check_finite=True
            # sim.check_finite(True)
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
        default="bathymetry", # initial smoothing: laplacian_rx0
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
        f = domain.plot(show_mesh=False, show_subdomains=True, tiling=domain.create_tiling())
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
        simstart = sim.load_restart(args.load_restart, decode_timedelta=False) # Old: , decode_timedelta=False) # decode_timedelta=False is to avoid problems with age_of_water

    simstart = datetime.datetime.strptime(args.start, "%Y-%m-%d %H:%M:%S")
    simstop = datetime.datetime.strptime(args.stop, "%Y-%m-%d %H:%M:%S")
    profile = setup if args.profile is not None else None
    run(
        sim,
        simstart,
        simstop,
        dryrun=args.dryrun,
        report=datetime.timedelta(days=interval_report),
        report_totals=datetime.timedelta(days=interval_report),
        profile=profile,
    )
