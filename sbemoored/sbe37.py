#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module sbemoored.sbe37"""


import os
import pathlib
import datetime
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import xarray as xr
import pandas as pd
import seabird
import gsw

import gvpy as gv


def proc(
    file,
    time_instrument,
    time_utc,
    data_out=None,
    figure_out=None,
    cal_time=None,
    show_plot=True,
    cut_time=None,
):
    """
    Combining SBE56 processing steps.

    Parameters
    ----------
    file : str or pathlib.Path
        Path to csv data file
    time_instrument : str or np.datetime64
        Instrument time at data download
    time_utc : str or np.datetime64
        UTC time at data download
    data_out : path object, optional
        Path to data output directory
    figure_out : path object, optional
        Path to figure output directory
    cal_time : np.datetime64 object, optional
        Time of post-deployment clock calibration
    show_plot : bool, optional
        Generate data plot. Default True.
    cut_time : np.datetime64, optional
        Cut time series short at this time. Time offset will be linearly scaled
        before it is applied.

    Returns
    -------
    t : xarray.DataArray
        DataArray with thermistor data
    """
    if type(file) != pathlib.PosixPath:
        file = pathlib.Path(file)
    # only read raw file if we haven't written the netcdf file yet
    filename = "{:s}.nc".format(file.stem)
    if data_out:
        savepath = data_out.joinpath(filename)
        if savepath.exists():
            print(
                "already processed\nreading netcdf file from\n{}".format(
                    savepath
                )
            )
            tds = xr.open_dataset(savepath)
            # Update time file stamp. This way, make still recognizes that
            # the file has been worked on.
            savepath.touch()
            savenc = False
        else:
            print("reading csv file")
            tds = read_sbe_cnv(file)
            savenc = True
    else:
        print("reading csv file")
        tds = read_sbe_cnv(file)
        savenc = False
    # cut short (if necessary) apply time drift
    tds = time_offset(
        tds, insttime=time_instrument, utctime=time_utc, cuttime=cut_time
    )
    # save to netcdf
    if savenc:
        save_nc(tds, data_out)
    # plot
    if show_plot:
        plot(tds, figure_out)

    return tds


def read_sbe_cnv(file, lat=0, lon=0):
    """
    Read Seabird SBE37 .cnv file, calculate derived variables
    and return as xarray.Dataset.

    Parameters
    ----------
    file : str
        Complete path to .cnv file
    lat : float, optional
        Latitude (used for gsw calculations). Defaults to zero.
    lon : float, optional
        Longitude (used for gsw calculations). Defaults to zero.

    Returns
    -------
    mc : xarray.Dataset
        Microcat data as Dataset with some metadata in the attributes.
    """
    # Read cnv file using Seabird package
    cnv = seabird.fCNV(file)

    # Parse the data structure generated by the Seabird package, including some
    # of the header information. Some older SBE37 data files do not contain
    # time as a variable, here we need to construct a time vector based on
    # start time and sampling frequency.
    if "timeJV2" in cnv.keys() or "timeJ" in cnv.keys():
        mc = parse_cnv_with_time(cnv)
    else:
        mc = parse_cnv_no_time(cnv)

    # Calculate oceanographic variables
    mc["SP"] = (["time"], gsw.SP_from_C(mc.c, mc.t, mc.p))
    if lat == 0 and lon == 0:
        print(
            "warning: absolute salinity, conservative temperature\n",
            "and density calculation may be inaccurate\n",
            "due to missing latitude/longitude",
        )
    mc["SA"] = (["time"], gsw.SA_from_SP(mc.SP, mc.p, lat=lat, lon=lon))
    mc["CT"] = (["time"], gsw.CT_from_t(mc.SA, mc.t, mc.p))
    mc["sg0"] = (["time"], gsw.sigma0(mc.SA, mc.CT))

    # Add attributes for some of the variables.
    attributes = {
        "p": dict(long_name="pressure", units="dbar"),
        "t": dict(long_name="in-situ temperature", units="°C"),
        "CT": dict(long_name="conservative temperature", units="°C"),
        "SA": dict(long_name="absolute salinity", units=r"kg/m$^3$"),
        "c": dict(long_name="conductivity", units="mS/cm"),
        "SP": dict(long_name="practical salinity", units=""),
        "sg0": dict(
            long_name=r"potential density $\sigma_0$", units=r"kg/m$^3$"
        ),
    }
    for k, att in attributes.items():
        if k in list(mc.variables.keys()):
            mc[k].attrs = att

    # Time offset attribute
    mc.attrs["time offset applied"] = 0

    return mc


def parse_cnv_no_time(cnv):
    # generate time vector
    entries = parse_header(cnv)
    tr = pd.period_range(
        start=entries["start_time"],
        freq="{}s".format(entries["interval"]),
        periods=entries["nvalues"],
    )
    mctime = tr.to_timestamp().to_numpy()
    # data vars
    dvars = {"pr": "p", "t090": "t"}
    mcdata = {}
    for k, di in dvars.items():
        if k in cnv.keys():
            mcdata[di] = (["time"], cnv[k])
    mc = xr.Dataset(data_vars=mcdata, coords={"time": mctime})
    mc.attrs["SN"] = entries["SN"]
    mc.attrs["file"] = cnv.attributes["filename"]
    mc.attrs["sbe_model"] = cnv.attributes["sbe_model"]
    mc.attrs["nvalues"] = entries["nvalues"]
    mc.attrs["start time str"] = entries["start_time_str"]
    mc.attrs["sampling interval"] = entries["interval"]

    # conductivity
    cvars = {"CNDC": "c"}
    cfac = 10 if entries["cnd_units"] == "S/m" else 1
    for k, di in cvars.items():
        if k in cnv.keys():
            # convert from S/m to mS/cm as this is needed for gsw.SP_from_C
            conductivity = cnv[k] * cfac
            mc[di] = (["time"], conductivity)
    return mc


def parse_cnv_with_time(cnv):
    entries = parse_header(cnv)
    # parse time
    if "timeJV2" in cnv.keys():
        mcyday = cnv["timeJV2"]
    elif "timeJ" in cnv.keys():
        mcyday = cnv["timeJ"]
    start_time_str_all = cnv.attributes["start_time"]
    start_time_str = start_time_str_all.split("[")[0]
    base_year = pd.to_datetime(start_time_str).year
    mctime = yday1_to_datetime64(base_year, mcyday)
    # let's make sure the first time stamp we generated matches the string in the cnv file
    assert pd.to_datetime(np.datetime64(mctime[0], "s")) == pd.to_datetime(
        start_time_str
    )

    # data vars
    dvars = {"prdM": "p", "tv290C": "t"}
    mcdata = {}
    for k, di in dvars.items():
        if k in cnv.keys():
            # print(di, ':', k)
            mcdata[di] = (["time"], cnv[k])
    mc = xr.Dataset(data_vars=mcdata, coords={"time": mctime})
    mc.attrs["SN"] = entries["SN"]
    mc.attrs["file"] = cnv.attributes["filename"]
    mc.attrs["sbe_model"] = cnv.attributes["sbe_model"]
    mc.attrs["nvalues"] = entries["nvalues"]
    mc.attrs["start time str"] = entries["start_time_str"]
    mc.attrs["sampling interval"] = entries["interval"]

    # conductivity
    cvars = {"cond0mS/cm": "c", "cond0S/m": "c"}
    for k, di in cvars.items():
        if k in cnv.keys():
            # convert from S/m to mS/cm as this is needed for gsw.SP_from_C
            if k == "cond0S/m":
                conductivity = cnv[k] * 10
            else:
                conductivity = cnv[k]
            mc[di] = (["time"], conductivity)
    return mc


def parse_header(cnv):
    hdr = cnv.raw_header()
    entries = dict(
        nvalues="# nvalues",
        start_time_str="# start_time",
        interval="# interval = ",
        cnd_units="conductivity",
    )
    for hi in hdr["descriptors"].split("\n"):
        for k, v in entries.items():
            if v in hi:
                xi = hi.find(" = ")
                entries[k] = hi[xi + 3 :]

    xi = entries["interval"].find(":")
    entries["interval"] = int(entries["interval"][xi + 2 :])
    entries["nvalues"] = int(entries["nvalues"])
    xi = entries["cnd_units"].find("[")
    entries["cnd_units"] = entries["cnd_units"][xi + 1 : xi + 4]
    if "[" in entries["start_time_str"]:
        entries["start_time_str"] = (
            entries["start_time_str"].split("[")[0].strip()
        )
    entries["start_time"] = pd.to_datetime(
        entries["start_time_str"]
    ).to_datetime64()

    # serial number
    entries2 = dict(
        temp_sn="* Temperature SN = ", cond_sn="* Conductivity SN = ",
    )
    sns = dict(temp_sn="", cond_sn="",)
    for hi in hdr["intro"].split("\n"):
        for k, v in entries2.items():
            if v in hi:
                xi = hi.find(" = ")
                sns[k] = hi[xi + 3 :]

    if sns["temp_sn"] == sns["cond_sn"]:
        entries["SN"] = sns["temp_sn"]

    return entries


def time_offset(tds, insttime, utctime, cuttime=None):
    """Calculate and apply time offset to SBE37 time series. Also cut time
    series short if necessary.

    Parameters
    ----------
    tds : xarray.DataArray
        SBE37 data structure
        
    insttime : np.datetime64
        Instrument time at end of time series
        
    utctime : np.datetime64
        UTC time at end of time series

    cuttime : np.datetime64, optional
        Cut time series short here. If provided, the time offset will be scaled
        linearly before it is applied.

    Returns
    -------
    tds : xarray.DataArray
        Data structure
    """

    if tds.attrs["time offset applied"]:
        print("time offset has already been applied")
    else:
        # cut time series short if necessary
        if cuttime is not None:
            starttime = tds.time[0].data
            overall_length = np.timedelta64(insttime - starttime, "s")
            cut_length = np.timedelta64(cuttime - starttime, "s")
            scale_factor = cut_length / overall_length
            print("scaling time offset by {:1.3f}".format(scale_factor))
            tds = tds.where(tds.time < cuttime, drop=True)
            tds.attrs["nvalues"] = len(tds.time)
        else:
            scale_factor = 1
        # calculate time offset
        offset = (
            np.timedelta64(insttime - utctime, "ms").astype("int")
            * scale_factor
        )
        tds.attrs["time drift in ms"] = offset
        # apply offset
        print(
            "applying time offset of {}ms".format(tds.attrs["time drift in ms"])
        )
        # generate linear time drift vector
        old_time = tds.time.copy()
        time_offset_linspace = np.linspace(
            0, tds.attrs["time drift in ms"], tds.attrs["nvalues"]
        )
        # convert to numpy timedelta64
        # this format can't handle non-integers, so we switch to nanoseconds
        time_offset = [
            np.timedelta64(int(np.round(ti * 1e6)), "ns")
            for ti in time_offset_linspace
        ]
        new_time = old_time - time_offset
        tds["time"] = new_time
        tds.attrs["time offset applied"] = 1

    return tds


def save_nc(tds, data_out):
    # save dataset
    filename = "{:s}.nc".format(tds.attrs["file"][:-4])
    savepath = data_out.joinpath(filename)
    print("Saving to {}".format(savepath))
    tds.to_netcdf(savepath)


def plot(tds, figure_out=None):
    # set up figure
    fig, ax = plt.subplots(
        nrows=3, ncols=1, figsize=(10, 9), sharex=True, constrained_layout=True
    )

    tds.p.plot(ax=ax[0])
    # plot a warning if time offset not applied
    if tds.attrs["time offset applied"] == 1:
        ax[0].text(
            0.05,
            0.9,
            "time offset of {:1.3f} seconds applied".format(
                tds.attrs["time drift in ms"] / 1000
            ),
            transform=ax[0].transAxes,
            backgroundcolor="w",
        )
    else:
        if tds.attrs["time drift in ms"] == 0:
            ax[0].text(
                0.05,
                0.9,
                "WARNING: time offset unknown",
                transform=ax[0].transAxes,
                color="red",
                backgroundcolor="w",
            )
        elif np.absolute(tds.attrs["time drift in ms"]) > 3.6e6:
            ax[0].text(
                0.05,
                0.9,
                "WARNING: time offset more than one hour, not applied",
                transform=ax[0].transAxes,
                color="red",
                backgroundcolor="w",
            )
        else:
            ax[0].text(
                0.05,
                0.9,
                "time offset not yet applied",
                transform=ax[0].transAxes,
                backgroundcolor="w",
            )
    tds.t.plot(ax=ax[1])
    tds.c.plot(ax=ax[2])

    for axi in ax:
        axi.grid()
        axi.set(title="SBE37 SN {}".format(tds.attrs["SN"]))
        axi.set(xlabel="")
        gv.plot.concise_date(axi)
    ax[0].invert_yaxis()

    if figure_out is not None or False:
        figurename = "{:s}.png".format(tds.attrs["file"][:-4])
        plt.savefig(figure_out.joinpath(figurename), facecolor="w", dpi=300)


def yday1_to_datetime64(baseyear, yday):
    """
    Convert year day (starting at yday 1) to numpy's datetime64 format.

    Parameters
    ----------
    baseyear : int
        Base year
    yday : float
        Year day

    Returns
    -------
    time : np.datetime64
        Time in numpy datetime64 format
    """
    base = datetime.datetime(baseyear, 1, 1, 0, 0, 0)
    time = [base + datetime.timedelta(days=ti) for ti in yday - 1]
    # convert to numpy datetime64
    time64 = np.array([np.datetime64(ti, "ms") for ti in time])
    return time64
