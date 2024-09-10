#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module sbemoored.sbe37"""

import os
import pathlib
import datetime
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import xarray as xr
import pandas as pd
import xmltodict
import gsw


def proc(
    file,
    time_instrument=None,
    time_utc=None,
    data_out=None,
    file_name=None,
    figure_out=None,
    cal_time=None,
    show_plot=True,
    cut_end=None,
    cut_beg=None,
    lat=None,
    lon=None,
    meta=None,
):
    """
    Combining SBE56 processing steps.

    Parameters
    ----------
    file : str or pathlib.Path
        Path to csv data file
    time_instrument : str or np.datetime64, optional
        Instrument time at data download
    time_utc : str or np.datetime64, optional
        UTC time at data download
    data_out : path object, optional
        Path to data output directory
    file_name : str, optional
        File name. Defaults to base name of raw file.
    figure_out : path object, optional
        Path to figure output directory
    cal_time : np.datetime64 object, optional
        Time of post-deployment clock calibration
    show_plot : bool, optional
        Generate data plot. Default True.
    cut_end : np.datetime64, optional
        Cut time series short at this time. Time offset will be linearly scaled
        before it is applied.
    cut_beg : np.datetime64, optional
        Cut time series at beginning.
    lat : float, optional
        Latitude.
    lon : float, optional
        Longitude.
    meta : dict, optional
        Meta data.

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
        if savepath.exists() and cut_end is None and cut_beg is None:
            print("already processed\nreading netcdf file from\n{}".format(savepath))
            tds = xr.open_dataset(savepath)
            # Update time file stamp. This way, make still recognizes that
            # the file has been worked on.
            savepath.touch()
            savenc = False
        else:
            print("reading csv file")
            tds = read_sbe_cnv(file, lat=lat, lon=lon)
            savenc = True
    else:
        print("reading csv file")
        tds = read_sbe_cnv(file, lat=lat, lon=lon)
        savenc = False
    # cut short (if necessary) apply time drift
    tds = time_offset(tds, insttime=time_instrument, utctime=time_utc, cuttime=cut_end)
    # cut at beginning
    if cut_beg is not None:
        mask = tds.time > cut_beg
        tds = tds.where(mask, drop=True)

    # read xmlcon and set attributes
    try:
        cfgp = read_xml_config(file)
        vars = dict(t="TemperatureSensor1", c="ConductivitySensor1", p="PressureSensor")
        for k, v in vars.items():
            tds[k].attrs["SN"] = int(cfgp[v].cal["SerialNumber"])
            tds[k].attrs["CalDate"] = cfgp[v].cal["CalibrationDate"]
    except:
        pass

    # add more attributes
    if isinstance(meta, dict):
        for k, v in meta.items():
            tds.attrs[k] = v
    if lon is not None:
        tds.attrs["lon"] = lon
    if lat is not None:
        tds.attrs["lat"] = lat

    # save to netcdf
    if savenc:
        save_nc(tds, data_out, file_name)
    # plot
    if show_plot:
        if file_name is not None:
            figure_name = f'{file_name[:file_name.find(".nc")]}.png'
        else:
            figure_name = None
        plot(tds, figure_out, figure_name)

    return tds


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
    elif insttime is None or utctime is None:
        print("no time readings provided; not applying any offset")
        tds.attrs["time offset applied"] = 0
        tds.attrs["time drift in ms"] = "N/A"
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
        # we also need to scale the time offset if the time series stops early
        else:
            scale_factor = 1
        # calculate time offset
        offset = np.timedelta64(insttime - utctime, "ms").astype("int") * scale_factor
        tds.attrs["time drift in ms"] = offset
        # apply offset
        print("applying time offset of {}ms".format(tds.attrs["time drift in ms"]))
        # generate linear time drift vector
        old_time = tds.time.copy()
        time_offset_linspace = np.linspace(
            0, tds.attrs["time drift in ms"], tds.attrs["nvalues"]
        )
        # convert to numpy timedelta64
        # this format can't handle non-integers, so we switch to nanoseconds
        time_offset = [
            np.timedelta64(int(np.round(ti * 1e6)), "ns") for ti in time_offset_linspace
        ]
        new_time = old_time - time_offset
        tds["time"] = new_time
        tds.attrs["time offset applied"] = 1

    return tds


def save_nc(tds, data_out, filename=None):
    # tds['time'] = gv.time.convert_units(tds.time, unit='ms')
    # save dataset
    if filename is None:
        filename = "{:s}.nc".format(tds.attrs["file"][:-4])
    savepath = data_out.joinpath(filename)
    print("Saving to {}".format(savepath))
    tds.to_netcdf(savepath)


def plot(tds, figure_out=None, figure_name=None):
    # set up figure
    fig, ax = plt.subplots(
        nrows=3, ncols=1, figsize=(10, 9), sharex=True, constrained_layout=True
    )

    tds.p.plot(ax=ax[0])
    # plot a warning if time offset not applied
    if tds.attrs["time drift in ms"] == "N/A":
        ax[0].text(
            0.05,
            1.05,
            "WARNING: no time offset provided",
            transform=ax[0].transAxes,
            color="red",
            backgroundcolor="w",
        )
    elif tds.attrs["time offset applied"] == 1:
        ax[0].text(
            0.05,
            1.05,
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
                1.05,
                "WARNING: time offset zero",
                transform=ax[0].transAxes,
                color="red",
                backgroundcolor="w",
            )
        elif np.absolute(tds.attrs["time drift in ms"]) > 3.6e6:
            ax[0].text(
                0.05,
                1.05,
                "WARNING: time offset more than one hour, not applied",
                transform=ax[0].transAxes,
                color="red",
                backgroundcolor="w",
            )
        else:
            ax[0].text(
                0.05,
                1.05,
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
        _concise_date(axi)
    ax[0].invert_yaxis()

    if figure_out is not None or False:
        if figure_name is None:
            figure_name = "{:s}.png".format(tds.attrs["file"][:-4])
        plt.savefig(figure_out.joinpath(figure_name), facecolor="w", dpi=300)


def clock_check(tds, plot=True, timedelta=200):
    """Check clock interval on SBE 37.

    Parameters
    ----------
    tds : xr.Dataset
        SBE 37 dataset with time variable.
    plot : bool, optional
        Plot time delta around time where clock changes. Default True.
    timedelta : int, optional
        Time deviation from median in ms. Default 200.

    Returns
    -------
    cutoff : np.datetime64
        First time stamp where clock is off. Returns None when the clock is fine.
    """

    dt = tds.time.diff(dim="time")
    dtm = dt.median(dim="time")
    cutoff = dt.where(
        np.absolute(dt - dtm) > np.timedelta64(timedelta, "ms"), drop=True
    ).time
    if cutoff.size > 0:
        cutoff = cutoff.data[0]
    else:
        cutoff = None
    print(cutoff)
    if plot and cutoff is not None:
        fig, ax = plt.subplots(
            constrained_layout=True,
            dpi=75,
            )
        t1 = cutoff - np.timedelta64(4, "h")
        t2 = cutoff + np.timedelta64(4, "h")
        tp = dt.sel(time=slice(t1, t2))
        tp.plot(ax=ax)
    return cutoff


def _find_xmlconfig(file):
    """Generate path to xml config file for current hex file.
    Config file needs to be in the same directory as the hex file."""
    name = file.stem
    # try upper case filename
    xmlfile = name.upper() + ".XMLCON"
    p = file.parent
    xmlfile = p.joinpath(xmlfile)
    # use os.listdir to find the actual case of the filename if the upper
    # case did not work.
    # if xmlfile.name not in os.listdir(os.path.dirname(xmlfile)):
    #     xmlfile = name.lower() + ".XMLCON"
    #     xmlfile = p.joinpath(xmlfile)
    return xmlfile


def read_xml_config(file):
    """Read xml config file."""
    xmlfile = _find_xmlconfig(file)
    try:
        with open(xmlfile) as fd:
            tmp = xmltodict.parse(fd.read())
    except OSError as e:
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), e.filename)
    tmp = tmp["SBE_InstrumentConfiguration"]
    tmp = tmp["Instrument"]
    sa = tmp["SensorArray"]["Sensor"]
    cfg = {}
    ti = 0
    ci = 0
    for si in sa:
        keys = si.keys()
        for k in list(keys):
            if "@" not in k and k != "NotInUse":
                if k == "TemperatureSensor":
                    ti += 1
                    kstr = "{}{}".format(k, ti)
                elif k == "ConductivitySensor":
                    ci += 1
                    kstr = "{}{}".format(k, ci)
                else:
                    kstr = k
                cfg[kstr] = si
                cfg[kstr]["cal"] = cfg[kstr][k]
                del cfg[kstr][k]
    cfgp = pd.DataFrame(cfg)
    cfgp = _xml_coeffs_to_float(cfgp)
    return cfgp


def _xml_coeffs_to_float(cfgp):
    # Convert calibration coefficients to floats.
    keep_strings = [
        "@SensorID",
        "SerialNumber",
        "CalibrationDate",
        "UseG_J",
    ]
    for k in cfgp.keys():
        for ki in cfgp[k]["cal"].keys():
            if isinstance(cfgp[k]["cal"][ki], str):
                if ki not in keep_strings:
                    cfgp[k]["cal"][ki] = float(cfgp[k]["cal"][ki])
            elif isinstance(cfgp[k]["cal"][ki], list):
                for i, li in enumerate(cfgp[k]["cal"][ki]):
                    for kli in li.keys():
                        cfgp[k]["cal"][ki][i][kli] = float(cfgp[k]["cal"][ki][i][kli])
        # We can't have None values in the xarray.Dataset later on
        # or otherwise it won't properly write to netcdf. Therefore,
        # convert any None items to 'N/A'
        for ki, v in cfgp[k]["cal"].items():
            if v is None:
                cfgp[k].cal[ki] = "N/A"
    return cfgp


def _concise_date(ax=None, minticks=3, maxticks=10, show_offset=True, **kwargs):
    """
    Better date ticks using matplotlib's ConciseDateFormatter.

    Parameters
    ----------
    ax : axis handle
        Handle to axis (optional).
    minticks : int
        Minimum number of ticks (optional, default 6).
    maxticks : int
        Maximum number of ticks (optional, default 10).
    show_offset : bool, optional
        Show offset string to the right (default True).

    Note
    ----
    Currently only works for x-axis

    See Also
    --------
    matplotlib.mdates.ConciseDateFormatter : For formatting options that
      can be used here.
    """
    if ax is None:
        ax = plt.gca()
    locator = mdates.AutoDateLocator(minticks=minticks, maxticks=maxticks)
    formatter = mdates.ConciseDateFormatter(locator, show_offset=show_offset, **kwargs)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(formatter)
    # remove axis label "time" if present
    if ax.get_xlabel() == "time":
        _ = ax.set_xlabel("")


def read_header(file):
    """Read header from SBE .cnv file

    Parameters
    ----------
    file : str or pathlib.Path
        .cnv file
    """
    header = []
    with open(file) as f:
        for s in f:
            if s.startswith("*END*"):
                header.append(s)
                break
            else:
                header.append(s)
    return header


def get_header_length(header):
    return len(header)


def parse_header(header):
    columns = []
    for h in header:
        if h.startswith("# name"):
            columns.append(h.split(" = ")[1].strip())
        elif h.startswith("# start_time"):
            tmp = h.split(" = ")[1]
            start_time_str = tmp.split(" [")[0]
        elif h.startswith("* Temperature SN"):
            sn = int(h.split(" = ")[1].strip())
        elif h.startswith("# nvalues"):
            nvalues = int(h.split(" = ")[1].strip())
        elif h.startswith("* sample interval"):
            sample_interval_str = h.split(" = ")[1]

    names = []
    units = []
    for c in columns:
        nam, u = c.split(": ")
        names.append(nam)
        units.append(u)

    # time
    start_time = pd.to_datetime(start_time_str).to_datetime64()
    base_year = pd.to_datetime(start_time_str).year

    return dict(
        sn=sn,
        nvalues=nvalues,
        names=names,
        units=units,
        start_time=start_time,
        start_time_str=start_time_str,
        base_year=base_year,
        sample_interval=sample_interval_str,
    )


def parse_time(ds):
    if "timeJV2" in ds:
        mcyday = ds["timeJV2"]
        mctime = yday1_to_datetime64(ds.base_year, mcyday.data)
    elif "timeJ" in ds:
        mcyday = ds["timeJ"]
        mctime = yday1_to_datetime64(ds.base_year, mcyday.data)
    else:
        interval = ds.sample_interval.split(" ")[0]
        tr = pd.period_range(
            start=ds.start_time,
            freq=f"{interval}s",
            periods=len(ds.index),
        )
        mctime = tr.to_timestamp().to_numpy()

    ds["time"] = (("index"), mctime)
    ds = ds.swap_dims(index="time")
    ds = ds.drop_vars(["index", "flag"])
    return ds


def parse_variables(ds):
    for k, v in ds.items():
        if k.startswith("pr"):
            ds = ds.rename({k: "p"})
            ds.p.attrs = dict(long_name="pressure", units="dbar")
        if k.startswith("t") and "90" in k:
            ds = ds.rename({k: "t"})
            ds.t.attrs = dict(long_name="temperature", units="°C")
        if k.startswith("c") and "S/" in k:
            ds = ds.rename({k: "c"})
            # convert from S/m to mS/cm as this is needed for gsw.SP_from_C
            if "S/m" in k:
                ds["c"] = ds.c * 10
            ds.c.attrs = dict(long_name="conductivity", units="mS/cm")
        if k.startswith("sal"):
            ds = ds.rename({k: "s"})
            ds.s.attrs = dict(long_name="salinity", units="psu")

    return ds


def gsw_calcs(ds, lat=None, lon=None):
    # Calculate oceanographic variables
    ds["SP"] = (["time"], gsw.SP_from_C(ds.c, ds.t, ds.p).data)
    if lat is None and lon is None:
        print(
            "warning: calculation of absolute salinity, conservative temperature, \n",
            "and density may be inaccurate due to missing latitude/longitude",
        )
        lat, lon = 0, 0
    ds["SA"] = (["time"], gsw.SA_from_SP(ds.SP, ds.p, lat=lat, lon=lon).data)
    ds["CT"] = (["time"], gsw.CT_from_t(ds.SA, ds.t, ds.p).data)
    ds["sg0"] = (["time"], gsw.sigma0(ds.SA, ds.CT).data)

    # Add attributes for some of the variables.
    attributes = {
        "CT": dict(long_name="conservative temperature", units="°C"),
        "SA": dict(long_name="absolute salinity", units=r"kg/m$^3$"),
        "SP": dict(long_name="practical salinity", units=""),
        "sg0": dict(long_name=r"potential density $\sigma_0$", units=r"kg/m$^3$"),
    }
    for k, att in attributes.items():
        if k in list(ds.variables.keys()):
            ds[k].attrs = att
    return ds


def read_sbe_cnv(file, lat=None, lon=None):
    """Read Seabird .cnv file, calculate derived variables and return as
    xarray.Dataset.

    Parameters
    ----------
    file : str or pathlib.Path
        Complete path to .cnv file
    lat : float, optional
        Latitude (used for gsw calculations).
    lon : float, optional
        Longitude (used for gsw calculations).

    Returns
    -------
    mc : xarray.Dataset
        Microcat data as Dataset with some metadata as attributes.
    """
    header = read_header(file)
    n = get_header_length(header)
    ph = parse_header(header)
    tmp = pd.read_csv(file, skiprows=n, names=ph["names"], sep=r"\s+")

    tmp = tmp.to_xarray()
    tmp.attrs["base_year"] = ph["base_year"]
    tmp.attrs["SN"] = ph["sn"]
    tmp.attrs["nvalues"] = ph["nvalues"]
    tmp.attrs["sample_interval"] = ph["sample_interval"]
    start_time = np.datetime_as_string(ph["start_time"], unit="s").replace("T", " ")
    tmp.attrs["start_time"] = start_time

    tmp.attrs["file"] = file.name

    tmp = parse_time(tmp)

    tmp = parse_variables(tmp)

    tmp = gsw_calcs(tmp, lat=lat, lon=lon)

    # Time offset attribute
    tmp.attrs["time offset applied"] = 0

    return tmp


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
    time64 = np.array([np.datetime64(ti, "ns") for ti in time])
    return time64
