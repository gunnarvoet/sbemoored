#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module sbemoored.sbe56"""


import os
import pathlib
import datetime
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import xarray as xr
import pandas as pd

import gvpy as gv


def proc(
    file,
    time_instrument=None,
    time_utc=None,
    data_out=None,
    figure_out=None,
    cal_time=None,
    show_plot=True,
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
    figure_out : path object, optional
        Path to figure output directory
    cal_time : np.datetime64 object, optional
        Time of post-deployment clock calibration
    show_plot : bool, optional
        Generate data plot. Default True.

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
            tds = xr.open_dataarray(savepath)
            # Update time file stamp. This way, make still recognizes that
            # the file has been worked on.
            savepath.touch()
            savenc = False
        else:
            print("reading csv file")
            tds = read_csv(file)
            savenc = True
    else:
        print("reading csv file")
        tds = read_csv(file)
        savenc = False
    # apply time drift
    if time_instrument is not None:
        tds = time_offset(tds, insttime=time_instrument, utctime=time_utc)
    else:
        tds.attrs["time drift in ms"] = np.nan
    # save to netcdf
    if savenc:
        save_nc(tds, data_out)
    # plot
    if show_plot:
        plot(tds, figure_out, cal_time)

    return tds


def read_csv_header(file):
    header = []
    out = dict()
    with open(file) as f:
        line = f.readline().strip()
        while line[0] == "%":
            if "Serial Number" in line:
                sni = line.find(" = ")
                out["SN"] = line[sni + 3 :]
            if "Instrument type" in line:
                sni = line.find(" = ")
                out["model"] = line[sni + 3 :]
            if "Firmware Version" in line:
                sni = line.find(" = ")
                out["firmware version"] = line[sni + 3 :]
            if "Calibration Date" in line:
                sni = line.find(" = ")
                out["calibration date"] = line[sni + 3 :]
            if "Source file" in line:
                parts = line.split("\\")
                out["source file"] = parts[-1]
            header.append(line.strip().replace("% ", ""))
            line = f.readline().strip()

    return out


def read_csv(file):
    tt = pd.read_csv(file, header=0, skiprows=range(11))
    timestr = tt.Date + " " + tt.Time
    timepd = pd.to_datetime(timestr.values)
    t = xr.DataArray(
        data=tt.Temperature.to_numpy(),
        coords=dict(time=timepd),
        dims="time",
        name="t",
    )

    header = read_csv_header(file)

    # calculate sampling period in s
    sampling_period = np.round(
        t.time[:100]
        .diff(dim="time")
        .median()
        .data.astype("timedelta64[ns]")
        .astype(int)
        / 1e9
    )
    t.attrs["units"] = "°C"
    t.attrs["long_name"] = "temperature"
    t.attrs["SN"] = header["SN"]
    t.attrs["model"] = header["model"]
    t.attrs["firmware version"] = header["firmware version"]
    t.attrs["file"] = header["source file"]
    t.attrs["calibration date"] = header["calibration date"]
    t.attrs["sample size"] = len(t)
    t.attrs["sampling period in s"] = sampling_period
    t.attrs["time offset applied"] = 0
    return t


def time_offset(tds, insttime, utctime, time_offset=None):
    """Calculate and apply time offset to SBE56 time series

    Parameters
    ----------
    tds : xarray.DataArray
        SBE56 data structure
    insttime : np.datetime64
        Instrument time at end of time series
    utctime : np.datetime64
        UTC time at end of time series
    time_offset : np.timedelta64, optional
        Time offset. If given, insttime and utctime are ignored.

    Returns
    -------
    tds : xarray.DataArray
    """

    if tds.attrs["time offset applied"]:
        print("time offset has already been applied")
    else:
        if time_offset is None:
            # calculate time offset
            offset = np.timedelta64(insttime - utctime, "ms").astype("int")
        else:
            offset = np.timedelta64(time_offset, "ms").astype("int")
        tds.attrs["time drift in ms"] = offset
        # apply offset
        print(
            "applying time offset of {}ms".format(tds.attrs["time drift in ms"])
        )
        # generate linear time drift vector
        old_time = tds.time.copy()
        time_offset_linspace = np.linspace(
            0, tds.attrs["time drift in ms"], tds.attrs["sample size"]
        )
        # convert to numpy timedelta64
        # this format can't handle non-integers, so we switch to nanoseconds
        time_offset_vec = [
            np.timedelta64(int(np.round(ti * 1e6)), "ns")
            for ti in time_offset_linspace
        ]
        new_time = old_time - time_offset_vec
        tds["time"] = new_time
        tds.attrs["time offset applied"] = 1

    return tds


def save_nc(tds, data_out):
    # save dataset
    filename = "{:s}.nc".format(tds.attrs["file"][:-4])
    savepath = data_out.joinpath(filename)
    print("Saving to {}".format(savepath))
    tds.to_netcdf(savepath)


def plot(solo, figure_out=None, cal_time=None):
    """Plot SBE56 time series.

    Parameters
    ----------
    solo : xarray.DataArray
        DataArray with thermistor data
    figure_out : path object, optional
        Path to figure output directory. Default None.
    cal_time : np.datetime64 or tuple
        Clock calibration(s) from warm water dip, either a single np.datetime64
        time or a tuple of two time stamps with pre- and post-deployment clock
        calibration.
    """
    # Generate a list of one or two time stamps and check format of time
    # stamp(s) provided.
    def time_ok(time):
        return True if type(time) == np.datetime64 else False

    if cal_time is not None:
        if time_ok(cal_time):
            cal = [cal_time]
        elif type(cal_time) == tuple:
            if time_ok(cal_time[0]) and time_ok(cal_time[1]):
                cal = [cal_time[0], cal_time[1]]
        else:
            raise TypeError(
                "Provide cal time(s) as single np.datetime64 or tuple thereof"
            )

    # Check if calibration time(s) are outside the time series.
    if cal_time is not None:
        for time in cal:
            if solo.time[-1] < time:
                print(f"clock cal time {time} is past end of time series")
            if solo.time[0] > time:
                print(f"clock cal time {time} is before start of time series")

    # Register how many cals to show.
    if cal_time is None:
        ncal = 0
        show_cal = False
    else:
        ncal = len(cal)
        show_cal = True

    # Set up figure.
    if show_cal:
        if ncal == 1:
            fig, [ax0, axcal] = plt.subplots(
                nrows=2, ncols=1, figsize=(10, 7), constrained_layout=True
            )
            axcal = [axcal]
        elif ncal == 2:
            fig = plt.figure(
                constrained_layout=True,
                figsize=(10, 7),
            )
            gs = fig.add_gridspec(2, 2)
            ax0 = fig.add_subplot(gs[0, :])
            ax1 = fig.add_subplot(gs[1, 0])
            ax2 = fig.add_subplot(gs[1, 1])
            axcal = [ax1, ax2]
    else:
        fig, ax0 = plt.subplots(
            nrows=1, ncols=1, figsize=(10, 4), constrained_layout=True
        )

    # Plot time series. Coarsen if it is too long to slow things down.
    if len(solo) > 1e5:
        coarsen_by = int(np.floor(60 / solo.attrs["sampling period in s"]))
        solo.coarsen(time=coarsen_by, boundary="trim").mean().plot(ax=ax0)
    else:
        solo.plot(ax=ax0)

    # Plot a warning if time offset not applied.
    if solo.attrs["time offset applied"] == 1:
        ax0.text(
            0.05,
            0.9,
            "time offset of {} seconds applied".format(
                solo.attrs["time drift in ms"] / 1000
            ),
            transform=ax0.transAxes,
            backgroundcolor="w",
        )
    else:
        if solo.attrs["time drift in ms"] == 0:
            ax0.text(
                0.05,
                0.9,
                "WARNING: time offset unknown",
                transform=ax0.transAxes,
                color="red",
                backgroundcolor="w",
            )
        elif np.absolute(solo.attrs["time drift in ms"]) > 3.6e6:
            ax0.text(
                0.05,
                0.9,
                "WARNING: time offset more than one hour, not applied",
                transform=ax0.transAxes,
                color="red",
                backgroundcolor="w",
            )
        else:
            ax0.text(
                0.05,
                0.9,
                "time offset not yet applied",
                transform=ax0.transAxes,
                backgroundcolor="w",
            )

    ax0.grid()
    ax0.set(title="SBE56 SN {}".format(solo.attrs["SN"]))
    ax0.set(xlabel="")
    gv.plot.concise_date(ax0)

    # Plot calibration.
    if show_cal:
        # determine time span of plot based on sampling period
        sampling_period = solo.attrs['sampling period in s']
        width = int(sampling_period * 20)
        for cal_time, axi in zip(cal, axcal):
            tmp = solo.sel(
                time=slice(
                    cal_time - np.timedelta64(width, "s"),
                    cal_time + np.timedelta64(width, "s"),
                )
            )
            if len(tmp) > 0:
                tmp.plot(ax=axi, marker=".")
                ylims = np.array(
                    [np.floor(tmp.min().data), np.ceil(tmp.max().data)]
                )
            else:
                ylims = np.array([1, 9])
            axi.plot(
                np.tile(cal_time, 2),
                ylims + np.array([1, -1]),
                linestyle="-",
                color="darkorchid",
                linewidth=1.5,
            )
            axi.annotate(
                "time calibration",
                (cal_time, ylims[0] + 0.5),
                xytext=(16, 8),
                textcoords="offset points",
                color="darkorchid",
                ha="left",
                backgroundcolor="w",
            )
            axi.set(
                xlim=[
                    cal_time - np.timedelta64(width, "s"),
                    cal_time + np.timedelta64(width, "s"),
                ],
                ylim=ylims,
                xlabel="",
            )
            axi.grid()
            gv.plot.concise_date(axi)
        if ncal == 2:
            axcal[-1].set(ylabel="")

    if figure_out is not None or False:
        figurename = "{:s}.png".format(solo.attrs["file"][:-4])
        plt.savefig(figure_out.joinpath(figurename), dpi=300)


def plot_old(solo, figure_out=None, cal_time=None):

    # check if cal_time is past end of time series
    if cal_time is not None:
        if solo.time[-1] < cal_time:
            show_cal = False
            print("clock cal time is past end of time series, not plotting")
        else:
            show_cal = True
    else:
        show_cal = False

    # set up figure
    if show_cal:
        fig, [ax0, ax1] = plt.subplots(
            nrows=2, ncols=1, figsize=(10, 7), constrained_layout=True
        )
    else:
        fig, ax0 = plt.subplots(
            nrows=1, ncols=1, figsize=(10, 4), constrained_layout=True
        )

    # plot time series. coarsen if it is too long to slow things down
    if len(solo) > 1e5:
        coarsen_by = int(np.floor(60 / solo.attrs["sampling period"]))
        solo.coarsen(time=coarsen_by, boundary="trim").mean().plot(ax=ax0)
    else:
        solo.plot(ax=ax0)
    # plot a warning if time offset not applied
    if solo.attrs["time offset applied"] == 1:
        ax0.text(
            0.05,
            0.9,
            "time offset of {} seconds applied".format(
                solo.attrs["time drift in ms"] / 1000
            ),
            transform=ax0.transAxes,
            backgroundcolor="w",
        )
    else:
        if solo.attrs["time drift in ms"] == 0:
            ax0.text(
                0.05,
                0.9,
                "WARNING: time offset unknown",
                transform=ax0.transAxes,
                color="red",
                backgroundcolor="w",
            )
        elif np.absolute(solo.attrs["time drift in ms"]) > 3.6e6:
            ax0.text(
                0.05,
                0.9,
                "WARNING: time offset more than one hour, not applied",
                transform=ax0.transAxes,
                color="red",
                backgroundcolor="w",
            )
        else:
            ax0.text(
                0.05,
                0.9,
                "time offset not yet applied",
                transform=ax0.transAxes,
                backgroundcolor="w",
            )

    ax0.grid()
    ax0.set(title="SBE56 SN {}".format(solo.attrs["SN"]))
    ax0.set(xlabel="")
    gv.plot.concise_date(ax0)

    # plot calibration
    if show_cal:
        tmp = solo.sel(
            time=slice(
                cal_time - np.timedelta64(60, "s"),
                cal_time + np.timedelta64(60, "s"),
            )
        )
        if len(tmp) > 0:
            tmp.plot(ax=ax1, marker=".")
            ylims = np.array(
                [np.floor(tmp.min().data), np.ceil(tmp.max().data)]
            )
        else:
            ylims = np.array([1, 9])
        ax1.plot(
            np.tile(cal_time, 2),
            ylims + np.array([1, -1]),
            linestyle="-",
            color="darkorchid",
            linewidth=1.5,
        )
        ax1.annotate(
            "time calibration",
            (cal_time, ylims[0] + 0.5),
            xytext=(8, 8),
            textcoords="offset points",
            color="darkorchid",
            ha="left",
            backgroundcolor="w",
        )
        ax1.set(
            xlim=[
                cal_time - np.timedelta64(60, "s"),
                cal_time + np.timedelta64(60, "s"),
            ],
            ylim=ylims,
            xlabel="",
        )
        ax1.grid()
        gv.plot.concise_date(ax1)

    if figure_out is not None or False:
        figurename = "{:s}.png".format(solo.attrs["file"][:-4])
        plt.savefig(figure_out.joinpath(figurename), facecolor="w", dpi=300)
