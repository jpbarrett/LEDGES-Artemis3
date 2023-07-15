#!/usr/bin/env python3
"""plotting Earth elevation"""

from __future__ import print_function
from builtins import range

import sys
import os
import argparse
from datetime import datetime
import math
import copy

from astroquery.jplhorizons import Horizons

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import FormatStrFormatter, ScalarFormatter
from matplotlib import gridspec

from calc_width import calc_width

def intrp_width(dB, df):
    dB_list =  list(df["suppression (dB)"])
    cone_widths = list(df["cone_width (deg)"])
    cone_err_upper = list(df["err_upper (deg)"])
    cone_err_lower = list(df["err_lower (deg)"])

    cw_val = np.interp(dB, dB_list, cone_widths)
    err_plus_val = np.interp(dB, dB_list, cone_err_upper)
    err_minus_val = np.interp(dB, dB_list, cone_err_upper)

    return cw_val, err_plus_val, err_minus_val


def plot_time_series(earth_table, sun_table, current_site, horizon_dB, antenna_dB, odir):
    """ plot as a time series """

    m2km = 1.0/1000.0
    MHz2KHz = 1000
    arcsec2Deg = 1.0/3600.0
    low_freq = 50*MHz2KHz
    high_freq = 100*MHz2KHz;
    current_site_height = max(0,current_site[3])*m2km ;
    total_dB = horizon_dB + antenna_dB

    freq_MHz = 50; #diffraction is worse at lower frequencies, so use low end of band, 50MHz
    freq_kHZ = freq_MHz*MHz2KHz

    #calculate cone_width, upper/lower error for this horizon_dB
    #cone_width, cone_err_upper, cone_err_lower = intrp_width(horizon_dB, dB_el_dataframe)
    print("SITE HEIGHT IN KM = ", current_site_height)
    cone_width, cone_err_upper, cone_err_lower = calc_width(freq_kHZ, current_site_height, horizon_dB)
    print(cone_width, cone_err_upper, cone_err_lower)

    #elevation limit
    elev_limit = (cone_width - 180.0)/2.0

    #get date string of start
    start_day_str = earth_table['datetime_str'][0]

    #measure time in JD since start
    x = list(earth_table['datetime_jd'])
    start_day = x[0]
    x -= start_day

    #earth and sun elevation
    y = list(earth_table['EL'])
    z = list(sun_table['EL'])

    #angle subtended by earth/sun at observation site in arcsec
    earth_ang_width = list(earth_table['ang_width'])
    sun_ang_width = list(sun_table['ang_width'])

    #earth locations
    yp = list()
    ym = list()

    #sun locations
    zp = list()
    zm = list()

    #elev limit locations
    lp = list()
    lm = list()

    for i in range(0, len(y)):
        yp.append(y[i] + (earth_ang_width[i]*arcsec2Deg)/2.0)
        ym.append(y[i] - (earth_ang_width[i]*arcsec2Deg)/2.0)
        lp.append(elev_limit + cone_err_upper)
        lm.append(elev_limit - cone_err_lower)

    for i in range(0, len(z)):
        zp.append(z[i] + (sun_ang_width[i]*arcsec2Deg)/2.0)
        zm.append(z[i] - (sun_ang_width[i]*arcsec2Deg)/2.0)


    #plot data
    fig = plt.figure(figsize=(11,8))
    plt.plot(x, yp, color='k')
    plt.plot(x, ym, color='k')
    plt.plot(x, zp, color='k')
    plt.plot(x, zm, color='k')

    #draw the necessary elevation limit with error
    lower_est = elev_limit - cone_err_lower
    upper_est = elev_limit + cone_err_upper

    # plt.axhline(y = elev_limit, color = 'r', linestyle = '-', label="Naive -70dB elevation @50MHz")
    plt.hlines(y = upper_est, xmin=min(x), xmax=max(x), color = 'k', linestyle = '-')
    plt.hlines(y = lower_est, xmin=min(x), xmax=max(x), color = 'k', linestyle = '-')

    ax  = plt.gca()
    ax.fill_between(x, ym, yp, label="Earth", color='g')
    ax.fill_between(x, zm, zp, label="Sun", color='b')
    ax.fill_between(x, lm, lp, label="Naive " + str(total_dB) + "dB elevation @50MHz", color='r')

    #some statistics about the Earth/Sun position
    total_duration = max(x) - min(x) #time duration in days
    earth_delta = upper_est - yp
    sun_delta = zp

    earth_frac_below = 0.0;
    sun_frac_below = 0.0
    both_frac_below = 0.0
    earth_longest_duration_below = 0.0
    sun_longest_duration_below = 0.0
    both_longest_duration_below = 0.0
    earth_duration_below = 0.0
    sun_duration_below = 0.0
    both_duration_below = 0.0

    for i in range(1,len(x)):
        if earth_delta[i] > 0.0:
            earth_frac_below += (x[i] - x[i-1])
            earth_duration_below += (x[i] - x[i-1])
            if earth_duration_below > earth_longest_duration_below:
                earth_longest_duration_below = earth_duration_below
        else:
            earth_duration_below = 0.0

        if sun_delta[i] < 0.0:
            sun_frac_below += (x[i] - x[i-1])
            sun_duration_below += (x[i] - x[i-1])
            if sun_duration_below > sun_longest_duration_below:
                sun_longest_duration_below = sun_duration_below
        else:
            sun_duration_below = 0.0

        if earth_delta[i] > 0.0 and sun_delta[i] < 0.0:
            both_frac_below += (x[i] - x[i-1])
            both_duration_below += (x[i] - x[i-1])
            if both_duration_below > both_longest_duration_below:
                both_longest_duration_below = both_duration_below
        else:
            both_duration_below = 0.0

    earth_frac_below /= total_duration
    sun_frac_below /= total_duration
    both_frac_below /= total_duration

    line1 = "Percent of time Earth fully below the elevation limit: " + str( round( earth_frac_below*100, 2) ) + "% \n"
    line2 = "Percent of time Sun fully below the horizon: "+ str( round( sun_frac_below*100, 2) ) + "% \n"
    line3 = "Percent of time both Sun and Earth fully below the limits: "+ str( round( both_frac_below*100, 2) ) + "% \n"
    line4 = "Earth longest period below the elevation limit (days): " + str(earth_longest_duration_below) + "\n"
    line5 = "Sun longest period below the horizon (days): " + str(sun_longest_duration_below) + "\n"
    line6 = "Both Sun and Earth longest period below limits (days): " + str(both_longest_duration_below) + "\n"
    info_text = line1+line2+line3+line4+line5+line6
    print(info_text)

    #add some information, and adjust plot
    plt.grid(True)
    plt.title('Earth/Sun elevation for site: ' + current_site[0] + "    (Lat, Lon, Height) = (" + str(current_site[1]) + "," + str(current_site[2]) + "," + str(current_site[3]) + ")" )
    plt.xlabel('Days since '+ start_day_str)
    plt.ylabel('Elevation (deg)')
    plt.legend()
    outfile = os.path.join( odir, current_site[0] + ".png")
    fig.tight_layout(rect=[0,0.16,1.0,1.0])
    fig.text(0.5, -0.01, info_text, horizontalalignment='center', wrap=True )

    plt.savefig(outfile)
    plt.close()

    site_data = { \
        'site_name' :current_site[0], \
        'earth_frac_below': round( earth_frac_below*100, 2), \
        'sun_frac_below' :  round( sun_frac_below*100, 2), \
        'earth_max_duration_below': earth_longest_duration_below, \
        'sun_max_duration_below': sun_longest_duration_below, \
        'both_max_duration_below' : both_longest_duration_below \
        }

    print(site_data)



    #next we plot the Earth elevation minus the elevaton limit vs the sun elevation 
    #plot data
    fig2 = plt.figure(figsize=(11,8))
    earth_delta = -1*earth_delta
    plt.plot(earth_delta, zp, color='k')
    plt.grid(True)
    plt.title('Earth Elevation relative to RFI elevation limit vs Sun Elevation' )
    plt.xlabel('Earth Elevation relative to RFI elevation limit (deg)')
    plt.ylabel('Sun Elevation (deg)')
    outfile = os.path.join( odir, current_site[0] + "_conops.png")


    ax  = plt.gca()
    earth_cutoff = 0.0
    sun_cutoff = 0.0
    ax.axvline(earth_cutoff,color = 'black',linestyle='dashed',lw=2)
    ax.axhline(sun_cutoff,color = 'black',linestyle='dashed',lw=2)
    if min(earth_delta) < earth_cutoff:
        ax.fill_between([min(earth_delta), earth_cutoff], min(zp), sun_cutoff, alpha=0.3, color='#00FF00', label="Observe/Warm")  #green '#1F98D0'blue
    if earth_cutoff <  max(earth_delta):
        ax.fill_between([earth_cutoff, max(earth_delta)], min(zp), sun_cutoff, alpha=0.3, color='#F9D307', label="Keep Warm")  # yellow
    if min(earth_delta) < earth_cutoff:
        ax.fill_between([min(earth_delta), earth_cutoff], sun_cutoff, max(zp), alpha=0.3, color='#F38D25', label="Charge/Observe" )  # orange
    if earth_cutoff < max(earth_delta):
        ax.fill_between([earth_cutoff, max(earth_delta)], sun_cutoff, max(zp), alpha=0.3, color='#DA383D', label="Charge/Uplink")  # red

    plt.legend()

    plt.savefig(outfile)
    plt.close()

    #calculated total integration time
    accumulated_integration_time = 0.0
    time_since_last_charge = 0.0
    int_time_array = list()
    for i in range(0, len(x)):
        int_time_array.append(0.0)

    for i in range(1,len(x)):
        if earth_delta[i] < 0.0 and zp[i] < 0.0 and time_since_last_charge < 2.0:
            accumulated_integration_time += (x[i] - x[i-1])
            time_since_last_charge += (x[i] - x[i-1])
        int_time_array[i] = accumulated_integration_time
        if zp[i] > 0.0:
            time_since_last_charge = 0.0;

    return site_data, x, int_time_array

    # print(site_data)
    # site_df = pd.DataFrame(site_data, columns=['site_name', 'earth_frac_below', 'sun_frac_below', 'earth_max_duration_below', "sun_max_duration_below", "both_max_duration_below"])
    # return site_df


#===== Main ==============================================

def main():
    """
    Plots the Earth/Sun elevation at the Artemis sites for some time period.
    Does not take any kind of local topography into account (assumes purely spherical Moon)
    Site positions are near but not exactly the center of the selected landing zones.
    """

    parser = argparse.ArgumentParser()
    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    #required.add_argument("-c", "--cone-width-data", dest='data_file', help="dB vs cone angle data table .csv file.", required=True)
    required.add_argument("-s", "--state-date", dest='start_date', help="start date of time period of interest in YYYY-MM-DD format", required=True)
    required.add_argument("-e", "--end-date", dest='end_date', help="end date of time period of interest in YYYY-MM-DD format", required=True)
    optional.add_argument("-o", "--output-dir", dest='output_directory', default="./", help="output directory, default is ./") #optional argument
    optional.add_argument("-p", "--period", dest='sampling_period', default="1", help="sampling period in days, default is 1") #optional argument
    optional.add_argument("-d", "--dB-horizon", dest="dB_horizon", default="-45", help="the required noise suppression by moon-earth geometry, default -45dB")
    optional.add_argument("-a", "--dB-antenna", dest="dB_antenna", default="-25", help="the required noise suppression by antenna orientation, default -25dBi")

    args = parser.parse_args()

    # data_file = os.path.abspath(args.data_file)

    #TODO validate date format
    start_date = args.start_date
    end_date = args.end_date

    odir = os.path.abspath(args.output_directory)
    period = args.sampling_period
    step = str(period) + 'd'
    dB_horizon = float(args.dB_horizon)
    dB_antenna = float(args.dB_antenna)


    print(start_date)
    print(end_date)
    print(odir)
    print(period)
    print(dB_horizon)
    print(dB_antenna)

    # #read data file generated with Neil Bassett's calc_width function
    # dB_el_df = pd.read_csv(data_file)
    # 
    # print(dB_el_df)

    #HEIGHTS extracted using https://quickmap.lroc.asu.edu/ 
    #error is probably 50-100m

    #all Artemis-3 site options
    #NAME, LAT, LON, HEIGHT (m)
    site_table = [
        ["Faustini_Rim_A",  -87.9403, 88.97176, -150.0],
        ["Peak_Near_Shackleton",  -88.80072, 123.12689, 1650.0],
        ["Connecting_Ridge", -89.44208, 219.38414, 1850.0],
        ["Connecting_Ridge_Extension", -89.00059, 258.81705, 1350.0 ],
        ["deGerlache_Rim_1", -88.66814, 289.25406, 1300.0 ],
        ["deGerlache_Rim_2", -88.27753, 295.47124, 200.0],
        ["deGerlache_Kocher_Massif", -85.78549, 243.59109, 2750.0],
        ["Haworth", -86.78431, 337.29546, 1650.0],
        ["Malapert_Massif", -85.97930, 0.06820, 4800.0],
        ["Leibnitz_Beta_Plateau", -85.42375, 31.97821, 6400.0],
        ["Nobile_Rim_1", -85.46393, 37.36865, 6050.0],
        ["Nobile_Rim_2", -83.97494, 58.85641, 1050.0],
        ["Amundsen_Rim", -84.21276, 69.68302, 450.0]
        ]
    # 
    # site_table = [
    #     ["Faustini_Rim_A",  -87.9403, 88.97176, -150.0]
    #     ]



    #JPL Horizons object ids
    sun_id = 10
    earth_id = 399
    moon_id = 301

    #constants
    m2km = 1.0/1000.0
    MHz2KHz = 1000.0
    freq_MHz = 50; #diffraction is worse at lower frequencies, so use low end of band, 50MHz
    freq_kHZ = freq_MHz*MHz2KHz
    height = 0.010; #assume and elevation of 10m (TODO use LRO digital elevation model)

    site_df = pd.DataFrame(columns=['site_name', 'earth_frac_below', 'sun_frac_below', 'earth_max_duration_below', "sun_max_duration_below", "both_max_duration_below"], dtype=object)

    time_ax = list()
    int_time_dict = {}
    for current_site in site_table:
        print(current_site[0], current_site[1], current_site[2], current_site[3])
        current_site_name = current_site[0]
        current_site_height = current_site[3]*m2km
        current_site_location = {'lon': current_site[2], 'lat': current_site[1], 'elevation': current_site_height, 'body': moon_id}
        earth_obj = Horizons(id=earth_id, location=current_site_location, epochs={'start':start_date, 'stop':end_date,'step':step}, )
        sun_obj = Horizons(id=sun_id, location=current_site_location, epochs={'start':start_date, 'stop':end_date,'step':step}, )
        earth_eph = earth_obj.ephemerides()
        sun_eph = sun_obj.ephemerides()
        print(earth_eph)
        tmp_df, x, int_time = plot_time_series( earth_eph, sun_eph, current_site, dB_horizon, dB_antenna, odir)
        time_ax = x
        int_time_dict[current_site_name] = int_time
        # site_df = site_df.append(tmp_df, ignore_index=True)

    site_df.to_csv(os.path.join(odir, "site_illumination.csv") )

    fig3 = plt.figure(figsize=(11,8))
    for site in int_time_dict:
        #next we plot the Earth elevation minus the elevaton limit vs the sun elevation 
        #plot data
        if int_time_dict[site][-1] > 0.0:
            plt.plot(time_ax, int_time_dict[site], label=site)

    plt.grid(True)
    plt.title('Integration time vs time' )
    plt.xlabel('Days since start')
    plt.ylabel('Max possible observation time.')
    plt.legend()
    outfile = os.path.join( odir, "site_integration.png")
    plt.savefig(outfile)
    plt.close()



#===== Official entry point ==============================

if __name__ == '__main__':
    main()
    sys.exit(0)
