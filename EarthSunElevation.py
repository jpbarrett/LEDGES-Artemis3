#!/usr/bin/env python3
"""Quick and dirty plotting of Earth/Sun elevation at Artemis-3 landing regions"""

import sys
import os
import argparse
import math
import json

#JPL horizons for ephemerides
from astroquery.jplhorizons import Horizons

import numpy as np
import matplotlib.pyplot as plt

LEDGES_DIR = os.environ.get('LEDGES_DIR')

#import Neil Bassett's RFI suppression calculator from local extern folder
if LEDGES_DIR is not None:
    from calc_width import calc_width

class LEDGESConstants(object):
    """minimal set of constants"""
    def __init__(self):
        self.m_to_km = 1.0/1000.0
        self.MHz_to_kHz = 1000
        self.arcsec_to_Deg = 1.0/3600.0
        self.target_freq_MHz = 50; #use low end of the band, since diffraction gets worse at lower frequencies
        self.target_freq_kHz = self.target_freq_MHz*self.MHz_to_kHz
        #JPL Horizons object ids
        self.sun_id = 10
        self.earth_id = 399
        self.moon_id = 301

def earth_sun_statistics(x, earth_delta, sun_delta):

    #some statistics about the Earth/Sun position
    total_duration = max(x) - min(x) #time duration in days

    cumm_earth_below = 0.0;
    cumm_sun_below = 0.0
    cumm_both_below = 0.0
    earth_longest_duration_below = 0.0
    sun_longest_duration_below = 0.0
    both_longest_duration_below = 0.0
    earth_duration_below = 0.0
    sun_duration_below = 0.0
    both_duration_below = 0.0

    for i in range(1,len(x)):
        if earth_delta[i] > 0.0:
            cumm_earth_below += (x[i] - x[i-1])
            earth_duration_below += (x[i] - x[i-1])
            if earth_duration_below > earth_longest_duration_below:
                earth_longest_duration_below = earth_duration_below
        else:
            earth_duration_below = 0.0

        if sun_delta[i] < 0.0:
            cumm_sun_below += (x[i] - x[i-1])
            sun_duration_below += (x[i] - x[i-1])
            if sun_duration_below > sun_longest_duration_below:
                sun_longest_duration_below = sun_duration_below
        else:
            sun_duration_below = 0.0

        if earth_delta[i] > 0.0 and sun_delta[i] < 0.0:
            cumm_both_below += (x[i] - x[i-1])
            both_duration_below += (x[i] - x[i-1])
            if both_duration_below > both_longest_duration_below:
                both_longest_duration_below = both_duration_below
        else:
            both_duration_below = 0.0
            
    results = dict()
    results["total_duration"] = total_duration
    results["earth_longest_duration_below"] = earth_longest_duration_below
    results["sun_longest_duration_below"] = sun_longest_duration_below
    results["both_longest_duration_below"] = both_longest_duration_below
    results["earth_duration_below"] = cumm_earth_below
    results["sun_duration_below"] = cumm_sun_below
    results["both_duration_below"] = cumm_both_below

    return results

def site_integration_time(x, earth_delta, sun_delta, battery_days):

    #crudely calculate total integration time
    accumulated_integration_time = 0.0
    time_since_last_charge = 0.0
    intg_time_array = list()
    for i in range(0, len(x)):
        intg_time_array.append(0.0)

    for i in range(1,len(x)):
        if earth_delta[i] < 0.0 and sun_delta[i] < 0.0 and time_since_last_charge < battery_days:
            accumulated_integration_time += (x[i] - x[i-1])
            time_since_last_charge += (x[i] - x[i-1])
        intg_time_array[i] = accumulated_integration_time
        if sun_delta[i] > 0.0: #assume time to charge batteries is nil
            time_since_last_charge = 0.0;

    results = dict()
    results["time"] = x 
    results["integration_time"] = intg_time_array

    return results


def plot_earth_sun_time_series(earth_table, sun_table, current_site, cone_info, total_dB, battery_days, odir):
    """ plot as a time series """

    lconst = LEDGESConstants()

    current_site_height = max(0,current_site[3])*lconst.m_to_km; #negative height not allowed

    cone_width = cone_info[0]
    cone_err_upper = cone_info[1]
    cone_err_lower = cone_info[2]
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
        yp.append(y[i] + (earth_ang_width[i]*lconst.arcsec_to_Deg)/2.0)
        ym.append(y[i] - (earth_ang_width[i]*lconst.arcsec_to_Deg)/2.0)
        lp.append(elev_limit + cone_err_upper)
        lm.append(elev_limit - cone_err_lower)

    for i in range(0, len(z)):
        zp.append(z[i] + (sun_ang_width[i]*lconst.arcsec_to_Deg)/2.0)
        zm.append(z[i] - (sun_ang_width[i]*lconst.arcsec_to_Deg)/2.0)

    #plot data
    fig = plt.figure(figsize=(11,8))
    plt.plot(x, yp, color='k')
    plt.plot(x, ym, color='k')
    plt.plot(x, zp, color='k')
    plt.plot(x, zm, color='k')

    #draw the necessary elevation limit with error
    lower_est = elev_limit - cone_err_lower
    upper_est = elev_limit + cone_err_upper

    plt.hlines(y = upper_est, xmin=min(x), xmax=max(x), color = 'r', linestyle = '-')
    plt.hlines(y = lower_est, xmin=min(x), xmax=max(x), color = 'r', linestyle = '-')
    plt.hlines(y = 0, xmin=min(x), xmax=max(x), color = 'k', linestyle = '-') #nominal horizon

    ax  = plt.gca()
    ax.fill_between(x, ym, yp, label="Earth", color='g')
    ax.fill_between(x, zm, zp, label="Sun", color='b')
    ax.fill_between(x, lm, lp, label="Naive " + str(total_dB) + "dB elevation @ " + str(lconst.target_freq_MHz) + " MHz", color='r')

    #site statistics
    earth_el_delta = upper_est - yp
    sun_el_delta = zp
    site_stats = earth_sun_statistics(x, earth_el_delta, sun_el_delta)
    integration_stats = site_integration_time(x, earth_el_delta, sun_el_delta, battery_days)

    earth_frac_below = site_stats["earth_duration_below"]/site_stats["total_duration"]
    sun_frac_below = site_stats["sun_duration_below"]/site_stats["total_duration"]
    both_frac_below = site_stats["both_duration_below"]/site_stats["total_duration"]

    line1 = "Percent of time Earth fully below the elevation limit: " + str( round( earth_frac_below*100, 2) ) + "% \n"
    line2 = "Percent of time Sun fully below the horizon: "+ str( round( sun_frac_below*100, 2) ) + "% \n"
    line3 = "Percent of time both Sun and Earth fully below the limits: "+ str( round( both_frac_below*100, 2) ) + "% \n"
    line4 = "Earth longest period below the elevation limit (days): " + str( site_stats["earth_longest_duration_below"]) + "\n"
    line5 = "Sun longest period below the horizon (days): " + str( site_stats["sun_longest_duration_below"] ) + "\n"
    line6 = "Both Sun and Earth longest period below limits (days): " + str( site_stats["both_longest_duration_below"] ) + "\n"
    info_text = line1+line2+line3+line4+line5+line6

    #add some information, and adjust plot
    plt.grid(True)
    plt.title('Earth/Sun elevation for site: ' + current_site[0] + "  (Lat, Lon, Height) = (" + str(current_site[1]) + "," + str(current_site[2]) + "," + str(current_site[3]) + ")" )
    plt.xlabel('Days since '+ start_day_str)
    plt.ylabel('Elevation (deg)')
    plt.legend()
    outfile = os.path.join( odir, current_site[0] + ".png")
    fig.tight_layout(rect=[0,0.16,1.0,1.0])
    fig.text(0.5, -0.01, info_text, horizontalalignment='center', wrap=True )

    plt.savefig(outfile)
    plt.close()
    
    return site_stats, integration_stats


def plot_integration_time(time_ax, intg_time_dict, battery_days):

    #plot the expected integration time vs time
    fig3 = plt.figure(figsize=(11,8))
    for site in intg_time_dict:
        if intg_time_dict[site][-1] > 0.0: #only plot sites which have non-zero integration time
            plt.plot(time_ax, intg_time_dict[site], label=site)
    
    plt.grid(True)
    plt.title('Total integration time vs time (battery_days = ' + str(battery_days) + ')' )
    plt.xlabel('Time since start (days)')
    plt.ylabel('Integration time (days).')
    plt.legend()
    outfile = os.path.join( odir, "site_integration.png")
    plt.savefig(outfile)
    plt.close()


#===== Main ==============================================

def main():
    """
    Plots the Earth/Sun elevation at the Artemis-3 sites for user specified time period.
    Does not take any kind of local topography into account other than site height (assumes purely spherical Moon)
    Site positions are near but not exactly the center of the selected landing zones.
    """
    
    if LEDGES_DIR is None:
        print("Environment not defined. Please run 'source ./setup_env.sh first'")
        return 1

    parser = argparse.ArgumentParser()
    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    #required.add_argument("-c", "--cone-width-data", dest='data_file', help="dB vs cone angle data table .csv file.", required=True)
    required.add_argument("-s", "--start-date", dest='start_date', help="start date of time period of interest in YYYY-MM-DD format", required=True)
    required.add_argument("-e", "--end-date", dest='end_date', help="end date of time period of interest in YYYY-MM-DD format", required=True)
    optional.add_argument("-b", "--battery-days", dest='battery_days', default=9999, help="time period (days) that batteries will support observing operations, default 9999 (no limit)") 
    optional.add_argument("-o", "--output-dir", dest='output_directory', default="./", help="output directory, default is ./") #optional argument
    optional.add_argument("-p", "--period", dest='sampling_period', default="1", help="sampling period in days, default is 1") #optional argument
    optional.add_argument("-d", "--dB-horizon", dest="dB_horizon", default="-45", help="the required noise suppression by moon-earth geometry, default -45 (dB)")
    optional.add_argument("-a", "--dB-antenna", dest="dB_antenna", default="-25", help="the required noise suppression by antenna orientation, default -25 (dBi)")
    optional.add_argument("-L", "--location-file", dest="location_file", default="./artemis3_locations.json", help="json file containing site location and elevation data")
    args = parser.parse_args()

    #TODO validate date format
    start_date = args.start_date
    end_date = args.end_date
    battery_days = args.battery_days
    odir = os.path.abspath(args.output_directory)
    period = args.sampling_period
    step = str(period) + 'd'
    dB_horizon = float(args.dB_horizon)
    dB_antenna = float(args.dB_antenna) #assume orientation of antenna gives this amount of suppression

    lconst = LEDGESConstants()

    site_list = []
    location_file = os.path.abspath(args.location_file)
    with open(location_file, "r") as locations:
        data = json.load(locations)
        if "sites" in data:
            site_list = data["sites"]

    time_ax = list()
    site_stats_dict = dict()
    intg_time_dict = dict()
    for current_site in site_list:
        print(current_site[0], current_site[1], current_site[2], current_site[3])
        current_site_name = current_site[0]
        current_site_height =  max(0,current_site[3])*lconst.m_to_km;  #negative height not allowed by calc_width function
        current_site_location = {'lon': current_site[2], 'lat': current_site[1], 'elevation': current_site_height, 'body': lconst.moon_id}
        
        #call Neil Bassett's calc_width, for cone_width and upper/lower error for this dB suppression
        cone_info = calc_width(lconst.target_freq_kHz, current_site_height, dB_horizon)

        #call JPL Horizons for earth/sun positions
        earth_obj = Horizons(id=lconst.earth_id, location=current_site_location, epochs={'start':start_date, 'stop':end_date,'step':step}, )
        sun_obj = Horizons(id=lconst.sun_id, location=current_site_location, epochs={'start':start_date, 'stop':end_date,'step':step}, )
        earth_eph = earth_obj.ephemerides()
        sun_eph = sun_obj.ephemerides()

        #TODO, add horizon mask calculated from LOLA DEM
        site_stats, intg_stats = plot_earth_sun_time_series(earth_eph, sun_eph, current_site, cone_info, dB_horizon+dB_antenna, battery_days, odir)

        time_ax = intg_stats["time"]
        intg_time_dict[current_site_name] = intg_stats
        site_stats_dict[current_site_name] = site_stats

    #dump site stats to file
    site_stats_file = os.path.join(odir,"site_illumination.json")
    print("dumping site stats to: ", site_stats_file)
    with open(site_stats_file, 'w') as ssf:
        json.dump(site_stats_dict, ssf)


#===== main entry point ==============================

if __name__ == '__main__':
    main()
    sys.exit(0)
