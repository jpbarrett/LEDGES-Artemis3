# LEDGES-Artemis3
Quick and dirty plotting of Earth/Sun elevation vs. time for LEDGES

Dependencies:
numpy
matplotlib
pandas
astroquery

From the main repository directory, run:
source ./setup_env.sh

Then to run the Earth, Sun elevation calculation and make the associated plots run:
./EarthSunElevation.py <options>

The options are given by ./EarthSunElevation.py -h and are:

usage: EarthSunElevation.py [-h] -s START_DATE -e END_DATE [-b BATTERY_DAYS] [-o OUTPUT_DIRECTORY] [-p SAMPLING_PERIOD] [-d DB_HORIZON] [-a DB_ANTENNA]

required arguments:
 -s START_DATE, --state-date START_DATE
                       start date of time period of interest in YYYY-MM-DD format
 -e END_DATE, --end-date END_DATE
                       end date of time period of interest in YYYY-MM-DD format

optional arguments:
 -b BATTERY_DAYS, --battery-days BATTERY_DAYS
                       time period (days) that batteries will support observing operations, default 9999 (no limit) 
 -o OUTPUT_DIRECTORY, --output-dir OUTPUT_DIRECTORY
                       output directory, default is ./
 -p SAMPLING_PERIOD, --period SAMPLING_PERIOD
                       sampling period in days, default is 1
 -d DB_HORIZON, --dB-horizon DB_HORIZON
                       the required noise suppression by moon-earth geometry, default -45dB
 -a DB_ANTENNA, --dB-antenna DB_ANTENNA
                       the required noise suppression by antenna orientation, default -25dBi
