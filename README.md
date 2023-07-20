# LEDGES-Artemis3
Quick and dirty plotting of Earth/Sun elevation vs. time for LEDGES

Dependencies:
* numpy 
* matplotlib
* astroquery 

From the main repository directory, run: source ./setup_env.sh \
Then to run the Earth, Sun elevation calculation and make the associated plots run: ./EarthSunElevation.py \<options\>

The options are given by ./EarthSunElevation.py -h and are as follows:

required arguments:
* \-s START_DATE, --start-date START_DATE
	+ start date of time period of interest in YYYY-MM-DD format
* \-e END_DATE, --end-date END_DATE
	+ end date of time period of interest in YYYY-MM-DD format

optional arguments:
*  -b BATTERY_DAYS, --battery-days BATTERY_DAYS
	+ time period (days) that batteries will support observing operations, default 9999 (no limit)
* \-o OUTPUT_DIRECTORY, --output-dir OUTPUT_DIRECTORY
	+ output directory, default is ./ 
* \-p SAMPLING_PERIOD, --period SAMPLING_PERIOD
	+ sampling period in days, default is 1 
* \-d DB_HORIZON, --dB-horizon DB_HORIZON
	+ the required noise suppression by moon-earth geometry, default -45 (dB) 
* \-a DB_ANTENNA, --dB-antenna DB_ANTENNA
	+ the required noise suppression by antenna orientation, default -25 (dBi) 
* \-L LOCATION_FILE, --location-file LOCATION_FILE
	+ json file containing site location and elevation data

Additional sites can be checked by creating a new json location file with the same format as artemis3_locations.json. Each site must be specified as a list entry containing: [name, lat, lon, and elevation], with coordinates in decimal degrees and elevation in meters. 
