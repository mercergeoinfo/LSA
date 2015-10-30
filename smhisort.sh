#!/bin/bash
# This script downloads an SMHI daily average temperature file and then sorts it into yearly files extracting only the date and temperature columns
#
# Get data file from SMHI
curl -C - -o data.csv http://opendata-download-metobs.smhi.se/api/version/1.0/parameter/2/station/178970/period/corrected-archive/data.csv
#
# Loop through data file and create subsets
for i in `seq 1995 2014`; do
    awk -v OFS="," -F";" 'BEGIN {print "Date,Temp"} $3 ~ /^'${i}'/ {print $3,$4}' data.csv > weather${i}.csv
done