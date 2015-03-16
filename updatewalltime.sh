#!/bin/sh
                          
for file in $(ls 131118*.jdf)
do
sed -e "s/#PBS -l walltime=1:00:00/#PBS -l walltime=6:00:00/ig" $file > /tmp/tempfile.tmp
mv /tmp/tempfile.tmp $file
done
