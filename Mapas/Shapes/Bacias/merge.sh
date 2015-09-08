#!/bin/bash

file= "merge.shp"

for i in $(ls *.shp)
do

      if [ -f "$file" ]
      then
           echo "creating final/merge.shp"
           ogr2ogr -f "ESRI Shapefile" -update -append $file $i -nln merge
      else
           echo "merging……"
      ogr2ogr -f "ESRI Shapefile" $file $i
fi

done
mkdir merged;

for f in *.shp;
do
    if [ -f merged/merged.shp ]
    then
        ogr2ogr -f "ESRI Shapefile" -update -append merged/merged.shp "$f" -nln Merged
    else
        ogr2ogr -f "ESRI Shapefile" merged/merged.shp "$f"
    fi;
done;


