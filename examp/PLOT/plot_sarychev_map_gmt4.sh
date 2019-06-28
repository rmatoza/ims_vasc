#!/bin/csh

# GMT4 Script

# clean temp files
rm grid_during.grd fnpix.cpt

# output postscript file
set psfile = sarychev_imsvasc.ps

gmtset BASEMAP_TYPE PLAIN
gmtset FRAME_PEN 2p
gmtset TICK_LENGTH -0.2c
gmtset ANOT_FONT_SIZE 20
gmtset HEADER_FONT_SIZE 25
gmtset ANOT_FONT_PRIMARY 1
gmtset LABEL_FONT_SIZE 14
gmtset PAPER_MEDIA Custom_14.5ix7.5i
gmtset BASEMAP_FRAME_RGB = +255/255/255
gmtset BASEMAP_FRAME_RGB = 0/0/0

# set range
set lonmin = 120.
set lonmax = 170.
set latmin = 30.
set latmax = 60.
set dl = 0.1

set XYRANGE = -R${lonmin}/${lonmax}/${latmin}/${latmax}
set PROJ = -JM145/45/6i

# set color palette
# get max value
set gridmax = `sort grid_during.txt -k3 | tail -n1 | awk '{print $3}'`
#set gridmax = 30000 # hardwire some max value instead
echo Max grid value number of pixels: $gridmax
makecpt -Cjet -T0/${gridmax}/1000 -Z -D > fnpix.cpt

# create GMT .grd file
awk '{print $1, $2, $3}' grid_during.txt | xyz2grd -: -V $XYRANGE -I0.1 -Ggrid_during.grd


# plot during
grdimage grid_during.grd $PROJ $XYRANGE -Cfnpix.cpt -P -V -K -Y3 >! $psfile
psbasemap -L165/34/34/250 -Bf5a10/f4a8:."During grid G@+d@+":NseW $PROJ -R -V -O -K >> $psfile
gmtset BASEMAP_FRAME_RGB = +0/0/0
psscale -D5i/-0.2i/2i/0.125ih -Cfnpix.cpt -Bf5e3a3e4:"Grid value (number of pixels)": --ANNOT_FONT_SIZE_PRIMARY=11p -O -K >> $psfile
gmtset BASEMAP_FRAME_RGB = +255/255/255
gmtset BASEMAP_FRAME_RGB = 0/0/0
psxy coastline-50m.dat $XYRANGE $PROJ -W1p,255/255/255  -m -O -K >> $psfile
cat ../IN/stationlist.txt | awk '{print $1, $2}' | psxy -: -Si0.6 -G255/0/0 -W1p -R -J -O -K >> $psfile
psxy sarychev.loc -: -St0.3 -W1.5p,255/255/0 -R -J -O -K >> $psfile
cat ../OUT/OUT_ICAT/icat* | awk '{print $7, $8}' | psxy -: -Sc0.3 -W1.5p,255/255/255 -R -J -O -K >> $psfile


# plot cleaned
# create GMT .grd file
awk '{print $1, $2, $4}' grid_proc.txt | xyz2grd -: -V $XYRANGE -I0.1 -Ggrid_clean.grd
grdimage grid_clean.grd $PROJ $XYRANGE -Cfnpix.cpt -P -V -O -K -X6.8i >> $psfile
psbasemap -L165/34/34/250 -Bf5a10/f4a8:."Cleaned grid G@+c@+":NseW $PROJ -R -V -O -K >> $psfile
gmtset BASEMAP_FRAME_RGB = +0/0/0
psscale -D5i/-0.2i/2i/0.125ih -Cfnpix.cpt -Bf5e3a3e4:"Grid value (number of pixels)": --ANNOT_FONT_SIZE_PRIMARY=11p -O -K >> $psfile
gmtset BASEMAP_FRAME_RGB = +255/255/255
gmtset BASEMAP_FRAME_RGB = 0/0/0
psxy coastline-50m.dat $XYRANGE $PROJ -W1p,255/255/255  -m -O -K >> $psfile
cat ../IN/stationlist.txt | awk '{print $1, $2}' | psxy -: -Si0.6 -G255/0/0 -W1p -R -J -O -K >> $psfile
psxy sarychev.loc -: -St0.3 -W1.5p,255/255/0 -R -J -O -K >> $psfile
cat ../OUT/OUT_ICAT/icat* | awk '{print $7, $8}' | psxy -: -Sc0.3 -W1.5p,255/255/255 -R -J -O >> $psfile


