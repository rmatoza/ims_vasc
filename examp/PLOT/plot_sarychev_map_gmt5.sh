#!/bin/csh

# GMT5 Script based on GMT4 version.
# Created by R. De Negri (rsd00@ucsb.edu) on Fri Jul  5 2019
# Last modification on Sun Jul  7 2019


# clean temp files
rm grid_during.grd fnpix.cpt

# output postscript file
set psfile = sarychev_imsvasc.ps

gmt gmtset MAP_FRAME_TYPE PLAIN
gmt gmtset MAP_FRAME_PEN 2p
gmt gmtset MAP_TICK_LENGTH -0.2c
gmt gmtset FONT_ANNOT_PRIMARY 20p,Helvetica-Bold,black
gmt gmtset FONT_TITLE 25p
gmt gmtset FONT_LABEL 14p
gmt gmtset PS_MEDIA Custom_14.5ix7.5i
gmt gmtset MAP_FRAME_PEN = 0/0/0
gmt gmtset MAP_TICK_PEN = 255/255/255

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
gmt makecpt -Cjet -T0/${gridmax}/1000 -Z -D > fnpix.cpt

# create GMT .grd file
awk '{print $1, $2, $3}' grid_during.txt | gmt xyz2grd -: -V $XYRANGE -I0.1 -Ggrid_during.grd

# plot during
gmt grdimage grid_during.grd $PROJ $XYRANGE -Cfnpix.cpt -E300 -P -V -K -Y3 >! $psfile
gmt psbasemap -BNseW+t"During grid G@+d@+" -Bxf5a10 -Byf4a8 $PROJ -R -V -O -K >> $psfile
gmt psbasemap -Lg165/34+c34+w250 --FONT_ANNOT_PRIMARY=255/255/255 $PROJ -R -V -O -K >> $psfile
gmt psscale -Dx4i/-0.2i+w2i/0.125i+h -Cfnpix.cpt -Bf5e3a3e4+l"Grid value (number of pixels)" --FONT_ANNOT_PRIMARY=11p --MAP_ANNOT_OFFSET_PRIMARY=0.4c --MAP_TICK_PEN=0/0/0 -O -K >> $psfile
gmt psxy coastline-50m.dat $XYRANGE $PROJ -W1p,255/255/255  -O -K >> $psfile
cat ../IN/stationlist.txt | awk '{print $1, $2}' | gmt psxy -: -Si0.6 -G255/0/0 -W1p -R -J -O -K >> $psfile
gmt psxy sarychev.loc -: -St0.3 -W1.5p,255/255/0 -R -J -O -K >> $psfile
cat ../OUT/OUT_ICAT/icat* | awk '{print $7, $8}' | gmt psxy -: -Sc0.3 -W1.5p,255/255/255 -R -J -O -K >> $psfile


# plot cleaned
# create GMT .grd file
awk '{print $1, $2, $4}' grid_proc.txt | gmt xyz2grd -: -V $XYRANGE -I0.1 -Ggrid_clean.grd
gmt grdimage grid_clean.grd $PROJ $XYRANGE -Cfnpix.cpt -E300 -P -V -O -K -X6.8i >> $psfile
gmt psbasemap -BNseW+t"Cleaned grid G@+c@+" -Bxf5a10 -Byf4a8 $PROJ -R -V -O -K >> $psfile
gmt psbasemap -Lg165/34+c34+w250 --FONT_ANNOT_PRIMARY=255/255/255 $PROJ -R -V -O -K >> $psfile
gmt psscale -Dx4i/-0.2i+w2i/0.125i+h -Cfnpix.cpt -Bf5e3a3e4+l"Grid value (number of pixels)" --FONT_ANNOT_PRIMARY=11p --MAP_ANNOT_OFFSET_PRIMARY=0.4c --MAP_TICK_PEN=0/0/0 -O -K >> $psfile
gmt psxy coastline-50m.dat $XYRANGE $PROJ -W1p,255/255/255  -O -K >> $psfile
cat ../IN/stationlist.txt | awk '{print $1, $2}' | gmt psxy -: -Si0.6 -G255/0/0 -W1p -R -J -O -K >> $psfile
gmt psxy sarychev.loc -: -St0.3 -W1.5p,255/255/0 -R -J -O -K >> $psfile
cat ../OUT/OUT_ICAT/icat* | awk '{print $7, $8}' | gmt psxy -: -Sc0.3 -W1.5p,255/255/255 -R -J -O >> $psfile

# Uncomment to create PDF from PS file
#ps2pdf $psfile
