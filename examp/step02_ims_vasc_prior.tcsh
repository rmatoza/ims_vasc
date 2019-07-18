#!/bin/tcsh

../bin/ims_vasc << DONE
IN/stationlist.txt
30. 60. 0.1
120. 170. 0.1
5.
0.33
5000.
0.1 5.0 2.0
500
0
IN/bullfilelist_prior.txt
2009 06 01 00 00 0.0
2009 06 11 00 00 0.0
OUT/OUT_PRIOR/out.pri.2009.0601.0611
DONE 

