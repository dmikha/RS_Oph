# version of pp to use
pp_version="pp_examples_Kelner"
#pp_version="pp_examples"

# observational data
array HESSII[5]
HESSII[1]="observations/sed_mono_090821_low.dat"
HESSII[2]="observations/sed_mono_100821_low.dat"
HESSII[3]="observations/sed_mono_110821_low.dat"
HESSII[4]="observations/sed_mono_120821_low.dat"
HESSII[5]="observations/sed_mono_130821_low_noshow.dat"

array HESS[5]
HESS[1]="observations/sed_stereo_090821_bias2p5.dat"
HESS[2]="observations/sed_stereo_100821_bias2p5.dat"
HESS[3]="observations/sed_stereo_110821_bias2p5.dat"
HESS[4]="observations/sed_stereo_120821_bias2p5.dat"
HESS[5]="observations/sed_stereo_130821_bias2p5.dat"

array LAT[5]
LAT[1]="observations/sed_LAT_090821_HC_REBIN.dat"
LAT[2]="observations/sed_LAT_100821_HC.dat"
LAT[3]="observations/sed_LAT_110821_HC.dat"
LAT[4]="observations/sed_LAT_120821_HC.dat"
LAT[5]="observations/sed_LAT_130821_HC_REBIN.dat"

#distance
d_to_earth=real(system("less myparameters.f90 | grep 'D_to_Earth =' | sed 's/^.*=//' | sed 's/_fp.*//'"))

factor = d_to_earth**2 * 1.1964951828635065e+38

#factor = factor / 3.


#extra functions
min(xx,yy) = (xx < yy) ? xx : yy
max(xx,yy) = (xx > yy) ? xx : yy
mcut(x,y)  = (x  > y ) ? 1  : 1/0



# plotting the data for day i
# HESS mono
# HESS stereo
# Fermi data
set macros
load 'add_data.cfg'
load 'multi_pannel.cfg'

# best plot range calculations
set term dumb
#electrons
plot \
   for [i=6:2:-1] "dat/syn_examples.dat" u (5.11e-4*$1):(8.187e-7*column(i)*$1**2/factor) w l ls i+1  notitle,\
   for [i=6:2:-1]  "dat/ic_examples.dat" u (5.11e-4*$1):(8.187e-7*column(i)*$1**2/factor) w l ls i+1  dt (5,5) notitle

data_max=GPVAL_DATA_Y_MAX

el_y_max=3.*data_max
el_y_min=data_max/(1.e5)

#protons
plot \
   for [i=6:2:-1] "dat/".pp_version.".dat" u (0.938e0*$1):(0.0015*column(i)*$1**2/factor) w l ls i+1  notitle

data_max=GPVAL_DATA_Y_MAX

pr_y_max=3.*data_max
pr_y_min=data_max/(1.e2)



set term postscript eps color colortext size 15cm,10cm  enhanced font "Helvetica" 16 linewidth 1
set bars small
set logscale xy
set output "ic_emission.eps"

if (!exists("MP_LEFT"))   MP_LEFT = .1
if (!exists("MP_RIGHT"))  MP_RIGHT = .95
if (!exists("MP_BOTTOM")) MP_BOTTOM = 0.1
if (!exists("MP_TOP"))    MP_TOP = 0.98
if (!exists("MP_GAP"))    MP_GAP = 0.

set multiplot layout 2,3 \
              margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP


set style line 1 lt 1 lc 'blue' lw 3 dt (5,5) 
set style line 2 lt 1 lc 'blue' lw 3
set style line 3 lt 1 lc 'red' lw 1 dt (5,5) 
set style line 4 lt 1 lc 'red' lw 1

set key samplen 6
set key left bottom
set yrange [1.e-13:5e-9]
set xrange [1.e-1:1.e4]

do for [i=1:5] {
unset title

if ( i == 1 ) {
@panel1
}
if ( i == 2 ) {
@panel2
}
if ( i == 3 ) {
@panel3
}
if ( i == 4 ) {
@panel4
}
if ( i == 5 ) {
@panel5
}


set label 1 "night ".i at 1.e2, 1.e-9
show label

plot "dat/ic_examples.dat" u (5.11e-4*$1):(8.187e-7*column(i+1)*$1**2/factor) w l ls 2 notitle,\
     @add_data
}
@panel6

set key spacing 2
set key center bottom
set key font 'Helvetica,20'
plot [1:10] [1:10] NaN ls 2 title "IC emission", @add_data_label

unset multiplot



set output "pp_emission.eps"
set multiplot layout 2,3 \
              margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP
set key samplen 6
set key left bottom
#set yrange [el_y_min*1.e2:el_y_max*3]
set yrange [1.e-13:5e-9]
set xrange [1.e-1:1.e4]

do for [i=1:5] {

if ( i == 1 ) {
@panel1
}
if ( i == 2 ) {
@panel2
}
if ( i == 3 ) {
@panel3
}
if ( i == 4 ) {
@panel4
}
if ( i == 5 ) {
@panel5
}

set label 1 "night ".i at 1.e2, 1.e-9
show label

if ( i == 0 ){
plot "pp_case_high/".pp_version.".dat" u (0.938e0*$1):(0.0015*column(i+1)*$1**2/factor) w l ls 2 lw 0.5 notitle "pp emission", @add_data
     }else{
plot "dat/".pp_version.".dat" u (0.938e0*$1):(0.0015*column(i+1)*$1**2/factor) w l ls 2 notitle,\
     @add_data			   
}
}
@panel6
set key spacing 2
set key center bottom 
set key font 'Helvetica,20'
plot [1:10] [1:10] NaN ls 2 title "p-p emission", @add_data_label

unset multiplot

