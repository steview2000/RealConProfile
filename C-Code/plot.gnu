#!/usr/bin/gnuplot -persist
#
#    
#    	G N U P L O T
#    	Version 4.6 patchlevel 4    last modified 2013-10-02 
#    	Build System: Linux x86_64
#    
#    	Copyright (C) 1986-1993, 1998, 2004, 2007-2013
#    	Thomas Williams, Colin Kelley and many others
#    
#    	gnuplot home:     http://www.gnuplot.info
#    	faq, bugs, etc:   type "help FAQ"
#    	immediate help:   type "help"  (plot window: hit 'h')
# set terminal x11  nopersist
# set output
unset clip points
set clip one
unset clip two
set bar 1.000000 front
set border 31 front linetype -1 linewidth 1.000
set timefmt z "%d/%m/%y,%H:%M"
set zdata 
set timefmt y "%d/%m/%y,%H:%M"
set ydata 
set timefmt x "%d/%m/%y,%H:%M"
set xdata 
set timefmt cb "%d/%m/%y,%H:%M"
set timefmt y2 "%d/%m/%y,%H:%M"
set y2data 
set timefmt x2 "%d/%m/%y,%H:%M"
set x2data 
set boxwidth
set style fill  empty border
set style rectangle back fc lt -3 fillstyle   solid 1.00 border lt -1
set style circle radius graph 0.02, first 0, 0 
set style ellipse size graph 0.05, 0.03, first 0 angle 0 units xy
set dummy x,y
set format x "% g"
set format y "% g"
set format x2 "% g"
set format y2 "% g"
set format z "% g"
set format cb "% g"
set format r "% g"
set angles radians
unset grid
set raxis
set key title ""
set key inside left top vertical Right noreverse enhanced autotitles nobox
set key noinvert samplen 4 spacing 1 width 0 height 0 
set key maxcolumns 0 maxrows 0
set key noopaque
unset label
unset arrow
set style increment default
unset style line
set style line 1  linetype 1 linecolor rgb "blue"  linewidth 1.000 pointtype 7 pointsize default pointinterval 0
set style line 2  linetype 1 linecolor rgb "red"  linewidth 1.000 pointtype 5 pointsize default pointinterval 0
set style line 3  linetype 1 linecolor rgb "green"  linewidth 1.000 pointtype 13 pointsize default pointinterval 0
set style line 10  linetype 1 linecolor rgb "black"  linewidth 1.000 pointtype 1 pointsize default pointinterval 0
unset style arrow
set style histogram clustered gap 2 title  offset character 0, 0, 0
unset logscale
set offsets 0, 0, 0, 0
set pointsize 1
set pointintervalbox 1
set encoding default
unset polar
unset parametric
unset decimalsign
set view 60, 30, 1, 1
set samples 100, 100
set isosamples 10, 10
set surface
unset contour
set clabel '%8.3g'
set mapping cartesian
set datafile separator "|"
unset hidden3d
set cntrparam order 4
set cntrparam linear
set cntrparam levels auto 5
set cntrparam points 5
set size ratio 0 1,1
set origin 0,0
set style data points
set style function lines
set xzeroaxis linetype -2 linewidth 1.000
set yzeroaxis linetype -2 linewidth 1.000
set zzeroaxis linetype -2 linewidth 1.000
set x2zeroaxis linetype -2 linewidth 1.000
set y2zeroaxis linetype -2 linewidth 1.000
set ticslevel 0.5
set mxtics default
set mytics default
set mztics default
set mx2tics default
set my2tics default
set mcbtics default
set xtics border in scale 1,0.5 mirror norotate  offset character 0, 0, 0 autojustify
set xtics autofreq  norangelimit
set ytics border in scale 1,0.5 nomirror norotate  offset character 0, 0, 0 autojustify
set ytics autofreq  norangelimit
set y2tics border in scale 1,0.5 nomirror norotate  offset character 0, 0, 0 autojustify
set y2tics autofreq  norangelimit
set ztics border in scale 1,0.5 nomirror norotate  offset character 0, 0, 0 autojustify
set ztics autofreq  norangelimit
set nox2tics
set cbtics border in scale 1,0.5 mirror norotate  offset character 0, 0, 0 autojustify
set cbtics autofreq  norangelimit
set rtics axis in scale 1,0.5 nomirror norotate  offset character 0, 0, 0 autojustify
set rtics autofreq  norangelimit
set title "" 
set title  offset character 0, 0, 0 font "" norotate
set timestamp bottom 
set timestamp "" 
set timestamp  offset character 0, 0, 0 font "" norotate
set rrange [ * : * ] noreverse nowriteback
set trange [ * : * ] noreverse nowriteback
set urange [ * : * ] noreverse nowriteback
set vrange [ * : * ] noreverse nowriteback
set xlabel "z [m]" 
set xlabel  offset character 0, 0, 0 font "" textcolor lt -1 norotate
set x2label "" 
set x2label  offset character 0, 0, 0 font "" textcolor lt -1 norotate
set xrange [ 0 : 2.2 ] noreverse nowriteback
set x2range [ * : * ] noreverse nowriteback
set ylabel "(T-T_{lin})/{/Symbol D}T " 
set ylabel  offset character 0, 0, 0 font "" textcolor lt -1 rotate by -270
set y2label "pressure [bar]" 
set y2label  offset character 0, 0, 0 font "" textcolor lt -1 rotate by -270
set yrange [ * : * ] noreverse nowriteback
set y2range [ * : * ] noreverse nowriteback
set zlabel "" 
set zlabel  offset character 0, 0, 0 font "" textcolor lt -1 norotate
set zrange [ * : * ] noreverse nowriteback
set cblabel "" 
set cblabel  offset character 0, 0, 0 font "" textcolor lt -1 rotate by -270
set cbrange [ * : * ] noreverse nowriteback
set zero 1e-08
set lmargin  -1
set bmargin  -1
set rmargin  -1
set tmargin  -1
set locale "en_US.utf8"
set pm3d explicit at s
set pm3d scansautomatic
set pm3d interpolate 1,1 flush begin noftriangles nohidden3d corners2color mean
set palette positive nops_allcF maxcolors 0 gamma 1.5 color model RGB 
set palette rgbformulae 7, 5, 15
set colorbox default
set colorbox vertical origin screen 0.9, 0.2, 0 size screen 0.05, 0.6, 0 front bdefault
set style boxplot candles range  1.50 outliers pt 7 separation 1 labels auto unsorted
set loadpath "/home/sweiss/.ownCloud/Templates/Gnuplot" 
set fontpath 
set psdir
set fit noerrorvariables
load '~/ownCloud/Templates/Gnuplot/beauty.gnu'
set border 31 front ls 91

f(x)=m*x+n;m=1;n=13;


set datafile separator "|"
set term post eps color enhanced size 4,3
set out 'Fig/6th_Iter.eps'
set multiplot layout 2,1
set yrange [0:0.2/27.]
set y2range [18.98:19.02]
set arrow 1 nohead from 1.1,0 to 1.1,0.2/27.
plot 'temp.csv' u 1:(($2-($1*27./2.2+13))/27.) title 'T' axis x1y1 w l ls 1,'pressure.csv' u 1:(1e-5*$2) title 'p' axis x1y2 w l ls 2

set ylabel "{/Symbol l} [J/K s m]"
set y2label "{/Symbol r} [kg/m^3]"
set yrange [0.0152:0.016]
set y2range [135:180]
set arrow 1 nohead from 1.1,0.0152 to 1.1,0.016
plot 'prop.csv' u 1:2 title '{/Symbol l}' axis x1y1 w l ls 1,'' u 1:3 title '{/Symbol r}' axis x1y2 w l ls 2

unset multiplot

set out

# The thermal expansion coefficient alpha:
set term post eps color enhanced size 4,2
set out 'Fig/alpha.eps'
alpha_m = 0.009189
set yrange [*:*]
set ylabel '({/Symbol a-a}_m)/{/Symbol a}_m)'
unset y2tics
set y2label ""
set arrow 1 nohead from 1.1,-0.3 to 1.1,0.6
plot 'prop.csv' u 1:(($4-alpha_m)/alpha_m) notitle w l ls 2
set out
set term x11
#    EOF
