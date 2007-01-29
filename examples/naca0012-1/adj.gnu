       reset
       set term png
       set out 'adj.png'

       set multiplot

       set size 0.5,0.5
       set origin 0.0,0.5
       set xran[-0.05:1.05]
       set nokey
#      set bmargin 1
       set title "Pressure/Friction coeffient" ,-1
#      set y2tics
       set ylabel '-Cp'
#      set y2label 'Cf'
       plot 'WALL.DAT' w lp lw 2 pt 6

       reset
       set size 0.5,0.5
       set origin 0.5,0.5
       set auto
       set nokey
       set logscale y
       set title "Residue" ,-1
       plot 'ADJ.RES' u 1:2 w l
       set nologscale y

       reset
       set size 0.5,0.5
       set origin 0.07,0.0
       set xran[-0.5:1.5]
       set yran[-1.0:1.0]
       set size ratio -1
       set noxtics
       set noytics
       set nokey
       set title "Mach number" ,-1
#      set bmargin 1
       plot 'ADJ.D' w l,'BD.DAT' w l lt 1 lw 2

       reset
       set size 0.5,0.5
       set origin 0.57,0.0
       set xran[-0.5:1.5]
       set yran[-1.0:1.0]
       set size ratio -1
       set noxtics
       set noytics
       set nokey
#      set bmargin 1
       set title "Pressure" ,-1
       plot 'ADJ.U' w l,'BD.DAT' w l lt 1 lw 2

       unset multiplot
#      set term x11
