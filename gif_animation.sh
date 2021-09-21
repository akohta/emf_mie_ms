#!/bin/bash

# settings
Delay=10
Loop=0

if [ $# -ne 1 ]; then
 echo "Usage : $0 image_directory"
 exit
fi

if [ ! -d $1 ]; then
 echo "$1 does not exist or is not a directory."
 exit
fi

convert -delay $Delay -loop $Loop $1/xy_Ex_*.png xy_Ex.gif &
convert -delay $Delay -loop $Loop $1/xy_Ey_*.png xy_Ey.gif &
convert -delay $Delay -loop $Loop $1/xy_Ez_*.png xy_Ez.gif &

convert -delay $Delay -loop $Loop $1/yz_Ex_*.png yz_Ex.gif &
convert -delay $Delay -loop $Loop $1/yz_Ey_*.png yz_Ey.gif &
convert -delay $Delay -loop $Loop $1/yz_Ez_*.png yz_Ez.gif &

convert -delay $Delay -loop $Loop $1/xz_Ex_*.png xz_Ex.gif &
convert -delay $Delay -loop $Loop $1/xz_Ey_*.png xz_Ey.gif &
convert -delay $Delay -loop $Loop $1/xz_Ez_*.png xz_Ez.gif &

wait

convert -delay $Delay -loop $Loop $1/xy_Hx_*.png xy_Hx.gif &
convert -delay $Delay -loop $Loop $1/xy_Hy_*.png xy_Hy.gif &
convert -delay $Delay -loop $Loop $1/xy_Hz_*.png xy_Hz.gif &

convert -delay $Delay -loop $Loop $1/yz_Hx_*.png yz_Hx.gif &
convert -delay $Delay -loop $Loop $1/yz_Hy_*.png yz_Hy.gif &
convert -delay $Delay -loop $Loop $1/yz_Hz_*.png yz_Hz.gif &

convert -delay $Delay -loop $Loop $1/xz_Hx_*.png xz_Hx.gif &
convert -delay $Delay -loop $Loop $1/xz_Hy_*.png xz_Hy.gif &
convert -delay $Delay -loop $Loop $1/xz_Hz_*.png xz_Hz.gif &

wait
