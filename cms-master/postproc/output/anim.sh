#!/bin/bash
#
name='2018_3'
#
anim_dir='figs_png_'${name}
echo ${anim_dir}
cp -r figs_png ${anim_dir}
rm figs_png/*.png
#
echo 'copy fin.'
#
#convert -layers optimize -delay 10 figs_png_${name}/*.png anim/anim_${name}.gif
convert -layers optimize -delay 10 figs_png_${name}/*.png anim/anim_${name}.gif
