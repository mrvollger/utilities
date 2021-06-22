#!/bin/bash



DIR="canu.assembly"
DIR="WH.assembly"
DIR=${1:-asm}
>&2 echo "Looking for canu dirs called: $DIR"



find . \( -type d -and \( \
	-path \*/$DIR/unitigging \
	-or -path \*/$DIR/canu-logs \
	-or -path \*/$DIR/unitigging.html.files \
	-or -path \*/$DIR/trimming.html.files \
	-or -path \*/$DIR/trimming \
	-or -path \*/$DIR/correction.html.files \
	-or -path \*/$DIR/correction \
	-or -path \*/$DIR/canu-scripts \
   	\) \) -prune -print #exec rm -rf "{}" \;


# xargs -n 200 -P 24 rm -r 

