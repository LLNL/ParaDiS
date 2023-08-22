#!/bin/bash

mins=15

if [ $1 ]; then
   mins=$1
fi

cmd="mxterm -g 300x80+60+10  1 1 $mins -q pbatch -A wbronze"

echo $cmd; $cmd
