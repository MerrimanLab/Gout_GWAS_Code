#! /bin/bash

# Just a script to rename X into 23:

sed 's/^X/23/g' $1 > ${1%.dat}.txt

