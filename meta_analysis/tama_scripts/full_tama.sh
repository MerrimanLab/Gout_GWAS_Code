#! /bin/bash

# Split the mantra.dat file into 5,000 SNPs chunk:

split -a 3 -l 25000 mantra.dat
FILE_NUM=$(ls -1 x* | wc -l)
echo $FILE_NUM

# Make a separate directory for each of the smaller files:
for i in $(seq 1 ${FILE_NUM}); do
mkdir mantra$i;
done

# Move each of the files into the corresponding directory:
array1=($(seq 1 ${FILE_NUM}))
array=(x{a..z}{a..z}{a..z})
for ((i=0;i<$FILE_NUM;i++)); do
mv ${array[$i]} mantra${array1[$i]}/mantra.dat;
done

# Copy the genome-wide `dmat.out` file and mantra.in file into all of the directories:
for i in $(seq 1 ${FILE_NUM}); do
cp dmat.out mantra$i/;
cp mantra.in mantra$i/;
done

# Make an `fname.in` file that specifies the names of the output files of the MANTRA program.
# It should be line-separated with the following:
# * Output file name
# * Parameter estimate output file name
# * Seed number

for i in $(seq 1 ${FILE_NUM}); do
echo output$i.out > mantra$i/fname.in;
echo output$i.beta.out >> mantra$i/fname.in;
echo 657089 >> mantra$i/fname.in;
done

# Edit the batch script and run array job

# Combine the outputs when finished:
# cat **/output{1..9258}.out > comb.out
# cat **/output{1..9258}.beta.out > comb.beta.out

