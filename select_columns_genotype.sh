#!/bin/sh

# I thank Dr. Gerardo Alves Fernandes Júnior for this challenge.

# Mandatory files (examples)
# data.txt
# 12000 0 2 1 0 5 0
# 12001 0 2 1 0 5 0
# 12002 0 2 1 0 5 0

# List of columns to be selected
# list.txt
# 2
# 3
# 5

# Text file
# ID.txt
# ID

# Rotate commands in order:

# informs you of the column total of the data
awk '{print NF}' data.txt | uniq

# creates the sequence of columns (column vector)
# seq 1 (total_colunas - 1) > seq
seq 1 6 > seq

# transpose the sequence (line vector
awk '{for (f=1;f<=NF;f++) col[f] = col[f]":"$f} END {for (f=1;f<=NF;f++) print col[f]}' seq | tr ':' ' ' >> temp1

# joins the ID file with the transposed sequence file
paste -d '' ID temp1 >> temp2

# adds the header to the data file
printf '0r temp2\nwq\n' | ed -s data.txt

# selects the columns from an external list
cols=($(sed '1!d;s/ /\n/g' data.txt | grep -nf list.txt | sed 's/:.*$//'))
cut -d ' ' -f 1$(printf ",%s" "${cols[@]}") data.txt >> temp.txt
cat temp.txt 
sed 1d temp.txt > output.txt
rm temp.txt temp1 temp2 seq.txt

# add header (if needed)
sed -e '1i\HEADER_HERE' data > data.txt

# Source: https://stackoverflow.com/questions/11098189/extract-columns-from-a-file-based-on-header-selected-from-another-file
# Source: https://unix.stackexchange.com/questions/25592/creating-a-sequence-of-numbers-one-per-line-in-a-file
# Source: https://stackoverflow.com/questions/1729824/an-efficient-way-to-transpose-a-file-in-bash

