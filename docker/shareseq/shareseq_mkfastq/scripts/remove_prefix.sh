#!/bin/bash

# if auto, compare file sizes and choose the one with larger total size
remove_prefix () {
for fastq in $@
do
    fastqnew=`cut -c 6-11 --complement <<< $fastq`
    mv $fastq $fastqnew
done
}

readarray fwd_arr < <(find $1 -maxdepth 1 -type f -name '__fwd_*.fastq.gz' -print0)
readarray rvs_arr < <(find $1 -maxdepth 1 -type f -name '__rvs_*.fastq.gz' -print0)

if [ ${#fwd_arr[@]} -gt 0  ] && [ ${#rvs_arr[@]} -eq 0 ]
then
for fastq in ${fwd_arr[@]}
do
    fastqnew=`cut -c 6-11 --complement <<< $fastq`
    mv $fastq $fastqnew
done
elif [ ${#fwd_arr[@]} -eq 0  ] && [ ${#rvs_arr[@]} -gt 0 ]
then
for fastq in ${rvs_arr[@]}
do
    fastqnew=`cut -c 6-11 --complement <<< $fastq`
    mv $fastq $fastqnew
done
elif [ ${#fwd_arr[@]} -gt 0 ] && [ ${#rvs_arr[@]} -gt 0  ]
then
fwd_size=`du -c ${fwd_arr[@]} | tail -n 1 | cut -f 1`
rvs_size=`du -c ${rvs_arr[@]} | tail -n 1 | cut -f 1`
if [ $fwd_size -gt $rvs_size ]
then
    remove_prefix ${fwd_arr[@]}
    rm -f ${rvs_arr[@]}
else
    remove_prefix ${rvs_arr[@]}
    rm -f ${fwd_arr[@]}
fi
fi