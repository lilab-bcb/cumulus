grep "^>" <(gunzip -c $1) | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt
cat $2 $1 > gentrome.fa.gz