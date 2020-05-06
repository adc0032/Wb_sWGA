#printf 'chr\n%.0s' {1..1166} >chr
#seq 1 1166 > num
#paste chr num
fasta1=$1
while read -r line
do
 set -- $line
 sed -i s/$1/$2/g $fasta1
done < scafmap.txt