# Code for subsampling of a fastq file (either no. reads or percentage of the total)

# This code takes a fastq file and outputs a subset containing k reads (k must be modified)
cat single.fastq | awk '{ printf("%s",$0); n++; if(n%4==0) {
printf("\n");} else { printf("\t");} }' |
awk -v k=10000 'BEGIN{srand(systime() + PROCINFO["pid"]);}{s=x++<k?x1:int(rand()*x);if(s<k)R[s]=$0}END{for(i in R)print R[i]}' |
awk -F"\t" '{print $1"\n"$2"\n"$3"\n"$4 > "single_sub.fastq"}'

# This code produces a percentage of the reads/fastq
cat single.fastq | awk '{ printf("%s",$0); n++; if(n%4==0) { printf("\n");} else { printf("\t");} }' |
awk 'BEGIN { srand(systime() + PROCINFO["pid"]); } NR % 4 == 0 { count++; } END { target = count * (percentage / 100); } { s = ++x <= target ? x : int(rand() * x); if (s <= target) { R[s] = $0; } } END { for (i in R) print R[i]; }' percentage="$your_percentage" |
awk -F"\t" '{print $1"\n"$2"\n"$3"\n"$4 > "single_sub.fastq"}

# N.B. I haven't actually tested this code yet - I will shortly
