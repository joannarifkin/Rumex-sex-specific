#Lep-MAP3
#remove some individuals...
#cut -f 7,45,55,58,87,110,113,119,125,151,167,176,184 --complement
#...
#final maps map1_for_la.txt, map2_for_la.txt and map3_for_la.txt
#informativemask = 1, 2 and 3
#final input data all_post.call_f.gz

#Lep-Anchor
#split into contigs
awk -f getgaps.awk hastate_28Sep2018_nbKXS.fasta >hastate_28Sep2018_nbKXS.gaps
awk '(NR==FNR){l[$1]=$2}(NR!=FNR){print $1 "\t" g[$1]+1 "\t" $2-1;g[$1]=$3}END{for (c in l) print c "\t" (g[c]+1) "\t" l[c]}' *.length *.gaps|awk '{++c[$1];print substr($1,1,index($1,";")-1)".c"c[$1]"\t1\t"$3-$2+1"\t"c[$1]"\tW\t"$0"\t+"}' >contigs.agp
awk -f la/makefasta.awk hastate_28Sep2018_nbKXS.fasta contigs.agp >contigs.fasta
#run Haplomerger + winmasker, create softmasked genome contigs_sm.fa.gz, chain all.chain.gz

#find and remove full haplotypes
zcat all.chain.gz|awk -f findFullHaplotypes.awk >fullHaplotypes50.txt
grep -v -F -w -f <(cut -f 2 fullHaplotypes50.txt) contigs.agp|sort -V >contigs_nh.agp
awk -f liftover.awk contigs.agp snps.txt >snps_c.txt

#reduce number of contigs
#see reduce_contigs

#make final contig fasta
awk -f la/makefasta.awk <(zcat contigs_sm.fa.gz) reduce_contigs/contigs2.agp|gzip >contigs2.fa.gz
#run Haplomerger on this, chain all2.chain.gz

#paf in in the scaffold coordinates, liftover of the paf to contigs2

#run Lep-Anchor normally
#see la/

#physical map

#snps_c.liftover2 contains the final contig names (after reduce contigs)

awk -f la/liftover.awk <(cat la/chr*.agp) snps_c.liftover2|sort -V >snps_c.liftover2.lg

for i in {1..5}
do
#echo "zcat all_post.call_f.gz|cut -f 7,45,55,58,87,110,113,119,125,151,167,176,184 --complement|java -cp lm3/ OrderMarkers2 data=- map=map_js.txt chromosome=$i evaluateOrder=<(grep LG$i snps_c.liftover2.lg|cut -f 3) improveOrder=0 phasingIterations=3 hyperPhaser=1 scale=0.05 2 0.5 >phys/order$i.txt 2>phys/order$i.err"
echo "zcat all_post.call_f.gz|awk -vs=151 -f zero.awk|java -cp lm3/ OrderMarkers2 data=- map=map_js.txt chromosome=$i evaluateOrder=<(grep LG$i snps_c.liftover2.lg|cut -f 3) improveOrder=0 phasingIterations=3 hyperPhaser=1 scale=0.05 2 0.5 useMorgan=1 >phys/order$i.txt 2>phys/order$i.err"
done >porder_commands.txt

source porder_commands.txt

for i in {1..5}
do
awk '(NR==FNR){m[$3]=$1"\t"$2}(NR!=FNR && /^[^#]/){print m[$1]"\t"$2"\t"$3}' snps_c.liftover2.lg phys/order$i.txt|cut -f 1-4 >phys/order$i.mapped
done

#remove some markers from map ends causing long gaps
