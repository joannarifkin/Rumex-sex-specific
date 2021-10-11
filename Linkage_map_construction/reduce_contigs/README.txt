#reduce number of short contigs...

zcat ../all.chain.gz|awk -f scalechain.awk|awk -vmaxD=4 -f withinchain.awk |gzip >within4.chain.gz


#take only contigs of 500bp or longer
rm test*.order
awk 'BEGIN{m=1}($3>=500){if ($6!=p) {if (n>=1000) {++m;n=0}; fn="test" m ".order"}; print $1"\t"$2"\t"$3 >fn; p=$6;++n}' ../contigs_nh.agp

rm test.eval*
for i in test*.order
do
        awk '{print $0"\t?\t1"}' $i >test.bed
        awk '{print $1"\t" int($3/2)"\t1\t1"}' $i >test.map
        zcat within4.chain.gz|java -cp ../bin/ PlaceAndOrientContigs map=test.map bed=test.bed evaluateAnchoring=$i chain=- improveAnchoring=1 numRuns=1 >>test.eval 2>>test.eval.err
done

awk '($3>=500)' ../contigs_nh.agp|cut -f 1-3 >full.order
awk '{print $0"\t?\t1"}' full.order >test.bed
awk '{print $1"\t" int($3/2)"\t1\t1"}' full.order >test.map
zcat within4.chain.gz|java -cp ../bin/ PlaceAndOrientContigs map=test.map bed=test.bed evaluateAnchoring=test.eval chain=- 2>full.err >full.eval

#calculate new contigs
awk -f makeagp2.awk full.eval|awk -vOFS="\t" '{if (!($1 in name)) name[$1]=$6; $1="P"name[$1];print}'|sort -V >contigs2.agp

#and their most likely orientation
awk 'BEGIN{map[0]="-";map[1]="+"}{l[$1]=$3;p[$1,$9]+=($3-$2+1)}END{for (i in l) print i"\t"l[i]"\t"map[p[i,"+"]>=p[i,"-"]]}' contigs2.agp |sort -V >contigs2.lo

#liftover to new contigs
cat ../fullHaplotypes50.txt >liftover.haplotypes
awk '($3<500){print 1"\t"$1"\t"$2"\t"$3}' ../contigs_nh.agp >>liftover.haplotypes
awk '(NR==FNR){l[$1]=$3}(NR!=FNR && /^[^#]/){if ($2>1) print "1\t" $1"\t1\t"($2-1); if ($3<l[$1]) print "1\t" $1"\t"($3+1)"\t"l[$1]}' ../contigs_nh.agp full.eval >>liftover.haplotypes

zcat ../all.chain.gz|java -cp ../bin/ LiftoverHaplotypes haplotypes=liftover.haplotypes map=<(awk '{print $0"\t"NR-1}' ../snps_c.txt) chain=- >../snps_c.liftover

awk -f ../liftover.awk contigs2.agp ../snps_c.liftover >snps_c.liftover2

awk -f lolinks.awk contigs2.lo >contigs2.paf &

