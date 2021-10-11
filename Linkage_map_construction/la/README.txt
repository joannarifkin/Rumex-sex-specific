cat ../contigs2.paf ../DTG-DNA-243_all_raw.lo.paf >all.paf
./lepanchor_wrapper2.sh -n 5 -m ../map3_for_la.txt -m ../map2_for_la.txt -m ../map1_for_la.txt -T 6 -t 5 -p all.paf -c ../all2.chain.gz -f ../contigs2.fa.gz
