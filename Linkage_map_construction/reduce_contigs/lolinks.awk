#awk -f lolinks.awk contigs2.lo >contigs2.paf
#contigs must be in the same order as in scaffolds...
BEGIN{
	alt["+"]="-"
	alt["-"]="+"
}
{
	for (j=1;j<=3;++j)
		data[NR,j]=$j
}
END{
	for (i = 1; i < NR; ++i)
		for (j = 1; j <= 4; ++j)
			if (i + j <= NR) {
				printLink(data[i,1], data[i,2], data[i,3], data[i+j,1], data[i+j,2], data[i+j,3],5-j)
				printLink(data[i,1], data[i,2], alt[data[i,3]], data[i+j,1], data[i+j,2], data[i+j,3],4-j)
				printLink(data[i,1], data[i,2], data[i,3], data[i+j,1], data[i+j,2], alt[data[i+j,3]],4-j)
				printLink(data[i,1], data[i,2], alt[data[i,3]], data[i+j,1], data[i+j,2], alt[data[i+j,3]],3-j)
			}
}
#print a link, c=contig, l=length, o=orientation, w=weight (integer)
function printLink(c1, l1, o1, c2, l2, o2, w     ,i) {
	if (substr(c1,1,index(c1,".c")) == substr(c2,1,index(c2, ".c")))
	for (i = 1; i <= w; ++i) {
        	if (o1 == "-")
			print "agp" ++n "\t200\t0\t100\t-\t" c1 "\t" l1 "\t0\t100\t60"
		else
			print "agp" ++n "\t200\t0\t100\t+\t" c1 "\t" l1 "\t" (l1-100) "\t" l1 "\t60"
		if (o2 == "-")
			print "agp" n "\t200\t100\t200\t-\t" c2 "\t" l2 "\t" (l2-100) "\t" l2 "\t60"
		else
			print "agp" n "\t200\t100\t200\t+\t" c2 "\t" l2 "\t0\t100\t60"
	}
}
