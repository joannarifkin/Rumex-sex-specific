#awk -vs=column1,column2,... -f zero.awk post.txt
#zeros genotype likelihoods (posteriors)
BEGIN{
	FS="\t"
	OFS="\t"
	n = split(s, sp, ",")
}
(NR<=7){print}
(NR>7){
	for (i = 1; i <= n; ++i)
		$(sp[i]) = "1 1 1 1 1 1 1 1 1 1"
	print
}
