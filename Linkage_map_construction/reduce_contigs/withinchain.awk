#take only chains between same scaffold... and maxD contigs apart
#zcat all.chain.gz|awk -f withinchain.awk >all_within.chain.gz
BEGIN{
	if (maxD == "")
		maxD = 3
}
{
        if ($1=="chain") {
		n1 = substr($3,1,index($3,".c"))
		n2 = substr($8,1,index($8,".c"))
		c1 = substr($3,index($3,".c")+2)+0
		c2 = substr($8,index($8,".c")+2)+0

                if (n1 == n2 && (c2-c1<=maxD && c2-c1>=-maxD))
                        keep = 1
                else
                        keep = 0
        }
        if (keep)
                print
}

