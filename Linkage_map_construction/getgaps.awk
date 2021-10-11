#awk -f getgaps.awk ref.fasta
#finds gaps consisting of "N" or "n"
{
        if (/>/) {
                if (name != $1) {
                        findNs(s)
                        s = ""
                }
                name = $1
        } else {
                s = s toupper($0)
        }
}
END{
        findNs(s)
}
function findNs(s   ,n,i,a,st,en,l) {
	n=split(s,a,"N")
	st = 0
	en = 0
	for (i=1;i<=n;++i) {
		l=length(a[i])
		if (l > 0) {
			if (st > 0)
				print substr(name, 2) "\t" st "\t" en
			st = en + l + 1
			en = st
		}
		else
			++en
	}
}


#substr is slow
#function findNs(s   ,n,i,pos){
#	n = length(s)
#        i = 0
#        while (1) {
#		match(s,/[Nn]+/)
#		if (RSTART > 0) {
#			print substr(name, 2) "\t" (i + RSTART)  "\t" (i + RSTART + RLENGTH - 1)
#			s = substr(s, RSTART + RLENGTH)
#			i += RSTART + RLENGTH - 1
#		} else
#			break
#        }
#}

