#calculate X2 values from a map
(/^[^#]/){
	if ($2 + $3 != prev) {
		printf $1 "\t" ++pos
		for (i = 7;i < NF; i+=3) {
			s1 = substr($i,1,length($i)/2)
			s2 = substr($i,length($i)/2+1)
			printf "\t" X3(s1,s2)
		}
		print ""
	}
	prev = $2 + $3
}
function mylog(a)
{
	if (a == 0)
		return 0
	else
		return log(a)
}
function X3(s1, s2  ,a1,a2,n,i,c,f0,f1,f2,l1,l0)
{
	split(s1,a1,"")
	n=split(s2,a2,"")
	for (i=1;i<=n;++i)
		++c[a1[i]+a2[i]]

#	print c[0],c[1],c[2]

	f0=c[0] / n
	f1=c[1] / n
	f2=c[2] / n

	l1 = c[0] * mylog(f0) + c[1] * mylog(f1) + c[2] * mylog(f2)
	l0 = (c[0] + c[2]) * mylog(0.25) + c[1] * mylog(0.5)
	return 2 * (l1 - l0)
}
