#calculate X2 values from a map
(/^[^#]/){
	if ($2 + $3 != prev) {
		printf $1 "\t" ++pos
		for (i = 7;i < NF; i+=3) {
			s1 = substr($i,1,length($i)/2)
			s2 = substr($i,length($i)/2+1)
			printf "\t" X(s1) "\t" X(s2)
		}
		print ""
	}
	prev = $2 + $3
}
function X(s   ,a,b,f,l1,l0)
{
	a = gsub(/0/,"",s)
	b = gsub(/1/,"",s)
	#print a,b
	f=a / (a + b)
	l1 = a * log(f) + b * log(1 - f)
	l0 = (a + b) * log(0.5)
	return 2 * (l1 - l0)
}
