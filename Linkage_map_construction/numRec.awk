#gets index of individual and number of recombinations...
/^[^#]/{
	str = ""
	str2 = ""
	f = 1
	si = 1
	for (i=7; i<=NF; i+=3) {
		str=str substr($i, 1, length($i)/2)
		str2=str2 substr($i, length($i)/2 + 1)
		for (j = 1; j <= length($i)/2; ++j)
			name[si++] = f "," j
		++f
	}
	n = length(str)
	for (i = 1; i <= n; ++i) {
		s = substr(str, i, 1)
		if (prev[i] != "" && s != prev[i]) {
			++r[i]
			pos[i, r[i]] = $1
		}
		prev[i] = s
	}

	for (i = 1; i <= n; ++i) {
		s = substr(str2, i, 1)
		if (prev2[i] != "" && s != prev2[i]) {
			++r2[i]
			pos2[i, r2[i]] = $1
		}
		prev2[i] = s
	}

}
END{
	for (i = 1; i <= n; ++i) {
		printf name[i] "\t" r[i]+0
		for (j = 1; j <= r[i]; ++j)
			printf "\t" pos[i, j]
		printf "\n" name[i] "\t" r2[i]+0 
		for (j = 1; j <= r2[i]; ++j)
			printf "\t" pos2[i, j]
		print ""
	}
}
