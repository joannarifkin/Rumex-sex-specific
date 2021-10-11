#scale score so that somewhat shorter alignments can be used...
BEGIN{
	if (maxScale == "")
		maxScale=2
	if (maxLength == "")
		maxLength = 1000
}
{
        if ($1=="chain") {
		#scale score so that somewhat shorter alignments can be used...
		if ($7-$6 < maxLength) {
			f = maxLength / ($7-$6)
			if (f > maxScale)
				f = maxScale
			$2 = int(f * $2)
		}
        }
        print
}
