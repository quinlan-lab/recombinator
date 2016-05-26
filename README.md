For mom and dad:

```
python recombinator.py \
	--ped K34175.ped \
	--vcf K34175.subset.snps.phased.vcf.gz \
	--families K34175 \
	| head -n 10
chrom	start	end	parent	family_id	same(1)_diff(2)	dad	mom	sib1	sib2	global_call_rate	global_depth_1_10_50_90
1	14932	14933	dad	K34175	1	G/A	G/G	G/G	G/G	1.00	41|43|51|61
1	14975	14976	dad	K34175	1	G/A	G/G	G/G	G/G	1.00	41|42|44|58
1	15189	15190	mom	K34175	2	G/G	G/A	G/G	G/A	1.00	32|37|52|68
1	15446	15447	dad	K34175	1	A/G	A/A	A/A	A/A	1.00	43|44|50|82
1	16256	16257	dad	K34175	1	G/C	G/G	G/C	G/C	1.00	34|34|36|52
1	16486	16487	mom	K34175	2	T/T	T/C	T/C	T/T	1.00	61|64|72|74
1	16496	16497	mom	K34175	1	A/A	A/G	A/G	A/G	1.00	58|59|66|77
1	17537	17538	mom	K34175	1	C/C	C/A	C/A	C/A	1.00	35|36|39|42
1	19917	19918	dad	K34175	2	G/C	G/G	G/G	G/C	1.00	36|38|45|49
```


Mark crossovers:

    python remove_noise.py --bg K34175.paternal.crossovers.bedgraph \
    | awk '{if ($4 =="1") {print $1"\t"$2"\t"$3"\t.\t.\t+"} else {print $1"\t"$2"\t"$3"\t.\t.\t-"}}' \
    | bedtools merge -s -i - \
    | bedtools groupby -g 1,4 -c 2,3 -o min,max \
    | awk '{if ($2=="+"){print $1"\t"$3"\t"$4"\t1"} else {print $1"\t"$3"\t"$4"\t2"}}' \
    | python report_crossovers.py --bg /dev/stdin
    1	8199262	8213868
    1	64090376	64096094
    1	165191701	165193662
