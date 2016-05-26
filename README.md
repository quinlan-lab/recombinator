Demo:

    python recombinator.py --ped K34175.ped \
         --vcf K34175.subset.snps.phased.vcf.gz \
         --families K34175 \
         --parent dad \
    > K34175.paternal.crossovers.bedgraph

    bgzip K34175.paternal.crossovers.bedgraph

    tabix -p bed K34175.paternal.crossovers.bedgraph.gz


For mom and dad:

    python recombinator.py \
        --ped K34175.ped \
        --vcf K34175.subset.snps.phased.vcf.gz \
        --families K34175 \
        --parent mom | head -n 10
    chrom	start	end	same(1)_diff(2)	dad	mom	sib1	sib2
	1	15189	15190	2	G/G	G/A	G/A	G/G
	1	16486	16487	2	T/T	T/C	T/T	T/C
	1	16496	16497	1	A/A	A/G	A/G	A/G
	1	17537	17538	1	C/C	C/A	C/A	C/A
	1	133159	133160	1	G/G	G/A	G/A	G/A
	1	135031	135032	1	G/G	G/A	G/A	G/A
	1	232959	232960	1	C/C	C/A	C/A	C/A
	1	252748	252749	2	C/C	C/A	C/C	C/A
	1	523752	523753	2	C/C	C/T	C/C	C/T

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
