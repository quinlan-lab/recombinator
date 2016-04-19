Demo:

    python recombinator.py --ped K34175.ped \
         --vcf K34175.subset.snps.phased.vcf.gz \
         --family K34175 \
         --parent dad \
    > K34175.paternal.crossovers.bedgraph

    bgzip K34175.paternal.crossovers.bedgraph

    tabix -p bed K34175.paternal.crossovers.bedgraph.gz


For mom:

    python recombinator.py \
        --ped K34175.ped \
        --vcf K34175.subset.snps.phased.vcf.gz \
        --family K34175 \
        --parent mom | head -n 10
    chrom	start	end	same(1)_diff(2)	dad_102622	mom_102621	template_102620	sib_102623
	1	15189	15190	2	G/G	G/A	G/A	G/G
	1	16486	16487	2	T/T	T/C	T/T	T/C
	1	16496	16497	1	A/A	A/G	A/G	A/G
	1	17537	17538	1	C/C	C/A	C/A	C/A
	1	133159	133160	1	G/G	G/A	G/A	G/A
	1	135031	135032	1	G/G	G/A	G/A	G/A
	1	232959	232960	1	C/C	C/A	C/A	C/A
	1	252748	252749	2	C/C	C/A	C/C	C/A
	1	523752	523753	2	C/C	C/T	C/C	C/T

For dad:

    python recombinator.py \
        --ped K34175.ped \
        --vcf K34175.subset.snps.phased.vcf.gz \
        --family K34175 \
        --parent dad | head -n 10
    chrom	start	end	same(1)_diff(2)	dad_102622	mom_102621	template_102620	sib_102623
	1	14932	14933	1	G/A	G/G	G/G	G/G
	1	15446	15447	1	A/G	A/A	A/A	A/A
	1	16256	16257	1	G/C	G/G	G/C	G/C
	1	19917	19918	2	G/C	G/G	G/C	G/G
	1	19918	19919	2	G/A	G/G	G/A	G/G
	1	19941	19942	1	G/C	G/G	G/C	G/C
	1	20155	20156	1	C/T	C/C	C/T	C/T
	1	20157	20158	1	A/C	A/A	A/A	A/A
	1	20165	20166	1	A/G	A/A	A/A	A/A