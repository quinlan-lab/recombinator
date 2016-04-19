Demo:

    python recombinator.py --ped K34175.ped \
         --vcf K34175.subset.snps.phased.vcf.gz \
         --family K34175 \
    > K34175.paternal.crossovers.bedgraph

    bgzip K34175.paternal.crossovers.bedgraph

    tabix -p bed K34175.paternal.crossovers.bedgraph.gz