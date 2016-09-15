recombinator
============

*recombinator* calls crossovers given a phased or unphased VCF with quartets.

when run like:

```
python recombinator.py \
    --ped $ped \
    --vcf $vcf \
    --region $chrom \
    --prefix xos
```

This will output a directory structure with crossovers per-family and then an aggregate file
with all families that looks like:
```
chrom	start	end	family_id	parent_id	informative-sites	informative-sites-r	change
6	5853504	5853541	fam1	ss443mom	385	3	0-1
6	5853571	5863131	fam1	ss443mom	3	86	1-0
6	6176624	6177027	fam1	ss443mom	86	3703	0-1
6	16398958	16402108	fam1	ss443mom	3703	1976	1-0
6	21955963	21959436	fam1	ss443mom	1976	2620	0-1
6	27163297	27163317	fam1	ss443mom	2620	9	1-0
6	27163370	27166175	fam1	ss443mom	9	1778	0-1
6	10120946	10122850	fam1	ss443dad	1973	5720	0-1
6	24851786	24852750	fam1	ss443dad	5720	3328	1-0
```

That's for unphased data. For phased data, it will also show the kid in
which the crossover is detected (not possible in unphased).

There will be files:
+ $prefix/$region.crossovers-unfiltered.bed
  - contains the only mimimally filtered crossovers that seem valid for all samples
+ $prefix/$region.crossovers.bed
  - contains a more filtered set (but may exclude small gene conversions) for all samples
+ $prefix/$region/ directory also contains the above to files per-family
+ $prefix/$region/$region.$family.png shows a plot of the crossovers like the example below.
+ $prefix/$region/$region.$family.$parent.bed.gz contains all of the informative sites
  in that parent. Any time the 'same' column changes, there is a putative crossover.

Example Plot
------------

First we show where there is a true crossover. The black points indicate informative variants.
The values of the black points indicate the crossover state where
the value is 1 if the alternate allele is at the same index (1st or 2nd) in both the parent
and the kid. The value is 0 if they are different. We can see that the genotype calls are extremely
good here, in many cases, we'd expect to see more single points outside of the blue blocks which
indicate the inferred state for a region. The red line indicates the location of a crossover.

![xo](https://cloud.githubusercontent.com/assets/1739/18555974/1697ac52-7b27-11e6-8faf-9659b2fd9c15.png "Clean Crossover")

Below, we show a case where there is no crossover for the entire chromosome. We do see a few noisy
points that are likely due to genotyping error.

![noxo](https://cloud.githubusercontent.com/assets/1739/18555978/19f7dd54-7b27-11e6-8da8-eb749b8093ce.png "No Crossovers")

