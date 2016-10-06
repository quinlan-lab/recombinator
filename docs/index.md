recombinator
============

*recombinator* calls crossovers, de novo mutations, given phased or unphased VCFs.

It will extract information, make plots, and report info about variants.


```
python -m recombinator \
    recombinator \
    --ped $ped \
    --vcf $vcf \
    --region $chrom \
    --prefix xos
```

This will output a directory structure with informative sites per-family and then an aggregate file
with all families that looks like:
```
chrom	start	end	family_id	parent_id	left-block	left-informative-sites	right-block	right-informative-sites	change
6	13362448	13371471	xx109	GOG012978	6:205799-13362449	5317	6:13371471-42449797	12572	0-1
6	42449796	42454760	xx109	GOG012978	6:13371471-42449797	12572	6:42454760-52769531	3372	1-0
6	52769530	52799232	xx109	GOG012978	6:42454760-52769531	3372	6:52799232-72183699	7443	0-1
6	72183698	72210615	xx109	GOG012978	6:52799232-72183699	7443	6:72210615-111961096	12339	1-0
6	111961095	111971518	xx109	GOG012978	6:72210615-111961096	12339	6:111971518-133578815	5639	0-1
6	133578814	133617587	xx109	GOG012978	6:111971518-133578815	5639	6:133617587-154285086	6058	1-0
6	154285085	154318283	xx109	GOG012978	6:133617587-154285086	6058	6:154318283-170919515	4753	0-1
6	110849230	110909829	xx109	GOG10572	6:7934314-110849231	36477	6:110909829-167867452	18916	1-0
6	167867451	167871017	xx109	GOG10572	6:110909829-167867452	18916	6:167871017-170919470	1267	0-1
6	11221929	11249151	x3750	GOG08928	6:217200-11221930	4411	6:11249151-14621353	1188	1-0
6	14621352	14626728	x3750	GOG08928	6:11249151-14621353	1188	6:14626728-138337903	54524	0-1
6	138337902	138352068	x3750	GOG08928	6:14626728-138337903	54524	6:138352068-170921794	16314	1-0
```

That's for unphased data. For phased data, it will also list the kid in
which the crossover is detected (not possible in unphased).

There will be files:
+ $prefix/$region.crossovers-unfiltered.bed
  - contains the only mimimally filtered crossovers that seem valid for all samples
+ $prefix/$region.crossovers.bed
  - contains a more filtered set (but may exclude small gene conversions) for all samples
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

Cohort Plots and Stats
======================

Once the recombinator script has run. The `cohort-plots` command
can be run as:
```
python -m recombinator cohort-plots --ped $ped --prefix $prefix crossovers.bed
```
and it will output a plot of crossover hot-spots separated by chromosome
and separating maternal from paternal. E.g.

![hotspots](https://cloud.githubusercontent.com/assets/1739/18610633/3b81d9a8-7cde-11e6-813c-9ff3286fce4d.png "hotspots")

Where we see a nice hotspot from fathers at about 65MB on chromosome 13.
The dashed lines indicate a z-score cutoff of >= 2.58, meaning values above (or below for males) that line
are above the 99.5% confidence of the data. We do not assume the data is normally distributed, but this makes
a reasonable cutoff for plotting.

A bed file is sent to stdout with:
+ the number of crossovers at each site
+ the parental sex of 'male', 'female', or 'both' where 'both is the sum of male and female crossovers.
+ the zscore of the number of crossovers.
It is expected that users will filter this to their own z-score cutoff and sex as needed.


It will also create an aggregate plot with the count of maternal and paternal crossovers:
![sex](https://cloud.githubusercontent.com/assets/1739/18605165/cb82df12-7c47-11e6-80da-0985482de14c.png "sex")

Enrichment
==========

After finding per-sample recombination sites, we will want to know if certain events are enriched in those crossovers. For example, we may expect that de novos are more likely to occur in a crossover
from the same sample.
We can evaluate this by comparing the observed number of overlaps between crossovers and de novos
in the same sample Vs the expected, which we derive by shuffling the sample-ids of the crossovers.
By shuffling the sample ids, we avoid bias due to genome composition at crossovers that would
come into play if we shuffled the spatial location of the events.

This can be done like this:
```
python -m recombinator enrichment $crossovers.all-samples.bed denovos.all-samples.bed --simulations 1000
```
To compare the sample-overlap in 1000 simulations to the single observed event.
This assumes that `crossovers.all-samples.bed` and `denovos.all-samples.bed` have a `sample_id`
column that indicates the sample in which the event occurred.

An example output looks like:

![dnxo](https://cloud.githubusercontent.com/assets/1739/18727794/a4b475b2-8007-11e6-9b77-85a4918359d1.png "DN::XO")

Where we see that, the observed de-novo::xo overlaps shown as the green line (57) is higher than
nearly all of the simulated values in blue.

While this example is for de novos, it can be used to compare crossovers to any data containg regions labelled by sample.

Phased vs Unphased
==================

In our experience, it's not simple to phase large cohorts with available methods. Beagle is quadratic in the
number of samples and shapeit does not give good results even with strict quality-filtering on input variants.

To get around this, the pipeline we recommend is to run `recombinator` on the unphased VCF. Then, run
`phased-from-unphased.sh` which will gather any pairs of variants between which there is a state-change.
It will combine those crossover-bounding variants with the (~2.5 million) omni-2.5 sites, output a with
only those sites, phase it, and run recombinator on the phased.

