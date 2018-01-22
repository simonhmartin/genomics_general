#### This is a collection of scripts for a range of genomic data processing and analysis.
#### Below are notes about some useful tools.

___

## Parsing VCF files

Most of my scripts use a processed `.vcf` format that I call `.geno`. This looks something like this:

```
#CHOM      POS      ind1      ind2      ind3
scaffold1  1        A/A       A/G       G|A
scaffold1  1        N/N       T/T       T|C
```

Missing data is denoted as `N`, and phased and unphased genotypes are shown conventionally with `|` and `/`.

The script `parseVCF.py` in the `VCF_processing` directory, will convert vcf to this format. It has various options for filtering based on read depth, genotype quality or any other flag in the `FORMAT` column of the vcf.

#### Example command:

```bash
python parseVCF.py -i input.vcf.gz --skipIndels --minQual 30 --gtf flag=DP min=5 | gzip > output.geno.gz
```
___

## Diversity and divergence analyses in sliding windows

The script `popgenWindows.py` computes some standard population genomic statistics in sliding windows:  *pi*,  *F*<sub>ST</sub> and *D<sub>XY</sub>*. It requires the script `genomics.py` to be present in the same directory, or in your Python path.

#### Example command
```bash
python popgenWindows.py -w 50000 -m 5000 -g input.geno.gz -o output.csv.gz -f phased -T 5 -p popA A1,A2,A3,A4 -p popB B1,B2,B3,B4,B6,B6 -p popC C1,C2
```
`python popgenWindows.py -h` Will print a full list of command arguments.

#### Notes

* Input is a `.geno` file as shown above. This can be gzipped (`.geno.gz`).
Output is a `.csv`. If you add `.gz` it will be gzipped.

* Genotypes can be incoded in various ways, as indicated by the `-f` flag.

`-f phased` means that genotypes include phase information (but they don't necessarily need to actually be phased):
```
scaffold1  1        A/A       G/G       G|A`
``` 

`-f pairs` is the same, but without the phase indicator:
```
scaffold1  1        AA       GG       AG
```

`-f haplo` means that genotypes are haploid bases:
```
scaffold1  1        A       G       G
```

`-f diplo` means that genotypes are single bases denoting the diploid genotype, using ambiguity codes for heterozygotes:
```
scaffold1  1        A       G       R
```

* The window can be defined based on genomic coordinates (`--windType coord`) or the number of sites (`--windType sites`). Windows will not cross contig/scaffold boundaries.

* The most common source of errors here involve the `-m` (`--minSites`) flag. If you are useing coordinate windows and have any sites with missing data, then `-m` must be set to a value smaller than the window size. If you have reduced representation data such as RADseq, you will need a much lower `-m` value (more like 1% of the window size or even less).

* If some samples are haploid and others are diploid, you can use one of the diploid formats, but indicate that certain samples are haploid by listing them after the `--haploid` flag. The script will force them to have haploid genotyps, and any apparently heterozygous genotype will be converted to `N`.

* The script can run on multiple cores (`-T` flag). Try different numbers, as using too many can slow the script down (due to the difficulty in sorting the outputs coming from the different cores).


___
## Compute ABBA-BABA statistics in sliding windows

The script `ABBABABAwindows.py` performs analyses described in [Martin et al. 2015, MBE](http://mbe.oxfordjournals.org/content/32/1/244.abstract?sid=a3d00925-b3fe-4214-b142-256739082832), compurting the *D* statistic and *f* estimators in windows across the genome. Like the script above, it requires `genomics.py`.

#### Example command

```bash
python ABBABABAwindows.py -g /zoo/disk1/shm45/vcf/set62/set62.chr21.DP5GQ30.AN100MAC1.diplo.gz -f diplo -o output.csv -w 100000 -m 100 -s 100000 -p P1 A1,A2,A3,A4 -p P2 B1,B2,B3,B4 -p P3  C1,C2,C3,C4 -p O D1,D2,D3,D4 -T 10 --minData 0.5
```
`python ABBABABAwindows.py -h` Will print a full list of command arguments.

#### Notes
* This script shares several command arguments with `popgenWindows.py`. And input formats are the same. Please see the notes for that script above.

* Four populations, with the names `P1`, `P2`, `P3` and `O` are requied.

####Output


|  Column Header  | Description |
|: --- :| --- | 
| `scaffold` | The scaffold the window is on (all windows are on a single scaffold) |
| `start` | window start position (inclusive) |
| `end` | window end position (NOTE, this can exceed the length of the scaffold) |
| `sites` | Number of genotypes sites in the input file in each window |
| `sitesUsed` | number of sites used to compute statistics (biallelic SNPs) |
| `ABBA` | Pseudo count of ABBA sites (including polymorphic sites) (See [Martin et al. 2015](https://doi.org/10.1093/molbev/msu269) Equation 2) |
| `BABA` | Pseudo count of BABA sites (including polymorphic sites) (See [Martin et al. 2015](https://doi.org/10.1093/molbev/msu269) Equation 3) |
| `D` | D statistic (see [Durand et al. 2011](https://doi.org/10.1093/molbev/msr048) Equation 2) |
| `fd` | fd admixture estimation (See [Martin et al. 2015](https://doi.org/10.1093/molbev/msu269) Equation 6) |
| `fdM` | Malinsky's modified fd to accomodate admixture between either P1 and P3 or P2 and P3 (See [Malinsky et al. 2015](https://doi.org/10.1126/science.aac9927) Supplementart Material Page 8) |

___
## Make trees for sliding windows

Two scripts in the `phylo/` directory will make trees in sliding windows: `phymlWindows.py` and `raxmlWindows.py`. As the names suggest they use [Phyml](http://www.atgc-montpellier.fr/phyml/) and [RAxML](http://sco.h-its.org/exelixis/web/software/raxml/), respectively.

#### Example command
```bash
python phyml_sliding_windows.py -T 10 -g input.phased.geno.gz --prefix output.phyml_bionj.w50 -w 50 --windType sites --model GTR 
```
`python phymlWindows.py -h` Will print a full list of command arguments.

#### Notes
* You need to have Phyml (or RAxML) installed on your machine. You can direct the script to the location of the executable. I recommend using an unthreaded version, since each window tree will run very quickly.

* The window can be defined based on genomic coordinates (`--windType coord`) or the number of sites (`--windType sites`). Windows will not cross contig/scaffold boundaries.

* Genotypes need to be in either the `phased`, `haplo` or `diplo` formats shown above (`diplo` format is not recommended, as heterozygous genotypes will be treated as single genotypes with ambiguous bases, which are ignored by Phyml.

* For the raxml script, you could also use `diplo` format, although I'm not sure whether the ambiguity codes will be used at all by RAxML. It is certainly better to use phased sequences if you can.

* If diploid genotypes are in the `phased` format, they will be split into haplotypes, and the suffixes '_A' and '_B' will be added to the sample names to distinguish the haplotypes.


* To make neighbour-joining trees, use the Phyml script, and set `--optimise n`, which tells it not to do any ML optimisation.










