# Scripts for parsing and manipulating VCF files 

## Contents

* [Parse VCF files](#parse-vcf-files)
* [Transfer VCF to a new coordinate system](#Transfer-VCF-to-a-new-coordinate-system)

___

## Parse VCF files

Most of my scripts use a processed `.vcf` format that I call `.geno`. This looks something like this:

```
#CHOM      POS      ind1      ind2      ind3
scaffold1  1        A/A       A/G       G|A
scaffold1  1        N/N       T/T       T|C
```

Missing data is denoted as `N`, and phased and unphased genotypes are shown conventionally with `|` and `/`.

NOTE if you want to use parallel processing, 

The script `parseVCF.py` will convert vcf to this format. It has various options for filtering based on read depth, genotype quality or any other flag in the `FORMAT` column of the vcf.

The script `parseVCFs.py` (note the 's' on the end) extends the usage of `parseVCF.py` by using [tabix](http://www.htslib.org/doc/tabix.html).
This allows **multi-threaded processing**, *and* **processing and merging of multiple VCFs**. 

#### Example commands:

**Base script (single threaded, and therefore slow)
**
```bash
python parseVCF.py -i input.vcf.gz --skipIndels --minQual 30 --gtf flag=DP min=5 | bgzip > output.geno.gz
```

**Wrapper script (multi threaded, and therefore faster)
**
```bash
python parseVCFs.py -i input1.vcf.gz -i input2.vcf.gz \
--skipIndels --minQual 30 --gtf flag=DP min=5 max=100 --threads 10 |
bgzip > output.geno.gz
```

#### Notes

* You can filter on any flag associated with the genotype, such as read depth (`DP`) and genotype likelihood (`PL`). Check the `FORMAT` column of your vcf file to see which flags are present. To add a genotype filter, add the argument `--gtf` followed by the flag title, minimum and/or maximum value, as in the example above.

* By default the script outputs the genotype for each sample at each site. You can also output a different field, such as `DP`, as long as it is present in the `FORMAT` column of the vcf. e.g. `--field DP`

* `minQual` refers to the `QUAL` column in the vcf, and not individual genotype qualities, for which you should use the `--gtf` argument.

* For parallel processing using tabix, the wrapper `parseVCFs.py` uses the contig lengths provided in the header of the VCF file. If you prefer to skip these, or if there is no header, you can provide a fasta index file using the `--fai` option.
___

## Transfer VCF to a new coordinate system

The script vcfChromTransfer.py takes as input VCF and [agp](https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/) file and outputs a VCF with chromosome names and positions transferred to the new coordinate system.

#### Example commands:

**Base script (single threaded, and therefore slow)
**
```bash
python vcfChromTransfer.py -v input.vcf.gz -a transfers.agp | bgzip > output.vcf.gz
```
#### Notes

* You will need to have [tabix](http://www.htslib.org/doc/tabix.html) installed and the input file must be indexed. You can do that with the command tabix `-p vcf input.vcf.gz`

* You can specify to only output one or a few chromosomes listed in the transfers file with the `--chroms` option. 

* If you don't have an agp. An alternative is to provide a text transfers file (`-t` option) with the columns:
1: New chromosome name
2: Start position on new chromosome
3: Snd position on new chromosome
4: Chromosome name in input file
5: Start position in input file
6: End position in input file
7: Strand of new chromosome relative to input file (+/-)
