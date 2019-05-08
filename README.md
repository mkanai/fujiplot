# :mount_fuji: Fuji plot :mount_fuji:
a circos representation of multiple GWAS results.

## Overview
Fuji plot (a circos representation of multiple GWAS results) was developed to efficiently summarize significant SNPs across dozens of traits ([Kanai, M. *et al*., *Nat. Genet.* 2018](http://dx.doi.org/10.1038/s41588-018-0047-6)). Technically, it is a wrapper script for [Circos](http://circos.ca/), generating input data and configuration files from GWAS results.

<p align="center"><img src="http://mkanai.github.io/img/Kanai2018_Fig1.svg" width="640px"></p>

Fuji plot is named after :mount_fuji: Mt. Fuji in Japan :mount_fuji:, whose ridge and skirts greatly resemble our figure!

<p align="center"><img src="https://i.imgur.com/QXlnKZk.jpg" width="480px"></p>

## Requirements
* [R](https://www.r-project.org/) with [dplyr](https://github.com/tidyverse/dplyr) and [stringr](https://github.com/tidyverse/stringr) packages.
  * To install these packages, type `install.packages("tidyverse")` in R console.
* [Circos](http://circos.ca/)
  * Please follow [their instructions](http://circos.ca/documentation/tutorials/configuration/installation/) to install Circos.
  * For the current version ([v0.69-6](http://circos.ca/distribution/circos-0.69-6.tgz)), you need to apply our [patch](https://gist.github.com/mkanai/be05f40f933112bfb70bb08076cdaa00) to avoid an infinite loop in text arrangements. Alternatively, you can download and use [our patched version](https://www.dropbox.com/s/z6jdwhj0o570fp8/circos-0.69-6-kanai.tgz?dl=0). I will report this issue soon.

## Usage
Please add `circos` to your `PATH` or set `CIRCOS_PATH` in `fujiplot.R` appropriately.

```{sh}
Rscript fujiplot.R input_example/input.txt input_example/traitlist.txt
```

## Inputs

Input files used in [Kanai, M. *et al.* (2018)](http://dx.doi.org/10.1038/s41588-018-0047-6) are provided under `./input_example`.

### Significant trait-associated loci (`input.txt`)

This file provides a concatenated list of all significant *trait-associated* loci (for each trait). It also defines a unique locus across multiple traits to highlight pleiotropic loci. The required fields are as follows:

* `LOCUS_ID`: Unique ID for each locus
* `CATEGORY`: Trait category
* `TRAIT`: Trait name
* `CHR`: Chromosome
* `BP`: Position (hg19)
* `MARKER`: Marker ID
* `GENE`: Gene symbol


### Trait list (`traitlist.txt`)

This file provides a list of traits and their categories. It defines an order of traits and a color of each category in a figure. The required fields are as follows:

* `CATEGORY`: Trait category
* `TRAIT`: Trait name
* `COLOR`: Category color

## Outputs

All output files will be created under `./output`.


### Circos outputs (`circos.{png,svg}`)
This is the raw output from Circos. Unfortunately, we could not figure out the appropriate way to make our [highlights](http://circos.ca/documentation/tutorials/highlights/) (radial lines) over [scatter plots](http://circos.ca/documentation/tutorials/2d_tracks/scatter_plots/). To this end, we edited the svg file using Adobe Illustrator and produced the final figure.

<p align="center"><img src="output/circos.svg" width="640px"></p>


### Optional: a bar plot (`barplot.pdf`)

This corresponds to Fig. 1a of [Kanai, M. *et al.* (2018)](http://dx.doi.org/10.1038/s41588-018-0047-6). Again, you need to do some Illustrator edits to make it pretty. You can suppress this output by setting `OUTPUT_BARPLOT = FALSE` in `fujiplot.R`.


## Citation
When using Fuji plot, please cite the following paper.

* Kanai, M., *et al*. Genetic analysis of quantitative traits in the Japanese population links cell types to complex human diseases. *Nat. Genet.* (2018) [doi:10.1038/s41588-018-0047-6](http://dx.doi.org/10.1038/s41588-018-0047-6)

## Contact
Masahiro Kanai (mkanai@g.harvard.edu)

http://mkanai.github.io/

## Links
* [JENGER](http://jenger.riken.jp/en/) (the lab website)
* [The BioBank Japan Project](https://biobankjp.org/english/index.html)
* [RIKEN Center for Integrative Medical Sciences](http://www.ims.riken.jp/english/)
* [National Bioscience Database Center Human Database](https://humandbs.biosciencedbc.jp/en/)
* [Circos](http://circos.ca/)
* [ClicO FS](http://codoncloud.com:3000/)
