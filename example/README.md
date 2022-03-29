# Quantification of a combinatorial CRISPR screen using pyCROQUET

The python Crispr Read to Oligo QUantification and Evaluation Tool or **pyCROQUET** command line application can be used to quantify guide abundance for single- and dual-guide CRISPR screens.

## Software installation

For this analysis we used the following software and versions:

* [SRA toolkit](https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.11.2/sratoolkit.2.11.2-ubuntu64.tar.gz) version 2.11.2
* [samtools](http://www.htslib.org/) version 1.15
* [pyCROQUET](https://github.com/cancerit/pycroquet) version 1.5.0 

The process used to install each of these tools prior to analysis can be found in [install_software_dependencies.ipynb](install_software_dependencies.ipynb). 

*Note: this is an example of how we installed the various applications - you should defer to the individual application instructions for installation on your system (which may differ).*

## Example dataset

The following paper contains experiments using several different libraries:

>Gonatopoulos-Pournatzis, T., Aregger, M., Brown, K.R. et al.   
>**Genetic interaction mapping and exon-resolution functional genomics with a hybrid Cas9–Cas12a platform.**.   
>*Nat Biotechnol 38, 638–648 (2020).*   
>https://doi.org/10.1038/s41587-020-0437-z

For the purpose of this example, we used pyCROQUET to quantify the digenic library which uses two guides to target paralogues. We selected this dataset as both the raw sequencing data (used for quantification) and raw counts (used for comparison) were available:

* **Raw counts and sample metadata** ([GEO series GSE144281](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144281)) 
* **Raw sequencing (FASTQ) and run metadata** ([SRA project PRJNA603290](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA603290&o=acc_s%3Aa))

To download the example dataset and required metadata, please see [GSE144281_paralogLibrary/metadata/download_and_prepare_example_metadata.ipynb](GSE144281_paralogLibrary/metadata/download_and_prepare_example_metadata.ipynb).

## Implementation

pyCROQUET has several 'modes' which tailor the analysis to the input data provided. For combinatorial CRISPR, we need the pyCROQUET `dual-guide` mode. 

While the script for quantifying the screen is not available with the publication, from the methods we can see that an alignment strategy using Bowtie allowed up to 3 mismatches prior to quantification.

> ***Dual-guide mapping and quantification***    
> *FASTQ files from paired-end sequencing were first processed to trim off flanking sequences up- and downstream of the guide sequence using a custom Perl script. Reads that did not contain the expected 3′ sequence, allowing up to two mismatches, were discarded. Preprocessed paired reads were then aligned to a FASTA file containing the library sequences using Bowtie (v.0.12.7) with the following parameters: -v 3 -l 18 -chunkmbs 256 –t <library_name>. The number of mapped read pairs for each dual-guide construct was then counted and merged, along with annotations, into a matrix.*

In the example, we run pyCROQUET four times, first using exact matching (no mismatches allowed) and then using increasingly less restrictive `rules` to allow 1, 2 or 3 mismatches (1M, 2M or 3M respectively). Finally, we do a simple comparison between the expected and observed counts using R.

Please see the following workbooks:

* Running pyCROQUET: [GSE144281_paralogLibrary/1_quantification_with_pyCROQUET.ipynb](GSE144281_paralogLibrary/1_quantification_with_pyCROQUET.ipynb)
* Comparing observed and expected counts [GSE144281_paralogLibrary/2_compare_expected_and_observed_counts.ipynb](GSE144281_paralogLibrary/2_compare_expected_and_observed_counts.ipynb)

## Questions

Please feel free to post any questions [here](https://github.com/cancerit/pycroquet/issues).
