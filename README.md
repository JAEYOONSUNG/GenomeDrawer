# GenomeDrawer
## Circos-based genome map visualization tool
![partitioning](https://github.com/user-attachments/assets/e72013fb-1b89-4085-8009-bde4b2d0120c)
![Untitled-3](https://github.com/user-attachments/assets/5c7f041d-9e7b-47e7-8605-82917d135192)




## Project Introduction

Welcome to our GitHub repository, where we're excited to share a series of workflows depicted genome map. This repository is composed of various scripts, each tailored to specific tasks within our broader research framework.

- Maintainers: Jae-Yoon Sung
- Contributors: Jae-Yoon Sung

## Key Features
User-Friendly Documentation: Detailed documentation is available to guide you through the installation, setup, and utilization of both the scripts and the database.

### Algorithms for analysis
- GC skew: Nucleotides divided by window 10,000 bp for GC skew calculation. GC skew is calculated as: GC skew was calculated using the formula (nC - nG)/(nC + nG), where nC and nG represent the counts of cytosine and guanine bases, respectively.

- COG analysis: COG analysis was performed using the eggNOG mapper. In cases where multiple COGs were assigned, only the first COG was used as the representative for mapping.


## Getting Started:
To begin using our resources, please follow the steps outlined in my instructions. 
Whether you're looking to integrate our scripts into your existing projects or explore our database for new insights, we've provided all the necessary instructions to get you started.

## Installation
### Requirements

The Genome_drawer is supported for macOS, Linux and Windows machines, which can provide an environment for using R.
It requires R version >=4.2.1 for release, and R version >=4.3 for devel.

The [EggNOG-mapper webserver](http://eggnog-mapper.embl.de), allows users to input sequences in FASTA format based on locus_tag identifiers and receive results in either XLSX or CSV format. Additionally, the standalone version available on GitHub is compatible with DNMB.

To download and install R, see the [R-project website](https://www.r-project.org/).

To download and install EggNOG-mapper, see the [EggNOG-mapper github](https://github.com/eggnogdb/eggnog-mapper).


#### Warning
The basic file for genomic analysis, known as a GenBank file, requires both sequence and annotation in full-format files such as gbff, gb, or gbk. Additionally, GenBank prefers a format based on the GeneMarkS2+ pipeline, and using a different annotation pipeline to obtain GenBank files may lead to errors.


## Anaylsis flow
## Prerequisites
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("Biostrings", "ComplexHeatmap", "Peptides"))

install.packages(c("qdap", "seqinr", "stringr", "stringi", "splitstackshape", "gtools", "ggplot2", "ggseqlogo", "circlize", "grid", "gridExtra","plyr", "dplyr", "tidyr", "readr", "reshape2", "data.table", "tibble", "openxlsx"))
```
- **Note:** If you encounter issues installing the qdap package, try installing it with the following command:
```r
install.packages("qdap", INSTALL_opts = "--no-multiarch")
```

   
- **Note:** Java and the rJava package must be installed and configured to enable .xlsx output using this package.
1. Install Java Development Kit (JDK):
        Download and install the appropriate JDK for your operating system from the Oracle website or OpenJDK.
2. Install rJava Package in R:
```r
install.packages("rJava")
library(rJava)
```

3. Set $JAVA_HOME Path:

You need to set the environment variable JAVA_HOME to point to the location of your JDK installation.
	
 •	On Windows:

1.	Install Java jdk (https://www.oracle.com/kr/java/technologies/downloads/)
2.	Check the “System Variables,” :
	•	Variable name: JAVA_HOME
	•	Variable value: The path to your JDK installation 

```bash
# Print the current value of the JAVA_HOME environment variable.
echo %JAVA_HOME%  #(e.g., C:\Program Files\Java\jdk-18)

# Set the JAVA_HOME environment variable to point to the Java Development Kit (JDK) installation.
# Replace [version] with your installed JDK version (e.g., jdk-18).
setx JAVA_HOME "C:\Program Files\Java\jdk[version]"

# Update the system PATH to include the bin directory of the JDK.
setx PATH "%JAVA_HOME%\bin;%path%"

# Check the installed Java version to confirm that the correct version is being used.
JAVA -version
```
3.	Restart R or RStudio.
	
 •	On macOS/Linux:
- If Xcode is not installed, you may encounter compiler issues during package installation. To resolve this, install Xcode from the App Store.
Add this line to your .bash_profile or .bashrc (depending on the shell):
```bash
# Navigate to your Java installation directory to check available Java versions
/Library/Java/JavaVirtualMachines/[my_java_folder]/Contents/Home # check my java list

# Open your .bash_profile (or .bashrc) file for editing
vi ~/.bash_profile # edit bash profile

# Press 'i' to enter insert mode in the vi editor
i # insert mode

# Add or update the JAVA_HOME environment variable with the path to your Java installation
export JAVA_HOME=/Library/Java/JavaVirtualMachines/[my_java_folder]/Contents/Home
# Add Java's bin directory to the system PATH variable so that Java commands can be run from the terminal
export PATH=${PATH}:$JAVA_HOME/bin

# Save the changes and exit the vi editor. ":wq!" means "write" (save) and "quit" (exit) forcefully
: # activate command line
wq! # save

# Apply the changes made to the .bash_profile or .bashrc immediately (without needing to restart the terminal)
source ~/.bash_profile  ## or ~/.bashrc #apply changes

# Verify that JAVA_HOME is set correctly by printing its value
echo $JAVA_HOME # validation
```

## Install GenomeDrawer R package
```r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
    
devtools::install_github("JAEYOONSUNG/GenomeDrawer")
```


## Run GenomeDrawer analysis
```r
setwd([GenBank and gff directory]) # Set the working directory to the location where your GenBank and gff files are stored.
library(GenomeDrawer)
run_drawer()
```
- Circular Representation:
 When completeness = TRUE, each contig is visualized as a separate circular map, representing it as an independent entity.
 This is particularly useful for analyzing genomes with multiple contigs that are not linearized or merged into a single chromosome.
- Contig-Length Based Linearization:
 When completeness = FALSE, all contigs are displayed within a single circular genome map, with each contig’s length proportionally calculated and represented as a segment within the map.
 This is ideal for visualizing fragmented assemblies in a comparative genome-wide context.
- **Note:** [Strain of interest].gff for merging data.

---
***Requirements***
**📄 Required GFF File Format**
Please ensure your GFF3 file follows the standard format, similar to those generated by NCBI RefSeq. Below is a minimal example with necessary headers and features:
```
##gff-version 3
#!gff-spec-version 1.21
#!processor NCBI annotwriter
#!genome-build ASM3037676v1
#!genome-build-accession NCBI_Assembly:GCF_030376765.1
#!annotation-date 07/01/2023 07:05:56
#!annotation-source NCBI RefSeq
##sequence-region NZ_CP128453.1 1 3650785
##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=1422
NZ_CP128453.1	RefSeq	region	1	3650785	.	+	.	ID=NZ_CP128453.1:1..3650785;Dbxref=taxon:1422;Is_circular=true;Name=ANONYMOUS;strain=EF60045
NZ_CP128453.1	RefSeq	gene	1	1350	.	+	.	ID=gene-QT235_RS00005;Name=dnaA;gene=dnaA;gene_biotype=protein_coding;locus_tag=QT235_RS00005
NZ_CP128453.1	Protein Homology	CDS	1	1350	.	+	0	ID=cds-WP_033027785.1;Parent=gene-QT235_RS00005;Name=WP_033027785.1;product=chromosomal replication initiator protein DnaA;protein_id=WP_033027785.1
```
- ##gff-version 3 must be declared at the top.
- Include genome metadata such as genome-build and annotation source.
- Feature entries should include region, gene, and CDS.
- Essential attributes: ID, Parent, Name, product, protein_id, etc.
---
**EggNOG-mapper** (Essential)
COG (Clusters of Orthologous Groups) annotations were performed using EggNOG-mapper[http://eggnog-mapper.embl.de], either via:
	•	Online version: http://eggnog-mapper.embl.de
	•	Local version: using the emapper.py script from the EggNOG-mapper package
```python
emapper.py --cpu 20 --mp_start_method forkserver --data_dir [eggnog_data directory] -o out --output_dir [eggnog_output] --temp_dir [eggnog_output] --override -m diamond --dmnd_ignore_warnings --dmnd_algo ctg -i [fasta] --evalue 0.001 --score 60 --pident 40 --query_cover 20 --subject_cover 20 --itype proteins --tax_scope auto --target_orthologs all --go_evidence non-electronic --pfam_realign none --report_orthologs --decorate_gff yes --excel

```
- **Note:** [Strain of interest].emapper.annotations.xlsx or [Strain of interest]emapper.annotations.csv output are used for merging data.


## Contributing
We welcome contributions from the community! If you have suggestions for improvements, additional scripts, or updates to the database, please see our contributing guidelines for more information on how to get involved.

## License
This project is released under MIT licence, which allows for both personal and commercial use, modification, and distribution of our work, provided that proper credit is given.

We hope our resources will prove invaluable to your research in systems biology. For any questions or feedback, please don't hesitate to reach out through our GitHub issues or contact section.

## Citation
If you use this piepline, please cite:
```
[DNMB] DNMB: A Strategic Blueprint for the Domestication of Geobacillus stearothermophilus as a Thermophilic Platform using the DNMB Suite.
             Jae-Yoon Sung, Mun Hoe Lee, Hyungbin Kim, Dariimaa Ganbat, Hyun-Woo Cho, Sang Jae Lee, Seong Bo Kim, and Dong-Woo Lee. 2024.
             XXX, XXX, https://doi.org/XXX
```
Please, cite also the underlying algorithm if it was used for the search step of DNMB:
```
[eggNOG-mapper v2] eggNOG-mapper v2: functional annotation, orthology assignments, and domain 
                   prediction at the metagenomic scale. Carlos P. Cantalapiedra, 
                   Ana Hernandez-Plaza, Ivica Letunic, Peer Bork, Jaime Huerta-Cepas. 2021.
                   Molecular Biology and Evolution, msab293, https://doi.org/10.1093/molbev/msab293

[Circos] Circos: an information aesthetic for comparative genomics.
	Martin Krzywinski  1 , Jacqueline Schein, Inanç Birol, Joseph Connors,
	Randy Gascoyne, Doug Horsman, Steven J Jones, Marco A Marra
	Genome research, 2009, 19.9: 1639-1645. https://doi: 10.1101/gr.092759.109
```

