Current development is awaiting a release use `minnow-velocity` for using the recent version of minnow 

# Minnow ( read level simulator for dscRNA-seq data)


Most analysis pipelines validate their results using known marker genes (which are not widely available for all types of analysis) and by using simulated data from gene-count-level simulators. Typically, the impact of using different read-alignment or UMI deduplication methods has not been widely explored. Assessments based on simulation tend to start at the level of assuming a simulated count matrix, ignoring the effect that different approaches for resolving UMI counts from the raw read data may produce. Here, we present minnow, a comprehensive sequence-level droplet-based single-cell RNA-seq (dscRNA-seq) experiment simulation framework.  Minnow accounts for important sequence-level characteristics of experimental scRNA-seq datasets and models effects such as PCR amplification,  CB (cellular barcodes) and UMI (Unique Molecule Identifiers) selection, and sequence fragmentation and sequencing. 

Minnow is a read level simulator for droplet based single cell RNA-seq data. Minnow simulates the reads by sampling sequences from the underlying de-Bruijn graph (using `--dbg`) of the reference transcriptome or alternatively just samples sequences from the reference transcriptome. As the `--dbg` option also enables other features of the software, it is useful to describe those.

<p align="center">
<img src="data/minnow_main_figure.001.jpeg">
</p>

# 'Updated tutorial to Minnow' (authors: Hirak Sarkar and Dongze He)




# Install the current active branch `minnow-velocity`
Minnow-velocity is currently available on Minnow's github repo under the minnow-velocity branch, for downloading, please
1. Open your Terminal and go to the folder you want to download minnow in.
2. Run the following codes:
``` shell
git clone --single-branch --branch minnow-velocity https://github.com/COMBINE-lab/minnow.git
cd minnow
mkdir build
cd build
cmake ..
make
```

**NOTICE: In the following steps, we assume that your working directory is in the `minnow/build` folder.**

# Step 0 -- Run Alevin for a dataset of the orgamism you will work on.
As Minnow need the BFH file dumpped from Alevin to make its probability files, we need to run alevin before simulation with `--dumpBfh` option specified.

### Step 0.1 -- build salmon index
Download the reference files from [GENCODE](https://www.gencodegenes.org), and build Salmon index.
```shell
mkdir data

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.pc_transcripts.fa.gz
gunzip data/gencode.vM25.pc_transcripts.fa.gz

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz
gunzip data/gencode.vM25.annotation.gtf.gz

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz
gunzip data/GRCm38.primary_assembly.genome.fa.gz

grep ">" GRCm38.primary_assembly.genome.fa | cut -d ">" -f 2 | cut -d " " -f 1 > GRCm38.primary_assembly.genome.chrnames.txt

cd ..

salmon index \
-t gencode.vM25.pc_transcripts.fa \
-i gencode.vM25.annotation.sidx --gencode -p 128 \
-d GRCm38.primary_assembly.genome.chrnames.txt
```

### Step 0.2 -- Run Alevin on a mouse dataset

1. Get data set from [Hermann et al., 2018](https://www.sciencedirect.com/science/article/pii/S2211124718316024?via%3Dihub). Here we use [bamtofastq](https://github.com/10XGenomics/bamtofastq) from 10X Genomics, please load before you run the code.

```shell
cd data
wget https://sra-pub-src-1.s3.amazonaws.com/SRR6459157/AdultMouse_Rep3_possorted_genome_bam.bam.1

mv AdultMouse_Rep3_possorted_genome_bam.bam.1 AdultMouse_Rep3_possorted_genome_bam.bam

bamtofastq --reads-per-fastq=500000000 AdultMouse_Rep3_possorted_genome_bam.bam FASTQtmp

mv FASTQtmp/Ad-Ms-Total-Sorted_20k_count_MissingLibrary_1_HK2GNBBXX/bamtofastq_S1_L006_I1_001.fastq.gz AdultMouseRep3_S1_L001_I1_001.fastq.gz

mv FASTQtmp/Ad-Ms-Total-Sorted_20k_count_MissingLibrary_1_HK2GNBBXX/bamtofastq_S1_L006_R1_001.fastq.gz AdultMouseRep3_S1_L001_R1_001.fastq.gz

mv FASTQtmp/Ad-Ms-Total-Sorted_20k_count_MissingLibrary_1_HK2GNBBXX/bamtofastq_S1_L006_R2_001.fastq.gz AdultMouseRep3_S1_L001_R2_001.fastq.gz
cd ..
```

2. Run Alevin to quantify the gene abundances based on the index generated above.

```shell
awk -F "\t" '$3 == "transcript" { print $9 }' data/gencode.vM25.annotation.gtf | tr -d ";\"" | awk '{print $4"\t"$2}' > data/gencode.vM25.annotation.tx2gene.tsv

salmon alevin -l ISR -i gencode.vM25.annotation.sidx \
-1 data/AdultMouseRep3_S1_L001_R1_001.fastq.gz \
-2 data/AdultMouseRep3_S1_L001_R2_001.fastq.gz \
-o alevin_out -p 36 \
--tgMap data/gencode.vM25.annotation.tx2gene.tsv \
--chromium \
--dumpFeatures --expectCells 1850 \
--dumpBfh
```

By specifying `--dumpBfh`, alevin will dump the BFH file in its result, which will be used in the following steps.


## Step 1 -- `minnow index`

The indexing phase of minnow is mostly concerned with creating de-Bruijn graph for the transcriptome, additionally to run minnow we need a transcript to gene mapping which can be obtained from gtf file. At the end, minnow index command will create a direcory with updated fasta file that has the following features:
- Any sequence in the fasta file with smaller length than the read length will be removed.
- PolyA tails will be clipped.
- `N`s will be replaced by random `A/T/G/C`.

## 1.1 De-Bruijn graph construction using `minnow index`

Create the de-Buijn graph on your **transcriptome file**  using the following command, note that **do not** pass your genome file into `minnow index`.

``` shell
src/minnow index -r data/gencode.vM25.pc_transcripts.fa -k 101 -f 20 --tmpdir tmp -p 10 -o minnow_index
```

There will be a `ref_k101_fixed.fa` in the `minnow index` folder, this fasta file is the updated fasta file we talked above. 


Followings are the options,

```
SYNOPSIS
        minnow index -r <ref_file>... -o <output_dir> [-f <filt_size>] [--tmpdir <twopaco_tmp_dir>] [-k <kmer_length>] [-p <threads>]

OPTIONS
        -r, --ref <ref_file>
                    path to the reference transcriptome fasta file

        -o, --output <output_dir>
                    directory where index is written

        -f, --filt-size <filt_size>
                    filter size to pass to TwoPaCo when building the reference dBG

        --tmpdir <twopaco_tmp_dir>
                    temporary work directory to pass to TwoPaCo when building the reference dBG

        -k, --klen <kmer_length>
                    length of the k-mer with which the dBG was built (default = 101)

        -p, --threads <threads>
                    total number of threads to use for building MPHF (default = 16)
```

While creating the index, `-k` option should be treated as `read length` for rest of the simulation.

## 1.2  Clean the transcriptome file header (optional)

If you just want to see the transcript name in your fasta file, please run the following option to clean the header in this file.
```shell
sed -e '/^>/ s/|.*//' data/gencode.vM25.pc_transcripts.fa > \
gencode.vM25.pc_transcripts_cleanHeader.fa
```


## Step 2 -- `minnow estimate`

You can find the `bfh.txt` file in the `alevin_out/alevin` folder, this file is what we need for `minnow estimate` to get geneProb and cellProb files, which specify the multimapping probability.

To run minnow estimate, please run 

```shell
src/minnow estimate -o minnow_estimate -r data/gencode.vM25.pc_transcripts.fa --g2t data/gencode.vM25.annotation.tx2gene.tsv --bfh alevin_out/alevin/bfh.txt
```
The `countProb.txt` file will be used in the following steps.




## Step 3 -- Prepare a read count matrix

Here we show how to use [splatter package](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1305-0) in R to simulate the read count matrix, you can [download this package from Bioconductor](https://www.bioconductor.org/packages/release/bioc/html/splatter.html).

1. create 100 genes with 150 cells data
Splatter is a package for the simulation of single-cell RNA sequencing count data, following codes are ran in R.

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("splatter")

library(splatter)

out_dir <- "data/test_splatter_data"
num_genes <- 100
num_cells <- 150
sim <- splatSimulate( 
	nGenes=num_genes, 
	batchCells=num_cells, 
	verbose = FALSE
)

dir.create(out_dir, showWarnings=FALSE)
write.table(colnames(sim), file= file.path(out_dir, "quants_mat_cols.txt"), quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(counts(sim), file= file.path(out_dir, "quants_mat.csv"), quote=FALSE, col.names=FALSE, row.names=FALSE, sep=",")  

```

2. Specify gene names for simulated data.
As we want to simulate date with real gene names, while splatter only returns meaningless gene names (for instance, gene1, gene2, etc), we need to randomly select gene names from our reference files such that the selected genes have at least one isoform with length longer than read length. Fortunately, if you have run `minnow estimate`, you will have a fixed reference file `ref_k101_fixed.fa` which contains only transcripts that meet this criteria.

```
grep ">" minnow_index/ref_k101_fixed.fa | awk -F "|" '{print $2}'|shuf -n 100 > data/test_splatter_data/quants_mat_rows.txt
```

If you are using Mac OS and `shuf` is not found, please download it from [Homebrew](https://brew.sh) by `brew install coreutils` or use a Linux machine.



## Step 4 -- `minnow simulate`

Now we can start our simulation process. There are two mode in minnow, here we give the example of ``--splatter-mode`.

We can run minnow simulate in `--splatter-mode` by running the following codes:

```shell
src/minnow simulate --splatter-mode \
--g2t  data/gencode.vM25.annotation.tx2gene.tsv \
--inputdir data/test_splatter_data \
--PCR 4 \
-r minnow_index/ref_k101_fixed.fa \
-e 0.01 \
-p 16 \
-o minnow_simulate \
--dbg \
--gfa minnow_index/dbg.gfa \
-w data/737K-august-2016.txt \
--countProb minnow_estimate/countProb.txt \
--custom \
--gencode 
```

Minnow generate two `fastq` files as the simulated read files, which can be analyzed by Alevin directly!


## (Optional) Step 5 -- Run Alevin on minnow simulated data

we will use alevin to re-estimate the read count matrix, all required files have already been downloaded in the `data` folder.

```
salmon alevin -l ISR -i gencode.vM25.annotation.sidx \
-1 minnow_simulate/hg_100_S1_L001_R1_001.fastq.gz \
-2 minnow_simulate/hg_100_S1_L001_R2_001.fastq.gz \
-o minnow_simulate/alevin_out -p 16 \
--tgMap data/gencode.vM25.annotation.tx2gene.tsv \
--chromium \
--dumpFeatures --expectCells 100 \
```

Now we are done with the process!

The truth read count matrix is in the `minnow_simulate/alevin` folder, and the re-estimated read count matrix is in the `minnow_simulate/alevin_out/alevin` folder, check it now!


### Things to be added
1. [x] Doublets
2. [ ] Empty-drops
3. [ ] Retained intron
4. [x] Clusters 
