# Create Reference for Count Tools

I'll use mm10 as an illustrative example.

## Make Reference with Cell Ranger

First, find reference data from 10X Genomics website for [mm10](https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build#mm10_3.0.0).

Sometimes, the built reference data are not provided, and you are only given instructions to build on your own. Pull the corresponding version of Cell Ranger docker image from [cumulusprod/cellranger](https://hub.docker.com/r/cumulusprod/cellranger/tags), and follow those instructions.

When finished, you'll have the following core files for building references for different count tools:

* ``fasta/genome.fa``
* ``genes/genes.gtf``
* ``star/``

Say the reference is in folder ``/projects/mm10`` on your local machine.

## Build Reference for StarSolo

Pull docker image from [quay.io/cumulus/starsolo](https://quay.io/repository/cumulus/starsolo?tab=tags) for building. In this example, I use tag ``2.7.6a``. Then start a docker container:

```
docker run -it --rm -v /projects/mm10:/ref -v /projects/starsolo-ref:/output quay.io/cumulus/starsolo:2.7.6a
```

Then inside the docker container, type:

```
STAR --runMode genomeGenerate --runThreadN 32 --genomeDir /output --genomeFastaFiles /ref/fasta/genome.fa --sjdbGTFfile /ref/genes/genes.gtf --genomeSAindexNbases 14 --genomeChrBinNbits 18
```

where ``32`` after ``--runThreadN`` option can be changed to adapt to number of available CPUs on your local machine.

When finished, exit docker container. Then type:

```
cd /projects
tar -czvf starsolo.tar.gz starsolo-ref
```

Finally, upload this zipped file to Google Bucket:

```
gsutil -m cp starsolo.tar.gz gs://regev-lab/resources/count_tools/mm10/
```

## Build Reference for Optimus

You don't need to build any new file for Optimus mm10 reference. On your local machine, copy those core files to Optimus reference folder:

```
mkdir /projects/optimus-ref
cp /projects/mm10/fasta/genome.fa /projects/optimus-ref
cp /projects/mm10/genes/genes.gtf /projects/optimus-ref
cp -r /projects/mm10/star projects/optimus-ref
```

Then you have to compress ``star`` folder into ``.tar`` format, which is required by Optimus:

```
cd /projects/optimus-ref
tar -cvf star.tar star
rm -r star
cd ..
```

After that, zip the whole reference folder:

```
tar -czvf optimus.tar.gz optimus-ref
```

And finally, upload this zipped file to Google Bucket:

```
gsutil -m cp optimus.tar.gz gs://regev-lab/resources/count_tools/mm10/
```

## Build Reference for BUSTools

Pull docker image from [shaleklab/kallisto-bustools](https://hub.docker.com/r/shaleklab/kallisto-bustools/tags) for building reference. In this example, I use tag ``0.24.4``. Then start a docker container:

```
docker run -it --rm -v /projects/mm10:/ref -v /projects/bustools-ref:/output shaleklab/kallisto-bustools:0.24.4
```

Then inside the docker container, type:

```
cd /output
kb ref -i transcriptome.idx -g transcripts_to_genes.txt -f1 cdna.fa /ref/fasta/genome.fa /ref/genes/genes.gtf
```

When finished, exit docker container. Then type:

```
cd /projects
tar -czvf bustools.tar.gz bustools-ref
```

Finally, upload this zipped file to Google Bucket:

```
gsutil -m cp bustools.tar.gz gs://regev-lab/resources/count_tools/mm10/
```

## Build Reference for Salmon Alevin

Pull docker image from [cumulusprod/alevin](https://hub.docker.com/r/cumulusprod/alevin/tags) for building reference. In this example, I use tag ``1.1``. Then start a docker container:

```
docker run -it --rm -v /projects/mm10:/ref -v /projects/bustools-ref:/bustools-ref -v /projects/alevin-ref:/output cumulusprod/alevin:1.1
```

Then inside the docker container, type:

```
cp /bustools-ref/cdna.fa /output
cd /output
gzip cdna.fa
salmon index -t cdna.fa.gz -p 32 -i salmon_index
rm cdna.fa.gz
```

The last Salmon index command is based on official tutorial at https://combine-lab.github.io/alevin-tutorial/2018/setting-up-resources/.

We also need a transcriptome to gene map file. This is also from Bustools reference, but we only need 2 columns from it. In python environment, type the following:

```python
import pandas as pd

df = pd.read_csv("/bustools-ref/transcripts_to_genes.txt", sep = '\t', header = None)
df_new = df[[0, 1]]
df_new.to_csv("/output/txp2gene.tsv", sep = '\t', header = False, index = False)
```

When finished, exit docker container. Then type:

```
cd /projects
tar -czvf alevin.tar.gz alevin-ref
```

Finally, upload this zipped file to Google Bucket:

```
gsutil -m cp alevin.tar.gz gs://regev-lab/resources/count_tools/mm10/
```

## Wrap Up

We are almost done. Next step is to update the index file, and add GS URL of mm10 reference to it.

The index file is at Cumulus GitHub repository: https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/workflows/count_tools/ref_index.tsv. Add a new line as follows in this file:

```
mm10	gs://regev-lab/resources/count_tools/mm10
```

where two elements are seperated by tab. Then upload the new index file to Google Bucket:

```
gsutil -m cp ref_index.tsv gs://regev-lab/resources/count_tools/
```

That's it!
