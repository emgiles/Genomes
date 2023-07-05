# GENOME ASSEMBLY AND ANNOTATION

### DATA PREPRATION
Create working directory with sequence files
Unzip fastq.gz subread files
### GENOME ASSEMBLY
Canu took several hours to run a machine with 16 cores
Run assembly using fastq files using Canu
```canu  -p project_prefix -d project_directory genomeSize=500m -pacbio pacbio_file.fastq```
.contigs.fasta is final filtered assembly results for diploid genomes
.unitigs.fasta is file that includes all alternative paths
.unassembled.fasta is a file of contigs with poor support
.report gives statistics of the raw data and final assembly including N50s and length
###### Check assembly with QUAST
```conda activate quast_purgedups_minimap```
```nohup quast -r /media/pablo/data1/emily/genomes/SG1_Dovetail_assembly/Scurria_scurra_uni1728-mb-hirise-3bs35_08-29-2020__hic_output.fasta -s -e --large -k -f -b VG2_canu_assembly.contigs.fasta & ```
Compare results to S. scurra assembly
Output is report.pdf

###### Check assembly with BUSCO
This example uses the mollusca database from ncbi. Also try with metazoa db.
```conda activate busco```
```nohup busco -o VG2_busco_mollusca -i /media/pablo/data1/emily/genomes/VG1_assembly/VG2_canu_assembly.contigs.fasta -l mollusca_odb10 -m genome & ```

### GENOME ASSEMBLY POST-PROCESSING
##### Remove duplications with purge_dups
See [PIPELINE INSTRUCTIONS in github](https://github.com/dfguan/purge_dups). 
Note, this software is difficult to install and run due to poor documentation online.
```conda activate quast_purgedups_minimap```

##### Check purged assembly with QUAST
```conda activate quast_purgedups_minimap```
```/home/pablo/anaconda2/envs/emily/bin/quast -r /media/pablo/data1/emily/genomes/SG1_Dovetail_assembly/Scurria_scurra_uni1728-mb-hirise-3bs35_08-29-2020__hic_output.fasta -s -e --large -k VG2_canu_purged_purged.fa```
comparing to S. scurra
output is report.pdf

##### Check purged assembly with BUSCO
```conda activate busco```
```nohup busco -o VG2_busco_mollusca -i /media/pablo/data1/emily/genomes/VG1_assembly/VG2_canu_purged_purged.fa -l mollusca_odb10 -m genome &```
using mollusca db

##### Scaffold with RAGTAG
```conda activate ragtag_blobtools_bwa```
```nohup ragtag.py scaffold /media/pablo/data1/emily/genomes/SG1_Dovetail_assembly/Scurria_scurra_uni1728-mb-hirise-3bs35_08-29-2020__hic_output.fasta /media/pablo/data1/emily/genomes/VG1_assembly/purge_dups_21oct2022/purged.fa &```
using S. scurra as a reference

###### Checking purged scaffolded assembly with QUAST
```conda activate quast_purgedups_minimap```
```nohup quast -r /media/pablo/data1/emily/genomes/SG1_Dovetail_assembly/Scurria_scurra_uni1728-mb-hirise-3bs35_08-29-2020__hic_output.fasta -e --large -k -f -b /media/pablo/data1/emily/genomes/VG1_assembly/RAGTAG4/ragtag_output/ragtag.scaffold.fasta &```

###### Checking purged scaffolded assembly with BUSCO
```conda activate busco```
```nohup busco -o VG2_busco_mollusca -i /media/pablo/data1/emily/genomes/VG1_assembly/VG2_canu_purged_purged_ragtag.scaffold.fasta -l mollusca_odb10 -m genome &```
(busco) pablo@pablo-Precision-5820-Tower:/media/pablo/data1/emily/genomes/VG1_assembly/QUAST_purged_woe_scaffolded/quast_results/results_2022_10_27_12_12_08/busco_stats$ 
```nohup busco -o VG2_busco_metazoa -i /media/pablo/data1/emily/genomes/VG1_assembly/VG2_canu_purged_purged_ragtag.scaffold.fasta -l metazoa_odb10 -m genome &```

### TRANSCRIPTOME AND PROTEIN LIST RECOVERY
##### Upload transcriptomes and check with busco
```conda activate busco```
```nohup busco -o SV_trinity_metazoa -i /media/pablo/data1/emily/genomes/VG1_assembly/transcriptomes/Ensamble_transcriptoma_SViridula.fasta -l metazoa_odb10 -m transcriptome &```
Should have high busco completeness. Duplications can be moderately high.
##### Check quality of raw reads
```conda activate fastqc```
```Fastqc 30j_1.fq.gz```
If the raw reads are of low quality be cautious about using transcriptome.
##### Align transcriptome to genome
HiSat and BWA-MEM are not recommended for alignments of transciptomes to genomes. Transcriptomes have multiple transcripts for the same genomic region. These types of alignments are best done with GMAP which can be installed via conda. 
Example command:
```gmap_build -d dbname  -D dbloc dbfasta```
GMAP recomends creating objects for -d and -D and dbfasta. This did not work so I ran with full paths.
###### Build GMAP database
```conda activate gmap```
```gmap_build -d viridula_masked -D /media/pablo/data1/emily/genomes/VG1_assembly/gmap/ /media/pablo/data1/Emily/genomes/VG1_assembly/VG2_canu_purged_purged_ragtag.scaffold.fasta.masked.nosemicolon  File was written to /media/pablo/data1/Emily/genomes/VG1_assembly/gmap//viridula_masked/viridula_masked.salcpchilddc```
###### Align transcripts with output to .sam
```gmap -D /media/pablo/data1/emily/genomes/VG1_assembly/gmap/ -d viridula_masked -B 5 -t 10 --input-buffer-size=1000000 --output-buffer-size=1000000 -f samse /media/pablo/data1/emily/genomes/VG1_assembly/transcriptomes/transcriptome2/Ensamble_s_viridula_31_y_143.fasta > viridula_masked.nosemicolon.transcriptome2.sam```
###### Use samtools (installed in another conda env) to convert from .sam to .bam
```conda activate samtools_blast```
```samtools view -S -b viridula_masked.nosemicolon.transcriptome2.sam > viridula_masked.nosemicolon.transcriptome2.bam```
###### Sort .bam
```samtools sort -t 10 -o viridula_masked.nosemicolon.transcriptome2_sorted.bam -T Gmap_temp viridula_masked.nosemicolon.transcriptome2.bam```
###### Check stats of alignment
```samtools flagstat viridula_masked.nosemicolon.transcriptome2_sorted.bam```
It is important that a high percentage of the transcripts map to the genome. In the case of viridula, 96.67% mapped.
##### Extract lists of proteins from Transcriptome
First, install TransDecoder and all dependencies.
###### Extract long open reading frames from transcriptome
```conda activate transdecoder```
```TransDecoder.LongOrfs -t Ensamble_s_viridula_31_y_143.fasta```
###### Predict the likely coding regions
```TransDecoder.Predict -t Ensamble_s_viridula_31_y_143.fasta```
Transdecoder finishes correctly if .bed .cds .gff3 and .pep files are created
###### Remove stop codons (*)
```conda activate sed```
```sed -e 's/\*//g' Ensamble_s_viridula_31_y_143.fasta.transdecoder.pep > Ensamble_s_viridula_31_y_143.fasta.transdecoder.faa```
###### Cluster and remove sequences that are 95% similar to an existing sequence
Install cd-hit via conda
```conda activate cd-hit```
```cd-hit -i Ensamble_s_viridula_31_y_143.fasta.transdecoder.faa -o Ensamble_s_viridula_31_y_143.fasta.transdecoder_CDHIT_95.faa -c 0.95```
###### remove special characters from headers
```conda activate sed```
```sed -e 's/:/_/g' Ensamble_s_viridula_31_y_143.fasta.transdecoder_CDHIT_95.faa | sed -e 's/=/_/g' | sed -e 's/:/_/g' >Ensamble_s_viridula_31_y_143.fasta.transdecoder_CDHIT_95_nospcharac.faa```
###### remove tildas from headers
```sed -e 's/~/_/g' Ensamble_s_viridula_31_y_143.fasta.transdecoder_CDHIT_95_nospcharac.faa > Ensamble_s_viridula_31_y_143.fasta.transdecoder_CDHIT_95_nospcharac2.faa``
### REPEAT MASKING
First install RepeatModeler and all associated software. Installation via conda didn't work, had to install manually. Installation takes a long time, many dependencies. See [this tutorial](https://darencard.net/blog/2022-10-13-install-repeat-modeler-masker/) for installation guide.

##### Build database for repeat modeling
```nohup /home/pablo/anaconda2/envs/repeatmodeler02/repeat-annotation/RepeatModeler-2.0.3/BuildDatabase -name viridula /media/pablo/data1/emily/genomes/VG1_assembly/VG2_canu_purged_purged_ragtag.scaffold.fasta```
This is quite fast, it is basically indexing the reference genome

##### Run RepeatModeler
```nohup /home/pablo/anaconda2/envs/repeatmodeler02/repeat-annotation/RepeatModeler-2.0.3/RepeatModeler -pa 8 -database viridula &```
This takes a long time... 72+hr to run
Important output files are viridula-families.fa and viridula-families.stk

##### Run RepeatClassifier
```nohup /home/pablo/anaconda2/envs/repeatmodeler02/repeat-annotation/RepeatModeler-2.0.3/RepeatClassifier -consensi viridula-families.fa -stockholm viridula-families.stk -engine ncbi &```
This takes 12-24hr.
Important output files are viridula-families-classified.stk and viridula-families.fa.classified

##### Run RepeatMasker
```nohup /home/pablo/anaconda2/envs/repeatmodeler02/repeat-annotation/RepeatMasker/RepeatMasker -pa 8 -lib /media/pablo/data1/emily/genomes/VG1_assembly/repeatmask04/viridula-families.fa.classified /media/pablo/data1/emily/genomes/VG1_assembly/VG2_canu_purged_purged_ragtag.scaffold.fasta -dir softmask_viridula -xsmall -poly &```
output files are placed in folder softmask_viridula. Most important output is the .tbl file that contains a summary of the masked elements and .fasta.masked file that is genome assembly softmasked.
This takes several hours to half a day.
copy softmasked assembly to working directory for postprocessing.

### PREPARE FOR ANNOTATION 
##### Remove all semicolons from sequence names
Install program sed previously using conda.
```conda activate sed```
```sed -e 's/;/_/g' VG2_canu_purged_purged_ragtag.scaffold.fasta.masked > VG2_canu_purged_purged_ragtag.scaffold.fasta.masked.nosemicolon```
very quick, no need for nohup
output file is file.semicolon

##### Acquire BRAKER via Singularity
BRAKER can also be installed manually, but installation of dependency AUGUSTUS is difficult. I decided to go for the Singularity container. Singularity must be installed previously.

Build BRAKER via Singularity in /home/usr
```singularity build braker3.sif docker://teambraker/braker3:latest```
Make BRAKER executable in /home/usr
```singularity exec braker3.sif braker.pl```
Install license key for GeneMark-ETP in home directory in /home/usr
```singularity exec braker3.sif print_braker3_setup.py```
Get test script 1 for BRAKER via Singularity in /home/usr
```singularity exec -B $PWD:$PWD braker3.sif cp /opt/BRAKER/example/singularity-tests/test1.sh .```
Get test script 2 for BRAKER via Singularity in /home/usr
```singularity exec -B $PWD:$PWD braker3.sif cp /opt/BRAKER/example/singularity-tests/test2.sh .```
Get test script 3 for BRAKER via Singularity in /home/usr
```singularity exec -B $PWD:$PWD braker3.sif cp /opt/BRAKER/example/singularity-tests/test3.sh .```
Export two bash environments to run test scripts
```export BRAKER_SIF=/home/pablo/braker3.sif```
Execute test script 1
```bash test1.sh```
Execute test script 2
```bash test2.sh```
Execute test script 3
```bash test3.sh```
### STRUCTURAL ANNOTATION i.e. PROTEIN ID
Was not able to solve errors trying to run braker from outside home directory. Sent a message to the author on github, waiting for fix.
```singularity exec braker3.sif braker.pl --species=Scurria_viridula05 --genome=VG2_canu_purged_purged_ragtag.scaffold.fasta.masked.nosemicolon --prot_seq=Ensamble_s_viridula_31_y_143.fasta.transdecoder_CDHIT_95_nospcharac2.faa --softmasking --workingdir=braker_14apr2023 --threads=15```
BRAKER finishes correctly if braker.aa, braker.gtf, augustus.hints.aa, and augustus.gtf are created.
###### Count number of proteins recovered
```grep -c ">" augustus.hints.aa```
###### Run busco on augustus.hints.aa
```conda activate busco```
```nohup busco -o VG2_busco_metazoa -i /media/pablo/disco1/pablo/scurria_genome_data/emily/viridula/braker_14apr2023/Augustus/augustus.hints.aa -l metazoa_odb10 -m proteins &```
Busco results in this case were C:96.8%, S:90.4%, D:6.4%
Can also run busco metazoa to verify. 
### FUNCTIONAL ANNOTATION
In brief, this step involves identifying functions and GO terms for the list of proteins recovered from BRAKER.If you do not already have a conda env for blast with blast installed, do so at this time.
##### SwissProt database recovery
I placed the databases in ~/miniconda3/envs/blast 
```wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz```
###### Unzip .fasta.gz
```gunzip uniprot_sprot.fasta.gz```
###### Index the database
```makeblastdb -in uniprot_sprot.fasta -out uniprot_sprot -dbtype prot```
##### Refseq database recovery
I was not entirely sure which Refseq database to use. There are many... for invertebrates, all, bacteria, etc.  I chose the one labeled as "complete"
```wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/complete/complete.1000.protein.faa.gz ```
###### unzip the database
```gunzip complete.1000.protein.faa.gz```
###### Index the database
```makeblastdb -in complete.1000.protein.faa -out refseq_complete_prot -dbtype prot```
##### Blast proteins against db
```nohup blastp -outfmt 6 -query /media/pablo/disco1/pablo/scurria_genome_data/emily/viridula/braker_14apr2023/Augustus/augustus.hints.aa -db uniprot_sprot -out /media/pablo/disco1/pablo/scurria_genome_data/emily/viridula/viridula_braker_14apr2023_augustushints.txt -num_threads 10 &```
Output format selected based on that suggested by [AHRD](https://github.com/groupschoof/AHRD#21-ahrd-example-usages)

 
