#!/usr/bin/env bash -l

#Alignment
 
for d in */
do
a="${d::-1}";

(cd "$d" && 

#1. Alignment to reference genome
echo $a Aalignment started!!
bwa mem -t 6 /mnt/9af4a3cb-07fd-489a-85f8-9e2253be3d37/Freyja/wastewater_Analysis_batch_2_03_11_2022/barcode84/index/nCoV-2019.reference.fasta fastq/FMR*.fastq.gz > "$a".sam 
echo $a alignment completed!


#2. Covert SAM into BAM
samtools view -@6 -b $a.sam > "$a".bam
echo $a bam file created!

#3.Sort the BAM file
samtools sort -@ 6 $a.bam -o "$a"_sorted.bam


echo $a sorting alignment completed!
eval "$(conda shell.bash hook)"
conda activate artic


## View Coverage 
samtools coverage -m --reference /home/fmr-nanopore/Downloads/fieldbioinformatics-1.2.1/primer_schemes/nCoV-2019/V1200/nCoV-2019.reference.fasta "$a"_sorted.bam > genome_coverage

echo $a Coverage completed!

#### iVar (always Check the primer version)
#Trimming Primers
ivar trim -e -b /home/fmr-nanopore/C-WAP-main/covidRefSequences/ARTICv4.bed -p trimmed_"$a" -i "$a"_sorted.bam -q 15 -m 20


echo $a Primer trimming completed!
#Sort the bam file
samtools sort -@ 6 trimmed_"$a".bam -o "$a"_resorted.bam

##### Variant calling
#mpileup
samtools mpileup -aa -A -d 10000 -B -Q 0 --reference /home/fmr-nanopore/Downloads/fieldbioinformatics-1.2.1/primer_schemes/nCoV-2019/V1200/nCoV-2019.reference.fasta -o "$a"_pile.up "$a"_resorted.bam

#Variant calling
cat "$a"_pile.up | ivar variants -p rawVarCalls -g /home/fmr-nanopore/C-WAP-main/covidRefSequences/covidGenomeAnnotation-NCBI.gff -r /home/fmr-nanopore/Downloads/fieldbioinformatics-1.2.1/primer_schemes/nCoV-2019/V1200/nCoV-2019.reference.fasta -m 10


bcftools mpileup -d 10000 -Ou -f /home/fmr-nanopore/Downloads/fieldbioinformatics-1.2.1/primer_schemes/nCoV-2019/V1200/nCoV-2019.reference.fasta "$a"_resorted.bam | bcftools call --ploidy 1 -mv -Oz -o "$a"_calls.vcf.gz

bcftools index "$a"_calls.vcf.gz 

echo $a vcf generation finished!
###Consensus Generation
cat /home/fmr-nanopore/Downloads/fieldbioinformatics-1.2.1/primer_schemes/nCoV-2019/V1200/nCoV-2019.reference.fasta | bcftools consensus "$a"_calls.vcf.gz > "$a"_consensus.fa
sed -i 's/MN908947.3/"$a"/g' "$a"_consensus.fa

echo $a consensus generated!


#Freyja
eval "$(conda shell.bash hook)"
conda activate freyja-env

freyja variants "$a"_resorted.bam --variants "$a"_freyja.variants.tsv --depths "$a"_freyja.depths.tsv --ref /home/fmr-nanopore/Downloads/fieldbioinformatics-1.2.1/primer_schemes/nCoV-2019/V1200/nCoV-2019.reference.fasta

freyja demix "$a"_freyja.variants.tsv "$a"_freyja.depths.tsv --output "$a"_freyja.demix --confirmedonly
mkdir demix
mv "$a"_freyja.demix demix/

freyja boot "$a"_freyja.variants.tsv "$a"_freyja.depths.tsv --nt 6 --nb 1000 --output_base "$a"_freyja_boot



## Plot generation
##Create demix folder with .demix file

freyja aggregate demix/ --output "$a"-aggregated.tsv

freyja plot --mincov 50 "$a"-aggregated.tsv --output "$a".pdf

freyja plot --lineages --mincov 50 "$a"-aggregated.tsv
echo $a freyja completed!

##################Annotation of variants
#Unzipping the file
gunzip "$a"_calls.vcf.gz 

#convert MN908947.3 to NC_045512.2 in VCF file
sed -i 's/MN908947.3/NC_045512.2/g' "$a"_calls.vcf

##Snpeff

java -Xmx8g -jar /home/fmr-nanopore/snpEff/snpEff.jar NC_045512.2 "$a"_calls.vcf > "$a"_ann_vcf
echo $a freyja completed!

echo $a snpeff completed
#Nanoplot
NanoPlot --prefix Nano --fastq fastq/*.fastq.gz -o "$a"_NanoPlot_result -t 8
echo NanoPlot completed for $a!!

);
done





