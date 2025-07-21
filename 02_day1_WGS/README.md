# 2025_KOGO_WGS

### 0. Type this code to enter your instance

```
ssh -p 7024 <계정번호>@[public_ip]
# ex) ssh -p 7024 edu04@59.26.46.104
```

### 1. BWA (Burrows-Wheeler Aligner) MEM (matching extension)

```
bwa mem -t 8 -M -Y \
 -R '@RG\tID:SRR11880780\tPL:ILLUMINA\tPU:SRR11880780\tSM:SRR11880780\tLB:SRR11880780' \
~/2025_KOGO_workshop/wgs/data/hg38/hg38.fa \
<(zcat ~/2025_KOGO_workshop/wgs/data/SRR11880780/SRR11880780_1.fastq.gz) \
<(zcat ~/2025_KOGO_workshop/wgs/data/SRR11880780/SRR11880780_2.fastq.gz) \
| samtools view -huS - \
| samtools sort -@ 2 -m 2G -o ~/2025_KOGO_workshop/wgs/results/SRR11880780.bam -O bam -T ~/2025_KOGO_workshop/wgs/results/SRR11880780.tmp
```

### 2. Markduplicates

```
java -jar /BiO/prog/picard/bin/picard.jar MarkDuplicates \
I=~/2025_KOGO_workshop/wgs/results/SRR11880780.bam \
O=~/2025_KOGO_workshop/wgs/results/dedup.SRR11880780.bam \
M=~/2025_KOGO_workshop/wgs/results/markdups_SRR11880780.txt \
ASSUME_SORT_ORDER=coordinate \
MAX_RECORDS_IN_RAM=2000000 \
COMPRESSION_LEVEL=1 \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=SILENT
```

```
samtools index dedup.SRR11880780.bam # required when sort order is not set to coordinate
```

### 3. Base Quality Score Recalibration (BQSR)

```
gatk BaseRecalibrator \
-R ~/2025_KOGO_workshop/wgs/data/hg38/hg38.fa \
-I ~/2025_KOGO_workshop/wgs/results/dedup.SRR11880780.bam \
--known-sites ~/2025_KOGO_workshop/wgs/data/hg38/dbsnp_146.hg38.vcf.gz \
--known-sites ~/2025_KOGO_workshop/wgs/data/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
--known-sites ~/2025_KOGO_workshop/wgs/data/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
-O ~/2025_KOGO_workshop/wgs/results/SRR11880780.recal.table
```

```
gatk ApplyBQSR \
-R ~/2025_KOGO_workshop/wgs/data/hg38/hg38.fa \
-I ~/2025_KOGO_workshop/wgs/results/dedup.SRR11880780.bam \
-O ~/2025_KOGO_workshop/wgs/results/SRR11880780.cram \
-bqsr ~/2025_KOGO_workshop/wgs/SRR11880780.recal.table
```

```
samtools index ~/2025_KOGO_workshop/wgs/results/SRR11880780.cram
```
