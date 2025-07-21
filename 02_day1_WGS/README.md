# 2025_KOGO_WGS

### 0. Type this code to enter your instance

```
ssh -p 7024 <계정번호>@[public_ip]
# ex) ssh -p 7024 edu04@59.26.46.104
```

### 1. BWA (Burrows-Wheeler Aligner) MEM (matching extension)

```
bwa mem -t 8 -M -Y -K 100000000 \
 -R '@RG\tID:SRR11880780\tPL:ILLUMINA\tPU:SRR11880780\tSM:SRR11880780\tLB:SRR11880780' \
~/2025_KOGO_workshop/wgs/data/hg38/hg38.fa \
<(zcat ~/2025_KOGO_workshop/wgs/data/SRR11880780/SRR11880780_1.fastq.gz) \
<(zcat ~/2025_KOGO_workshop/wgs/data/SRR11880780/SRR11880780_2.fastq.gz) \
| samtools view -huS - \
| samtools sort -@ 2 -m 2G -o SRR11880780.bam -O bam -T SRR11880780.tmp
```

### 2. Markduplicates

```
java -jar /BiO/prog/picard/bin/picard.jar MarkDuplicates \
I=SRR622461_sub.bam \
O=dedup.SRR622461_sub.bam \
M=markdups_SRR622461_sub.txt \
ASSUME_SORT_ORDER=queryname \
MAX_RECORDS_IN_RAM=2000000 \
COMPRESSION_LEVEL=1 \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=SILENT \
PROGRAM_RECORD_ID=null \
ADD_PG_TAG_TO_READS=false \
READ_NAME_REGEX=null
```
```
samtools index dedup.SRR622461_sub.bam
```

### 3. Base Quality Score Recalibration (BQSR)

```
gatk BaseRecalibrator \
-R /BiO/home/edu{i}/2025_KOGO_workshop/wgs/data/hg38/hg38.fa \
-I dedup.SRR622461_sub.bam \
--known-sites /BiO/home/edu{i}/2025_KOGO_workshop/wgs/data/hg38/dbsnp_146.hg38.vcf.gz \
--known-sites /BiO/home/edu{i}/2025_KOGO_workshop/wgs/data/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
--known-sites /BiO/home/edu{i}/2025_KOGO_workshop/wgs/data/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
-O SRR622461_sub.recal.table

# edu 번호를 본인에 맞게 바꾸어주셔야합니다
```

```
gatk ApplyBQSR \
-R /BiO/home/edu{i}/2025_KOGO_workshop/wgs/data/hg38/hg38.fa \
-I dedup.SRR622461_sub.bam \
-O SRR622461_sub.cram \
-bqsr ./SRR622461_sub.recal.table
```

```
samtools index SRR622461_sub.cram
```
