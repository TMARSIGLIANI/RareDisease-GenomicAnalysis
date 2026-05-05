#!/bin/bash

# Variant calling pipeline for rare Mendelian disease diagnosis in family trios.
# Covers alignment, QC, variant calling, annotation and coverage tracks on chr20.
# Usage: bash raredisease.sh
# Dependencies: bowtie2, samtools, fastqc, qualimap, freebayes, bcftools, vep, bedtools, multiqc

for trio in trio_1 trio_2 trio_3 trio_4 trio_5; do
  echo "=== WORKING ON $trio ==="
  cd $trio

  # Paired-end alignment to chr20. Read groups added to tag each member in the VCF.
  # Piped directly into samtools to avoid intermediate SAM files.
  echo "--- Aligning reads ---"
  bowtie2 -1 HG00445.targets_R1.fq.gz -2 HG00445.targets_R2.fq.gz -x /home/BCG2026_exam/chr20 --rg-id "child" --rg "SM:child" | samtools view -Sb | samtools sort -o ${trio}_child.bam
  bowtie2 -1 HG00446.targets_R1.fq.gz -2 HG00446.targets_R2.fq.gz -x /home/BCG2026_exam/chr20 --rg-id "father" --rg "SM:father" | samtools view -Sb | samtools sort -o ${trio}_father.bam
  bowtie2 -1 HG00447.targets_R1.fq.gz -2 HG00447.targets_R2.fq.gz -x /home/BCG2026_exam/chr20 --rg-id "mother" --rg "SM:mother" | samtools view -Sb | samtools sort -o ${trio}_mother.bam

  echo "--- Indexing BAMs ---"
  samtools index ${trio}_child.bam
  samtools index ${trio}_father.bam
  samtools index ${trio}_mother.bam

  # FastQC for per-file QC, Qualimap for exome-specific coverage metrics.
  # Both aggregated by MultiQC at the end.
  echo "--- Running FastQC ---"
  fastqc *.bam

  echo "--- Running Qualimap ---"
  qualimap bamqc -bam ${trio}_child.bam --feature-file ../chr20_ILMN_Exome_2.0_Plus_Panel.hg38_padded.bed -outdir ${trio}_child
  qualimap bamqc -bam ${trio}_father.bam --feature-file ../chr20_ILMN_Exome_2.0_Plus_Panel.hg38_padded.bed -outdir ${trio}_father
  qualimap bamqc -bam ${trio}_mother.bam --feature-file ../chr20_ILMN_Exome_2.0_Plus_Panel.hg38_padded.bed -outdir ${trio}_mother

  # Joint haplotype-based variant calling on the trio.
  # Filters: min mapping quality 20, min 5 alt reads, min base quality 10, min coverage 10x.
  echo "--- Variant calling ---"
  freebayes -f /home/BCG2026_exam/chr20.fa -m 20 -C 5 -Q 10 --min-coverage 10 ${trio}_child.bam ${trio}_father.bam ${trio}_mother.bam > ${trio}.vcf

  bgzip ${trio}.vcf
  bcftools index ${trio}.vcf.gz

  # Genotype filter based on expected inheritance pattern (RR=hom ref, RA=het, AA=hom alt).
  # trio_1/4/5: autosomal recessive | trio_2: maternal dominant | trio_3: de novo
  if [ $trio == "trio_1" ] || [ $trio == "trio_4" ] || [ $trio == "trio_5" ]; then
    FILTER='GT[0]="AA" && GT[1]="RA" && GT[2]="RA"'
  elif [ $trio == "trio_2" ]; then
    FILTER='GT[0]="RA" && GT[1]="RR" && GT[2]="RA"'
  elif [ $trio == "trio_3" ]; then
    FILTER='GT[0]="RA" && GT[1]="RR" && GT[2]="RR"'
  fi

  # Filter by target regions, samples, genotype pattern and call quality.
  echo "--- Filtering variants ---"
  bcftools view -R ../chr20_ILMN_Exome_2.0_Plus_Panel.hg38_padded.bed ${trio}.vcf.gz | bcftools view -S ../samples.txt | bcftools view -i "${FILTER}" | bcftools filter -i 'QUAL>20' -Ov -o ${trio}.cand.vcf

  # Annotate with functional impact, allele frequencies (1000G, gnomAD) and
  # pathogenicity predictions (SIFT, PolyPhen).
  echo "--- Annotating with VEP ---"
  vep -i ${trio}.cand.vcf -o ${trio}.vep_annotated.vcf --vcf --cache --offline --assembly GRCh38 --dir_cache /data/vep_cache --use_given_ref --mane --pick_allele --af --af_1kg --af_gnomade --max_af --sift b --polyphen b --no_fasta

  # Keep HIGH/MODERATE impact variants with MAX_AF < 0.01% or absent from population databases.
  echo "--- Filtering VEP output ---"
  filter_vep -i ${trio}.vep_annotated.vcf -o ${trio}.vep_filtered.vcf --filter "(IMPACT is HIGH or IMPACT is MODERATE) and (not MAX_AF or MAX_AF < 0.0001)"

  # BedGraph coverage tracks for IGV visualization, capped at 100x.
  echo "--- Computing coverage tracks ---"
  bedtools genomecov -ibam ${trio}_father.bam -bg -trackline -trackopts 'name="father"' -max 100 > ${trio}_fatherCov.bg
  bedtools genomecov -ibam ${trio}_mother.bam -bg -trackline -trackopts 'name="mother"' -max 100 > ${trio}_motherCov.bg
  bedtools genomecov -ibam ${trio}_child.bam -bg -trackline -trackopts 'name="child"' -max 100 > ${trio}_childCov.bg

  cd ..
done

# Aggregate all QC reports into a single interactive HTML report.
echo "--- Running MultiQC ---"
multiqc trio_*/

echo "==== END ===="