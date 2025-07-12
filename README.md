# WES Variant Calling Pipeline

This repository contains a complete pipeline for variant calling from **Whole Exome Sequencing (WES)** data using tools like `FastQC`, `fastp`, `BWA`, `Picard` and `GATK`.
This is part of a project that involved variant calling from WES of families with primary open angle glaucoma.
Variant calling was done parallely on split bam files based on chromosomes and then merged together for speed and memory efficiency. If your aim does not allow that please stick to regular HaplotypeCaller or HaplotypeCallerSpark.
Change the arguments specified to allocate resources as per your needs. 

---

## 📁 Project Structure

```
.
├── data/                   # Raw sequencing data (not tracked)
├── output/                 # Output directory (not tracked)
├── tools/                 # Tools and binaries
├── Ref/                   # Reference files (e.g., hg38, dbs)
├── variant_calling.sh     # Main pipeline script
├── data.txt               # Sample names (e.g., SRR32421545)
├── .gitignore
└── README.md
```



## ⚙️ Tools Used

- `FastQC` – Quality check
- `fastp` – Trimming and adapter removal
- `BWA-MEM` – Read alignment
- `Picard` – Mark duplicates and read group assignment
- `GATK` – Variant calling and BQSR

---

## 🔁 Pipeline Steps

1. **Input:** Paired-end FASTQ files (NovaSeq, Hybrid Selection) or binary SRA file
2. **Quality Control:** `FastQC`, `fastp`
3. **Alignment:** `bwa mem` → `SAM` → `BAM` → sorted BAM
4. **Mark Duplicates:** `Picard`
5. **Base Quality Recalibration (BQSR):** `GATK BaseRecalibrator` and `ApplyBQSR`
6. **Variant Calling:** `GATK HaplotypeCaller` (parallel per chromosome)

---

## 🚀 How to Run

```bash
bash variant_calling.sh
Make sure to update data.txt with your SRA run accessions or sample names.

📦 Requirements
Linux with bash
Java ≥ 8
Python 3

Tools installed in tools/ directory and added to $PATH

📝 Notes
This pipeline is optimized for human WES data.

Large files like .fastq, .bam, and .vcf.gz are excluded via .gitignore.

For known sites in BQSR, download:

1000G_phase1.snps.high_confidence.hg38.vcf.gz

Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

📧 Contact
YadhuKrishnan
Email: yadhuryk@gmail.com
GitHub: @krishnanyadhu
