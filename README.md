# RareDisease-GenomicAnalysis
A genomics workflow to detect rare Mendelian diseases in simulated family trio data.

## Why this project
This is part of the final exam for the Genomics and Transcriptomics course — a core module of the Master's in Bioinformatics for Computational Genomics at UniMi and PoliMi. The idea behind the course is pretty straightforward: learn how real genomic data gets handled in practice, from raw sequencing reads all the way to something biologically meaningful, working mostly in Unix. This repository was created for purely educational purposes and is intended to demonstrate the functionality of the pipeline rather than provide a production-ready tool. As such, applying it to different datasets would require adapting the pipeline variables, file paths, and directory structure to match the new input. It is not designed to work out of the box with arbitrary data.

## What I did
I worked on simulated exome sequencing data from 5 family trios — two parents and a child — and tried to reach a diagnosis for each of them. The core challenge is leveraging the trio structure to filter out noise and zero in on the rare variant actually responsible for the condition. The analysis is scoped to chromosome 20, GRCh38 reference, to keep things feasible within the project. The workflow goes through the usual steps of a variant calling pipeline: read alignment, sorting and indexing with SAMtools, target region handling with BEDTools to make sure we're only looking at the exome, and variant calling with FreeBayes, that works particularly well in this kind of family-based setting. Candidate variants were then annotated and prioritized using Ensembl VEP, focusing on high-impact mutations with low allele frequency. The pipeline also integrates a comprehensive quality control assessment via FastQC, Qualimap, and MultiQC to ensure data reliability. 

## Student
Tiziano Marsigliani
