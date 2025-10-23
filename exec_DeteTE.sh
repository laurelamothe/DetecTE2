#!/bin/bash
sbatch -o . -J "" DetecTE_master.sh --fasta test.fasta --database db_test.fasta --assembler TRINITY --outdir . --dataset_name ""