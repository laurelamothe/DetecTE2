# DetecTE2

**A pipeline for the detection and annotation of transposable elements in transcriptomes**

DetecTE is a pipeline enabling the extraction of active transposable element (TE) sequences in transcriptomes without requiring a reference genome. It also clusters TE sequences into families to estimates their diversity. 

The pipeline integrates multiple stages, including:
- TE detection using reciprocal BLAST (`tblastn` and `blastx`) between transcripts and a protein TE database
- Superfamily assignment (e.g. Copia, Gypsy, Bel)
- Clustering into transposable element families

## ✨ Features

- **TE detection** using `blastx` and `tblastn`
- **Superfamily assignment** based on the majority annotation of BLAST hits
- **Clustering** of elements by families by a `blastn` of sequences against themselves
- **Parallel execution** across multiple transcriptomes

## 🧰 Requirements

- Python version 3.12
- BLAST+ suite version 2.14.0
- Execution with a **SLURM** job scheduler

## 📦 Installation

Clone the repository:

```bash
git clone https://github.com/laurelamothe/DetecTE2
cd DetecTE2
```

## 🚀 Usage
### 1) Format the protein Transposable elements database
The sequences title should be formatted as follows:<br>
`>sequence_name#CLASS/GROUP/SUPERFAMILLY sequence description`<br>
In particular, the pipeline will use the CLASS, GROUP and SUPERFAMILLY information to perform the annotation and the clustering.
Example database:
```pgsql
>DIRS_element#LTR/DIRS/DIRS  this is a DIRS sequence
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
>ISL2EU_element#DNA/PIF/ISL2EU this is a ISL2EU sequence
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
>PiggyBac_element#DNA/PIGGYBAC/PiggyBac this is a PiggyBac sequence
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
...
```
We recommand using a curated database enriched with sequences from the studied species.

### 2) Update the list of superfamilies
The file `superfamilly_list.txt` should be edited to contain every annotation included in the TE database. 
The annotation should follow the format of the database sequence titles:
```pgsql
LTR/DIRS/DIRS
DNA/PIF/ISL2EU
DNA/PIGGYBAC/PiggyBac
...
```

### 3) Run the pipeline
The pipeline architecture depends on the SLURM job scheduler. Ensure it is installed before running DetecTE2.
```pgsql
sbatch DetecTE_main.sh --fasta  /path/to/transcriptome/transcripts.fasta \               
                       --database /path/to/database/db.fasta \ 
                       --assembler ASSEMBLER \
                       --outdir /path/to/output_directory/ \
                       --prefix prefix 
...
```
*Notes*
- To detect TE sequences in several transcriptomes, paste the path to each file. The sequences will be pooled into families independently of the transcriptome in which they were detected. 

- The assembler names should be either SPADES or TRINITY. For multiple transcriptome inputs, enter the assembler name for each transcriptome, in the same order. 

- The prefix is specific to the transcriptome input. For multiple transcriptome inputs, indicate the corresponding number of prefixes.

You can run the pipeline with the example dataset: 
```pgsql
sbatch -o ./report.out DetecTE_main.sh --fasta ./test.fasta --database ./db_test.fasta --assembler TRINITY --outdir ./results --prefix test
...
```
