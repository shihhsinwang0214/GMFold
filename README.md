# GMFOLD: Graph matching for DNA-aptamer secondary structure classification and machine learning interpretabiity

This repository contains the code and experiments to replicate the results presented in the paper titled "GMFOLD: Graph matching for DNA-aptamer secondary structure classification and machine learning interpretability"

**Paolo Climaco**^[1], **Noelle Mitchell**^[2], **Matthew Tyler**^[3], **Kyungae Yang**[4],  **Anne Andrews**^[2,5,6]. and **Andrea L. Bertozzi**^[3,5]

[1] Institut für Numerische Simulation, Universität Bonn, Germany
[2] Department of Chemistry and BioChemistry, University of California, Los Angeles, 90095, CA, USA\
[3] Department of Mathematics, University of California, Los Angeles, 90095, CA, USA\
[4] Department of Medicine, Columbia University Irving Medical Center, New York, NY, 10032, USA\
[5] California NanoSystems Institute, University of California, Los Angeles, 90095, CA, USA\
[6] Departments of Psychiatry and Biobehavioral Sciences and Bioengineering, Semel Institute for Neuroscience and Human Behavior, and Hatos Center for Neuropharmacology , Los Angeles, 90095, CA, USA


Contact climaco@ins.uni-bonn for questions about code and data.

## Directory Structure
```plaintext

├── data/                           # Folder containing experimental data analyzed in our experiments
|    ├── dna.csv                        # CSV file with DNA sequences data from Seqfold github repository
|    ├── fold_published.csv             # CSV file with sequences and associated secondary structures and energies computed with GMfold
|    ├── Raw_files/                     # Subfolder containing raw, unprocessed experimental data
|    ├── published_clean_files/         # Subfolder with cleaned and processed data ready for analysis
|
├── notebooks/                      # Folder containing Jupyter notebooks for large scale data analysis 
|    ├── comparison.ipynb               # Notebook comparing different DNA folding algorithms: MGfold, Seqfold 2_0 and Seqfold
|    ├── folding_constraint.txt         # Text file outlining constraints required to get mfold/UNAfold energies
|    ├── fold_data.ipynb                # Notebook for processing and folding DNA sequences at scale
|    ├── read_excels_files.ipynb        # Notebook for reading and processing Excel files with experimental data
|    ├── Seq_vs_Seq2.ipynb              # Notebook comparing speed of Seqfold and  Seqfold 2.0 for DNA folding
|    ├── similarity_search.ipynb        # Notebook for similarity search of DNA secondary structures
|    ├── topic_modeling.ipynb           # Notebook applying topic modeling techniques 
|
├── src/                           # Folder containing source code for DNA folding
|    ├── dna.py                         # Script containing DNA enthalpy and entropy parameters foe energy computations
|    ├── gm_energy_functions.py         # Script containing energy functions used in MGfold
|    ├── GMfold.py                      # Main script for DNA folding using graph matching methods: Implementation of MGfold
|    ├── graph_matching.py              # Script for performing graph matching on DNA structures
|    ├── seqfold2_0.py                  # Script implementing the SeqFold 2.0 algorithm for DNA sequence folding
|    ├── seqfold_accessed_07_24/        # Folder with version of SeqFold (version 0.7.17) on which we based our code (Accesed on July 2024)
|           ├── fold_accessed_07_24     # Script implementing Seqfold (version 0.7.17)
            ├── LICENCE_seqfold         # Text file including seqfold licence agreement
|    ├── Types.py                       # Script defining custom data types for DNA folding operations
|    ├── utils.py                       # Utility script with helper functions used across DNA folding workflows

└── environment.yml            # python environment with packages required to run code.
└── README.md                   # Project README file.

```