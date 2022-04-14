# Visualizing the Human Genome
### UC Berkeley - W209 Data Visualization - Final Project
### Authors:
- Eric Ellestad
- Whitney Fee
- Angel Ortiz Nuñez


## Division of Responsibilities:
- Whitney Fee:
  - Tableau Visualizations
  - CSS and Website Formatting
  - Definitions

- Angel Ortiz Nuñez
  - Altair Charts - Gene Chromosome
  - Altair Charts - Protein Coding Composition by Chromosome & Genome
  - Altair Interactivity Across Charts
  - Altair Gene Function Charts

- Eric Ellestad:
  - Protein viewer
  - Bootstrap Formatting
  - Gene Expression Charts
  - Flask Integration



## Introduction

[The Human Genome Project](https://www.genome.gov/human-genome-project) was completed in 2003 which ushered in the genomics era by sequencing an entire human genome for the very first time. In the following two decades, advances in modern high-throughput genetic sequencing technologies have made DNA sequencing faster, cheaper, and widely available. This has lead to a proliferation of genomic "big data" and the availability of large-scale public databases of detailed human genomic data.

The human genome consists of 23 chromosome pairs that include around 20,000 genes and over 6 billion base pairs in total DNA length. Given that the scale and complexity of genomic data can make these topics difficult to grasp intuitively, our project is an attempt to create an interactive visualization tool of the human genome that stitches together the macromolecular structure with the building blocks all within the context of the central dogma of molecular biology: DNA gets transcribed into RNA which gets translated into Proteins that drive many bodily functions.


Folder Structure:
```
.
├── Procfile
├── README.md
├── app.py
├── data
│   ├── chromosomes
│   ├── genes
│   └── genome_genes
├── requirements.txt
├── static
│   ├── charts
│   ├── css
│   └── js
└── templates
└── venv
    ├── bin
    ├── lib
    └── pyvenv.cfg
```

Launch Virtual Machine:
`source venv/bin/activate`

Code to Launch Flask App from the Command Line:
`python app.py`
