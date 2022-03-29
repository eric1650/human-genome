# Visualizing the Human Genome
### By Eric Ellestad, Whitney Fee, and Angel Ortiz Nuñez
### UC Berkeley - W209 Data Visualization - Final Project

The human genome consists of 23 chromosome pairs that include 20,000 protein-coding genes and over 6 billion base pairs in total DNA length. The Human Genome Project kicked off the genomics era and modern high-throughput genetic sequencing technologies have resulted in large scale publicly available databases of detailed human genetic information.

Our project is attempting to create an intuitive and interactive visualization tool of the human genome that shows the chromosomal location and subcomponent organization of each gene. Protein-coding genes will include additional visualizations regarding their RNA transcripts and information regarding the resulting post-translation protein such as its functional role in the body and a 3D visualization of the folded structure of the protein. The fully-zoomed out overview will reveal the organization and scale of the human genome while the zoomed-in details will reveal information about specific genes and their biological function.


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

Requirements:
- Flask
- Altair

Code to Launch Flask App from the Command Line:
`python app.py`
