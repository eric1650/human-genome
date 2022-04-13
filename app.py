##############################
# Flask Implementation
# Visualizing the Human Genome
# Updated 03/25/2022
##############################

from flask import Flask, render_template, jsonify
import altair as alt
import pandas as pd
import json
import zipfile


# Create Flask App Instance
app = Flask(__name__)

##############################
# Genome Data
##############################

genome = pd.read_csv('data/genome_genes/genome.csv.zip')
gene_names = list(genome.gene_name.unique())
gene_names.sort()

chr_composition = pd.read_csv('data/gene_composition/chr_composition.csv')
genome_composition = pd.read_csv('data/gene_composition/genome_composition.csv')

##############################
# Flask routes
##############################

# Render index.html homepage
@app.route('/')
def index(gene_names=gene_names):
    return render_template('index.html', gene_names=gene_names)


# Render molstar.html viewer of selected gene
@app.route('/molstar/<gene_name>')
def molstar(gene_name, genome=genome):
    gene = genome[genome.gene_name == gene_name]
    gene = gene[gene.feature == 'gene'].reset_index()
    uniprot_id = gene.loc[0, 'uniprot_id']
    return render_template('molstar.html', uniprot_id=uniprot_id)


# Retrieve Protein Info for a Specific Gene
@app.route('/protein_info/<gene_name>')
def protein_info(gene_name, genome=genome):
    gene = genome[genome.gene_name == gene_name]
    gene = gene[gene.feature == 'gene'].reset_index()

    chromosome = gene.loc[0, 'chromosome']
    gene_start = str(gene.loc[0, 'start'])
    gene_end = str(gene.loc[0, 'end'])
    protein_name = gene.loc[0, 'protein_name']
    protein_function = gene.loc[0, 'protein_function']
    aa_seq = gene.loc[0, 'aa_sequence']
    aa_length = len(aa_seq)

    return jsonify(
    chromosome=chromosome,
    gene_start=gene_start,
    gene_end=gene_end,
    name=protein_name,
    function=protein_function,
    aa_seq=aa_seq,
    aa_length=aa_length
    )


# Render definitions.html for part 1
@app.route('/definitions')
def definitions():
    return render_template('definitions.html')


##############################
# Altair routes
##############################

# Render altair pie chart of gene composition by genome
@app.route('/chart/gene_composition_genome/')
def gene_composition_genome(genome_composition=genome_composition):

    pie_chart = alt.Chart(genome_composition, title= 'Percentage of DNA in Genome that are Genes').mark_arc().encode(
        theta= alt.Theta('Percent of Genome:Q'),
        color=alt.Color('Protein Coding:N', title= 'Gene Label'),
        tooltip=['Protein Coding', 'Percent of Genome']
    )

    return pie_chart.to_json()

# Render altair bar chart of gene composition by chromosome
@app.route('/chart/gene_composition_chr/')
def gene_composition_chr(chr_composition=chr_composition):

    protein_coding_bar_bp = alt.Chart(chr_composition, title='Amount of DNA within Each Chromosome that are in Genes').mark_bar().encode(
        x = alt.X('Chromosome:N', sort=None),
        y = alt.Y('Base Pair Length:Q', title="Total Amount of DNA Base Pairs"),
        color = 'Protein Coding',
        tooltip = alt.Tooltip(['Chromosome','Base Pair Length','Protein Coding'])
    ).properties(
        width = 750
    )

    # protein_coding_bar_pct = alt.Chart(chr_composition, title='Percent of DNA in Each Chromosome that is in a Gene').mark_bar().encode(
    #     x = alt.X('Chromosome:N', sort=None),
    #     y = alt.Y('Percent of Chromosome:Q', scale=alt.Scale(domain=[0,100])),
    #     color = 'Protein Coding',
    #     tooltip = alt.Tooltip(['Chromosome','Percent of Chromosome','Protein Coding'])
    # ).properties(
    #     width = 750
    # )

    # chart = alt.vconcat(protein_coding_bar_bp & protein_coding_bar_pct)

    return protein_coding_bar_bp.to_json()


# Render altair chart of gene components
@app.route('/chart/gene_location/<gene_name>')
def gene_location(gene_name, genome=genome, chr_composition=chr_composition):

    gene = genome[genome.gene_name == gene_name].reset_index()
    chromosome = gene.loc[0,'chromosome']
    chr_length = chr_composition[chr_composition.Chromosome == chromosome]['Base Pair Length'].sum()
    domain = [0, chr_length]

    chart = alt.Chart(gene).mark_point(size=150, filled=True).encode(
        x = alt.X('start:Q', title=f"DNA Base Pair Position on {chromosome}", scale=alt.Scale(domain=domain)),
        y = alt.Y('chromosome', title=None),
        color = alt.Color('gene_name', title=None),
        tooltip = alt.Tooltip(['gene_name', 'chromosome' ,'start', 'end'])
    ).properties(width=750, title=f"Location of {gene_name} on {chromosome}")

    return chart.to_json()


# Render altair chart of gene components
@app.route('/chart/gene_components/<gene_name>')
def gene_components(gene_name):

    # gene = genome[genome.gene_name == gene_name]
    path = 'data/genes/' + gene_name + '.csv.zip'
    gene = pd.read_csv(path)
    chromosome = list(gene[gene.feature == 'gene']['chromosome'])[0]

    domain_min = gene.start.min()
    domain_max = gene.end.max()
    domain = [domain_min, domain_max]

    bar = alt.Chart(gene[gene.feature.isin(['gene','exon','UTR','CDS'])]).mark_bar().encode(
        x=alt.X('start:Q', scale=alt.Scale(domain=domain), title=f"DNA Base Pair Position on {chromosome}"),
        x2='end:Q',
        y=alt.Y('feature:N', sort=None, title='Gene Component'),
        color = alt.Color('feature:N', title='Gene Component'),
        tooltip = alt.Tooltip(['gene_name','feature','start','end'])
    ).properties(
        width = 750,
        title=f'Components of Gene {gene_name}'
    ).interactive()

    hist_features = alt.Chart(gene[gene.feature.isin(['exon','UTR','CDS'])]).mark_bar().encode(
        x = alt.X('count(feature)', title=f'Count of Each Component Type in Gene {gene_name}'),
        y = alt.Y('feature:N', sort='-x', title='Gene Component'),
        color = alt.Color('feature:N', title='Gene Component'),
        tooltip = alt.Tooltip(['gene_name','feature','count(feature)'])
    ).properties(width=750)

    gene_layout = bar & hist_features

    return gene_layout.to_json()


# render altair chart of gene transcription
@app.route('/chart/gene_expression/<gene_name>')
def gene_expression(gene_name):

    # import gene
    path = 'data/genes/' + gene_name + '.csv.zip'
    gene = pd.read_csv(path)
    chromosome = list(gene[gene.feature == 'gene']['chromosome'])[0]

    # get exons only
    transcripts = gene[gene.feature != 'gene']
    exons = transcripts[transcripts.feature == 'exon']

    # get spliced components
    spliced = transcripts[transcripts.feature.isin(['UTR','CDS'])]

    # get translated components
    translated = spliced[spliced.feature.isin(['CDS'])]

    # set chart domain
    domain_min = gene.start.min()
    domain_max = gene.end.max()
    domain = [domain_min, domain_max]

    # build transcription chart
    transcription = alt.Chart(exons).mark_bar().encode(
        x = alt.X('start', scale=alt.Scale(domain=domain), title=f"DNA Base Pair Position on {chromosome}"),
        x2 = 'end',
        y = alt.Y('chromosome:N', title=None),
        color = 'feature',
        tooltip = alt.Tooltip(['gene_name','feature','start','end'])
    ).properties(
        width = 750,
        title = f"1. Transcription: Exons are the portions of a gene that get transcribed into RNA"
    )

    # build splicing chart
    splicing = alt.Chart(spliced).mark_bar().encode(
        x = alt.X('start', scale=alt.Scale(domain=domain), title=f"DNA Base Pair Position on {chromosome}"),
        x2 = 'end',
        y = alt.Y('chromosome:N', title=None),
        color = 'feature',
        tooltip = alt.Tooltip(['gene_name','feature','start','end'])
    ).properties(
        width = 750,
        title = f"2. Splicing: Untranslated Regions (UTR) of Exons are spliced out leaving only Coding DNA Sequences (CDS)"
    )

    translation = alt.Chart(translated).mark_bar().encode(
        x = alt.X('start', scale=alt.Scale(domain=domain), title=f"DNA Base Pair Position on {chromosome}"),
        x2 = 'end',
        y = alt.Y('chromosome:N', title=None),
        color = 'feature',
        tooltip = alt.Tooltip(['gene_name','feature','start','end'])
    ).properties(
        width = 750,
        title = f"3. Translation: The remaining CDS regions get translated into an Amino Acid Chain"
    )

    chart = alt.vconcat(transcription & splicing & translation)

    return chart.to_json()



# render altair chart of protein composition by amino acid
@app.route('/chart/protein_composition/<gene_name>')
def protein_composition(gene_name):

    # gene = genome[genome.gene_name == gene_name]
    path = 'data/genes/' + gene_name + '.csv.zip'
    gene = pd.read_csv(path)

    gene = gene[gene.feature == 'gene'].reset_index()
    protein_name = gene.loc[0, 'protein_name']
    aa_seq = gene.loc[0, 'aa_sequence']

    aa_dict = {
        'Cysteine (Cys)': 'C', 'Aspartic Acid (Asp)': 'D', 'Serine (Ser)': 'S',
        'Glutamine (Gln)': 'Q', 'Lysine (Lys)': 'K', 'Tryptophan (Trp)': 'W',
        'Asparagine (Asn)': 'N', 'Proline (Pro)': 'P', 'Threonine (Thr)': 'T',
        'Phenylalanine (Phe)': 'F', 'Alanine (Ala)': 'A', 'Glycine (Gly)': 'G',
        'Isoleucine (Ile)': 'I', 'Leucine (Leu)': 'L', 'Histidine (His)': 'H',
        'Arginine (Arg)': 'R', 'Methionine (Met)': 'M', 'Valine (Val)': 'V',
        'Glutamic Acid (Glu)': 'E', 'Tyrosine (Tyr)': 'Y'
        }

    aa_count = pd.DataFrame(columns = ['Amino_Acid','Count'])
    for key, value in aa_dict.items():
        aa_count.loc[value] = [key, 0]

    for aa in aa_seq:
        aa_count.loc[aa,'Count'] += 1

    amino_acid_composition = alt.Chart(aa_count).mark_bar().encode(
        x = alt.X('Count', title='Count of Amino Acid in Protein'),
        y = alt.Y('Amino_Acid:N', sort='-x', title='Amino Acid'),
        tooltip = alt.Tooltip(['Amino_Acid','Count']),
        color = 'Amino_Acid'
    ).properties(
        width = 750,
        title=f'Count of each Amino Acid in Protein from Gene: {gene_name}'
    )

    return amino_acid_composition.to_json()

##############################
# Run Flask
##############################

# Run Flask App
if __name__ == '__main__':
    app.run(debug=True)
