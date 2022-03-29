##############################
# Flask Implementation
# Visualizing the Human Genome
# Updated 03/25/2022
##############################

from flask import Flask, render_template
import altair as alt
import pandas as pd
import json
import zipfile


# Create Flask App Instance
app = Flask(__name__)

##############################
# Genome Data
##############################

# direct_url="https://drive.google.com/uc?export=download&confirm=9iBg&id=1-8H7V0kzGlx6vPpRM4koQ7vr_KtIXgKn"
# genome = pd.read_csv(direct_url)

# chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
# for chromosome in chromosomes:
#     path = 'data/chromosomes/' + chromosome + '.csv.zip'
#     if chromosome == chromosomes[0]:
#         genome = pd.read_csv(path)
#     else:
#         data = pd.read_csv(path)
#         genome = pd.concat([genome, data], ignore_index=True)
genome = pd.read_csv('data/genome_genes/genome.csv.zip')
gene_names = list(genome.gene_name.unique())
gene_names.sort()

##############################
# Flask routes
##############################

# Render index.html homepage
@app.route('/')
def index(gene_names=gene_names):
    return render_template('index.html', gene_names=gene_names)

# Render molstar.html viewer of selected gene
@app.route('/molstar/<gene_name>')
def molstar(genome=genome, gene_name='TP53'):
    gene = genome[genome.gene_name == gene_name]
    gene = gene[gene.feature == 'gene'].reset_index()
    uniprot_id = gene.loc[0, 'uniprot_id']
    return render_template('molstar.html', uniprot_id=uniprot_id)

# Render definitions.html for part 1
@app.route('/definitions_part_1')
def definitions_part_1():
    return render_template('definitions_part_1.html')

# Render definitions.html for part 2
@app.route('/definitions_part_2')
def definitions_part_2():
    return render_template('definitions_part_2.html')

# Render definitions.html for part 3
@app.route('/definitions_part_3')
def definitions_part_3():
    return render_template('definitions_part_3.html')

# Render definitions.html for part 4
@app.route('/definitions_part_4')
def definitions_part_4():
    return render_template('definitions_part_4.html')

# Render definitions.html for parts 5, 6, 7
@app.route('/definitions_part_5_6_7')
def definitions_part_5_6_7():
    return render_template('definitions_part_5_6_7.html')


##############################
# Altair routes
##############################

# Render altair chart of gene overview
@app.route('/chart/gene_overview/')
def gene_overview():
    chart = json.load('static/charts/gene_overview.vl.json')
    return chart

# Render altair chart of gene components
@app.route('/chart/gene_components/<gene_name>')
def gene_components(gene_name='TP53'):

    # gene = genome[genome.gene_name == gene_name]
    path = 'data/genes/' + gene_name + '.csv.zip'
    gene = pd.read_csv(path)

    transcripts = gene[gene.feature != 'gene']

    domain_min = gene.start.min()
    domain_max = gene.end.max()
    domain = [domain_min, domain_max]

    bar = alt.Chart(transcripts).mark_bar().encode(
        x=alt.X('start:Q', scale=alt.Scale(domain=domain)),
        x2='end:Q',
        y='feature:N',
        color = 'transcript_type:N',
        tooltip = alt.Tooltip(['gene_name','feature','transcript_type','exon_id','start','end'])
    ).properties(
        width = 750,
        title=f'Components of Gene {gene_name}'
    ).interactive()

    hist_features = alt.Chart(transcripts).mark_bar().encode(
        y = 'count(feature)',
        x = alt.X('feature:N', sort='-y')
    ).properties(width = 350)

    hist_transcript_type = alt.Chart(transcripts).mark_bar().encode(
        y = 'count(transcript_type)',
        x = alt.X('transcript_type:N', sort='-y')
    ).properties(width = 350)

    gene_layout = bar & (hist_features | hist_transcript_type)

    return gene_layout.to_json()


# render altair chart of gene transcription
@app.route('/chart/transcription/<gene_name>')
def transcription(gene_name='TP53'):

    # gene = genome[genome.gene_name == gene_name]
    path = 'data/genes/' + gene_name + '.csv.zip'
    gene = pd.read_csv(path)

    transcripts = gene[gene.feature != 'gene']
    exons = transcripts[transcripts.feature == 'exon']
    exons = exons[exons.transcript_type != 'nonsense_mediated_decay']

    domain_min = gene.start.min()
    domain_max = gene.end.max()
    domain = [domain_min, domain_max]

    transcription = alt.Chart(exons).mark_bar().encode(
        x = alt.X('start', scale=alt.Scale(domain=domain), title='Location of Gene Features'),
        x2 = 'end',
        y = alt.Y('seqname', title='Chromosome'),
        color = 'feature',
        tooltip = alt.Tooltip(['gene_name','feature','transcript_type','exon_id','start','end'])
    ).properties(
        width = 750,
        title = f"Transcribed Exons of Gene {gene_name}"
    ).interactive()

    return transcription.to_json()


# render altair chart of gene splicing
@app.route('/chart/splicing/<gene_name>')
def splicing(gene_name='TP53'):

    # gene = genome[genome.gene_name == gene_name]
    path = 'data/genes/' + gene_name + '.csv.zip'
    gene = pd.read_csv(path)

    transcripts = gene[gene.feature != 'gene']
    spliced = transcripts[transcripts.feature.isin(['UTR','CDS','start_codon','stop_codon'])]
    spliced = spliced[spliced.transcript_type != 'nonsense_mediated_decay']

    domain_min = gene.start.min()
    domain_max = gene.end.max()
    domain = [domain_min, domain_max]

    splicing = alt.Chart(spliced).mark_bar().encode(
        x = alt.X('start', scale=alt.Scale(domain=domain), title='Location of Gene Features'),
        x2 = 'end',
        y = alt.Y('seqname', title='Chromosome'),
        color = 'feature',
        tooltip = alt.Tooltip(['gene_name','feature','transcript_type','exon_id','start','end'])
    ).properties(
        width = 750,
        title = f"Spliced mRNA Components of Gene {gene_name}"
    ).interactive()

    return splicing.to_json()


# render altair chart of gene translation
@app.route('/chart/translation/<gene_name>')
def translation(gene_name='TP53'):

    # gene = genome[genome.gene_name == gene_name]
    path = 'data/genes/' + gene_name + '.csv.zip'
    gene = pd.read_csv(path)

    transcripts = gene[gene.feature != 'gene']
    spliced = transcripts[transcripts.feature.isin(['UTR','CDS','start_codon','stop_codon'])]
    spliced = spliced[spliced.transcript_type != 'nonsense_mediated_decay']
    translated = spliced[spliced.feature.isin(['CDS','start_codon','stop_codon'])]

    domain_min = gene.start.min()
    domain_max = gene.end.max()
    domain = [domain_min, domain_max]

    translation = alt.Chart(translated).mark_bar().encode(
        x = alt.X('start', scale=alt.Scale(domain=domain), title='Location of Gene Features'),
        x2 = 'end',
        y = alt.Y('seqname', title='Chromosome'),
        color = 'feature',
        tooltip = alt.Tooltip(['gene_name','feature','transcript_type','exon_id','start','end'])
    ).properties(
        width = 750,
        title = f"Translated CDS Regions of Gene {gene_name}"
    ).interactive()

    return translation.to_json()


# render altair chart of protein composition by amino acid
@app.route('/chart/protein_composition/<gene_name>')
def protein_composition(gene_name='TP53'):

    # gene = genome[genome.gene_name == gene_name]
    path = 'data/genes/' + gene_name + '.csv.zip'
    gene = pd.read_csv(path)

    gene = gene[gene.feature == 'gene'].reset_index()
    protein_name = gene.loc[0, 'protein_name']
    aa_seq = gene.loc[0, 'aa_sequence']

    aa_dict = {'Cys': 'C', 'Asp': 'D', 'Ser': 'S', 'Gln': 'Q', 'Lys': 'K',
        'Trp': 'W', 'Asn': 'N', 'Pro': 'P', 'Thr': 'T', 'Phe': 'F', 'Ala': 'A',
        'Gly': 'G', 'Ile': 'I', 'Leu': 'L', 'His': 'H', 'Arg': 'R', 'Met': 'M',
        'Val': 'V', 'Glu': 'E', 'Tyr': 'Y'}

    aa_count = pd.DataFrame(columns = ['Amino_Acid','Count'])
    for key, value in aa_dict.items():
        aa_count.loc[value] = [key, 0]

    for aa in aa_seq:
        aa_count.loc[aa,'Count'] += 1

    amino_acid_composition = alt.Chart(aa_count).mark_bar().encode(
        x = 'Count',
        y = alt.Y('Amino_Acid:N', sort='-x'),
        tooltip = alt.Tooltip(['Amino_Acid','Count'])
    ).properties(
        width = 750,
        title=f'Distribution of Amino Acids in Protein: {protein_name}'
    )

    return amino_acid_composition.to_json()

##############################
# Run Flask
##############################

# Run Flask App
if __name__ == '__main__':
    app.run(debug=True)
