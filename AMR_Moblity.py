#This script was made by Alan Risal (geneCodeAlan)

#Pipeline involves interpreting the results from ABRicate and mobileelementfinder csv outputs 
#It is specifically for these two forms of csv outputs but can be adjusted 

#The code is looking for colocations---that is proximal instances between AMR genes and mobile elements
#By parsing and cross checking the values of AMR gene indices with IS segments, AMobilityReader is used to create a numerical 
#representation of the risk of mobility of a specific gene. 
#Obviously plasmid contigs are assumed mobile however it is necccessary to observe proximal relationships within chromosomal 
#structures to predict mobility of seemingly dormant AMR genes. 

#By coupling two powerful bioinformatics tools we can help make sense of mobility, additionally other factors can be worked into this calculation
#for example: TE's ICE's Integrons etc. Many more to be added to the sum of weights predicting mobility 

#Added pipeline features -- 
# 1. Visualization return results in MatPlotLib using Biopython to format the prokaryotic gene accordingly 
# 2. TE identifiers -- Can stick with CARD database -- but should also run a TE screener to find Transposase genes
#   a) Using transposase genetic values -- we identify transposase proximity to the other genes and label it
#   b) Since IS and MGE Finder already labels TE's we don't need to refile it, just need to map transposase's 


#Resources and imports used 
# -- ABRicate --db card for tabular result of AMR genes
# -- mobilegenefinder -- good job at locating IS -- want to try with a different one 
""""
IS - insertion sequence
cn - composite transposon (named after length of intermediary region and defining ISs, e.g., cn_12127_IS26)
Tn - (unit) transposon
MITE - miniature inverted repeat
ICE - integrative conjugative element
CIME - cis-mobilizable element
IME - integrative mobilizable element

"""
# -- integron finder -- Galaxy local 



import csv
import tkinter as tk
from tkinter import filedialog

#Graphing packs
from Bio.Graphics import GenomeDiagram
from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import os

from reportlab.pdfgen import canvas
from PyPDF2 import PdfReader, PdfWriter

def add_title_to_pdf(input_pdf, output_pdf, title):
    # Create a temporary PDF with the title
    temp_title_pdf = "temp_title.pdf"
    c = canvas.Canvas(temp_title_pdf, pagesize=(20*cm, 20*cm))
    c.setFont("Helvetica-Bold", 24)
    c.drawCentredString(10*cm, 19*cm, title)
    c.save()

    # Merge the title PDF and the map PDF
    reader1 = PdfReader(temp_title_pdf)
    reader2 = PdfReader(input_pdf)
    writer = PdfWriter()

    # Overlay title on the first page
    page = reader2.pages[0]
    page.merge_page(reader1.pages[0])
    writer.add_page(page)

    # Add any additional pages from the original map PDF
    for i in range(1, len(reader2.pages)):
        writer.add_page(reader2.pages[i])

    with open(output_pdf, "wb") as f_out:
        writer.write(f_out)

    # Clean up temp file
    os.remove(temp_title_pdf)

def load_amr_genes(tabular_path):

    """
    Load AMR genes from ABRicate tabular file (tab-delimited, this is just how ABRicate returns its files)
    It is not csv and instead in .tab so we need to adjust a little
    
    """
    amr_genes = []
    with open(tabular_path, 'r', newline='') as f:
        # Find the header line (starts with # in the Abricate result, #File)
        while True:
            header_line = f.readline()
            if not header_line:
                raise ValueError("No header line found in AMR file.")
            if header_line.startswith('#'):
                break
        # Remove the leading '#' and any whitespace, then split for header (tab-delimited)
        header = header_line.lstrip('#').strip().split('\t')
        print("Parsed header fields:", header)
        reader = csv.DictReader(f, fieldnames=header, delimiter='\t')
        for i, row in enumerate(reader):
            print(f"Row {i}")  #print each row
            try:
                #Load each gene to confirm parsing
                start = int(row['START'].strip())
                end = int(row['END'].strip())
                amr_genes.append({
                    'gene': row['GENE'],
                    'start': min(start, end),
                    'end': max(start, end),
                    'product': row['PRODUCT'],
                    'resistance': row['RESISTANCE']
                })
                print(f"Loaded AMR gene: {row['GENE']} ({min(start, end)}-{max(start, end)})")
            except (ValueError, KeyError) as e:
                print(f"Skipping AMR row due to error: {e} | Row: {row}")
                continue
    return amr_genes

def load_is_elements(csv_path):
    """Load IS elements from mobile element finder CSV"""
    is_elements = []
    with open(csv_path, 'r') as f:
        # Read and filter out comment lines
        lines = [line for line in f if not line.startswith('#')]
        reader = csv.DictReader(lines)
        for row in reader:
            if row.get('type') in ['insertion sequence', 'composite transposon', 'ice']:
                try:
                    start = int(row['start'])
                    end = int(row['end'])
                    is_elements.append({
                        'name': row['name'],
                        'start': min(start, end),
                        'end': max(start, end),
                        'type': row['type']
                    })
                    print(f"Loaded IS element: {row['name']} ({min(start, end)}-{max(start, end)})")
                except (ValueError, KeyError) as e:
                    print(f"Skipping IS row due to error: {e} | Row: {row}")
                    continue
    return is_elements

def check_proximity(amr_gene, is_element, max_distance=15000): #adjust threshold if you would like! 
    """Check if AMR gene overlaps or is within proximity threshold of IS/MGE"""
    # Overlap: any base in common of the indices
    if (amr_gene['start'] <= is_element['end'] and amr_gene['end'] >= is_element['start']):
        if (is_element['start'] <= amr_gene['start'] and amr_gene['end'] <= is_element['end']):
            print("  -> gene_within_IS")
            return 'gene_within_IS', 0
        else:
            print("  -> overlap")
            return 'overlap', 0

    distance_before = amr_gene['start'] - is_element['end']
    distance_after = is_element['start'] - amr_gene['end']
    min_distance = min(distance_before, distance_after, key=abs)
    if abs(min_distance) <= max_distance:
        print(f"  -> within_15kb ({abs(min_distance)} bp)")
        return 'within_15kb', abs(min_distance)

    return None, None

def find_proximities(amr_genes, is_elements):
    """Find all AMR genes in proximity to IS elements"""
    results = []
    for gene in amr_genes:
        for is_elem in is_elements:
            relation, distance = check_proximity(gene, is_elem)
            if relation:
                results.append({
                    'amr_gene': gene['gene'],
                    'amr_start': gene['start'],
                    'amr_end': gene['end'],
                    'is_name': is_elem['name'],
                    'is_type': is_elem['type'],
                    'is_start': is_elem['start'],
                    'is_end': is_elem['end'],
                    'relationship': relation,
                    'distance': distance,
                    'product': gene['product'],
                    'resistance': gene['resistance']
                })
    return results

def plot_circular_map(amr_genes, is_elements, seq_length, job_title, output_file="chromosome_map.pdf"):
    gd_diagram = GenomeDiagram.Diagram("AMR and IS elements")
    amr_track = gd_diagram.new_track(1, name="AMR Genes", greytrack=False)
    amr_set = amr_track.new_set()
    is_track = gd_diagram.new_track(2, name="IS Elements", greytrack=False)
    is_set = is_track.new_set()

    # Add AMR genes
    for gene in amr_genes:
        feature = SeqFeature(
            FeatureLocation(gene['start'], gene['end'], strand=1)
        )
        amr_set.add_feature(
            feature,
            color=colors.red, label=True, name=gene['gene']
        )

    # Add IS elements
    for elem in is_elements:
        feature = SeqFeature(
            FeatureLocation(elem['start'], elem['end'], strand=1)
        )
        is_set.add_feature(
            feature,
            color=colors.blue, label=True, name=elem['name']
        )

    gd_diagram.draw(
    format="circular",
    circular=True,
    pagesize=(20*cm, 20*cm),
    start=0,
    end=seq_length,
    )

    # Ensure output directory exists
    output_dir = "chromosome_maps"
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, f"{job_title}_chromosome_map.pdf")
    gd_diagram.write(output_path, "PDF")
    add_title_to_pdf(output_path, output_path, str(job_title))

    print(f"Map saved to {job_title}_chromosome_map.pdf")

def extract_proximal_features(proximities):
    # Use sets to avoid duplicates
    amr_seen = set()
    is_seen = set()
    amr_genes = []
    is_elements = []
    for p in proximities:
        amr_key = (p['amr_gene'], p['amr_start'], p['amr_end'])
        is_key = (p['is_name'], p['is_start'], p['is_end'])
        if amr_key not in amr_seen:
            amr_genes.append({
                'gene': p['amr_gene'],
                'start': p['amr_start'],
                'end': p['amr_end'],
                'product': p['product'],
                'resistance': p['resistance']
            })
            amr_seen.add(amr_key)
        if is_key not in is_seen:
            is_elements.append({
                'name': p['is_name'],
                'start': p['is_start'],
                'end': p['is_end'],
                'type': p['is_type']
            })
            is_seen.add(is_key)
    return amr_genes, is_elements

def main():
    root = tk.Tk()
    root.withdraw()

    job_title = input("Enter job Title of Run, Ensure Uniqueness: ")
    print("Please select the ABRicate .tabular file (AMR gene results) in the file dialog.")
    amr_tabular_path = filedialog.askopenfilename(
        title="Select ABRicate .tabular file",
        filetypes=[("Tabular files", "*.tab"), ("All files", "*.*")]
    )
    if not amr_tabular_path:
        print("No ABRicate file selected. Exiting.")
        return

    is_csv_path = filedialog.askopenfilename(
        title="Select MobileElementFinder CSV file",
        filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
    )
    if not is_csv_path:
        print("No IS element CSV file selected. Exiting.")
        return
    print("Please select the full FASTA sequence")
    ref_seq = filedialog.askopenfilename(
        title="Select full FASTA file",
        filetypes=[("FASTA files", "*.fna"),("All files", ".")]
    )
    if not ref_seq:
        print("No IS element CSV file selected. Exiting.")
        return

    # Get sequence length from FASTA
    with open(ref_seq, "r") as handle:
        record = next(SeqIO.parse(handle, "fasta"))
        length = len(record.seq)

    # Load data
    amr_genes = load_amr_genes(amr_tabular_path)
    is_elements = load_is_elements(is_csv_path)
    
    # Find proximities
    proximities = find_proximities(amr_genes, is_elements)
    
    # Output results
    if proximities:
        print("AMR genes in proximity to IS elements:")
        print("{:<15} {:<15} {:<15} {:<15} {:<15} {:<20} {:<10}".format(
            'AMR Gene', 'IS Element', 'Relationship', 'Distance(bp)', 
            'AMR Location', 'Product', 'Resistance'))
        for p in proximities:
            amr_loc = f"{p['amr_start']}-{p['amr_end']}"
            print("{:<15} {:<15} {:<15} {:<15} {:<15} {:<20} {:<10}".format(
                p['amr_gene'],
                p['is_name'],
                p['relationship'],
                p['distance'],
                amr_loc,
                p['product'][:17] + '...' if len(p['product']) > 20 else p['product'],
                p['resistance'][:10] + '...' if len(p['resistance']) > 13 else p['resistance']
            ))
    else:
        print("No AMR genes found in proximity to IS elements")
    
    output_choice = input("Return (A)ll annotated AMR/IS or (P)roximal only? [A/P]: ").strip().upper()
    if output_choice == "A":
        job_title += "_All"
        plot_circular_map(amr_genes, is_elements, length, job_title)
        pass
    elif output_choice == "P":
        # Only plot proximal AMR/IS features
        job_title += "_Proximal"
        prox_amr, prox_is = extract_proximal_features(proximities)
        plot_circular_map(prox_amr, prox_is, length, job_title)
        pass
    else: 
        print("Invalid Key Value")
        pass


if __name__ == "__main__":
    main()
