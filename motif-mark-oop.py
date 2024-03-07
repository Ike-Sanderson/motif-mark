#!/usr/bin/env python
from __future__ import annotations

import cairo
import re
import argparse

parser = argparse.ArgumentParser(description="A program to input a fasta file and a motif file and output an image of genes, exons, and motifs")
parser.add_argument("-f", type=str, help="filename1.fasta", required=True)
parser.add_argument("-m", type=str, help="filename2.txt", required=True)

args = parser.parse_args()
# Global variables
fasta_input = args.f
motifs = args.m
INTERMEDIATE_FASTA = 'intermediate.fa'
MARGIN_X = 10

motif_dict: dict[str, str] = {} #create dictionary with ambiguous motif as key: all possible pyrim. replacement combos as values
#eg 'YGCY': '[CTU]GC[CTU]'  
# motif_color_dict = {} # create dictionary that will hold motif as key and pycairo RGB as value

def process_fasta(file:str):
    '''This function generates an intermediate one-line fasta file.'''
    with open(fasta_input, "r") as fin, open(INTERMEDIATE_FASTA, 'w') as fout:
        first_time = True
        for line in fin:
            line = line.strip()
            if line.startswith('>'):  # If it's a header line
                if not first_time:
                    fout.write("\n")
                fout.write(line)
                fout.write("\n")
                first_time = False
            else:
                fout.write(line)                

def process_intermediate_fasta(file: str) -> list[Gene]:
    my_gene_list = []
    my_exon_list = []
    with open(file) as fin:
        gene_num = 1
        y = 50
        while True:     # this savvy while true loop works such that every time through is one iteration of a gene
            header = fin.readline()
            seq = fin.readline()
            if header == "":
                # EOF
                break
            
            ### start gene processing ###
            match = re.match(r'^>(\S+)', header)  # Extract gene name using regex; name will be whatever comes after >
            gene_name = match.group(1)
            bp_length = len(seq)
            gene = Gene(gene_name, gene_num, bp_length)
            gene.draw_gene(context,gene_num)
            my_gene_list.append(gene)
            gene_num += 1
            
            ### start exon processing ###
            exon_case_matches = [(m.start(0), m.end(0)) for m in re.finditer('[A-Z+]+', seq)]
            exon_start = exon_case_matches[0][0] 
            exon_end = exon_case_matches[0][1]
            exon_length = exon_end - exon_start + 1
            exon = Exon(gene_num)
            exon.draw_exon(context, exon_start, exon_length, (gene_num - 1))
            my_exon_list.append(exon)
            
            ### start motif processing ###
            fasta_seq_as_upper = seq.upper()
            # e.g., ambiguous_motif = YYY
            for ambiguous_motif, regex_motif in motif_dict.items():
                motif_start_indices = []   
                motif_matches = re.finditer(regex_motif, fasta_seq_as_upper) 
                for match in motif_matches:
                    motif_start = match.start()
                    motif_end = match.end()
                    motif_size = motif_end - motif_start
                    motif = Motif()
                    motif.draw_motif(motif_start, motif_size, (gene_num - 1), ambiguous_motif)
            ### end processing
            y += 100 # advance y value for the next gene
    return my_gene_list

def read_in_motifs(motifs):
    with open(motifs,'r') as fin2:
        for line in fin2:
            line = line.strip()
            line = line.upper() #turning all motif bases into upper case (i will use upper case sequences later)
            replacement_motifs = re.sub(r'Y', '[CTU]', line) # replace Y with any of three valid pyrimidine bases: 'C', 'T', and 'U'
            motif_dict[line] = replacement_motifs
    return motif_dict


class Gene:
    def __init__(self, name: str, gene_num: int, length):
        ## Data ##
        self.name = name
        self.gene_num = gene_num
        self.length = length

    ## Methods ##
    def __repr__(self) -> str:
        return f"Gene({self.name}, {self.gene_num}, {self.length})"
    
    def draw_gene(self, context, gene_num):
        self.context = context
        # Set line width, font size, and color 
        context.set_line_width(4)
        context.set_font_size(30)
        context.select_font_face("Arial",cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        context.move_to(MARGIN_X, (gene_num * 100) + 45)  # Draw gene name slightly below the line
        context.set_source_rgb(0, 0, 0)
        context.show_text(self.name)        
        context.stroke()
        # Draw gene in black
        context.set_source_rgb(0, 0, 0)  # Black color for base pair lines
        context.move_to(MARGIN_X, (gene_num * 100))
        context.line_to(MARGIN_X + self.length, (gene_num * 100))
        context.stroke()


class Exon:
    def __init__(self, gene_num: int):
        self.gene_num = gene_num

    #def build_exon(self, context, exon_start, exon_length, y) -> Exon:
    def draw_exon(self, context, exon_start, exon_length, gene_num):
        self.context = context
        context.set_line_width(24)
        context.set_source_rgb(0.4, 0.9, 0.4)  # green
        context.move_to(MARGIN_X + exon_start, (gene_num * 100))
        context.line_to(MARGIN_X + (exon_start + exon_length), (gene_num * 100))
        context.stroke()

class Motif:
    def __init__(self):
        pass
    ## Methods ##
    def draw_motif(self,motif_start, motif_size, gene_num, ambiguous_motif):
        self.context = context
        rectangle_height = 24
        context.set_line_width(1)
        if ambiguous_motif == "YGCY":
            context.set_source_rgb(1, 0, 0)  # red
        elif ambiguous_motif == "GCAUG":
            context.set_source_rgb(1, 0, 1)  # pink
        elif ambiguous_motif == "CATAG":
            context.set_source_rgb(0, 0, 1)  # blue
        elif ambiguous_motif == "YYYYYYYYYY":
            context.set_source_rgb(0, 0, 0)  # black
        else:
            context.set_source_rgb(0, 0, 0)
        context.move_to(MARGIN_X + motif_start, 0)
        context.rectangle((MARGIN_X + motif_start), (gene_num * 100) - 12, motif_size, rectangle_height) 
        context.fill()

## Main Code ##
read_in_motifs(motifs)
# Surface and context
surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, (771+ 50), 600) # after data exploration, using length of max gene length as width of context#
context = cairo.Context(surface)

# Set background color to white
context.set_source_rgb(1,1,1) # white
context.rectangle(0,0,(771+ 50), 600) #using length of max gene length as width of context#
context.fill()

#process_fasta("/Users/ikesanderson/bioinfo/Bi625/oop_motif/Figure_1.fasta") #this calls the function to read the fasta file
process_fasta(fasta_input) #this calls the function to read the fasta file
genes = process_intermediate_fasta(INTERMEDIATE_FASTA) # makes intermediate file in which all processing will occur
context.set_source_rgb(0, 0, 0) # setting color to black here as base color resets it from white for subsequent draw commands

##legend
#legend font size, type
context.set_font_size(24)
context.select_font_face("Arial",cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
#legend labels and colors
legend_labels = ["ygcy", "GCAUG", "catag", "YYYYYYYYYY"]
legend_colors = [(1, 0, 0), (1, 0, 1), (0, 0, 1), (0, 0, 0)]

#set starting position
x = 50
y = 500
# draw legend
for label, color in zip(legend_labels, legend_colors):
    # set color for legend text
    context.set_source_rgb(*color)
    # draw legend text
    context.move_to(x,y)
    context.show_text(label)
    y += 30 #move to next line

#draw title
context.set_source_rgb(0, 0, 0)
context.set_font_size(35)
context.select_font_face("Arial",cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
context.move_to(203,50)
context.show_text("Motif Mark - Ike Sanderson")

# make output file name = file name of input, saved as .png
# eg Figure_1.fasta -> Figure_1.png
name = fasta_input.split(".")
output_name = f'{name[0]}.png'
surface.write_to_png(output_name) # Save the surface to a PNG file

#################################
