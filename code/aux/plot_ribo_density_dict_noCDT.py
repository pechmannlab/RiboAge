#!/usr/bin/env python

## Script without need for codon dwell times (CDT) to plot density profiles.
## Modified by Pedro Bordignon
## April 2018

## ----------------------------------------
## scikit-ribo
## ----------------------------------------
## a module for visualization riboseq data
## ----------------------------------------
## author: Han Fang
## contact: hanfang.cshl@gmail.com
## website: hanfang.github.io
## date: 10/28/2016
## ----------------------------------------

from __future__ import print_function, division
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import sys
import pandas as pd
import seaborn as sns
import numpy as np
# import pybedtools as pbt
import argparse
import multiprocessing
import json
import random





AA_cDNA_dict = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
"TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
"TAT":"Y", "TAC":"Y", "TAA":"STOP", "TAG":"STOP",
"TGT":"C", "TGC":"C", "TGA":"STOP", "TGG":"W",
"CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
"CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
"CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
"CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
"ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
"ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
"AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
"AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
"GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
"GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
"GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
"GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"}


def returnCoverageOnGene(geneName):
    # check if the gene exist
    if geneName not in set(riboDf["gene"]):
        sys.stderr.write("[error]\tthe database does not have gene: " + str(geneName) + "\n")
        return
    # retrieve the data for this gene, construct numpy array
    geneDf = riboDf[riboDf["gene"] == geneName]

    riboCnt = np.array(geneDf["ribosome_count"])
    codonNts = np.array(geneDf["codon"])
    codonAA = np.array(geneDf["AA"])

    genes_dict[geneName] = [riboCnt.tolist(), codonNts.tolist(), codonAA.tolist()]


## ----------------------------------------
## main
## ----------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="input data frame, required")
    # parser.add_argument("-c", help="codon dwell time, required")
    parser.add_argument("-g", help="gene of interest, or type [all] - automatically plot all genes, required")
    parser.add_argument("-o", help="output path, required")
    ## check if there is any argument
    if len(sys.argv) <= 1:
        parser.print_usage()
        sys.exit(1)
    else:
        args = parser.parse_args()
    ## process the file if the input files exist
    if (args.i != None and args.g != None and args.o != None):
        sys.stderr.write("[status]\tprocessing the input file: " + args.i + "\n")
        # genes_dict = {}
        manager = multiprocessing.Manager()
        genes_dict = manager.dict()
        df_fn = args.i
        # codon_fn = args.c
        gene = args.g
        output = args.o
        # create folders
        # cmd = "mkdir -p " + output + "/plots"
        # os.system(cmd)
        # all genes or one gene
        # fig = figures(df_fn, codon_fn, output)
        # fig.loadDat()
        if gene != "all":
            sys.stderr.write("[execute]\tplotting ribosome coverage for gene: " + str(gene) + "\n")
            fig.plotCoverageOnGene(gene)
        else:
            sys.stderr.write("[execute]\tcomputing ribosome coverage for each gene" + "\n")
            # fig.plotAllGenes()
            # returnAllGenes(df_fn, codon_fn, output)
            # codonDT = pd.read_table(codon_fn, header=0)
            riboDf = pd.read_table(df_fn,  header=0)
            # riboDf = pd.merge(riboDf, codonDT, how="left", on="codon")
            riboDf['AA'] = riboDf['codon'].map(AA_cDNA_dict)

            # loop over all genes and plot
            geneNames = set(riboDf["gene"])
            print(len(geneNames), "genes to parse")
            # cpus = 16 if multiprocessing.cpu_count() > 16 else max(1, multiprocessing.cpu_count() - 1)
            pool = multiprocessing.Pool(6)
            # [val for i, val in enumerate(random.sample(geneNames, 10))]
            pool.map_async(returnCoverageOnGene, geneNames)
            pool.close()
            pool.join()
            # cpt = 0
            # for i in geneNames:
            #     print("[status]\trunning for "+i+" gene.\t"+str(cpt))
            #     returnCoverageOnGene(riboDf, i)
            #     cpt += 1

            os.system("mkdir -p " + output)
            with open(output + '/ALL_genes_profile_dict.json', 'w') as gd:
            	json.dump(genes_dict.copy(), gd)

            sys.stderr.write("[status]\tFinished returning all gene profiles." + "\n")
        ## end
        sys.stderr.write("[status]\tPlotting module finished" + "\n")
    else:
        sys.stderr.write("[error]\tmissing argument" + "\n")
        parser.print_usage()
