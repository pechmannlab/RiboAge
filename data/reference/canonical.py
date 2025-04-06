import os, sys
import numpy as np
import pandas as pd
import gzip as gz
import glob as glob
import multiprocessing as mp

from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC







#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def extract_uniprot_ids(uniprotIN):
    """
    match Uniprot ID to Ensembl canonical transcript ID's
    from Uniprot annotation
    """

    uniprot = pd.read_csv(uniprotIN, header=0, index_col=False, sep='\t')
    print(uniprot)


    """
    resultDF = pd.DataFrame(columns=['uniprot', 'canonical'])
    for i in range(len(uniprot)):
        current_entry = uniprot.iloc[i]
        current_uni = current_entry['Entry']
        current_ensembl = str(current_entry['Ensembl'])
        current_altern = str(current_entry['Alternative products (isoforms)'])

        if current_ensembl[0:4] == "ENST": 
            if "Sequence=Displayed;" in current_altern:
                current_canonical = current_altern.split("Sequence=Displayed;")[0].split("IsoId=")[-1].replace(';','').split()[0]

                current_ensembl = current_ensembl.split(';')
                for j in range(len(current_ensembl)):
                    j_ensembl = current_ensembl[j]
                    if '[' in j_ensembl:
                        j_ensembl = j_ensembl.split('[')
                        current_key = j_ensembl[1].replace(']', '').split()[0]
                        #print(np.array(current_key.split()[0]), np.array(current_canonical.split())[0])
                        if np.array(current_key.split()[0]) == np.array(current_canonical.split()[0]):        # wrangle weird string formatting
                            canonical_transcript = j_ensembl[0]
                            break

            else:
                list_ensembl = current_ensembl.strip(';').split(';')
                canonical_transcript = list_ensembl[0]

            resultDF.loc[len(resultDF)] = (current_uni, canonical_transcript.split('.')[0])

    print(resultDF)
    #resultDF.to_csv("../data/processed/uniprot_canonical.txt", header=True, index=False, sep='\t')

    return resultDF
    """



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def cdsseq_dict(fileIN):
    """
    parse dictionary of sequences
    """

    seq = SeqIO.parse(fileIN, "fasta")

    seq_dict = {}
    for s in seq:
        current_name = s.name
        if '|' in current_name:
            current_name = current_name.split('|')[1]
        else:
            current_name = current_name.split(' ')[0]


        current_seq = s.seq
        seq_dict[current_name] = current_seq

    return seq_dict














if __name__ == '__main__':



    uniprot = pd.read_csv('uniprotkb_taxonomy_name_elegans_AND_mod_2024_10_05.tsv', header=0, index_col=False, sep='\t')
    #print(uniprot[['Entry', 'Gene Names']])
    #print(uniprot.columns)


    protein_seq = cdsseq_dict('Caenorhabditis_elegans.fa')
    #print(protein_seq)


    rna_seq = cdsseq_dict('Caenorhabditis_elegans.WBcel235.cds.all.fa')
    #print(rna_seq)
    print(list(rna_seq.keys())[0:10])

  
    """  
    canonicalDF = pd.DataFrame(columns=['uniprot', 'gene', 'canonicalID'])
    output = open('Caenorhabditis_elegans.WBcel235.cds.canonical.fa', 'w')
    list_j = []
    for i in list(protein_seq.keys()) :
        current_uni = list(uniprot[uniprot['Entry']==i]['Gene Names'])
        if len(current_uni) > 0:
            current_uni = str(current_uni[0]).split()
            current_id = current_uni[-1]
            #print(i, current_id)
            current_len = int(len(protein_seq[i]))

            list_current = []
            for j in list(rna_seq.keys()):
                if current_id == j[0:len(current_id)]:
                    list_current.append(j)

            for j in list_current:
                isoform_len = int(len(rna_seq[j])/3. - 1)
                print(i, current_id, current_len, j, isoform_len)

                if current_len == isoform_len and j[0:len(current_id)] not in list_j:
                    print(i, current_id, current_len, j, isoform_len)
                    canonicalDF.loc[len(canonicalDF)] = (i, current_id, j)

                    output.write(">" + j + "\n")
                    output.write( str(rna_seq[j]) + '\n')
                    list_j.append(j[0:len(current_id)])

    output.close()
    canonicalDF.to_csv("cele_canonical.txt", header=True, index=False, sep='\t')
    """


    list_rRNA = ['T09B4.23', 'F31C3.7', 'F31C3.11', 'F31C3.9', 'F31C3.8', 'MTCE.7', 
                 'MTCE.33', 'ZK218.16', 'ZK218.17', 'ZK218.18', 'ZK218.19', 'ZK218.12', 
                 'ZK218.20', 'Y102A5D.5', 'Y102A5D.6', 'Y102A5D.7', 'Y102A5D.8', 
                 'Y102A5D.9', 'Y102A5D.10', 'Y102A5D.11', 'Y102A5D.12', 'T27C5.18']

   
    tmp_tRNA = pd.read_csv("tRNA.txt", header=None, index_col=False, sep='\t')
    tmp_tRNA = tmp_tRNA[8] 
    list_tRNA = [x.split(';')[1].replace('transcript_id', '').replace(' ', '').replace('"','') for x in tmp_tRNA]
 
    #"""
    nc_seq = cdsseq_dict('c_elegans.PRJNA13758.WS292.ncRNA_transcripts.fa')

    output = open('Cele_rRNA.fasta', 'w')
    for i in set(list_rRNA + list_tRNA):
        if i in nc_seq.keys():
            current_seq = nc_seq[i]
            print(current_seq)
            current_seq = str(current_seq)
            current_seq = current_seq.replace("U", "T")
            output.write(">" + i + "\n")
            output.write( current_seq + '\n')
        else:
            print(i)
    output.close()
    #"""
