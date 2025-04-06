import os, sys
import numpy as np
import pandas as pd
import time
import itertools as it
import pickle
import subprocess




#--------------------------------------------------------------------------------------------
def parse_abundance(species, anno):
    """
    Parse and average kallisto abundance estimates
    """

    # load and merge individual files
    list_df = []
    for i in list(anno['Run']):
        abundance = pd.read_csv('../data/kallisto/kallisto_'+species+'/kallisto_'+i+'/abundance.tsv', header=0, index_col=False, sep='\t')
        tmpdf = pd.DataFrame(data=np.array(abundance['tpm']), columns=[i], index=list(abundance['target_id']) )
        list_df.append(tmpdf)
    abundanceDF = pd.concat(list_df, axis=1)

    # betweenArrays quantile normalization
    data = np.array(abundanceDF)
    #mm = np.mean(data, axis=1)
    mm = np.min(data, axis=1) > 1
    data = data[mm,:]
    newIndex = abundanceDF.index[mm]

    np.savetxt("tmp.all", np.around(data, 5), fmt='%.5f'  )
    cmd_bowtie = 'R --vanilla < qn.R > log'
    output = subprocess.run(cmd_bowtie, shell=True)
    data_q = np.loadtxt("tmp.all.N")
    os.remove("tmp.all")
    os.remove("tmp.all.N")

    abundanceDF = pd.DataFrame(data=data_q, columns=abundanceDF.columns, index=newIndex)

    # average over replicas
    res = np.zeros(( len(abundanceDF), 3 )) * np.nan
    for ix, i in enumerate(sorted(list(set(anno['AGE'])))):
        current_idx = list(anno[anno['AGE']==i]['Run'])
        current_df = abundanceDF[current_idx]
        res[:,ix] = np.nanmean(current_df, 1)

    resDF = pd.DataFrame(data=np.around(res, 3), index=abundanceDF.index, columns=['young', 'middle', 'old'])
    resDF.to_csv('../data/processed/'+species+'_abundance.txt', header=True, index=True, sep='\t')

    return resDF





#--------------------------------------------------------------------------------------------
def codon_demand(seqs, abundance, species):
    """
    Computes codon usage as codon frequences weighted by transcript abundance
    """
    result = np.zeros(( 61, 3 ))
    for ax, age in enumerate(['young', 'middle', 'old']):
        res = np.zeros(( len(codons_nonstop) ))
        for i in list(seqs.keys()):
            current_seq = seqs[i]
            current_counts = np.zeros(( len(codons_nonstop) ))
            if i in list(abundance.index):
                current_abund = abundance.loc[i][age]
                for jx, j in enumerate(codons_nonstop):
                    current_counts[jx] = np.sum(np.array(current_seq) == j)

                current_demand = current_abund * current_counts
                res += current_demand
    
        demand = res/np.sum(res)
        #demand = res 

        result[:,ax] = demand

    resultDF = pd.DataFrame(data=np.around(result, 4), columns=['young', 'middle', 'old'], index=codons_nonstop)
    resultDF.to_csv('../data/processed/'+species+'_codonusage.txt', header=True, index=True, sep='\t')

    return resultDF






if __name__ == '__main__':


    codons_nonstop = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT', 
    		  'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT',
    		  'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT',
    		  'TAC', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT',]
   
    sequence_codons = pickle.load(open('../data/processed/scer/scer_scikit_codonseq.pkl', 'rb'))
    sequence_cele = pickle.load(open('../data/processed/cele/cele_scikit_codonseq.pkl', 'rb'))
  


    scer_anno = pd.read_csv("../data/processed/scer/scer_anno.txt", header=0, index_col=False, sep='\t')
    cele_anno = pd.read_csv("../data/processed/cele/cele_anno.txt", header=0, index_col=False, sep='\t')

    scer_abundance = parse_abundance('scer', scer_anno)
    cele_abundance = parse_abundance('cele', cele_anno)

    res = codon_demand(sequence_codons, scer_abundance, 'scer')
    res = codon_demand(sequence_cele, cele_abundance, 'cele')

    scer_tAI = pd.read_csv("../data/yeast.ws.gcn", header=0, index_col=None, sep='\t')
    scer_CSC = pd.read_csv("../data/CSC/scer_CSC.txt", header=0, index_col=0, sep='\t')

    # bring different codon indices together in DF
    nTEdf = pd.DataFrame(columns=['codon', 'tAI', 'CSC', 'CU_young', 'CU_middle', 'CU_old'])
    for i in res.index:
        current_cu = np.array( res.loc[i] )
        current_tai = list(scer_tAI[scer_tAI['Codon']==i]['tAI'])[0]
        current_csc = scer_CSC.loc[i]['Coller']
        nTEdf.loc[len(nTEdf)] = (i, current_tai, current_csc, current_cu[0], current_cu[1], current_cu[2])
    nTEdf.to_csv("../data/processed/optimality_cu.txt", header=True, index=False, sep='\t')