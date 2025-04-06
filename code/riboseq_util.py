import os, sys
import numpy as np
import pandas as pd
import pickle
import itertools as it

from sklearn.metrics import roc_auc_score
from scipy import stats


# utility functions and postanalysis of Riboseq modelling


#--------------------------------------------------------------------------------------------
# parse abundance data of translation genes
scer_transl = pd.read_csv("../data/reference/scer_translation.txt", header=0, index_col=False, sep='\t')
scer_abund = pd.read_csv("../data/processed/scer_abundance.txt", header=0, index_col=0, sep='\t')
abundDF = pd.DataFrame(columns=['ORF', 'exp_young', 'exp_middle', 'exp_old', 'category', 'complex'])
for i in range(len(scer_transl)):
    current_orf = scer_transl.iloc[i]['Scer']
    current_category = scer_transl.iloc[i]['category']
    current_complex = scer_transl.iloc[i]['complex']
    if current_orf in scer_abund.index:
        current_abundance = np.around(np.array(scer_abund.loc[current_orf]), 3)
        abundDF.loc[len(abundDF)] = (current_orf, current_abundance[0], current_abundance[1], current_abundance[2], current_category, current_complex)
abundDF.to_csv("../data/processed/scer_translation_abundance.txt", header=True, index=False, sep='\t')



cele_transl = pd.read_csv("../data/reference/cele_translation.txt", header=0, index_col=False, sep='\t')
cele_abund = pd.read_csv("../data/processed/cele_abundance.txt", header=0, index_col=0, sep='\t')
abundDF = pd.DataFrame(columns=['ORF', 'exp_young', 'exp_middle', 'exp_old', 'category', 'complex'])
for i in range(len(cele_transl)):
    current_orf = cele_transl.iloc[i]['Name']
    current_category = cele_transl.iloc[i]['category']
    current_complex = cele_transl.iloc[i]['complex']
    if current_orf in cele_abund.index:
        current_abundance = np.around(np.array(cele_abund.loc[current_orf]), 3)
        abundDF.loc[len(abundDF)] = (current_orf, current_abundance[0], current_abundance[1], current_abundance[2], current_category, current_complex)
abundDF.to_csv("../data/processed/cele_translation_abundance.txt", header=True, index=False, sep='\t')



#--------------------------------------------------------------------------------------------
def arrival_probs(ACdf, tRNAdf):
    """
    Compute tRNA arrival probabilities for all codons. 
    Codons without cognate tRNA yield 0. 
    ACdf: codon-anticodon-AA dataframe
    tRNAdf: input tRNA counts and abundances
    """

    # constants
    l = 1.5e-8					# effective length [ref]
    D = 8.42e-11				# diffusion coefficient [ref]
    V = 4.2e-17					# cell volume [ref]

    # diffusion
    tau = l**2 / (6 * D)		# transition time between locations [ref]
    n_loc = V / (l**3)			# approx number of cell locations [ref]

    lmbda = np.zeros(( len(ACdf) )) * np.nan
    for cx, codon in enumerate(list(ACdf['codon'])):
        current_anticodon = list(ACdf[ACdf['codon']==codon]['anticodon'])[0]
        if current_anticodon in list(tRNAdf['Anticodon']):
            current_tRNA_count = list(tRNAdf[tRNAdf['Anticodon']==current_anticodon]['Count'])[0]
            current_tRNA_abundance = list(tRNAdf[tRNAdf['Anticodon']==current_anticodon]['Abundance'])[0]
        else:
            current_tRNA_count = 0
            current_tRNA_abundance = 0       
        lmbda[cx] = (current_tRNA_abundance / n_loc) / tau
        #lmbda[cx] = (current_tRNA_count / n_loc) / tau

    pa = lmbda / np.sum(lmbda)

    pa_n = np.zeros(( len(pa) ))
    pa_n[pa>0] = 1 - pa[pa>0]

    resultDF = ACdf.copy()
    resultDF['pa'] = pa
    resultDF['pa_n'] = pa_n
    resultDF['lambda'] = lmbda

    return resultDF


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def compute_pb(pos_weights, nt_energies):
 
    x = list(it.product((0,1,2,3),(0,1,2,3),(0,1,2,3)))
    y = np.array(list(it.product(x,x)))
    z = ["".join(a) for a in list(np.array(nt)[np.array(x)]) ]

    ENG = np.reshape( np.sum( pos_weights * nt_energies[y[:,0,:], y[:,1,:][:,::-1] ], axis=1), (64, 64) )
  
    delta_A = np.exp(-ENG)  # stability of codon/anticodon interaction
    delta_A = np.transpose(delta_A) # transpose 2x to normlize along axis=1
    delta_A /= np.sum(delta_A, axis=0)
    pbMAT = np.transpose(delta_A)

    E2 = pd.DataFrame(data=pbMAT, columns=z, index=z)

    return E2




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def compute_roc(Emat):
    """
    bla
    """
    rocDF = pd.DataFrame(columns=['pb', 'label'])
    for i in list(anticodonDF['codon']):
        current_codon_aa = list(anticodonDF[anticodonDF['codon']==i]['AA'])[0]
        for j in list(trnaDF['Anticodon']):
    
            if j in list(anticodonDF['anticodon']):
                current_ac_aa = list(anticodonDF[anticodonDF['anticodon']==j]['AA'])[0]
                current_ac_codon = list(anticodonDF[anticodonDF['anticodon']==j]['codon'])[0]
                current_energy = Emat.loc[i][j]
                #print(i, j, current_codon_aa, current_ac_aa, current_energy)
                if current_codon_aa == current_ac_aa:
                    rocDF.loc[len(rocDF)] = (current_energy, 1)
                else:
                    rocDF.loc[len(rocDF)] = (current_energy, 0)

    auc = roc_features(np.array(rocDF['pb']), np.array(rocDF['label']))

    return rocDF, auc



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def roc_features(x, y):
    """
    get ROC AUC for select input features
    sum of feature counts as simple model of total number of putative binding sites
    """

    auc = roc_auc_score(y, x)

    return np.around(auc, 4)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def extract_rocDF(preOUT, anticodonDF, trnaDF, N=20):
    """
    Compute ROC for anticodon-codon binding probabilities from model parameters
    """

    #load model parameters
    E = pickle.load(open('../data/processed/'+preOUT+'_E.pkl', 'rb'))
    P = np.loadtxt('../data/processed/'+preOUT+'_P.txt')
    CC = np.loadtxt('../data/processed/'+preOUT+'_CC.txt')

    # sort by model performance
    CC_order = np.argsort(CC)
    CC_order = CC_order[::-1]

    idx_top = CC_order[0:N]
    resmat = np.zeros(( N, 16 )) * np.nan
    list_df = []
    for ix, i in enumerate(idx_top):			# [0:1]
        resmat[ix,:] =  np.reshape(E[:,:,i], (1,np.prod(np.shape(E[:,:,i]))) ) 

        current_P = P[:,i]
        current_E = E[:,:,i]

        E2 = compute_pb(current_P, current_E)

        rocDF, auc = compute_roc(E2)
        rocDF.columns = ['pb_'+str(ix), 'label_'+str(ix)]
        list_df.append(rocDF)
 
    rocDF = pd.concat(list_df, axis=1)
    sel_columns = ['label_0'] + ['pb_'+str(x) for x in range(N)]
    rocDF = rocDF[sel_columns]
    rocDF.to_csv('../data/processed/'+preOUT+'_ROC.txt', header=True, index=False, sep='\t')


    labels = np.array([['AA', 'AC', 'AG', 'AT'],
				   ['CA', 'CC', 'CG', 'CT'],
	               ['GA', 'GC', 'GG', 'GT'],
	               ['TA', 'TC', 'TG', 'TT']])
    labels = np.reshape(labels, (1,np.prod(np.shape(labels)))).flatten()
    engDF = pd.DataFrame(data=np.around(resmat, 5), columns=labels)
    engDF.to_csv('../data/processed/'+preOUT+'_E20.txt', header=True, index=False, sep='\t')

    P2 = P[:,idx_top]
    np.savetxt('../data/processed/'+preOUT+'_P20.txt', P2)


    return rocDF



anticodonDF = pd.read_csv("../data/tRNA/anticodon2.txt", header=0, index_col=None, sep='\t')
trnaDF = pd.read_csv('../data/tRNA/yeast_tRNA_count.csv', delimiter=",", header=0)

rocDF = extract_rocDF('scer2_basic_young', anticodonDF, trnaDF, N=20)
rocDF = extract_rocDF('scer2_basic_middle', anticodonDF, trnaDF, N=20)
rocDF = extract_rocDF('scer2_basic_old', anticodonDF, trnaDF, N=20)

rocDF = extract_rocDF('scer2_cu_young', anticodonDF, trnaDF, N=20)
rocDF = extract_rocDF('scer2_cu_middle', anticodonDF, trnaDF, N=20)
rocDF = extract_rocDF('scer2_cu_old', anticodonDF, trnaDF, N=20)


anticodonDF = pd.read_csv("../data/tRNA/anticodon2.txt", header=0, index_col=None, sep='\t')
trnaDF = pd.read_csv('../data/tRNA/cele_tRNA_count.txt', delimiter=" ", header=0)

rocDF = extract_rocDF('cele2_basic_young', anticodonDF, trnaDF, N=20)
rocDF = extract_rocDF('cele2_basic_middle', anticodonDF, trnaDF, N=20)
rocDF = extract_rocDF('cele2_basic_old', anticodonDF, trnaDF, N=20)

rocDF = extract_rocDF('cele2_cu_young', anticodonDF, trnaDF, N=20)
rocDF = extract_rocDF('cele2_cu_middle', anticodonDF, trnaDF, N=20)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
codons_nonstop = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT', 
    		  'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT',
    		  'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT',
    		  'TAC', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT',]

x = np.ones(( 20 )) * 0.4
s = np.ones(( 41 )) * 0.6




def mean_dwelltime(x, s, cu):
    """
    fit tRNA wobble coefficient x and competition factor s from riboseq average RDs
    indexing: one file is by sorted codons, one by sorted amino acids
    abundances are experimental yeast tRNA abundances
    wobble parameters x < 0.5, i.e. the major codon has to have a weight > 0.5
    """

    tau = 4.45e-7
    n = 1.24e7

    result = np.zeros(( 61 ))

    result[39] = (tau * n)/(  212721 * (1-x[0]) * s[0])				# GCT 
    result[37] = (tau * n)/(  212721 * x[0] * s[0])					# GCC 
    result[36] = (tau * n)/(  384682 * (1-x[1]) * s[1])				# GCA 
    result[38] = (tau * n)/(  384682 * x[1] * s[1])  				# GCG 
    result[27] = (tau * n)/(  428095 * (1-x[2]-x[3])* s[2])			# CGT 
    result[25] = (tau * n)/(  428095 * x[2] * s[2])					# CGC 
    result[24] = (tau * n)/(  428095 * x[3] * s[2])					# CGA 
    result[26] = (tau * n)/(  60792 * s[3])							# CGG 
    result[8]  = (tau * n)/(  683001 * s[4])						# AGA 
    result[10] = (tau * n)/(  106890 * s[5])						# AGG 
    result[3]  = (tau * n)/(  378101 * x[4] * s[6])					# AAT 
    result[1]  = (tau * n)/(  378101 * (1-x[4]) * s[6]) 			# AAC 
    result[35] = (tau * n)/(  1789916 * x[5] * s[7])				# GAT 
    result[33] = (tau * n)/(  1789916 * (1-x[5]) * s[7])			# GAC
    result[56] = (tau * n)/(  106137 * x[6] * s[8])					# TGT 
    result[54] = (tau * n)/(  106137 * (1-x[6]) * s[8])				# TGC 
    result[32] = (tau * n)/(  1618235 * s[9])						# GAA 
    result[34] = (tau * n)/(  393379 * s[10])						# GAG 
    result[16] = (tau * n)/(  436225 * s[11])						# CAA 
    result[18] = (tau * n)/(  89093 * s[12])						# CAG 
    result[43] = (tau * n)/(  826342 * x[7] * s[13])				# GGT 
    result[41] = (tau * n)/(  826342 * (1-x[7]) * s[13])			# GGC 
    result[40] = (tau * n)/(  287325 * s[14])						# GGA 
    result[42] = (tau * n)/(  129338 * s[15])						# GGG 
    result[19] = (tau * n)/(  727396 * x[8] * s[16])				# CAT 
    result[17] = (tau * n)/(  727396 * (1-x[8]) * s[16])			# CAC 
    result[15] = (tau * n)/(  916236 * (1-x[9]) * s[17])			# ATT 
    result[13] = (tau * n)/(  916236 * x[9] * s[17])				# ATC 
    result[12] = (tau * n)/(  54352 * s[18])						# ATA 
    result[57] = (tau * n)/(  271592 * s[19])						# TTA 
    result[59] = (tau * n)/(  709768 * s[20])						# TTG 
    result[31] = (tau * n)/(  63503 * x[10] * s[21])				# CTT 
    result[29] = (tau * n)/(  63503 * (1-x[10]) * s[21])			# CTC 
    result[28] = (tau * n)/(  94608 * (1-x[11]) * s[22])			# CTA 
    result[30] = (tau * n)/(  94608 * x[11] * s[22])				# CTG 
    result[0]  = (tau * n)/(  222273 * s[23])						# AAA 
    result[2]  = (tau * n)/(  397111 * s[24])						# AAG 
    result[14] = (tau * n)/(  266917 * s[25])						# ATG 
    result[60] = (tau * n)/(  104391 * x[12] * s[26])				# TTT 
    result[58] = (tau * n)/(  104391 * (1-x[12]) * s[26])			# TTC 
    result[23] = (tau * n)/(  57720 * (1-x[13]) * s[27])			# CCT 
    result[21] = (tau * n)/(  57720 * x[13] * s[27])				# CCC 
    result[20] = (tau * n)/(  286531 * (1-x[14]) * s[28])			# CCA 
    result[22] = (tau * n)/(  286531 * x[14] * s[28])				# CCG 
    result[53] = (tau * n)/(  792870 * (1-x[15]) * s[29])			# TCT 
    result[51] = (tau * n)/(  792870 * x[15] * s[29])				# TCC 
    result[50] = (tau * n)/(  86913 * s[30])						# TCA 
    result[52] = (tau * n)/(  34793 * s[31])						# TCG 
    result[11] = (tau * n)/(  249948 * x[16] * s[32])				# AGT 
    result[9]  = (tau * n)/(  249948 * (1-x[16]) * s[32])			# AGC 
    result[7]  = (tau * n)/(  381835 * (1-x[17]) * s[33])			# ACT 
    result[5]  = (tau * n)/(  381835 * x[17] * s[33])				# ACC 
    result[4]  = (tau * n)/(  105862 * s[34])						# ACA 
    result[6]  = (tau * n)/(  113104 * s[35])						# ACG 
    result[55] = (tau * n)/(  298500 * s[36])						# TGG 
    result[49] = (tau * n)/(  343867 * x[18] * s[37])				# TAT 
    result[48] = (tau * n)/(  343867 * (1-x[18]) * s[37])			# TAC 
    result[47] = (tau * n)/(  806173 * (1-x[19]) * s[38])			# GTT 
    result[45] = (tau * n)/(  806173 * x[19] * s[38])				# GTC 
    result[44] = (tau * n)/(  145923 * s[39])						# GTA 
    result[46] = (tau * n)/(  118162 * s[40])						# GTG 


    return result



xx = np.ones(( 20 )) * 0.1
dt = mean_dwelltime(xx, s, scer_codonusage['young'])
print(dt)

optDF  = pd.read_csv("../data/processed/optimality_cu.txt", header=0, index_col=None, sep='\t')
 
ribodensities = pd.read_csv("../data/processed/scer_ribodensities.txt", header=0, index_col=None, sep='\t')
ribodensities.index = list(ribodensities['codon'])
ribodensities = ribodensities.loc[codons_nonstop]

RD = ribodensities.sort_values(by=['codon'])
true = RD['RD_young']

best_cc = np.corrcoef(dt, true)[0,1]
print(best_cc)

for i in range(1000000):
    idx = np.random.choice(np.arange(20))
    proposal_xx = np.copy(xx)
    proposal_xx[idx] = 0.5 * np.random.rand()

    dt = mean_dwelltime(proposal_xx, s, scer_codonusage['young'])
    cc = np.corrcoef(dt, true)[0,1]

    if cc > best_cc:
        best_cc = np.copy(cc)
        print(best_cc)
        xx = np.copy(proposal_xx)
        np.savetxt("tmp.xx", xx)


dt = mean_dwelltime(xx, s, scer_codonusage['young'])
print(xx)
print(s)
print(dt)
np.savetxt("tmp.dt", dt)
