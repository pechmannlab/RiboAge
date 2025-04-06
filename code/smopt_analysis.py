import os, sys
import numpy as np
import pandas as pd
import pickle
import subprocess


# Assorted functions and analyses to extract information from SMoPT simulation output
# SP (2025)


codons_nonstop = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT', 
    		  'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT',
    		  'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT',
    		  'TAC', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT',]
   
sequence_codons = pickle.load(open('../data/processed/scer/scer_scikit_codonseq.pkl', 'rb'))




#--------------------------------------------------------------------------------------------
def parse_smopt_profiles(preOUT):
    """
    Parse output data from SMoPT simulations
    """

    dwelltimes = pd.read_csv("../data/SMoPT/smopt_"+preOUT+"_dwelltimes.txt", header=0, index_col=False, sep='\t')
    collisions = pd.read_csv("../data/SMoPT/smopt_"+preOUT+"_collisions.txt", header=0, index_col=False, sep='\t')

    c_dict = {}
    p_dict = {}
    t_dict = {}

    for i in sorted(set(collisions['ORF'])):
        L = len(sequence_codons[i])
        current_profile = np.zeros(( L )) * np.nan
        current_collision = np.zeros(( L )) * np.nan
        current_timeprof = np.zeros(( L )) * np.nan

        current_df = collisions[collisions['ORF']==i]
        for j in range(len(current_df)):
            current_position = current_df.iloc[j]['position']-1
            if current_position >=0 and current_position < L:
                current_collision[current_position] += 1

        current_dt = dwelltimes[dwelltimes['ORF']==i]
        for j in range(len(current_dt)):
            current_pos = current_dt.iloc[j]['position']-1
            if current_pos >= 0 and current_pos < L:
                current_profile[current_pos] += 1

        for j in range(L):
            current_dwell = np.array(current_dt[current_dt['position']==(j+1)]['dwelltime'])
            if len(current_dwell) > 0:
                current_timeprof[j] = np.nanmean(current_dwell)
        #print(current_timeprof)

        print(i)
        c_dict[i] = current_collision
        p_dict[i] = current_profile
        t_dict[i] = current_timeprof

    pickle.dump(c_dict, open('../data/SMoPT/profiles_'+preOUT+'_collisions.pkl', 'wb'))
    pickle.dump(p_dict, open('../data/SMoPT/profiles_'+preOUT+'_elongation.pkl', 'wb'))
    pickle.dump(t_dict, open('../data/SMoPT/profiles_'+preOUT+'_dwelltimes.pkl', 'wb'))



parse_smopt_profiles('young2')
parse_smopt_profiles('middle2')
parse_smopt_profiles('old2')






#--------------------------------------------------------------------------------------------
# computes dwell times from simulation output
# adds synthetic collision noise (drawn from expdist) 


list_arrays = []
list_seq = []

dwelltimes_young = pickle.load(open('../data/SMoPT/profiles_young2_dwelltimes.pkl', 'rb'))
dwelltimes_middle = pickle.load(open('../data/SMoPT/profiles_middle2_dwelltimes.pkl', 'rb'))
dwelltimes_old = pickle.load(open('../data/SMoPT/profiles_old2_dwelltimes.pkl', 'rb'))

scer_codonusage = pd.read_csv("../data/processed/scer_codonusage.txt", header=0, index_col=0, sep='\t')


list_allorfs = sorted(set(list(dwelltimes_young.keys())+list(dwelltimes_middle.keys())+list(dwelltimes_old.keys())))

for i in list_allorfs:
    if i in list(dwelltimes_young.keys()) and i in list(dwelltimes_middle.keys()) and i in list(dwelltimes_old.keys()):
        current_seq = sequence_codons[i]
        L = len(current_seq)
        current_young = np.array(dwelltimes_young[i])
        current_middle = np.array(dwelltimes_middle[i])
        current_old = np.array(dwelltimes_old[i])
        current_dwell = np.transpose(np.array([current_young, current_middle, current_old]))

        list_arrays.append(current_dwell)
        list_seq.append(current_seq)
        #print(i)

all_arrays = np.concatenate(list_arrays)
print(np.shape(all_arrays))

np.savetxt("../data/SMoPT/simRD_all.txt", all_arrays)
all_seq = np.concatenate(list_seq)
dwellt = np.zeros(( 61, 3 ))


for jx, j in enumerate(codons_nonstop):
    idx = np.array(all_seq) == j
    current_data = all_arrays[idx,:]
    dwellt[jx, :] = np.nanmedian(current_data, axis=0)

dwellDF = pd.DataFrame(data=dwellt, columns=['young', 'middle', 'old'])
dwellDF.to_csv("../data/SMoPT/codon_dwelltimes.txt", header=True, index=False, sep='\t')


dwellr = np.zeros(( 61 ))
dwell_hypo = np.zeros(( 61, 50 ))
data_normal = all_arrays[:, 0]
all_sample = np.zeros(( len(data_normal), 50))
ICdf = pd.DataFrame(columns=['pct', 'IC'])
for t in range(2, 50):
    for jx, j in enumerate(codons_nonstop):
        codon_coefficient = scer_codonusage.loc[j]['old']/scer_codonusage.loc[j]['young']
        idx = np.array(all_seq) == j
        current_data = data_normal[idx]
       
        #codon_coefficient = 1		# without codon-usage factor

        rr = np.zeros((100))		# randomization actually not needed as very stable
        for r in range(100):
            current_sample = np.copy(current_data)
            N_rand = int(codon_coefficient*(t/100)*len(current_data))
            idx_rand = np.random.choice(np.arange(len(current_data)), N_rand, replace=False )

            current_sample[idx_rand] = np.random.exponential( scale=1, size=len(idx_rand) ) 

            rr[r] = np.median(current_sample) 
        dwell_hypo[jx, t] = np.around(np.mean(rr), 4) #np.nanmedian(current_data) 
        all_sample[idx,t] = current_sample

    breaks = np.arange(5000)/100
    histo_sample = all_sample[:,t]
    histo_sample = histo_sample[histo_sample>0]		# no zeros (shouldn't be there anyways, unless 'rounding errors')
    h, hx = np.histogram(histo_sample, bins=breaks, density=True)
    h = h/100

    S = 0
    for jx, j in enumerate(codons_nonstop):
        current_dwell = dwell_hypo[jx, t]
        current_x = np.where(current_dwell > breaks)[0][-1]
        current_p = h[current_x]

        S -= current_p * np.log(current_p)

    print(t, S)
    ICdf.loc[len(ICdf)] = (t, S)

  

ICdf.to_csv("../data/SMoPT/IC_hypo.txt", header=True, index=False, sep='\t')
np.savetxt('../data/SMoPT/dwell_hypo.txt', dwell_hypo)

#ICdf.to_csv("../data/SMoPT/IC_hypo_CU1.txt", header=True, index=False, sep='\t')
#np.savetxt('../data/SMoPT/dwell_hypo_CU1.txt', dwell_hypo)








#--------------------------------------------------------------------
def pval_adjust(pval):
    """
    correct for multiple testing
    Benjamini & Hochberg FDR method (R funct p.adjust)
    """
    pval = np.array(pval)
    if np.size(pval) > 1:
        padj = np.nan * np.ones( (np.size(pval) ) )
        nn = np.isnan(pval)
        pval = pval[~nn]
        n = len(pval)
        i = np.array( range(n)[::-1]) + 1
        o = np.array( sorted(range(n), key=lambda k: pval[k])[::-1] )
        ro = np.array( sorted(o, key=lambda k: o[k]) )
        adj = np.minimum.accumulate( float(n)/np.array(i) * pval[o] )
        for i in range(len(adj)):
            adj[i] = min(1, adj[i])
        padj[~nn] = adj[ro]
    else:
        padj = pval

    return padj



from scipy.stats import fisher_exact
#oddsr, p = fisher_exact(table, alternative='two-sided')

young_collisions = pickle.load(open('../data/SMoPT/profiles_young2_collisions.pkl' ,'rb'))
young_elongation = pickle.load(open('../data/SMoPT/profiles_young2_elongation.pkl', 'rb'))
young_dwelltimes = pickle.load(open('../data/SMoPT/profiles_young2_dwelltimes.pkl', 'rb'))

old_collisions = pickle.load(open('../data/SMoPT/profiles_old2_collisions.pkl' ,'rb'))
old_elongation = pickle.load(open('../data/SMoPT/profiles_old2_elongation.pkl', 'rb'))
old_dwelltimes = pickle.load(open('../data/SMoPT/profiles_old2_dwelltimes.pkl', 'rb'))


YAL003Wdf = pd.DataFrame({'collision_young': young_collisions['YAL033W'],
						  'elongation_young': young_elongation['YAL033W'],
						  'collision_old':old_collisions['YAL033W'],
						  'elongation_old':old_elongation['YAL033W'],
						  'dwelltime_young': young_dwelltimes['YAL033W'],
						  'dwelltime_old':  old_dwelltimes['YAL033W']})

YAL003Wdf.to_csv("../data/SMoPT/YAL003W.txt", header=True, index=False, sep='\t')




YBR170Cdf = pd.DataFrame({'collision_young': young_collisions['YBR170C'],
						  'elongation_young': young_elongation['YBR170C'],
						  'collision_old':old_collisions['YBR170C'],
						  'elongation_old':old_elongation['YBR170C'],
						  'dwelltime_young': young_dwelltimes['YBR170C'],
						  'dwelltime_old':  old_dwelltimes['YBR170C']})

YBR170Cdf.to_csv("../data/SMoPT/YBR170C.txt", header=True, index=False, sep='\t')




def sites_collision(collprofile):
    """
    upstream codons linked to ribosome collisions
    """

    theta_coll = 0          # min height of collision peak
    theta_stretch = 1		# gap between two collision regions 

    coll_idx = np.where(collprofile > theta_coll)[0]
    #coll_ends = np.where( (coll_idx[1:] - coll_idx[0:(len(coll_idx)-1)]) > theta_stretch )[0] + 1

    return coll_idx


# extracts CU at collision sites and upstream codons
cy = 0
ey = 0
co = 0
eo = 0
list_arrays = []
list_seq = []
list_table = []
resDF = pd.DataFrame(columns=['ORF', 'OR', 'pval'])
counts_young_site = np.zeros((61))
counts_young_upstr = np.zeros((61))
counts_old_site = np.zeros((61))
counts_old_upstr = np.zeros((61))
counts_young_rand = np.zeros(( 61, 1000 ))
counts_old_rand = np.zeros(( 61, 1000))
collDF = pd.DataFrame(columns=['n_young', 'n_old', 'L'])
for i in list(young_collisions.keys()):
    if i in list(old_collisions.keys()):

        current_seq = sequence_codons[i]
        L = len(current_seq)

        current_c_y = young_collisions[i]
        current_e_y = young_elongation[i]
        current_c_o = old_collisions[i]
        current_e_o = old_elongation[i]

        cy += np.sum(current_c_y)
        ey += np.sum(current_e_y)
        co += np.sum(current_c_o)
        eo += np.sum(current_e_o)

        current_table = np.array([[np.sum(current_c_o), np.sum(current_e_o)],[np.sum(current_c_y), np.sum(current_e_y)]])
        OR, pval = fisher_exact(current_table, alternative='two-sided')
        resDF.loc[len(resDF)] = (i, OR, pval)

        coll_young = np.zeros(( L ))
        coll_young[current_e_y > 0] = current_c_y[current_e_y > 0] / current_e_y[current_e_y > 0]

        coll_old = np.zeros(( L ))
        coll_old[current_e_o > 0] = current_c_o[current_e_o > 0] / current_e_o[current_e_o > 0]

        current_diff = coll_old - coll_young
        list_arrays.append(current_diff)
        list_seq.append(current_seq)

        gene_table = np.zeros(( L, 7 )) * np.nan
        for j in range(L):
            gene_table[j,0:4] = np.array([ current_c_o[j], current_e_o[j], current_c_y[j], current_e_y[j] ])
            
            if current_e_y[j] > 0 and current_e_o[j] > 0:		# has to be coverage
                if current_c_y[j] > 0 or current_c_o[j] > 0:
                    gene_table[j,4], gene_table[j,5] = fisher_exact(np.reshape(gene_table[j,0:4], (2,2)), alternative='two-sided')
        
        if np.sum(np.isnan(gene_table[:,5])) < np.shape(gene_table)[0]: 
            gene_table[:,6] = pval_adjust(gene_table[:,5])
        list_table.append(gene_table)

        coll_young_pos = sites_collision(coll_young)
        for j in coll_young_pos:
            if j + 10 < L:
                coll_site = current_seq[j]
                coll_upstr = current_seq[j+10]		# +10
                counts_young_site[codons_nonstop.index(coll_site)] += 1
                counts_young_upstr[codons_nonstop.index(coll_upstr)] += 1

        coll_old_pos = sites_collision(coll_old)
        for j in coll_old_pos:
            if j + 10 < L:
                coll_site = current_seq[j]
                coll_upstr = current_seq[j+10]		# +10
                counts_old_site[codons_nonstop.index(coll_site)] += 1
                counts_old_upstr[codons_nonstop.index(coll_upstr)] += 1

        collDF.loc[len(collDF)] = (len(coll_young_pos), len(coll_old_pos), L)

        for r in range(1000):
            rand_young_pos = np.random.choice( np.arange(L-1), len(coll_young_pos), replace=False)
            rand_old_pos = np.random.choice( np.arange(L-1), len(coll_old_pos), replace=False)
            for j in rand_young_pos:
                rand_site = current_seq[j]
                counts_young_rand[codons_nonstop.index(rand_site), r] += 1
            for j in rand_old_pos:
                rand_site = current_seq[j]
                counts_old_rand[codons_nonstop.index(rand_site), r] += 1


resDF['p_adj'] = pval_adjust(resDF['pval'])
resDF.to_csv("../data/SMoPT/collisions_pergene.txt", header=True, index=False, sep='\t')


c_sites = np.sum(np.array(collDF), axis=0)
print(c_sites)
c_sites = np.array([ c_sites[0]/c_sites[2], c_sites[1]/c_sites[2] ]) * 100
c_colls = np.array([ cy/ey, co/eo ]) * 100
print(cy, ey, co, eo)
collDF = pd.DataFrame({'age': ['young', 'old'], 'pct_collisions': list(c_colls), 'pct_sites': list(c_sites)}) 
collDF.to_csv("../data/SMoPT/collisions_counts.txt", header=True, index=False, sep='\t')


ccDF = pd.DataFrame({'codon': codons_nonstop, 'young_site': counts_young_site, 'young_upstr': counts_young_upstr, 
					 'old_site': counts_old_site, 'old_upstr': counts_old_upstr})
ccDF.to_csv("../data/SMoPT/collisions_codoncounts.txt", header=True, index=False, sep='\t')


np.savetxt("../data/SMoPT/rand_young.txt", counts_young_rand)
np.savetxt("../data/SMoPT/rand_old.txt", counts_old_rand)


all_arrays = np.concatenate(list_arrays)
all_seq = np.concatenate(list_seq)
all_tables = np.concatenate(list_table)

diff = np.zeros(( 61 ))
aaDF = pd.DataFrame(columns=['Codon', 'CO', 'EO', 'CY', 'EY', 'OR', 'pval'])
for jx, j in enumerate(codons_nonstop):
    idx = np.array(all_seq) == j
    score = all_arrays[idx]
    diff[jx] = np.nanmedian(score[score!=0])

    current_counts = np.array(np.sum(all_tables[idx,:][:,0:4], axis=0), dtype=int)
    current_OR, current_p = fisher_exact(np.reshape(current_counts, (2,2)), alternative='two-sided')

    aaDF.loc[len(aaDF)] = (j, current_counts[0], current_counts[1], current_counts[2], current_counts[3], current_OR, current_p)

aaDF['p_adj'] = pval_adjust(aaDF['pval'])
np.savetxt('tmp.diff', diff)
aaDF.to_csv('tmp.aadf', header=True, index=False, sep='\t')


