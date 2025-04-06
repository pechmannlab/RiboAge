import os, sys
import numpy as np
import pandas as pd
import itertools as it
import glob
import json
import csv
import pickle
import gzip as gz
import subprocess




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def filter_normalize_data(speciesIN, threshold_ribo, threshold_mm):
    """
    final data QC
    min length: 100aa
    clip 7aa at beginning and end for biases, equiv to ~20nt of riboseq biases
    min coverage: 1
    """

    ribo_dict = pickle.load( open("../data/processed/"+speciesIN+"/"+speciesIN+"_scikit_mat_pre.pkl", 'rb') )
    filterDF = pd.read_csv("../data/processed/"+speciesIN+"/"+speciesIN+"_filter.txt", header=0, index_col=False, sep='\t')
    anno = pd.read_csv("../data/processed/"+speciesIN+"/"+speciesIN+"_anno.txt", header=0, index_col=False, sep='\t')

    list_age = sorted(list(set(anno['AGE'])))
    idx_young = np.where(anno['AGE']==list_age[0])[0]
    idx_middle = np.where(anno['AGE']==list_age[1])[0]
    idx_old = np.where(anno['AGE']==list_age[2])[0]

    final_ribo = {}
    counter = 0
    for i in ribo_dict.keys():
    	if i in list(filterDF['ORF']):
            current_ribo = ribo_dict[i]     
            current_el = list(filterDF[filterDF['ORF']==i]['EL'])[0]
            current_mm = list(filterDF[filterDF['ORF']==i]['MM'])[0]
            counter += 1

            avrg_all = np.mean(current_ribo, axis=1)
            avrg_age = np.array([np.mean(np.sum(current_ribo[idx_young,:], 0)), np.mean(np.sum(current_ribo[idx_middle,:], 0)), np.mean(np.sum(current_ribo[idx_old,:], 0))])
            avrg_null = np.sum(current_ribo == 0, axis=1) < np.shape(current_ribo)[1]/2.

            if np.all(avrg_age >= 1) and current_mm <= threshold_mm:
                for j in range(current_ribo.shape[0]):
                    current_average = np.nanmean(current_ribo[j,:])
                    if current_average >= threshold_ribo:				# min coverage
                        current_ribo[j,:] /= current_average			# normalize to 1 mean
                    else:
                        current_ribo[j,:] = np.nan

                final_ribo[i] = current_ribo

    print(len(final_ribo))

    # save output
    pickle.dump(final_ribo, open('../data/processed/'+speciesIN+'/'+speciesIN+'_scikit_mat.pkl', 'wb'))
   


#--------------------------------------------------------------------------------------------
def load_data(speciesIN):
    """
    load required base data
    """

    scikit_data = pickle.load( open("../data/processed/"+speciesIN+"/"+speciesIN+"_scikit_mat.pkl", 'rb') )
    mm_consensus = pickle.load( open("../data/processed/"+speciesIN+"/"+speciesIN+"_mm_consensus.pkl", 'rb') )

    # there seems to be a problem with including/excluding stop codons in pysam with different gtfs/settings required for scer vs. cele
    dt2 = {}
    mm2 = {}
    for i in scikit_data.keys():
        current_data = scikit_data[i]
        if i in list(mm_consensus.keys()):
            current_mm = (np.array(mm_consensus[i], dtype=bool) == False) 
            if np.shape(current_data)[1] == (len(current_mm)):
                dt2[i] = current_data
                mm2[i] = current_mm
            elif np.shape(current_data)[1] == (len(current_mm) -1) :
                current_mm = current_mm[0:np.shape(current_data)[1]]
                dt2[i] = current_data
                mm2[i] = current_mm
    
    return dt2, mm2 		#scikit_data, mm_consensus



#--------------------------------------------------------------------------------------------
def correlation_check(data, mm, annotation):
    """
    Computes correlation between replicas within and between input groups (ages)
    data: dictionary of data matrices
    annotation: rows in matrices for each gene
    """
    list_age = sorted(list(set(annotation['AGE'])))

    idx_young = np.where(annotation['AGE']==list_age[0])[0]
    ndx_young = np.where(annotation['AGE']!=list_age[0])[0]
    idx_middle = np.where(annotation['AGE']==list_age[1])[0]
    ndx_middle = np.where(annotation['AGE']!=list_age[1])[0]
    idx_old = np.where(annotation['AGE']==list_age[2])[0]
    ndx_old = np.where(annotation['AGE']!=list_age[2])[0]
  
    resultDF = pd.DataFrame(columns=['self', 'non', 'young', 'nyoung', 'middle', 'nmiddle', 'old', 'nold', 'young_middle', 'young_old', 'middle_old'])

    for gene in data.keys():

        current_data = data[gene]
        current_mm = mm[gene]
        current_data = current_data[:,current_mm]

        a_self = []
        a_non  = []

        # within young vs young to not young
        l_young = []
        for i in it.combinations(idx_young,2):
            current_corr = np.corrcoef(current_data[i[0],:], current_data[i[1],:])[0,1]
            l_young.append(current_corr)
            a_self.append(current_corr)
        if np.any(np.isnan(np.array(l_young))==False):
            corr_young = np.around(np.nanmean(np.array(l_young)), 3)
        else:
            corr_young = np.nan

        n_young = []
        for i in idx_young:
            for j in ndx_young:
                current_corr = np.corrcoef(current_data[i,:], current_data[j,:])[0,1]
                n_young.append(current_corr)
                a_non.append(current_corr)
        if np.any(np.isnan(np.array(n_young))==False):
            corr_notyoung = np.around(np.nanmean(np.array(n_young)), 3)
        else:
            corr_notyoung = np.nan

        # within middle & not middle
        l_middle = []
        for i in it.combinations(idx_middle,2):
            current_corr = np.corrcoef(current_data[i[0],:], current_data[i[1],:])[0,1]
            l_middle.append(current_corr)
            a_self.append(current_corr)
        if np.any(np.isnan(np.array(l_middle))==False):
            corr_middle = np.around(np.nanmean(np.array(l_middle)), 3)
        else:
            corr_middle = np.nan

        n_middle = []
        for i in idx_middle:
            for j in ndx_middle:
                current_corr = np.corrcoef(current_data[i,:], current_data[j,:])[0,1]
                n_middle.append(current_corr)
                a_non.append(current_corr)
        if np.any(np.isnan(np.array(n_middle))==False):
            corr_notmiddle = np.around(np.nanmean(np.array(n_middle)), 3)
        else:
            corr_notmiddle = np.nan

        # within old & not old
        l_old = []
        for i in it.combinations(idx_old,2):
            current_corr = np.corrcoef(current_data[i[0],:], current_data[i[1],:])[0,1]
            l_old.append(current_corr)
            a_self.append(current_corr)
            #if np.isnan(current_corr) == False:
            #    np.savetxt("TMP/cele_middle_"+gene, current_data[[i[0], i[1]],:])                
        if np.any(np.isnan(np.array(l_old))==False):
            corr_old = np.around(np.nanmean(np.array(l_old)), 3)
        else:
            corr_old = np.nan

        n_old = []
        for i in idx_old:
            for j in ndx_old:
                current_corr = np.corrcoef(current_data[i,:], current_data[j,:])[0,1]
                n_old.append(current_corr)
                a_non.append(current_corr)
        if np.any(np.isnan(np.array(n_old))==False):
            corr_notold = np.around(np.nanmean(np.array(n_old)), 3)
        else:
            corr_notold = np.nan

        if np.any(np.isnan(np.array(a_self))==False):
            corr_self = np.around(np.nanmean(np.array(a_self)), 3)
        else:
            corr_self = np.nan

        if np.any(np.isnan(np.array(a_non))==False):
            corr_non  = np.around(np.nanmean(np.array(a_non)), 3)
        else:
            corr_non = np.nan

        l_young_middle = []
        for i in idx_young:
            for j in idx_middle:
                current_corr = np.corrcoef(current_data[i,:], current_data[j,:])[0,1]
                l_young_middle.append(current_corr)
        if np.any(np.isnan(np.array(l_young_middle))==False):
            corr_young_middle = np.around(np.nanmean(np.array(l_young_middle)), 3)
        else:
            corr_young_middle = np.nan

        l_young_old = []
        for i in idx_young:
            for j in idx_old:
                current_corr = np.corrcoef(current_data[i,:], current_data[j,:])[0,1]
                l_young_old.append(current_corr)
        if np.any(np.isnan(np.array(l_young_old))==False):
            corr_young_old = np.around(np.nanmean(np.array(l_young_old)), 3)
        else:
            corr_young_old = np.nan

        l_middle_old = []
        for i in idx_middle:
            for j in idx_old:
                    current_corr = np.corrcoef(current_data[i,:], current_data[j,:])[0,1]
                    l_middle_old.append(current_corr)
        if np.any(np.isnan(np.array(l_middle_old))==False):
            corr_middle_old = np.around(np.nanmean(np.array(l_middle_old)), 3)
        else:
            corr_middle_old

        chck = np.array([corr_self, corr_non, corr_young, corr_notyoung, corr_middle, corr_notmiddle, corr_old, corr_notold, corr_young_middle, corr_young_old, corr_middle_old])
        #if np.all( np.isnan(chck) == False):
        #    if np.all( chck > 0 ):
        resultDF.loc[len(resultDF)] = (corr_self, corr_non, corr_young, corr_notyoung, corr_middle, corr_notmiddle, corr_old, corr_notold, corr_young_middle, corr_young_old, corr_middle_old)


    return resultDF



#--------------------------------------------------------------------------------------------
def consensus_var(data, mm, codon_seq, annotation, preOUT):
    """
    Computes normalized average profiles by age and extracts codon ribosome densities
    data: dictionary of data matrices for each gene
    annotation: rows in matrices for each gene
    """
    list_age = sorted(list(set(annotation['AGE'])))
    list_idx  = [np.where(annotation['AGE']==x)[0] for x in list_age]
    list_columns = ['young', 'middle', 'old']


    idx_young = np.where(annotation['AGE']==list_age[0])[0]
    idx_middle = np.where(annotation['AGE']==list_age[1])[0]
    idx_old = np.where(annotation['AGE']==list_age[2])[0]
   

    trim = 7 	# omit first and last 7 codons (~20nt) (tho no change if not)

    naivemean_dict = {}

    list_arrays = []
    list_seq = []
    list_orf = []
    for gene in list(data.keys()):
        current_data = data[gene]
        current_data2 = np.zeros(( len(list_age), np.shape(current_data)[1] )) * np.nan
        for j in range(len(list_age)):
            current_data2[j,:] =  np.nanmean( current_data[np.where(annotation['AGE']==list_age[j])[0], :] , axis=0 ) 
            
        current_seq = np.array( codon_seq[gene] )
        current_data_trim = current_data2[:,trim:-trim]
        current_seq_trim = current_seq[trim:-trim]
        current_orf_trim = np.repeat(gene, len(current_seq_trim))

        list_arrays.append(current_data_trim)
        list_seq.append(current_seq_trim)
        list_orf.append(current_orf_trim)

    all_arrays = np.concatenate(list_arrays, axis=1)
    all_seq = np.concatenate(list_seq)
    all_orf = np.concatenate(list_orf)

    print(np.shape(all_arrays))
    print(np.shape(all_seq))

    # !! use R for quantile normalization - change if/when there is time !! 
    np.savetxt("tmp.all", np.around(np.transpose(all_arrays), 4), fmt='%.4f'  )
    cmd_bowtie = 'R --vanilla < qn.R > log'
    output = subprocess.run(cmd_bowtie, shell=True)
    all_arrays_q = np.loadtxt("tmp.all.N")
    all_arrays_q = np.transpose(all_arrays_q)
    os.remove("tmp.all")
    os.remove("tmp.all.N")
    print(np.shape(all_arrays_q))


    meanDF = pd.DataFrame(columns=['codon', 'aa', 'RD_young', 'RD_middle', 'RD_old'])
    for i in codons_nonstop:
        current_idx = all_seq == i
        current_data = all_arrays_q[:,current_idx]
        current_aa = gencode.get(i)

        list_a = []
        list_a2 = []
        for j in range(len(list_age)):
            current_age =  current_data[j, :]                    
            list_a.append(np.around( np.nanmedian(current_age), 3) )

        meanDF.loc[len(meanDF)] = [i] + [current_aa] + list_a #+ list_a2 + list_dd
    meanDF.to_csv("../data/processed/"+preOUT+"_ribodensities.txt", header=True, index=False, sep='\t')


    pairDF = pd.DataFrame(columns=['codon', 'RD_young', 'RD_middle', 'RD_old'])
    for gene in list(data.keys()):
        current_data = all_arrays_q[:,all_orf==gene]
        current_seq = all_seq[all_orf==gene]
     

       
        for pos in range( len(current_seq) - 2 ):      # +2 so that can also extract codon pairs
            current_codon_pair = current_seq[pos] + current_seq[pos+1]
            list_pos = np.zeros(( len(list_idx) ))*np.nan
            for i in range(len(list_idx)):
                current_age = current_data[i,:]
                list_pos[i] = np.around( np.nanmean( np.array([current_age[pos], current_age[pos+1]]) ), 3)
            list_pos = [str(x) for x in list_pos]
            if os.path.exists("TMP_codonpair/"+current_codon_pair+".txt.gz"):
                output = gz.open("TMP_codonpair/"+current_codon_pair+".txt.gz", "at")
            else:
                output = gz.open("TMP_codonpair/"+current_codon_pair+".txt.gz", "wt")
                output.writelines( "\t".join(list_columns)+"\n" )
            output.writelines( "\t".join(list_pos)+"\n" )
            output.close()


    files_codonpair = glob.glob("TMP_codonpair/*.txt.gz")
    for i in files_codonpair:		# [0:1]
        current_data = pd.read_csv(i, header=0, index_col=None, sep='\t', compression='gzip')
        current_data = np.array(current_data)
        current_codonpair = i.split('/')[1].split('.')[0]

        ### !!! dimensions of imported data !!! ###
        list_a = np.zeros(( np.shape(current_data)[1] ))*np.nan
        for j in range(np.shape(current_data)[1]):
            list_a[j] = np.around( np.nanmedian(current_data[:,j]), 3) 
        list_a = [str(x) for x in list_a]
        pairDF.loc[len(pairDF)] = [current_codonpair] + list_a
        os.remove(i)
    pairDF.to_csv("../data/processed/"+preOUT+"_pairdensities.txt", header=True, index=False, sep='\t', na_rep="NA")
   
    return all_arrays_q	 




#--------------------------------------------------------------------------------------------
def get_rand_mc(dataIN):
    """
    get MC from randomized profile
    """
    current_matrix = np.copy(dataIN)
    current_matrix = shuffle_data(current_matrix)

    if np.any( np.isnan(current_matrix)==False ):
        current_mc = np.nanmean(current_matrix, 0)
    else:
        current_mc = np.zeros(( np.shape(current_matrix)[1] ))*np.nan

    return current_mc

#--------------------------------------------------------------------------------------------
def shuffle_data(orf_mat):
    """
    return matrix of randomized profiles
    """
    N = orf_mat.shape[0]
    rand_mat = np.copy(orf_mat)
    for i in range(N):
        local_state = np.random.RandomState(seed=None)
        rand_mat[i,:] = local_state.choice(rand_mat[i,:], len(rand_mat[i,:]), replace=False)

    return rand_mat

#--------------------------------------------------------------------------------------------
def compute_theta(dataIN, mmIN, anno, preOUT):
    """
    get thresholds from randomized consensus
    """

    list_idx  = [np.where(anno['AGE']==x)[0] for x in sorted(list(set(anno['AGE'])))] 

    list_orfs = list( dataIN )
    mc_dict = {}
    theta_df = pd.DataFrame(columns=['ORF', 'p3_10_young', 'p3_90_young', 'p3_10_middle', 'p3_90_middle', 'p3_10_old', 'p3_90_old'])

    for ix, orf in enumerate(list_orfs):
        current_data = dataIN[orf]

        theta_lo = np.zeros(( len(list_idx) )) * np.nan
        theta_hi = np.zeros(( len(list_idx) )) * np.nan
        for j in range(len(list_idx)):
            current_age = current_data[list_idx[j],:]

            #print(ix, orf, current_data.shape[1], len(current_mm))
            if np.any(np.isnan(current_age)==False):
                current_mc = np.nanmean(current_age, axis=0)
                mc_dict[orf] = current_mc

                max_iter = 100        
                pool = mp.Pool(processes=20)
                output = pool.map(get_rand_mc, [current_age for iteration in range(max_iter)])
                output = np.array(output)
                pool.close()
        
                output3 = np.zeros(( output.shape[0], output.shape[1]-2 ))
                for rand_experiment in range(output3.shape[0]):
                    for position in range(output3.shape[1]-2):     #to get kmers of length 3
                        output3[rand_experiment, position] = np.mean(output[rand_experiment, position:position+3])
  
                theta_lo[j] = np.around( np.percentile(output3[np.isnan(output3)==False], 10), 5)           
                theta_hi[j] = np.around( np.percentile(output3[np.isnan(output3)==False], 90), 5)
    
        theta_df.loc[len(theta_df)] = [orf, theta_lo[0], theta_hi[0], theta_lo[1], theta_hi[1], theta_lo[2], theta_hi[2]]
        print([orf, theta_lo[0], theta_hi[0], theta_lo[1], theta_hi[1], theta_lo[2], theta_hi[2]])

    theta_df.to_csv("../data/processed/theta_"+preOUT+".txt", header=True, index=False, sep='\t')
    print(theta_df)


#--------------------------------------------------------------------------------------------
def count_kmer(gene_list, codon_seqs, R, kmer_size=3):
    """
    loop through yeast transcriptome to get general expected redundancies of codon kmers
    this will generate a smallish list of candidate kmers within seconds, so speeds up search

    codon_seqs: dictionary of codon sequence as array of codons
    R: number of occurrances in genome
    kmer_size: length of kmer in number of codons

    MM: 'yes', if multi-mapping, exclude mm regions and count; 'no' count all
    """

    kmer = kmer_size
    MM = 'no'

    list_seqfile = list( codon_seqs.keys() )
    kmer_dict = {}

    for orf in gene_list:
        if orf in list_seqfile:
            current_seq = np.array(codon_seqs[orf])

            for pos in range(len(current_seq) - (kmer + 1) ):
                if MM == 'yes' and orf in list( mm_consensus.keys() ):
                    current_mm = mm_consensus[orf]
                    if np.all(current_mm[pos:(pos+kmer)]):              # check that no kmer position is MM
                        current_kmer = "".join( current_seq[pos:pos+kmer])
                        if current_kmer in kmer_dict.keys():
                            kmer_dict[current_kmer] += 1
                        else:
                            kmer_dict[current_kmer] = 1

                elif MM == 'no':
                    current_kmer = "".join( current_seq[pos:pos+kmer])
                    if current_kmer in kmer_dict.keys():
                        kmer_dict[current_kmer] += 1
                    else:
                        kmer_dict[current_kmer] = 1

    new_dict = {}
    list_redundant = []
    for k in kmer_dict.keys():
        if kmer_dict[k] > R:
            if k not in list_redundant:
        	    list_redundant.append(k)
  
    return list_redundant

#--------------------------------------------------------------------------------------------
def extract_kmers(data_mc, codon_seqs, anno, preOUT):
    """
    extracts all kmers of length 3 codons as per following criteria:
    - redundancy of at least 20 occurrences in set of analyzed genes
    - using pre-computed per-gene theta-thresholds ('theta.txt')
    - at least 10 uccurrences below and 10 above gene-specific thresholds
    """

    theta = pd.read_csv("../data/processed/theta_"+preOUT+".txt", header=0, index_col=False, sep='\t')

    trim = 7			# omit first and last 7 codons / ~20 nt per gene due to known biases of translation initiation and termination
    kmer = 3

    list_idx  = [np.where(anno['AGE']==x)[0] for x in sorted(list(set(anno['AGE'])))] 

    #list_orfs = list( data_mc.keys() )
    list_orfs = list( theta['ORF'] )
    list_candidates = count_kmer(list_orfs, codon_seqs, 20, kmer)

    print(len(list_orfs), len(list_candidates))

    f_kmer = open("../data/processed/kmer_"+preOUT+"_all.txt", "wt")	# write/read more memory efficient!
    f_kmer.writelines( "\t".join(['ORF', 'kmer', 'position', 'young', 'middle', 'old'])+"\n" )

    for ix, orf in enumerate(list_orfs):
        print(ix, orf)

        current_thetas = np.array(theta[theta['ORF']==orf][['p3_10_young', 'p3_90_young', 'p3_10_middle', 'p3_90_middle', 'p3_10_old', 'p3_90_old']]).flatten()
        if np.all( np.isnan(current_thetas)==False ):

            current_consensus = data_mc[orf]
            current_sequence = codon_seqs[orf]

            for pos in range( trim, len(current_sequence) - (trim + kmer) ):	# omit first/last 7/20 positions and allow for kmer length
                current_kmer = current_sequence[pos:pos+kmer]
                current_kmer2 = "".join(current_kmer)

                if current_kmer2 in list_candidates:								# omit multi-mapping positions
                    current_status = np.zeros(( len(list_idx) ), dtype=int) 
                    for j in range(len(list_idx)):
                        current_score = np.nanmean(current_consensus[list_idx[j],:][:,pos:pos+kmer], 0)
                        current_theta_lo =  current_thetas[(j*2)]
                        current_theta_hi =  current_thetas[(j*2)+1]    

                        if np.nanmean( current_score) < current_theta_lo:		# <=
                            current_status[j] = -1
                           
                        elif np.nanmean( current_score) > current_theta_hi:
                            current_status[j] = 1
                    
 
                    list_pos = [str(x) for x in [orf]+[current_kmer2]+[pos]+list(current_status) ]
                    f_kmer.writelines( "\t".join(list_pos)+"\n" )
    f_kmer.close()
    kmer_df = pd.read_csv("../data/processed/kmer_"+preOUT+"_all.txt", header=0, index_col=None, sep='\t')

  

    # filter hits
    def filter_kmertable(kmertable, theta_min=10):
        """
        subfunction to only keep kmers with 10 cases of each lo/hi
        theta_min = 10 each (below / above)
        """

        list_kmers = list(set(kmertable['kmer']))
        list_keep = np.array([])

        for kmer in list_kmers:
            current_df = kmertable[kmertable['kmer']==kmer]
            if len(current_df) >= (2.*theta_min) :
                current_scores = np.array(kmertable[['young', 'middle', 'old']])
                current_lo = np.sum(current_scores == -1, axis=0)
                current_hi = np.sum(current_scores ==  1, axis=0)
                if np.any(current_lo >= theta_min) and np.any(current_hi >= theta_min):
                    current_idx = np.where( np.array( kmertable['kmer']==kmer) )[0]
                    list_keep = np.concatenate((list_keep, current_idx))

        Fout_result = '../data/processed/kmer_'+preOUT+'_filtered.txt'
        result_df = kmertable.iloc[list_keep]
        result_df.to_csv(Fout_result, header=True, index=False, sep='\t')

        return result_df

    #df_filtered = filter_kmertable(kmer_df, 10)


    def filter_kmertable_nonDT(kmertable):
        """
        subfunction to only keep kmers with 15/5 cases of each lo/hi
        col_id: 'theta0595' or 'theta1090'
        hi: > 16 hi and < 4 lo
        lo: > 16 lo and < 4 hi
        """

        theta_min = 10
        theta_hi = 16
        theta_lo = 4

        # only kmers below/above gene-specific thresholds

        list_kmers = list(set(kmertable['kmer']))
        list_keep_hi = np.array([])
        list_keep_lo = np.array([])

        for kmer in list_kmers:
            current_df = kmertable[kmertable['kmer']==kmer]
            if len(current_df) >= (2.*theta_min) :
                current_scores = np.array(current_df[['young', 'middle', 'old']])
                current_lo = np.sum(current_scores == -1, axis=0)
                current_hi = np.sum(current_scores ==  1, axis=0)

                if np.any( (current_hi >= theta_hi) * (current_lo <= theta_lo) ):
                    current_idx = np.where( np.array( kmertable['kmer']==kmer) )[0]
                    list_keep_hi = np.concatenate((list_keep_hi, current_idx))

                elif np.any( (current_hi <= theta_lo) * (current_lo >= theta_hi) ):
                    current_idx = np.where( np.array( kmertable['kmer']==kmer) )[0]
                    list_keep_lo = np.concatenate((list_keep_lo, current_idx))


        Fout_result = '../data/processed/kmer_'+preOUT+'_filtered_nonDT+.txt'
        result_df_nonDT_hi = kmertable.iloc[list_keep_hi]
        result_df_nonDT_hi.to_csv(Fout_result, header=True, index=False, sep='\t')

        Fout_result = '../data/processed/kmer_'+preOUT+'_filtered_nonDT-.txt'
        result_df_nonDT_lo = kmertable.iloc[list_keep_lo]
        result_df_nonDT_lo.to_csv(Fout_result, header=True, index=False, sep='\t')

        return result_df_nonDT_hi, result_df_nonDT_lo 

    #df_filtered_nonDT_hi, df_filtered_nonDT_lo = filter_kmertable_nonDT(kmer_df)

    
    # compile count table
    kmer_counts = pd.DataFrame(columns=['kmer', 'count', 'lo_young', 'hi_young', 'class_young', 'lo_middle', 'hi_middle', 'class_middle', 'lo_old', 'hi_old', 'class_old'])
    list_col = ['young', 'middle', 'old']
    theta_min = 10
    theta_hi = 16
    theta_lo = 4
    for kmer in list(set(kmer_df['kmer'])):
        current_lo = np.zeros(( len(list_idx) ), dtype=int) 
        current_hi = np.zeros(( len(list_idx) ), dtype=int)
        current_class = np.repeat('none', len(list_idx))

        current_df = kmer_df[kmer_df['kmer']==kmer]
        current_scores = np.array(current_df[['young', 'middle', 'old']])
        current_lo = np.sum(current_scores == -1, axis=0)
        current_hi = np.sum(current_scores ==  1, axis=0)

        for j in range(len(list_idx)):
  
            if np.any(current_lo[j] >= theta_min) and np.any(current_hi[j] >= theta_min):
                current_class[j] = 'DT' 
            elif (current_hi[j] >= theta_hi) and (current_lo[j] <= theta_lo) :
                current_class[j] = "nDT+"
            elif (current_hi[j] <= theta_lo) and (current_lo[j] >= theta_hi) :
                current_class[j] = "nDT-"
            else:
                current_class[j] = 'all'
        
            #current_lo[j] = np.sum( np.array(current_df[list_col[j]]) == -1 )
            #current_hi[j] = np.sum( np.array(current_df[list_col[j]]) ==  1 )
 
        kmer_counts.loc[len(kmer_counts)] = (kmer, len(current_df), current_lo[0], current_hi[0], current_class[0], current_lo[1], current_hi[1], current_class[1], current_lo[2], current_hi[2], current_class[2])
        kmer_counts.to_csv('../data/processed/kmer_'+preOUT+'_counts.txt', header=True, index=False, sep='\t')
  
    return kmer_df #, df_filtered, df_filtered_nonDT_hi, df_filtered_nonDT_lo



#--------------------------------------------------------------------------------------------
def kmer_frequencies(kmertable_all, kmertable_filtered, kmertable_nonDT_hi, kmertable_nonDT_lo, data_mm, codon_seqs, preOUT):
    """
    codon frequencies in kmers
    background: mRNA sequences of set of 'good' genes (minus multi-mapping positions)
    redundant: all kmers with 20x redundancy in 'good' genes, i.e. 'all' from extract kmer function
    kmer_filtered: hits at 10%/90% thresholds
    """

    def codon_bgfreq(codon_seqs, data_mm):
        """
        get codon background frequencies from mRNA seqs
        seqs: dictionary of yeast mRNA sequences
        data_mm: dictionary of multi-mapping boolean
        """
        codon_counts = np.zeros(( len(codons_nonstop) ))
        list_orfs = list( data_mm.keys() )

        for ix, orf in enumerate(list_orfs):
            if orf in list(codon_seqs.keys()):
                current_seq = codon_seqs[orf]
                current_mm = data_mm[orf]

                for pos in range( np.min( np.array([len(current_mm), len(current_seq)]) ) ):
                    if current_mm[pos] and (current_seq[pos] in codons_nonstop):
                        current_index = codons_nonstop.index(current_seq[pos])
                        codon_counts[current_index] += 1
        codon_counts = np.around( codon_counts / np.sum(codon_counts), 5)

        return codon_counts


    def codonfreqs_kmerdf(kmertable):
        """
        get codon frequencies from kmertable
        """      
        codon_counts_kmer = np.zeros(( len(codons_nonstop) ))
        for kmer in kmertable['kmer']:
            current_kmer_codons = [ kmer[(i*3):((i*3)+3)] for i in range(3) ] # ! hard coded for length L=3
            for codon in current_kmer_codons:
                current_index = codons_nonstop.index(codon)
                codon_counts_kmer[current_index] += 1        
        codon_counts_kmer /= np.sum(codon_counts_kmer)

        return np.around(codon_counts_kmer, 5)

    #kmertable_threshold = kmertable_all[kmertable_all['threshold']==1]
    kmertable_all2      = kmertable_all[kmertable_all['threshold']==0]


    cc_bg = codon_bgfreq(codon_seqs, data_mm)
    cc_all  = codonfreqs_kmerdf(kmertable_all2)			# without hits
    cc_theta = codonfreqs_kmerdf(kmertable_filtered)
    cc_nDT_hi = codonfreqs_kmerdf(kmertable_nonDT_hi)   # min 16 max 4 at 1090
    cc_nDT_lo = codonfreqs_kmerdf(kmertable_nonDT_lo)   # min 16 max 4 at 1090

    output = pd.DataFrame({'codon': list(codons_nonstop), 
                            'kmer_theta': list(cc_theta), 
                            'redundant': list(cc_all), 
                            'background': list(cc_bg),
                            'nDThi': list(cc_nDT_hi),
                            'nDTlo': list(cc_nDT_lo) } )  
    output.to_csv("../data/processed/kmer_"+preOUT+"_frequencies.txt", header=True, index=False, sep='\t')

    return output




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def histogram_difference(data1, data2):
    """
    computes weighted difference between two histograms
    preset to [0,3] in 299 bins / 300 breaks
    #    r1, p1 = mannwhitneyu(data_i, data_j, use_continuity=False, alternative=None)
    #    r2, p2 = wilcoxon(data_i, data_j, zero_method='wilcox', correction=False)
    """

    bins = np.arange(1000)/10
    #print(len(data1), len(data2))

    histo1, bin_edges = np.histogram(data1, bins=bins, density=True)
    histo2, bin_edges = np.histogram(data2, bins=bins, density=True)

    D = np.around( np.sum( np.sqrt(np.square(histo1 - histo2)) * bins[1:] ), 3)

    return D









if __name__ == '__main__':

    sequence_codons = pickle.load(open('../data/processed/scer/scer_scikit_codonseq.pkl', 'rb'))
    sequence_cele = pickle.load(open('../data/processed/cele/cele_scikit_codonseq.pkl', 'rb'))


    # define codon list as global var
    codons_nonstop = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT', 
    		  'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT',
    		  'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT',
    		  'TAC', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT',]
   

    gencode = pd.DataFrame(
        {'codon':['TTT', 'TTC', 'TTA', 'TTG', 'TCT', 'TCC', 'TCA', 'TCG', 'TAT', 'TAC', 'TAA', 'TAG', 'TGT', 'TGC', 'TGA', 'TGG',
        'CTT', 'CTC', 'CTA', 'CTG', 'CCT', 'CCC', 'CCA', 'CCG', 'CAT', 'CAC', 'CAA', 'CAG', 'CGT', 'CGC', 'CGA', 'CGG',
        'ATT', 'ATC', 'ATA', 'ATG', 'ACT', 'ACC', 'ACA', 'ACG', 'AAT', 'AAC', 'AAA', 'AAG', 'AGT', 'AGC', 'AGA', 'AGG',
        'GTT', 'GTC', 'GTA', 'GTG', 'GCT', 'GCC', 'GCA', 'GCG', 'GAT', 'GAC', 'GAA', 'GAG', 'GGT', 'GGC', 'GGA', 'GGG'],
        'aa':['F', 'F', 'L', 'L', 'S', 'S', 'S', 'S', 'Y', 'Y', '*', '*', 'C', 'C', '*', 'W',
        'L', 'L', 'L', 'L', 'P', 'P', 'P', 'P', 'H', 'H', 'Q', 'Q', 'R', 'R', 'R', 'R',
        'I', 'I', 'I', 'M', 'T', 'T', 'T', 'T', 'N', 'N', 'K', 'K', 'S', 'S', 'R', 'R',
        'V', 'V', 'V', 'V', 'A', 'A', 'A', 'A', 'D', 'D', 'E', 'E', 'G', 'G', 'G', 'G'] })
    gencode = dict(gencode.values)



    tAI = [0.431, 0.616, 1, 0.27, 0.246, 0.488, 0.14, 0.677, 0.677, 0.123, 0.278, 0.054, 0.123, 0.576, 0.616, 0.8, 0.554, 0.431, 0.239, 0.189, 0.616, 0.089, 0.197, 0.123, 0, 0.266, 0.062, 0.369, 0.185, 0.062, 0.059, 0.027, 0.862, 0.985, 0.399, 0.433, 0.308, 0.488, 0.099, 0.677, 0.185, 0.985, 0.182, 0.433, 0.123, 0.621, 0.163, 0.862, 0.493, 0.216, 0.185, 0.488, 0.121, 0.677, 0.246, 0.369, 0.108, 0.431, 0.616, 0.754, 0.27]

    inhibitory = ['CGACCG', 'CGAGCG', 'CGAATA', 'CTCCCG', 'CGACGA', 'CGACGG', 'CGACTG', 'GTACCG', 'GTGCGA', 'GTACGA', 'CTGCCG', 'CTGCGA', 'ATACGA', 'ATACGG', 'AGGCGG', 'CTGATA', 'AGGCGA']
    good_inhibitory = [ "CGAATA", "GTACCG", "GTGCGA", "GTACGA", "CTGCCG", "CTGCGA", "ATACGA", "ATACGG", "AGGCGG", "CTGATA", "AGGCGA"]
        

    ## thresholds for filtering
    theta_ribo = 1    # minimum of 1 read per position on average
    theta_mm   = 0.3  # maximum of 30% multi-mapping

    filter_normalize_data('scer', theta_ribo, theta_mm)
    filter_normalize_data('cele', theta_ribo, theta_mm)
    # - final numbers 
    # scer all: 6283
    # scer filtered: 4306
    # cele all canonical: 14021
    # cele filtered: 4363


    scer_data, scer_mm = load_data('scer')
    cele_data, cele_mm = load_data('cele')
    print(len(scer_data), len(cele_data))

    scer_anno = pd.read_csv("../data/processed/scer/scer_anno.txt", header=0, index_col=False, sep='\t')
    cele_anno = pd.read_csv("../data/processed/cele/cele_anno.txt", header=0, index_col=False, sep='\t')


    # - correlations within and between replicas
    corrDF = correlation_check(scer_data, scer_mm, scer_anno)
    corrDF.to_csv("../data/processed/scer_corr.txt", header=True, index=False, sep='\t', na_rep='NA')

    corrDF = correlation_check(cele_data, cele_mm, cele_anno)
    corrDF.to_csv("../data/processed/cele_corr.txt", header=True, index=False, sep='\t', na_rep='NA')


    # - ribosome densities (RDs)
    scer_allrd = consensus_var(scer_data, scer_mm, sequence_codons, scer_anno, 'scer')
    cele_allrd = consensus_var(cele_data, cele_mm, sequence_cele, cele_anno, 'cele')
 
    scer_allrd = pd.DataFrame(np.transpose(np.around(scer_allrd, 3)) )
    scer_allrd.to_csv("../data/processed/scer_allRD.txt", header=False, index=False, sep='\t', na_rep="NA")
    cele_allrd = pd.DataFrame(np.transpose(np.around(cele_allrd, 3)) )
    cele_allrd.to_csv("../data/processed/cele_allRD.txt", header=False, index=False, sep='\t', na_rep="NA")


    # compute DT sequences
    compute_theta(scer_data, scer_mm, scer_anno, 'scer')
    kmer_all = extract_kmers(scer_data, sequence_codons, scer_anno, 'scer')

    compute_theta(cele_data, cele_mm, cele_anno, 'cele')
    kmer_all = extract_kmers(cele_data, sequence_cele, cele_anno, 'cele')


    # save individual examples
    data_YBR212W = scer_data['YBR212W']
    np.savetxt("../data/processed/riboseq_YBR212W.txt", data_YBR212W)
    data_YMR070W = scer_data['YMR070W']
    np.savetxt("../data/processed/riboseq_YMR070W.txt", data_YMR070W)