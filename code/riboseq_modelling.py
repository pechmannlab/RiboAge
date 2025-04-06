import os, sys
import numpy as np
import pandas as pd
import time
import itertools as it
import multiprocessing as mp
import pickle
from scipy import stats



#--------------------------------------------------------------------------------------------
def arrival_probs(ACdf, tRNAdf):
    """
    Compute tRNA arrival probabilities for all codons. 
    Codons without cognate tRNA yield 0. 
    ACdf: codon-anticodon-AA dataframe
    tRNAdf: input tRNA counts and abundances
    """

    # constants from Landerer et al. (MBE 2024) [actually not needed as they cancel out]
    l = 1.5e-8					# effective length
    D = 8.42e-11				# diffusion coefficient
    V = 4.2e-17					# cell volume 

    # diffusion
    tau = l**2 / (6 * D)		# transition time between locations
    n_loc = V / (l**3)			# approx number of cell locations

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



#--------------------------------------------------------------------------------------------
def binding_probs(ACdf, E2):
    """
    Computes tRNA-codon binding probability given matrix of all binding probabilities
    """
    
    nt = ['A', 'C', 'G', 'T']

    list_codons = list(ACdf['codon'])
    list_tRNAs = list(ACdf[ACdf['lambda']>0]['anticodon'])

    tRNAdf = ACdf[ACdf['lambda']>0]
    tRNAdf.index = tRNAdf['anticodon']

    pa =  tRNAdf['pa']
    pb =  E2.loc[list_codons][list_tRNAs]

    pa_n = 1 - pa
    #start_time = time.time()

    resultDF = pd.DataFrame(columns=['codon', 'AA', 'p_decoding'])
    for cx, codon in enumerate(list_codons):
        current_AA = list(ACdf[ACdf['codon']==codon]['AA2'])[0]
        current_syn = list(tRNAdf[tRNAdf['AA2']==current_AA]['codon'])

        p_decoding = np.zeros(( len(current_syn) ))*np.nan
        for sx, c in enumerate(current_syn):
            current_AC = list(ACdf[ACdf['codon']==c]['anticodon'])[0]
            if current_AC in list(tRNAdf['anticodon']):
                codon_pb = pb.loc[codon]

                current_pb = pb.loc[codon, current_AC]  
                current_pa = tRNAdf.loc[current_AC, 'pa']
                current_pa_n = 1-current_pa

                current_ar = tRNAdf.loc[current_AC, 'lambda']
                other_ar = np.array(tRNAdf['lambda'])[np.array(tRNAdf['anticodon'])!=current_AC]
                other_pb = np.array(codon_pb)[np.array(tRNAdf['anticodon'])!=current_AC]

                current_pb_n = current_pb * np.prod( np.power( (1-other_pb), (other_ar)/current_ar )  )

                p_decoding[sx] =  (current_pa * current_pb) + (current_pa_n * current_pb_n)
         
        resultDF.loc[len(resultDF)] = (codon, current_AA, np.sum(p_decoding) )

    return resultDF



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def compute_pb(pos_weights, nt_energies):
    """
    Computes binding probabilities from positional weights and nt interaction energies
    """

    x = list(it.product((0,1,2,3),(0,1,2,3),(0,1,2,3)))
    y = np.array(list(it.product(x,x)))
    z = ["".join(a) for a in list(np.array(nt)[np.array(x)]) ]

    ENG = np.reshape( np.sum( pos_weights * nt_energies[y[:,0,:], y[:,1,:][:,::-1] ], axis=1), (64, 64) )
  
    delta_A = np.exp(-ENG)  # stability of codon/anticodon interaction
    delta_A = np.transpose(delta_A) # transpose to normalize along axis=1
    delta_A /= np.sum(delta_A, axis=0)
    pbMAT = np.transpose(delta_A)

    E2 = pd.DataFrame(data=pbMAT, columns=z, index=z)
   
    return E2



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def compute_reltRNA(ACdf, tRNAdf, EE):
    """
    normalized tRNA abundances
    split codon usage into tRNAs based on relative binding energies
    codon usage has to be included in ACdf
    """

    #tRNAdf = tRNAdf[tRNAdf['Neigh_Nucl']!="N"]			# !! only for Scer !! 

    normDF = ACdf.copy()
    cognate = np.zeros(( len(normDF) ), dtype=int)
    for ix, i in enumerate(normDF['anticodon']):
        if i in list(tRNAdf['Anticodon']):
            cognate[ix] = 1
    normDF['tRNA'] = list(cognate)

    rel_usage = np.zeros(( len(normDF) )) 
    list_ac = list(normDF['anticodon'])
    for i in list(normDF['codon']):
        current_codon_aa = list(normDF[normDF['codon']==i]['AA2'])[0]
        current_ac = normDF[normDF['AA2']==current_codon_aa]
        current_ac = list( current_ac[current_ac['tRNA']==1]['anticodon'] )
        ac_idx = [list_ac.index(i) for i in list(current_ac) ]

        current_cu = list(normDF[normDF['codon']==i]['CU'])[0]

        current_pb = np.array( EE.loc[i][current_ac] )
        current_rel = current_pb/np.sum(current_pb)

        rel_usage[ac_idx] += (current_cu * current_rel) * 100

    normDF['rCU'] = rel_usage 

    list_idx = []
    for i in list(tRNAdf['Anticodon']):
        current_idx = list_ac.index(i)
        list_idx.append(current_idx)

    rCU = normDF.iloc[list_idx]['rCU']
    relAbundance = np.around( np.array(tRNAdf['Abundance']) / np.array(rCU), 3)
    tRNAdf2 = pd.DataFrame({'Anticodon': list(tRNAdf['Anticodon']), 'Count': list(tRNAdf['Count']), 'Abundance': list(relAbundance)} )

    return tRNAdf2



#--------------------------------------------------------------------------------------------
def optimize_weights(ACdf, tRNAdf, Emat, Pmat, CR, seed=0):
    """
    : 1000 steps all params, full E and P sampled
    : 1000 steps gradient down, E and P subsampled in 2 steps
    : 1000 steps single params, single params of E and P subsampled in 2 steps
    """
 
    np.random.seed(seed=seed)

    probs_arrival = arrival_probs(ACdf, tRNAdf, )
    probs_decoding = binding_probs(probs_arrival, Pmat, Emat)
    DR = probs_decoding.sort_values(by=['codon'])['p_decoding'] 
    best_CC = np.corrcoef( CR, DR )[0,1]
    #best_CC = stats.spearmanr( CR, DR ).correlation
    best_E = np.copy(Emat)
    best_P = np.copy(Pmat)
    best_DR = np.copy(DR)


    p_min = 0.2
    p_max = 1   
    iter_phase = 1000 
    subsample_gradient = np.array([np.around( p_min + (p_max - p_min) * ((iter_phase - i)/iter_phase) , 3) for i in np.arange(iter_phase)]) #**2		# linear or quadratic 
    subsample = np.concatenate((np.repeat(1, iter_phase), subsample_gradient, np.repeat(0, iter_phase) ) )

    iter_max = len(subsample)

    for i in range(iter_max):
        subsample_prob = subsample[i]
        if np.max(probs_decoding['p_decoding']) < 0.0001:	# edge case all to zero then stop
            break

        # optimizing energy matrix
        if subsample_prob == 1:							# full sampling of E and P
            proposal_E = np.random.normal(best_E, 0.5)
            proposal_P = np.ones((3))
            proposal_P[0:2] = np.random.normal(best_P[0:2], 0.5)
            proposal_P[proposal_P < 1] = np.abs(proposal_P[proposal_P < 1]) + 1 	# > 1

        elif subsample_prob < 1 and subsample_prob > 0:	# subsampling of E 
            idx_E = np.reshape( np.random.binomial(1, subsample_prob, size=16), (4,4) )
            proposal_E = np.copy(best_E)
            proposal_E[idx_E] = np.random.normal(best_E[idx_E], 0.5)

        elif subsample_prob == 0:						# sampling one entry of E
            idx_E = np.random.choice(np.arange(4), 2, replace=True)
            proposal_E = np.copy(best_E)
            proposal_E[idx_E] = np.random.normal(best_E[idx_E], 0.5)

        probs_decoding = binding_probs(probs_arrival, proposal_P, proposal_E)
        DR = probs_decoding.sort_values(by=['codon'])['p_decoding'] 
        CC = np.corrcoef( CR, DR )[0,1]
        if CC > best_CC:
            best_E = np.copy(proposal_E)
            best_P = np.copy(proposal_P)
            best_CC = np.copy(CC)
            best_DR = np.copy(DR)


        # optimizing position matrix P
        if subsample_prob < 1 and subsample_prob > 0:	# subsampling of P
            idx_P = np.random.choice(np.arange(2), 1)	# pick one
            proposal_P = np.copy(best_P)
            if np.random.random(1) < 0.9: 		# one new value
                proposal_P[idx_P] = np.random.normal(best_P[idx_P], 0.5)      
            else:								# flip values
                proposal_P[[0,1]] = best_P[[1,0]]
            proposal_P[proposal_P < 1] = np.abs(proposal_P[proposal_P < 1]) + 1 	

        elif subsample_prob == 0:
            idx_P = np.random.choice(np.arange(2), 1)
            proposal_P = np.copy(best_P)
            if np.random.random(1) < 0.95: 		# one new value
                proposal_P[idx_P] = np.random.normal(best_P[idx_P], 0.5)
            else:								# flip values
                proposal_P[[0,1]] = best_P[[1,0]]
            proposal_P[proposal_P < 1] = np.abs(proposal_P[proposal_P < 1]) + 1

        probs_decoding = binding_probs(probs_arrival, proposal_P, proposal_E)
        DR = probs_decoding.sort_values(by=['codon'])['p_decoding'] 
        CC = np.corrcoef( CR, DR )[0,1]
        if CC > best_CC:
            best_E = np.copy(proposal_E)
            best_P = np.copy(proposal_P)
            best_CC = np.copy(CC)
            best_DR = np.copy(DR)


    print(seed, "final best:", best_CC)
    #np.savetxt("tmp.cc", np.array([CR, best_DR]))

    return best_E, best_P, best_CC, best_DR



#--------------------------------------------------------------------------------------------
def optimize_weights_tRNA(ACdf, tRNAdf, Emat, Pmat, CR, seed=0):
    """
    Optimize tRNA abundances
    FC constraint has to be changed in code below
    : 1000 steps all params
    : 1000 steps gradient down
    : 1000 steps single params
    """
 
    tRNAdf2 = tRNAdf.loc[tRNAdf['Abundance'] > 0]
 

    np.random.seed(seed=seed)

    probs_arrival = arrival_probs(ACdf, tRNAdf2)
    probs_decoding = binding_probs(probs_arrival, Pmat, Emat)
    DR = probs_decoding.sort_values(by=['codon'])['p_decoding'] 
    best_CC = np.corrcoef( CR, DR )[0,1]
    best_DR = np.copy(DR)
    best_tRNA = tRNAdf2.copy()

    p_min = 0.2
    p_max = 1   
    iter_phase = 1000 
    subsample_gradient = np.array([np.around( p_min + (p_max - p_min) * ((iter_phase - i)/iter_phase) , 3) for i in np.arange(iter_phase)]) #**2		# linear or quadratic 
    subsample = np.concatenate((np.repeat(1, iter_phase), subsample_gradient, np.repeat(0, iter_phase) ) )
  

    if np.max(tRNAdf['Abundance']) > 1000:
        sample_sd = 10000		# scer / abundance
    else:
        sample_sd = 1			# cele / GCN

    constraint = 2 #1.5	#3. #2.  		# CHANGE FC CONSTRAINT HERE 
    iter_max = len(subsample)
    for i in range(iter_max):
        if np.max(probs_decoding['p_decoding']) < 0.0001:	# edge case all to zero then stop
            break

        subsample_prob = subsample[i]
        proposal_tRNA = best_tRNA.copy()
        proposal_abundance = np.array(proposal_tRNA['Abundance'])

        if subsample_prob == 1:							# full sampling 
            proposal_abundance = np.random.normal(np.array(best_tRNA['Abundance']), sample_sd) 
            if np.all(proposal_abundance > 0.01): 			# 0 # only positive tRNA abundances! 
                current_constraint = proposal_abundance / np.array(tRNAdf2['Abundance'])
                check_constraint = (current_constraint < constraint) * (current_constraint > (1./constraint) )		# within fold-change constraint
                if np.all(check_constraint):	# constrained fold-change to WT
                    proposal_tRNA['Abundance'] = proposal_abundance

        elif subsample_prob < 1 and subsample_prob > 0:	# subsampling 
            idx_E =  np.random.binomial(1, subsample_prob, size=len(proposal_tRNA))
            proposal_abundance[idx_E] = np.random.normal(np.array(best_tRNA['Abundance'])[idx_E], sample_sd)
            if np.all(proposal_abundance > 0.01):
                current_constraint = proposal_abundance / np.array(tRNAdf2['Abundance'])
                check_constraint = (current_constraint < constraint) * (current_constraint > (1./constraint) )
                if np.all(check_constraint):
                    proposal_tRNA['Abundance'] = proposal_abundance

        elif subsample_prob == 0:						# sampling one entry
            idx_E = np.random.choice(np.arange(len(proposal_tRNA)), replace=False)
            proposal_abundance[idx_E] = np.random.normal(np.array(best_tRNA['Abundance'])[idx_E], sample_sd)
            if np.all(proposal_abundance > 0.01):
                current_constraint = proposal_abundance / np.array(tRNAdf2['Abundance'])
                check_constraint = (current_constraint < constraint) * (current_constraint > (1./constraint) )
                if np.all(check_constraint):
                    proposal_tRNA['Abundance'] = proposal_abundance


        probs_arrival = arrival_probs(ACdf, proposal_tRNA)
        probs_decoding = binding_probs(probs_arrival, Pmat, Emat)
        DR = probs_decoding.sort_values(by=['codon'])['p_decoding'] 
        CC = np.corrcoef( CR, DR )[0,1]

        if CC > best_CC:
            best_CC = np.copy(CC)
            best_DR = np.copy(DR)
            best_tRNA = proposal_tRNA.copy()

    print(seed, "final best:", best_CC)

    return best_tRNA, best_CC, best_DR





#--------------------------------------------------------------------------------------------
def optimize_weights_normalized(ACdf, tRNAdf, Emat, Pmat, CR, seed=0):
    """
    Optimize energy weights for normalized tRNA abundances
    : 1000 steps all params, full E and P sampled
    : 1000 steps gradient down, E and P subsampled in 2 steps
    : 1000 steps single params, single params of E and P subsampled in 2 steps
    """
 
    np.random.seed(seed=seed)

    E2 = compute_pb(Pmat, Emat)
    norm_tRNAdf = compute_reltRNA(ACdf, tRNAdf, E2)
    probs_arrival = arrival_probs(ACdf, norm_tRNAdf)
    probs_decoding = binding_probs(probs_arrival, E2) #Pmat, Emat)
    DR = probs_decoding.sort_values(by=['codon'])['p_decoding'] 

    best_CC = np.corrcoef( CR, DR )[0,1]
    best_E = np.copy(Emat)
    best_P = np.copy(Pmat)
    best_DR = np.copy(DR)

    p_min = 0.2
    p_max = 1   
    iter_phase = 1000 
    subsample_gradient = np.array([np.around( p_min + (p_max - p_min) * ((iter_phase - i)/iter_phase) , 3) for i in np.arange(iter_phase)]) #**2		# linear or quadratic 
    subsample = np.concatenate((np.repeat(1, iter_phase), subsample_gradient, np.repeat(0, iter_phase) ) )

    iter_max = len(subsample)

    for i in range(iter_max):
        subsample_prob = subsample[i]
        if np.max(probs_decoding['p_decoding']) < 0.0001:	# edge case all to zero then stop
            break

        # optimizing energy matrix
        if subsample_prob == 1:							# full sampling of E and P
            proposal_E = np.random.normal(best_E, 0.5)
            proposal_P = np.ones((3))
            proposal_P[0:2] = np.random.normal(best_P[0:2], 0.5)
            proposal_P[proposal_P < 1] = np.abs(proposal_P[proposal_P < 1]) + 1 	# > 1

        elif subsample_prob < 1 and subsample_prob > 0:	# subsampling of E 
            idx_E = np.reshape( np.random.binomial(1, subsample_prob, size=16), (4,4) )
            proposal_E = np.copy(best_E)
            proposal_E[idx_E] = np.random.normal(best_E[idx_E], 0.5)

        elif subsample_prob == 0:						# sampling one entry of E
            idx_E = np.random.choice(np.arange(4), 2, replace=True)
            proposal_E = np.copy(best_E)
            proposal_E[idx_E] = np.random.normal(best_E[idx_E], 0.5)


        E2 = compute_pb(proposal_P, proposal_E)
        norm_tRNAdf = compute_reltRNA(ACdf, tRNAdf, E2)
        probs_arrival = arrival_probs(ACdf, norm_tRNAdf)
        probs_decoding = binding_probs(probs_arrival, E2) #proposal_P, proposal_E)
        DR = probs_decoding.sort_values(by=['codon'])['p_decoding'] 
        CC = np.corrcoef( CR, DR )[0,1]
        if CC > best_CC:
            #print(CC)
            best_E = np.copy(proposal_E)
            best_P = np.copy(proposal_P)
            best_CC = np.copy(CC)
            best_DR = np.copy(DR)


        # optimizing position matrix P
        if subsample_prob < 1 and subsample_prob > 0:	# subsampling of P
            idx_P = np.random.choice(np.arange(2), 1)	# pick one
            proposal_P = np.copy(best_P)
            if np.random.random(1) < 0.9: 		# one new value
                proposal_P[idx_P] = np.random.normal(best_P[idx_P], 0.5)      
            else:								# flip values
                proposal_P[[0,1]] = best_P[[1,0]]
            proposal_P[proposal_P < 1] = np.abs(proposal_P[proposal_P < 1]) + 1 	

        elif subsample_prob == 0:
            idx_P = np.random.choice(np.arange(2), 1)
            proposal_P = np.copy(best_P)
            if np.random.random(1) < 0.95: 		# one new value
                proposal_P[idx_P] = np.random.normal(best_P[idx_P], 0.5)
            else:								# flip values
                proposal_P[[0,1]] = best_P[[1,0]]
            proposal_P[proposal_P < 1] = np.abs(proposal_P[proposal_P < 1]) + 1

        E2 = compute_pb(proposal_P, proposal_E)
        norm_tRNAdf = compute_reltRNA(ACdf, tRNAdf, E2)
        probs_arrival = arrival_probs(ACdf, norm_tRNAdf)
        probs_decoding = binding_probs(probs_arrival, E2) #proposal_P, proposal_E)
        DR = probs_decoding.sort_values(by=['codon'])['p_decoding'] 
        CC = np.corrcoef( CR, DR )[0,1]
        if CC > best_CC:
            #print(CC)
            best_E = np.copy(proposal_E)
            best_P = np.copy(proposal_P)
            best_CC = np.copy(CC)
            best_DR = np.copy(DR)


    print(seed, "final best:", best_CC)
    #np.savetxt("tmp.cc", np.array([CR, best_DR]))

    return best_E, best_P, best_CC, best_DR



#--------------------------------------------------------------------------------------------
def optimize_weights_parallel(packaged_input):
    """
    import/export for parallelization
    """

    ACdf = packaged_input[0]
    tRNAdf = packaged_input[1]
    Emat = packaged_input[2]
    Pmat = packaged_input[3]
    CR = packaged_input[4]
    seeed = packaged_input[5]

    optE, optP, optCC, optDR = optimize_weights(ACdf, tRNAdf, Emat, Pmat, CR, seeed)

    return optE, optP, optCC, optDR


#--------------------------------------------------------------------------------------------
def optimize_weights_tRNA_parallel(packaged_input):
    """
    import/export for parallelization
    """

    ACdf = packaged_input[0]
    tRNAdf = packaged_input[1]
    Emat = packaged_input[2]
    Pmat = packaged_input[3]
    CR = packaged_input[4]
    seeed = packaged_input[5]

    opt_tRNA, optCC, optDR = optimize_weights_tRNA(ACdf, tRNAdf, Emat, Pmat, CR, seeed)

    return opt_tRNA, optCC, optDR


#--------------------------------------------------------------------------------------------
def optimize_weights_normalized_parallel(packaged_input):
    """
    import/export for parallelization
    """

    ACdf = packaged_input[0]
    tRNAdf = packaged_input[1]
    Emat = packaged_input[2]
    Pmat = packaged_input[3]
    CR = packaged_input[4]
    seeed = packaged_input[5]

    optE, optP, optCC, optDR = optimize_weights_normalized(ACdf, tRNAdf, Emat, Pmat, CR, seeed)

    return optE, optP, optCC, optDR





#--------------------------------------------------------------------------------------------
def mc_optimization(ACdf, tRNAdf, matEn, matPos, trueCR, iterN, preOUT, modus='energy'):
    """
    ACdf, tRNAdf to compute p_arrival: computed probabilities of tRNA arrivals
    matEn: matrix of nucleotide effective interaction energies
    matPos: matrix of position penalties
    trueCR: 'codon translation rates' from riboseq as 'ground truth'
    iterN: number of runs
    preOUT: prefix for output file names
    """
 
    result_E = np.zeros(( 4, 4, iterN )) * np.nan
    result_P = np.zeros(( 3, iterN )) * np.nan
    result_CC = np.zeros(( iterN )) * np.nan
    result_DR = np.zeros(( 61, iterN )) * np.nan 
    result_tRNA = []


    if modus == 'energy':               # refit anticodon-codon energies
        pool = mp.Pool(processes=22)
        pool_set = np.arange(iterN)
        if len(np.shape(trueCR)) == 1:
            pool_inputs = [(ACdf, tRNAdf, matEn, matPos, trueCR, s) for s in pool_set]
        elif len(np.shape(trueCR)) > 1: 
            pool_inputs = [(ACdf, tRNAdf, matEn, matPos, trueCR[:,s], s) for s in pool_set]
        result = pool.map(optimize_weights_parallel, pool_inputs )  
        for rx, res in enumerate(result):
            result_E[:,:,rx] = res[0]
            result_P[:,rx] = res[1]
            result_CC[rx] = res[2]
            result_DR[:,rx] = res[3]
        pool.close()
        pool.join()

        pickle.dump(result_E, open('../data/processed/'+preOUT+'_E.pkl', 'wb'))
        np.savetxt('../data/processed/'+preOUT+'_P.txt', np.around(result_P, 5), fmt='%.5f')
        np.savetxt('../data/processed/'+preOUT+'_CC.txt', np.around(result_CC, 5) , fmt='%.5f')
        np.savetxt('../data/processed/'+preOUT+'_DR.txt', np.around(result_DR, 5), fmt='%.5f')

        print(result_CC)
        print( "best found:", np.around(np.max(result_CC), 5) )


    elif modus == 'tRNA':           # refit/optimize tRNA abundances
        print('now in tRNA modus')
        pool = mp.Pool(processes=22)
        pool_set = np.arange(iterN)
        pool_inputs = [(ACdf, tRNAdf, matEn[:,:,s], matPos[:,s], trueCR, s) for s in pool_set]
        result = pool.map(optimize_weights_tRNA_parallel, pool_inputs )  
        for rx, res in enumerate(result):
            result_tRNA.append( res[0]['Abundance'] )
            result_CC[rx] = res[1]
            result_DR[:,rx] = res[2]
        pool.close()
        pool.join()

        res_tRNA = pd.concat(result_tRNA, axis=1)
        res_tRNA.to_csv('../data/processed/'+preOUT+'_tRNA.txt', sep='\t')
        np.savetxt('../data/processed/'+preOUT+'_CC.txt', np.around(result_CC, 5) , fmt='%.5f')
        np.savetxt('../data/processed/'+preOUT+'_DR.txt', np.around(result_DR, 5), fmt='%.5f')

        print(result_CC)
        print( "best found:", np.around(np.max(result_CC), 5) )


    elif modus == 'normalized':     # refit anticodon-codon energies for normalized tRNA abundances
        print('now in normalized modus')

        pool = mp.Pool(processes=22)
        pool_set = np.arange(iterN)
        if len(np.shape(trueCR)) == 1:
            pool_inputs = [(ACdf, tRNAdf, matEn, matPos, trueCR, s) for s in pool_set]
        elif len(np.shape(trueCR)) > 1: 
            pool_inputs = [(ACdf, tRNAdf, matEn, matPos, trueCR[:,s], s) for s in pool_set]
        result = pool.map(optimize_weights_normalized_parallel, pool_inputs )  
        for rx, res in enumerate(result):
            result_E[:,:,rx] = res[0]
            result_P[:,rx] = res[1]
            result_CC[rx] = res[2]
            result_DR[:,rx] = res[3]
        pool.close()
        pool.join()

        pickle.dump(result_E, open('../data/processed/'+preOUT+'_E.pkl', 'wb'))
        np.savetxt('../data/processed/'+preOUT+'_P.txt', np.around(result_P, 5), fmt='%.5f')
        np.savetxt('../data/processed/'+preOUT+'_CC.txt', np.around(result_CC, 5) , fmt='%.5f')
        np.savetxt('../data/processed/'+preOUT+'_DR.txt', np.around(result_DR, 5), fmt='%.5f')

        print(result_CC)
        print( "best found:", np.around(np.max(result_CC), 5) )





#--------------------------------------------------------------------------------------------
def optinp_tRNA(preOUT):
    """
    Extract parameter sets from best performing models
    """

    top = 10
    mpl = 20

    E = pickle.load(open('../data/processed/'+preOUT+'_E.pkl', 'rb'))
    P = np.loadtxt('../data/processed/'+preOUT+'_P.txt')
    CC = np.loadtxt('../data/processed/'+preOUT+'_CC.txt')

    CC_order = np.argsort(CC)
    CC_order = CC_order[::-1]

    newE = np.zeros(( np.shape(E) )) * np.nan
    newP = np.zeros(( np.shape(P) )) * np.nan

    idx = 0
    for i in range(top):
        current_E = E[:,:,i]
        current_P = P[:,i]
        for j in range(mpl):
            newE[:,:,idx] = np.copy(current_E)
            newP[:,idx] = np.copy(current_P)
            idx += 1

    return newE, newP







if __name__ == '__main__':


    nt = ['A', 'C', 'G', 'T']

    # define codon list as global var
    codons_nonstop = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT', 
    		  'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT',
    		  'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT',
    		  'TAC', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT',]
   


    # indices as global vars to speed up
    x = list(it.product((0,1,2,3),(0,1,2,3),(0,1,2,3)))
    y = np.array(list(it.product(x,x)))
    z = ["".join(a) for a in list(np.array(nt)[np.array(x)]) ]



    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    #########
    # YEAST #
    #########

    ribodensities = pd.read_csv("../data/processed/scer_ribodensities.txt", header=0, index_col=None, sep='\t')
    ribodensities.index = list(ribodensities['codon'])
    ribodensities = ribodensities.loc[codons_nonstop]

    anticodonDF = pd.read_csv("../data/tRNA/anticodon2.txt", header=0, index_col=None, sep='\t')
    trnaDF = pd.read_csv('../data/tRNA/yeast_tRNA_count.csv', delimiter=",", header=0)
    
    
    RD = ribodensities.sort_values(by=['codon'])
    CR = 1/RD['RD_young']
    CR_middle = 1/RD['RD_middle']
    CR_old = 1/RD['RD_old']


    EEmat = np.array([[-0.1, -0.1, -0.1, -2],
    			   [-0.1, -0.1, -2, -0.1],
    			   [-0.1, -2, -0.1, -0.1], 
    			   [-2, -0.1, -0.1, -0.1]])

    posPen = np.array([3, 3, 1])


    
    
    # 1. optimizing binding energies ('basic')
    mc_optimization(anticodonDF, trnaDF, EEmat, posPen, CR, 200, "scer2_basic_young")
    mc_optimization(anticodonDF, trnaDF, EEmat, posPen, CR_middle, 200, "scer2_basic_middle")
    mc_optimization(anticodonDF, trnaDF, EEmat, posPen, CR_old, 200, "scer2_basic_old")

    # 2. optimizing tRNA abundances: run for each FC constraint
    top10_E, top10_P = optinp_tRNA("scer_basic_young")    
    mc_optimization(anticodonDF, trnaDF, top10_E, top10_P, CR, 200, "scer2_tRNA_2_young", modus='tRNA')
    mc_optimization(anticodonDF, trnaDF, top10_E, top10_P, CR_middle, 200, "scer2_tRNA_2_middle", modus='tRNA')
    mc_optimization(anticodonDF, trnaDF, top10_E, top10_P, CR_old, 200, "scer2_tRNA_2_old", modus='tRNA')

    # 3. optimizing energies with normlized CU 
    scer_codonusage = pd.read_csv("../data/processed/scer_codonusage.txt", header=0, index_col=0, sep='\t')
    acDF_young = anticodonDF.copy()
    acDF_young['CU'] = np.array(scer_codonusage['young'])
    acDF_middle = anticodonDF.copy()
    acDF_middle['CU'] = list(scer_codonusage['middle'])
    acDF_old = anticodonDF.copy()
    acDF_old['CU'] = list(scer_codonusage['old'])

    mc_optimization(acDF_young, trnaDF, EEmat, posPen, CR, 200, "scer3_cu_young", modus='normalized')
    mc_optimization(acDF_middle, trnaDF, EEmat, posPen, CR_middle, 200, "scer3_cu_middle", modus='normalized')
    mc_optimization(acDF_old, trnaDF, EEmat, posPen, CR_old, 200, "scer3_cu_old", modus='normalized')



    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ########
    # CELE #
    ########

    ribodensities = pd.read_csv("../data/processed/cele_ribodensities.txt", header=0, index_col=None, sep='\t')
    ribodensities.index = list(ribodensities['codon'])
    ribodensities = ribodensities.loc[codons_nonstop]

    anticodonDF = pd.read_csv("../data/tRNA/anticodon2.txt", header=0, index_col=None, sep='\t')
    trnaDF = pd.read_csv('../data/tRNA/cele_tRNA_count.txt', delimiter=" ", header=0)
       
    RD = ribodensities.sort_values(by=['codon'])
    CR = 1/RD['RD_young']
    CR_middle = 1/RD['RD_middle']
    CR_old = 1/RD['RD_old']

    #np.savetxt("tmp.true", CR_old)


    EEmat = np.array([[-0.1, -0.1, -0.1, -2],
    			   [-0.1, -0.1, -2, -0.1],
    			   [-0.1, -2, -0.1, -0.1], 
    			   [-2, -0.1, -0.1, -0.1]])

    posPen = np.array([3, 3, 1])



    # 1. energy opt
    mc_optimization(anticodonDF, trnaDF, EEmat, posPen, CR, 200, "cele2_basic_young")
    mc_optimization(anticodonDF, trnaDF, EEmat, posPen, CR_middle, 200, "cele2_basic_middle")
    mc_optimization(anticodonDF, trnaDF, EEmat, posPen, CR_old, 200, "cele2_basic_old")

    # 2. tRNA opt
    top10_E, top10_P = optinp_tRNA("cele2_basic_young")    
    mc_optimization(anticodonDF, trnaDF, top10_E, top10_P, CR, 200, "cele2_tRNA_2_young", modus='tRNA')
    mc_optimization(anticodonDF, trnaDF, top10_E, top10_P, CR_middle, 200, "cele2_tRNA_2_middle", modus='tRNA')
    mc_optimization(anticodonDF, trnaDF, top10_E, top10_P, CR_old, 200, "cele2_tRNA_2_old", modus='tRNA')

    #3. codon usage energy opt 
    cele_codonusage = pd.read_csv("../data/processed/cele_codonusage.txt", header=0, index_col=0, sep='\t')
    acDF_young = anticodonDF.copy()
    acDF_young['CU'] = np.array(cele_codonusage['young'])
    acDF_middle = anticodonDF.copy()
    acDF_middle['CU'] = list(cele_codonusage['middle'])
    acDF_old = anticodonDF.copy()
    acDF_old['CU'] = list(cele_codonusage['old'])

    mc_optimization(acDF_young, trnaDF, EEmat, posPen, CR, 200, "cele3_cu_young", modus='normalized')
    mc_optimization(acDF_middle, trnaDF, EEmat, posPen, CR_middle, 200, "cele3_cu_middle", modus='normalized')
    mc_optimization(acDF_old, trnaDF, EEmat, posPen, CR_old, 200, "cele3_cu_old", modus='normalized')
