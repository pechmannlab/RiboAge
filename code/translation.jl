using CSV;
using DataFrames;


# Stochastic Model of Protein Translation (SMoTP): implementation of Shah et al (Cell 2013)
# (c) S Pechmann | sebastian@pechmannlab.net  (2025)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# parameters
Rf = Int64(2e5);				# Initial number of free ribosomes
tot_ribo = 2e5;           		# Total ribosomes
tot_tRNA = 3.3e6;               # Total tRNAs

Nr = 1.56e6;					# Number available discrete ribosome positions
tau_r = 5e-4;					# characteristic time ribosome

Nt = 1.24e7;					# Number available descrete tRNA positions
tau_t = 4.45e-7;				# chracteristic time tRNA
s = 7.78e-4;					# tRNA competition coefficient


codons_nonstop = ["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"];
codons_stop = ["TAA", "TAG", "TGA"];


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load mRNA sequences into dict of codon seqs
infile = open("../data/SMoPT/codonseqs.txt", "r")
seqs = Dict()
for line in eachline(infile)
    current_line = split(line, "\t")
    if length(current_line) == 2
        current_orf = current_line[1]
        current_seq = current_line[2]
        if length(current_seq)/3 <= 5000						# max length of 5000 codons (longer than longest gene in data)
        	codon_seq = []
        	for i in 1:Int(floor(length(current_seq)/3))		# don't need to floor
            	current_start = ((i-1)*3)+1
            	current_end = ((i-1)*3)+3
        		current_codon = current_seq[current_start:current_end]
        		push!(codon_seq, current_codon)
        	end
        	seqs[current_orf] = codon_seq
        end
    end
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#load abundance data and scale to total expression of 60000 mRNA copies
initiation = CSV.read("../data/SMoPT/abundance.txt", header=1, delim='\t')
abundance = CSV.read("../data/processed/scer_abundance.txt", header=1, delim='\t')

sel_init = []							# only abundance/genes for which init prob avail
for i in 1:(size(abundance)[1])
	current_name = abundance[i,1]
	if current_name in initiation[:,:ORF]
		push!(sel_init, i)
	end
end
abundance = abundance[sel_init, :]

# scale to total expression of simulation system, parameterized as in Shah et al. 
for i in 2:4
	abundance[:,i] = Int64.(round.(abundance[:,i] /sum(abundance[:,i])*60000))
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# initiate dict of matrices with ribosome positions
g = Dict()
N = size(abundance)[1]
for i in 1:N
	current_name = abundance[i,1]
	current_exp = abundance[i,2]			# change here for young (2), middle (3), old (4) !! 
	if current_name in keys(seqs)
		if current_name in initiation[:,:ORF]
			current_seq = seqs[current_name]
			g[current_name] = zeros(Int64, current_exp, length(current_seq) )
		end
 	end
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DF with initiation probs based on sequence specific value and number of free mRNAs
list_orf = []
list_pinit = []
list_copies = []
for i in keys(g)
	current_data = g[i]
    current_init = initiation[findall(initiation[:,:ORF] .== i),:IniProb][1]
	current_fi = sum( sum(current_data[:,1:11], dims=2) .== 0 )
    push!(list_orf, i)
    push!(list_pinit, current_init)
    push!(list_copies, current_fi)
end
init = DataFrame(ORF=list_orf, p_init=list_pinit, n_mRNA=list_copies)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load tRNA file and parse data and indices
trna = CSV.read("../data/SMoPT/translationDF.txt", header=1, delim=' '); 

list_tRNA = sort(unique(trna[:tRNA]))
tRNA_counts = zeros(Int64, length(list_tRNA))
tRNAmap = DataFrame(Codon=[], tRNA=[], index=[], weight=[])
tRNAdict = Dict()			# dict of indices

for i in 1:length(codons_nonstop)
	current_codon = codons_nonstop[i]
	current_tRNA = trna[findall(trna[:,:Codon] .== current_codon),:tRNA][1]
    current_coef = trna[findall(trna[:,:Codon] .== current_codon),:Coef][1]
   	current_idx = findall(list_tRNA .== current_tRNA)[1]
    push!(tRNAmap, [current_codon, current_tRNA, current_idx, current_coef] )

    tRNAdict[current_codon] = current_idx
    if current_codon in list_tRNA
    	tRNA_counts[current_idx] = trna[findall(trna[:,:Codon] .== current_codon),:GCN][1]
    end
end

# scale tRNA counts to experimental vals of total tRNA expression
tRNA_counts =  Int.(round.(tRNA_counts/sum(tRNA_counts)*tot_tRNA))
tRNA_idx = Int.(tRNAmap[:index])
tRNA_weight = Float64.(tRNAmap[:weight])


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# initialize dictionary with stacks for ribosome positions
# ribosome_positions: dict of lists for each codon with ribosomes waiting on that codon
# ribosome_counts: summarizing vector with counts of ribosomes on each codon
# collisions: dict of waiting times of blocked ribosomes to be added after resolve of collision

ribosome_positions = Dict()
for i in codons_nonstop
	ribosome_positions[i] = []
end

collisions = Dict()
ribosome_counts = zeros(Int64, length(codons_nonstop) )


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Main simulation 
preOUT = "young2"				# output prefix for logfiles, CHANGE HERE

n_e_times = zeros(Int64, 61);
count_initiation = 0;
count_elongation = 0;
count_collisions = 0;
count_terminations = 0;

itr = 0;					# iteration counter
t = 0;						# start at t=0
t_max = 150;				# max simulation time
t_reporting = 100;			# burn-in time after which output is written

file_dwell = open("smopt_"*preOUT*"_dwelltimes.txt", "w");
write(file_dwell, "time", "\t", "ORF", "\t", "Codon", "\t", "position", "\t", "dwelltime", "\n");

file_collision = open("smopt_"*preOUT*"_collisions.txt", "w")
write(file_collision, "time", "\t", "ORF", "\t", "Codon", "\t", "position", "\n");

file_log = open("smopt_"*preOUT*"_log.txt", "w");
write(file_log, "iter", "\t", "time", "\t", "Rf", "\t", "fraction_Rf", "\t", "initation", "\t", "elongation", "\t", "collision", "\t", "terminantion", "\n")

while t < t_max

	global itr += 1
 
	# Compute initiation and elongation rates
	r_init = (init[:p_init] .* init[:n_mRNA]) * (Rf / (tau_r * Nr));
	r_elong = (ribosome_counts .* tRNA_counts[tRNA_idx] .* tRNA_weight ) * (s / (tau_t * Nt))

	# Inverse rate
	inv_rate = 1/(sum(r_init) + sum(r_elong))
	global t += inv_rate

	P = [sum(r_init), sum(r_elong)]
	P /= sum(P)
	rr1 = rand()


	if itr % 100000 == 0		# write to logfile every 100000 steps
		write(file_log, string(itr), "\t", string(round(t, digits=5)), "\t", string(Rf), "\t", string(round(Rf/tot_ribo, digits=3)), "\t", string(count_initiation), "\t", string(count_elongation), "\t", string(count_collisions), "\t", string(count_terminations), "\n")
	end

	if itr % 500000 == 0		# print update to screen every 500000 steps
		println([itr, round(t, digits=2), P, Rf, Rf/tot_ribo, count_initiation, count_elongation, count_collisions, count_terminations])
	end



	if rr1 < cumsum(P)[1]
		#println("initiation")
		global count_initiation += 1
        global Rf -= 1											# minus 1 free ribosome

		# Pick gene based on initiable mRNAs and initiation rate
		p_init = r_init / sum(r_init)		# prob of inititation for each gene
		rr = rand()
		gene_choice = findall( rr .< cumsum(p_init) )[1];
		current_name = init[:ORF][gene_choice]

        # draw one mRNA copy
		current_gene = g[current_name]
		current_copy = rand(findall(vec(sum(current_gene[:,1:10], dims=2) .== 0)))

		# occupy position 1 by a new ribosome
		g[current_name][current_copy,1] = 1
		current_codon = seqs[current_name][1]

		current_position = 1
		current_collision = current_gene[current_copy, (current_position+1):minimum([current_position+11, length(seqs[current_name]) ] ) ]
		if all(current_collision .== 0)				# add to stack of elongation competent ribosomes
			ribosome_counts[findall(codons_nonstop .== current_codon)[1] ] += 1	# occupied codons
			push!(ribosome_positions[current_codon], [current_name, current_copy, 1, round(t, digits=5)] )	# stack of ribosome positions
		else										# put on hold
			current_key = current_name*"_"*string(current_copy)*"_"*string(current_position)
			collisions[current_key] = round(t, digits=5)
			global count_collisions += 1
			if t > t_reporting
				write(file_collision, string(round(t, digits=7)), "\t", current_name, "\t", current_codon, "\t", string(current_position+1), "\n");
			end
		end

		# update number of initiable mRNAs
		init[init[:ORF].==current_name, :n_mRNA] = sum( sum(current_gene[:,1:10], dims=2) .== 0 )


	else
		#println("elongation")
		# pick a codon for elongation
		p_elong = r_elong / sum(r_elong)
		rr = rand()
		codon_choice = findall( rr .< cumsum(p_elong) )[1] 
		current_codon_elongation = codons_nonstop[codon_choice]

		# pick random ribosome sitting on codon of choice
		random_ribosome = rand(1:length(ribosome_positions[current_codon_elongation])) 
		current_name, current_copy, current_position, previous_time = ribosome_positions[current_codon_elongation][random_ribosome]
		current_gene = g[current_name]
		# delete current ribosome from stack
		deleteat!(ribosome_positions[current_codon_elongation], random_ribosome)
		ribosome_counts[findall(vec(codons_nonstop .== current_codon_elongation))[1]] -= 1
		g[current_name][current_copy,current_position] = 0

		if current_position + 1 < length(seqs[current_name]) 
			global count_elongation += 1
			g[current_name][current_copy,current_position+1] = 1
			next_codon = seqs[current_name][current_position+1]
			n_e_times[findall(vec(codons_nonstop .== current_codon_elongation))[1]] += 1
			delta_time = t - previous_time
			if t > t_reporting
				write(file_dwell, string(round(t, digits=5)), "\t", current_name, "\t", current_codon_elongation, "\t", string(current_position), "\t", string(round(delta_time, digits=7)), "\n");
			end
		
			current_collision = current_gene[current_copy, (current_position+2):minimum([current_position+11, length(seqs[current_name]) ] ) ]		# + 2 because ribosome_positions already updated
			if all(current_collision .== 0)				# add to stack of elongation competent ribosomes
				ribosome_counts[findall(codons_nonstop .== next_codon)[1] ] += 1	# occupied codons
				push!(ribosome_positions[next_codon], [current_name, current_copy, current_position+1, round(t, digits=5)] )	# stack of ribosome positions
			else										# collision! put on hold
				current_key = current_name*"_"*string(current_copy)*"_"*string(current_position+1)
				collisions[current_key] = round(t, digits=5)
				global count_collisions += 1
				if t > t_reporting
					write(file_collision, string(round(t, digits=7)), "\t", current_name, "\t", next_codon, "\t", string(current_position+1), "\n");
				end
			end

			if current_position > 11		# re-add unblocked ribosome to list of elongating
				past_collision = current_gene[current_copy, current_position-10]
				if past_collision == 1		
					past_codon = seqs[current_name][current_position-10]
					past_key = current_name*"_"*string(current_copy)*"_"*string(current_position-10)	
					#orig_time = get(collisions, past_key, round(t, digits=5) )
					orig_time = collisions[past_key]
					delete!(collisions, past_key)
					push!(ribosome_positions[past_codon], [current_name, current_copy, current_position-10, orig_time])
					ribosome_counts[findall(vec(codons_nonstop .== past_codon))[1]] += 1
				end
			end

			tRNA_counts[tRNAdict[current_codon_elongation]] -= 1		# block current tRNA 
			if current_position > 1 # position 0 doesn't release previous tRNA
				previous_codon = seqs[current_name][current_position-1]
				tRNA_counts[tRNAdict[previous_codon]] += 1	# free up previous tRNA 
			end


		else		# termination at end of sequence 
			global count_terminations += 1
			global Rf += 1
			g[current_name][current_copy,current_position] = 0
			if current_position > 1 # position 1 doesn't release previous tRNA
				previous_codon =  seqs[current_name][current_position-1]
				tRNA_counts[tRNAdict[previous_codon]] += 1		# free up previous tRNA
			end

			# if termination unblocks previously stalled ribosme, re-add to stack
			if current_position > 11
				past_collision = current_gene[current_copy, current_position-10]
				if past_collision == 1		# re-add unblocked ribosome to list of elongating
					past_codon = seqs[current_name][current_position-10]
					past_key = current_name*"_"*string(current_copy)*"_"*string(current_position-10)	
					#orig_time = get(collisions, past_key, round(t, digits=5) )
					orig_time = collisions[past_key]
					delete!(collisions, past_key)
					push!(ribosome_positions[past_codon], [current_name, current_copy, current_position-10, orig_time])
					ribosome_counts[findall(vec(codons_nonstop .== past_codon))[1]] += 1
				end
			end
		end

		# update initiable mRNA copies for current gene
		init[init[:ORF].==current_name, :n_mRNA] = sum( sum(current_gene[:,1:10], dims=2) .== 0 )
	end

end


close(file_dwell)
close(file_collision)
close(file_log)

println("initiation:",count_initiation)
println("elongation:", count_elongation)
println("collision:", count_collisions)
println("termination:",count_terminations)

println(t)
println(n_e_times)
