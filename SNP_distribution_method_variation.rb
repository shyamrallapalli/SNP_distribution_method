#encoding: utf-8
require_relative 'lib/reform_ratio'
require_relative 'lib/write_it'
require_relative 'lib/stuff'
require_relative 'lib/mutation'
require_relative 'lib/SDM'
require_relative 'lib/snp_dist'
require_relative 'lib/plot'

require 'pp'
require 'benchmark'
require 'csv'


if ARGV.empty?
	puts "Please specify a (1) dataset, a (2) name for the output folder, a (3) threshold to discard the contigs which a ratio below it and a (4) factor to calculate the ratio (1, 0.1, 0.01...) (5) kind of cross: back or out"
else 
	dataset = ARGV[0] 
	file = ARGV[1]
	threshold = ARGV[2].to_i #degree of filtering:  100, 50, 10, 5
	adjust = ARGV[3]
	cross = ARGV[4]
	puts "Looking for SNPs in #{dataset}"
	puts "Output will be in #{dataset}/#{file}"
	puts "A factor of #{adjust} will be used to calculate the ratio"
	puts "#{cross}-cross"
	if threshold > 0
		puts "Filtering step on"
	elsif threshold == 0
		puts "Filtering step off"
	else 
		puts "Not valid filtering value, plese specify 0 to skip filtering and 1 to allow it"
		exit
	end 
	#puts "All the contigs with a ratio lower than #{threshold} will not be considered in the analysis"
end 

###############################################################
######Files
loc = "arabidopsis_datasets/#{dataset}"
file = "#{loc}/#{file}"
vcf_file = "#{loc}/snps.vcf"
fasta_file = "#{loc}/frags.fasta"
fasta_shuffle = "#{loc}/frags_shuffled.fasta"

hm_list = WriteIt.file_to_ints_array("#{loc}/hm_snps.txt") # Get SNP distributions
ht_list = WriteIt.file_to_ints_array("#{loc}/ht_snps.txt")
##############################################################
##############################################################

#Create arrays of ids of those fragments that contain  SNPs in the vcf file and hashes with the id and position eh
snp_data, hm, ht, frag_pos = Stuff.snps_in_vcf(vcf_file)

frag_pos_hm = frag_pos[:hom]
frag_pos_ht = frag_pos[:het]

#Create hashes for fragments ids and SNP position for the correct genome - reference 

dic_pos_hm =  Stuff.dic_id_pos(hm, hm_list)
dic_pos_ht =  Stuff.dic_id_pos(ht, ht_list)
 ########



##Create dictionaries with the id of the fragment as the key and the NUMBER of SNPs as value
dic_hm = Stuff.create_hash_number(hm)
dic_ht = Stuff.create_hash_number(ht)

##Create array with ordered fragments (fromf fasta_file) and from shuffled fragments (fasta_shuffle)
frags = ReformRatio.fasta_array(fasta_file)
frags_shuffled = ReformRatio.fasta_array(fasta_shuffle)

##From the previous array take ids and lengths and put them in 2 separate new arrays
ids_ok, lengths_ok, id_len_ok = ReformRatio.fasta_id_n_lengths(frags)
ids, lengths, id_len = ReformRatio.fasta_id_n_lengths(frags_shuffled)


genome_length = ReformRatio.genome_length(fasta_file)

##Define snps in hashes (fragment id as key and snp density as value). Create also lists 

##Assign the number of SNPs to each fragment in the shuffled list (hash)
##If a fragment does not have SNPs, the value assigned will be 0.
ok_hm, snps_hm = Stuff.define_snps(ids_ok, dic_hm)
ok_ht, snps_ht = Stuff.define_snps(ids_ok, dic_ht)
shuf_hm, shuf_snps_hm = Stuff.define_snps(ids, dic_hm)
shuf_ht, shuf_snps_ht = Stuff.define_snps(ids, dic_ht)

dic_ratios, ratios, ids_short, dic_ratios_inv  = Stuff.important_ratios(snps_hm, snps_ht, ids_ok, threshold, adjust) 
dic_ratios_shuf, ratios_shuf, ids_short_shuf, dic_ratios_inv_shuf = Stuff.important_ratios(shuf_snps_hm, shuf_snps_ht, ids, threshold, adjust) 

s_hm, s_snps_hm = Stuff.define_snps(ids_short, dic_hm)
s_ht, s_snps_ht = Stuff.define_snps(ids_short, dic_ht)

##if the ratio is not higher than a given threshold, eliminate the ids from the ids array
##and then eliminate those SNPs positions from the files, create new files.  
shuf_short_ids = Stuff.important_ids(ids_short, ids)
hm_sh = Stuff.important_pos(ids_short, dic_pos_hm)
ht_sh = Stuff.important_pos(ids_short, dic_pos_ht)


contigs_discarded = ids.length - ids_short.length
puts "#{contigs_discarded} contigs out of #{ids.length} discarded"


 shuf_hm, shu_snps_hm = Stuff.define_snps(shuf_short_ids, dic_hm)

#Define SNPs per fragment in the shuffled fasta array and then normalise the value of SNP density per fragment length

dic_shuf_hm_norm = Stuff.normalise_by_length(lengths, shuf_hm)

##Iteration: look for the minimum value in the array of values, that will be 0 (fragments without SNPs) and put the fragments 
#with this value in a list. Then, the list is cut by half and each half is added to a new array (right, that will be used 
#to reconstruct the right side of the distribution, and left, for the left side)
puts "\n"


perm_hm, perm_ratio, mut, hyp_positions = SDM.calling_SDM(dic_shuf_hm_norm, dic_ratios_inv_shuf, cross, dic_pos_hm)


puts "Hypothetical positions carrying the causal mutation #{hyp_positions}"



#Define SNPs in the recently ordered array of fragments.
dic_or_hm, snps_hm_or = Stuff.define_snps(perm_hm, dic_hm)
dic_or_ht, snps_ht_or = Stuff.define_snps(perm_hm, dic_ht)


###Calculate ratios and delete those equal to or lower than 1 so only the important contigs remain.
#dic_ratios, ratios = Stuff.important_ratios(snps_hm, snps_ht, ids_ok)
dic_expected_ratios, expected_ratios, exp_ids_short, exp_inv_ratios = Stuff.important_ratios(snps_hm_or, snps_ht_or, perm_hm, threshold, adjust)

#Take IDs, lenght and sequence from the shuffled fasta file and add them to the permutation array 

fasta_perm = Stuff.create_perm_fasta(perm_hm, frags_shuffled, ids)


#Create new fasta file with the ordered elements
File.open("#{loc}/frags_ordered_thres#{threshold}.fasta", "w+") do |f|
  fasta_perm.each { |element| f.puts(element) }
end

fasta_ordered = "arabidopsis_datasets/#{dataset}/frags_ordered_thres#{threshold}.fasta"
frags_ordered = ReformRatio.fasta_array(fasta_ordered)
ids_or, lengths_or, id_len_or = ReformRatio.fasta_id_n_lengths(frags_ordered)


###Calculate size of the group of fragments that have a high hm/ht ratio
contig_size = (genome_length/ids_ok.length).to_f
center = contig_size*(perm_hm.length) 
puts "The length of the group of contigs that have a high hm/ht ratio is #{center.to_i} bp"
puts "______________________"
 

Dir.mkdir("#{file}")

#Create arrays with the lists of SNP positions in the new ordered file.
het_snps, hom_snps = ReformRatio.perm_pos(frags_ordered, snp_data)
WriteIt::write_txt("#{file}/perm_hm", hom_snps) 
WriteIt::write_txt("#{file}/perm_ht", het_snps)
WriteIt::write_txt("#{file}/hm_snps_short", hm_sh) # save the SNP distributions for the best permutation in the generation
WriteIt::write_txt("#{file}/ht_snps_short", ht_sh)



Mutation::density_plots(contig_size, ratios, expected_ratios, hom_snps, het_snps, center, file, mut, frag_pos_hm) 



