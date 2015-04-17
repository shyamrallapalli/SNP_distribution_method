#encoding: utf-8
require_relative 'lib/reform_ratio'
require_relative 'lib/write_it'
require_relative 'lib/stuff'
require_relative 'lib/mutation'
require_relative 'lib/SDM'
require_relative 'lib/snp_dist'

require 'pp'
require 'benchmark'
require 'csv'

dataset = ARGV[0] 
perm = ARGV[1]
threshold = ARGV[2]
adjust = ARGV[3]

# puts "A factor of #{adjust} will be used to calculate the ratio"
# puts "All the contigs with a ratio lower than #{adjust} will not be considered in the analysis"

######Files
vcf_file = "arabidopsis_datasets/#{dataset}/snps.vcf"
fasta_file = "arabidopsis_datasets/#{dataset}/frags.fasta"
fasta_shuffle = "arabidopsis_datasets/#{dataset}/frags_shuffled.fasta"

hm_list = WriteIt.file_to_ints_array("arabidopsis_datasets/#{dataset}/hm_snps.txt") # Get SNP distributions
ht_list = WriteIt.file_to_ints_array("arabidopsis_datasets/#{dataset}/ht_snps.txt")


#Create list of fragment ids that contain  SNPs from vcf file and hash of hashes (
#frag_pos) that contains type of SNP as hard key, fragment ids as key and array of SNP position as value. 
snp_data, hm, ht, dic_pos_hm, dic_pos_ht = Stuff.snps_in_vcf(vcf_file)
dic_pos_hm =  Stuff.dic_id_pos(hm, hm_list)
dic_pos_ht =  Stuff.dic_id_pos(ht, ht_list)

#Create hashes for fragments ids and SNP position
 
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

#ratios
dic_ratios, ratios, ids_short, ratios_by_len, dic_ratios_by_len = Stuff.important_ratios(snps_hm, snps_ht, ids_ok, id_len_ok, threshold, adjust)

ids_short.flatten!

s_hm, s_snps_hm = Stuff.define_snps(ids_short, dic_hm)
s_ht, s_snps_ht = Stuff.define_snps(ids_short, dic_ht)

##if the ratio is not higher than a given threshold, eliminate the ids from the ids array
##and then eliminate those SNPs positions from the files, create new files.  
shuf_short_ids = Stuff.important_ids(ids_short, ids)
hm_sh = Stuff.important_pos(ids_short, dic_pos_hm)
ht_sh = Stuff.important_pos(ids_short, dic_pos_ht)


shuf_hm, shu_snps_hm = Stuff.define_snps(shuf_short_ids, dic_hm)


#Define SNPs per fragment in the shuffled fasta array and then normalise the value of SNP density per fragment length

dic_shuf_hm_norm = Stuff.normalise_by_length(lengths, shuf_hm)

##Iteration: look for the minimum value in the array of values, that will be 0 (fragments without SNPs) and put the fragments 
#with this value in a list. Then, the list is cut by half and each half is added to a new array (right, that will be used 
#to reconstruct the right side of the distribution, and left, for the left side)

perm_hm = SDM.sorting(dic_shuf_hm_norm)

##Measuree time of SDM. Eventually add time needed for the remaining steps until we define the mutation
 Benchmark.bm do |b|
    b.report {10.times do ; perm_hm = SDM.sorting(dic_shuf_hm_norm);  end}
end

#Define SNPs in the recently ordered array of fragments.
dic_or_hm, snps_hm_or = Stuff.define_snps(perm_hm, dic_hm)
dic_or_ht, snps_ht_or = Stuff.define_snps(perm_hm, dic_ht)

###Calculate ratios and delete those equal to or lower than 1 so only the important contigs remain.
#dic_ratios, ratios = Stuff.important_ratios(snps_hm, snps_ht, ids_ok)
dic_expected_ratios, expected_ratios, exp_rat_by_len = Stuff.important_ratios(snps_hm_or, snps_ht_or, perm_hm, id_len, threshold, adjust)

#Take IDs, lenght and sequence from the shuffled fasta file and add them to the permutation array 

fasta_perm = Stuff.create_perm_fasta(perm_hm, frags_shuffled, ids)

#Create new fasta file with the ordered elements
File.open("arabidopsis_datasets/#{dataset}/frags_ordered#{threshold}.fasta", "w+") do |f|
  fasta_perm.each { |element| f.puts(element) }
end

fasta_ordered = "arabidopsis_datasets/#{dataset}/frags_ordered#{threshold}.fasta"
frags_ordered = ReformRatio.fasta_array(fasta_ordered)

#Create arrays with the lists of SNP positions in the new ordered file.
het_snps, hom_snps = ReformRatio.perm_pos(frags_ordered, snp_data)

###Calculate size of the group of fragments that have a high hm/ht ratio
contig_size = (genome_length/ids_ok.length).to_f
center = contig_size*(perm_hm.length)
puts "The length of the group of contigs that have a high hm/ht ratio is #{center.to_i} bp"



causal, candidate, percent = Mutation.define(hm_sh, ht_sh, hom_snps, het_snps, genome_length, ratios, expected_ratios)


Dir.mkdir("arabidopsis_datasets/#{dataset}/#{perm}")
Dir.chdir("arabidopsis_datasets/#{dataset}/#{perm}") do
	WriteIt::write_txt("perm_hm", hom_snps) # save the SNP distributions for the best permutation in the generation
	WriteIt::write_txt("perm_ht", het_snps)
  WriteIt::write_txt("hm_snps_short", hm_sh) # save the SNP distributions for the best permutation in the generation
  WriteIt::write_txt("ht_snps_short", ht_sh)
	File.open("mutation.txt", "w+") do |f|
		f.puts "The length of the group of contigs that form the peak of the distribution is #{center.to_i} bp"
		f.puts "Location of causal mutation in correctly ordered genome: #{causal}"
		f.puts "Candidate SNP position in permutation: #{candidate}"
		f.puts "Shift #{percent} %"
	end
end

distribution_plots = Mutation.distribution_plot(center, ratios, expected_ratios, dataset, perm)

csv = "arabidopsis_datasets/#{dataset}/ratio_positions#{threshold}_#{adjust}.csv"
Stuff::csv_pos_ratio(csv, dic_pos_hm, dic_ratios)
