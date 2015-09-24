#encoding: utf-8

require_relative 'lib/write_it'
require_relative 'lib/stuff'
require_relative 'lib/mutation'
require_relative 'lib/SDM'
require_relative 'lib/snp_dist'
require_relative 'lib/plot'
require_relative 'lib/ratio_filtering'
require_relative 'lib/output'
require_relative 'lib/reform_ratio'

require 'pp'
require 'benchmark'
require 'csv'
require 'yaml'
require 'fileutils'

if ARGV.length == 6
	dataset = ARGV[0]
	output_folder = ARGV[1]
	threshold = ARGV[2].to_i #degree of filtering:  100, 50, 10, 5
	adjust = ARGV[3]
	cross = ARGV[4]
	log_folder = ARGV[5]
	puts "Looking for SNPs in #{dataset}"
	puts "Output will be in #{dataset}/#{output_folder}"
	if threshold > 0
		puts "Filtering step on: #{threshold}% selected"
	elsif threshold == 0
		puts "Filtering step off. "
	else
		puts "Not valid filtering value, plese specify 0 to skip filtering or a positive integer to allow it"
		exit
	end
	puts "A factor of #{adjust} will be used to calculate the ratio"
else
	puts "Please specify a (1) dataset, a (2) name for the output folder, a (3) threshold to discard the contigs which a ratio below it and a (4) factor to calculate the ratio (1, 0.1, 0.01...) (5) kind of cross: back or out"
end

########Files#####################################################
loc = "#{dataset}"
output_folder  = "#{loc}/#{output_folder}"
vcf_file = "#{loc}/snps.vcf"
fasta_file = "#{loc}/frags.fasta"
fasta_shuffle = "#{loc}/frags_shuffled.fasta"
logfold  = "#{loc}/#{log_folder}"
FileUtils.mkdir_p "#{loc}/#{log_folder}"

hm_list = WriteIt.file_to_ints_array("#{loc}/hm_snps.txt") #create arrays for SNP densities
ht_list = WriteIt.file_to_ints_array("#{loc}/ht_snps.txt")
##############################################################
##############################################################

warn "I am here at step: 1\n"

####[1] Open VCF file
snp_data, hm, ht, frag_pos = Stuff.snps_in_vcf(vcf_file)
File.open("#{logfold}/1_1_snp_data.yml", "w") do |file|
  file.write snp_data.to_yaml
end
File.open("#{logfold}/1_2_frag_pos.yml", "w") do |file|
  file.write frag_pos.to_yaml
end
WriteIt::write_txt("#{logfold}/1_3_hm_positions", hm)
WriteIt::write_txt("#{logfold}/1_4_ht_positions", ht)

warn "I am here at step: 2\n"

##Hash with the fragments ids and the homozygous and heterozygous SNP positions
frag_pos_hm = frag_pos[:hom]
frag_pos_ht = frag_pos[:het]



##Hashes with fragments ids and SNP positions for the correctly ordered genome
dic_pos_hm =  Stuff.dic_id_pos(hm, hm_list)
dic_pos_ht =  Stuff.dic_id_pos(ht, ht_list)
File.open("#{logfold}/2_1_dic_pos_hm.yml", "w") do |file|
  file.write dic_pos_hm.to_yaml
end
File.open("#{logfold}/2_2_dic_pos_ht.yml", "w") do |file|
  file.write dic_pos_ht.to_yaml
end


##Create hashes with the id of the fragment as the key and the NUMBER of SNPs as value
dic_hm = Stuff.create_hash_number(hm)
dic_ht = Stuff.create_hash_number(ht)




####[2] Open FASTA files containing the ordered and unordered contigs
##Create array with ordered fragments (from fasta_file) and from shuffled fragments (fasta_shuffle)
frags = Stuff.fasta_array(fasta_file)
frags_shuffled = Stuff.fasta_array(fasta_shuffle)

warn "I am here at step: 3\n"

##From the previous array take ids and lengths and put them in 2 separate new arrays
ids_ok, lengths_ok, id_len_ok = Stuff.fasta_id_n_lengths(frags)
ids, lengths, id_len = Stuff.fasta_id_n_lengths(frags_shuffled)

genome_length = Stuff.genome_length(fasta_file)

average_contig = genome_length/ids.length

warn "I am here at step: 4\n"
warn "#{genome_length}\t#{average_contig}\n"

##Assign the number of SNPs to each fragment in the shuffled list (hash)
##If a fragment does not have SNPs, the value assigned will be 0.
ok_hm, snps_hm = Stuff.define_snps(ids_ok, dic_hm)
File.open("#{logfold}/4_1_ok_hm.yml", "w") do |file|
  file.write ok_hm.to_yaml
end
WriteIt::write_txt("#{logfold}/4_2_snps_hm", snps_hm)

ok_ht, snps_ht = Stuff.define_snps(ids_ok, dic_ht)
shuf_hm, shuf_snps_hm = Stuff.define_snps(ids, dic_hm)
shuf_ht, shuf_snps_ht = Stuff.define_snps(ids, dic_ht)


warn "I am here at step: 5\n"

dic_ratios, ratios, ids_short, dic_ratios_inv  = Ratio_filtering.important_ratios(snps_hm, snps_ht, ids_ok, threshold, adjust)
File.open("#{logfold}/5_1_dic_ratios.yml", "w") do |file|
  file.write dic_ratios.to_yaml
end
File.open("#{logfold}/5_2_dic_ratios_inv.yml", "w") do |file|
  file.write dic_ratios_inv.to_yaml
end
WriteIt::write_txt("#{logfold}/5_3_ratios", ratios)
WriteIt::write_txt("#{logfold}/5_4_ids_short", ids_short)

dic_ratios_shuf, ratios_shuf, ids_short_shuf, dic_ratios_inv_shuf = Ratio_filtering.important_ratios(shuf_snps_hm, shuf_snps_ht, ids, threshold, adjust)

##Redefine the arrays of SNPs after discarding the contigs that fell below the threshold provided. We refered to them as the "important contigs" and the SNPs on those are the "important positions"

warn "I am here at step: 6\n"

s_hm, s_snps_hm = Stuff.define_snps(ids_short, dic_hm)

shuf_short_ids = Ratio_filtering.important_ids(ids_short, ids)
hm_sh = Ratio_filtering.important_pos(ids_short, dic_pos_hm)
ht_sh = Ratio_filtering.important_pos(ids_short, dic_pos_ht)
WriteIt::write_txt("#{logfold}/6_1_shuf_short_ids", shuf_short_ids)
WriteIt::write_txt("#{logfold}/6_2_hm_sh", hm_sh)
WriteIt::write_txt("#{logfold}/6_3_ht_sh", ht_sh)

##Calculate how many contigs were discarded
shuf_hm, shu_snps_hm = Stuff.define_snps(shuf_short_ids, dic_hm)
File.open("#{logfold}/6_4_shuf_hm.yml", "w") do |file|
  file.write shuf_hm.to_yaml
end

warn "I am here at step: 7\n"

####[4] SDM
## Calculate scores (number of homozygous SNPs in each contig divided by fragment length)

dic_shuf_hm_norm = SDM.normalise_by_length(lengths, shuf_hm)
File.open("#{logfold}/7_1_dic_shuf_hm_norm.yml", "w") do |file|
  file.write dic_shuf_hm_norm.to_yaml
end

##Iteration: look for the minimum value in the array of values, that will be 0 (fragments without SNPs) and put the fragments
#with this value in a list. Then, the list is cut by half and each half is added to a new array (right, that will be used
#to reconstruct the right side of the distribution, and left, for the left side)


warn "I am here at step: 8\n"

perm_hm, perm_ratio, mut, hyp_positions = SDM.calling_SDM(dic_shuf_hm_norm, dic_ratios_inv_shuf, dic_pos_hm, cross, average_contig)

WriteIt::write_txt("#{logfold}/8_1_perm_hm", perm_hm)
WriteIt::write_txt("#{logfold}/8_2_perm_ratio", perm_ratio)
WriteIt::write_txt("#{logfold}/8_3_mut", mut)
WriteIt::write_txt("#{logfold}/8_4_hyp_positions", hyp_positions)

puts "Hypothetical positions carrying the causal mutation #{hyp_positions}"

# #Define SNPs in the r ordered array of fragments.
dic_or_hm, snps_hm_or = Stuff.define_snps(perm_hm, dic_hm)
dic_or_ht, snps_ht_or = Stuff.define_snps(perm_hm, dic_ht)


warn "I am here at step: 9\n"

#thres = 0
#pp thres
#Calculate ratios in the contig permutation obtained from SDM
dic_expected_ratios, expected_ratios, exp_ids_short, exp_inv_ratios = Ratio_filtering.important_ratios(snps_hm_or, snps_ht_or, perm_hm, threshold, adjust)

File.open("#{logfold}/9_1_dic_expected_ratios.yml", "w") do |file|
  file.write dic_expected_ratios.to_yaml
end
File.open("#{logfold}/9_2_exp_inv_ratios.yml", "w") do |file|
  file.write exp_inv_ratios.to_yaml
end
WriteIt::write_txt("#{logfold}/9_3_ratios", expected_ratios)
WriteIt::write_txt("#{logfold}/9_4_ids_short", exp_ids_short)

####[5] Outputs

#Create FASTA file for the contig permutation obtained from SDM
fasta_perm = Output.create_perm_fasta(perm_hm, frags_shuffled, ids)

File.open("#{loc}/frags_ordered_thres#{threshold}.fasta", "w+") do |f|
  fasta_perm.each { |element| f.puts(element) }
end

fasta_ordered = "#{loc}/frags_ordered_thres#{threshold}.fasta"
frags_ordered = Stuff.fasta_array(fasta_ordered)
ids_or, lengths_or, id_len_or = Stuff.fasta_id_n_lengths(frags_ordered)

contig_size = (genome_length/ids_ok.length).to_f
center = contig_size*(perm_hm.length)
puts "The length of the group of contigs that have a high Hom/het ratio is #{center.to_i} bp"
puts "______________________"

#Make Output directory
Dir.mkdir("#{output_folder}")

# #Create arrays with the  SNP positions in the new ordered file.
het_snps, hom_snps = ReformRatio.perm_pos(frags_ordered, snp_data)


WriteIt::write_txt("#{output_folder}/hyp_positions", hyp_positions)
WriteIt::write_txt("#{output_folder}/perm_hm", hom_snps)
WriteIt::write_txt("#{output_folder}/perm_ht", het_snps)
WriteIt::write_txt("#{output_folder}/hm_snps_short", hm_sh) # save the SNP distributions for the best permutation in the generation
WriteIt::write_txt("#{output_folder}/ht_snps_short", ht_sh)

####[6] Plots

##Plot expected vs SDM ratios, QQplots

Mutation::density_plots(contig_size, ratios, expected_ratios, hom_snps, het_snps, center, output_folder, mut, frag_pos_hm)
