# encoding: utf-8

require_relative 'lib/file_rw'
require_relative 'lib/mutation'
require_relative 'lib/ratio_filtering'
require_relative 'lib/SDM'
require_relative 'lib/vcf'

require 'pp'
require 'benchmark'
require 'csv'
require 'yaml'
require 'fileutils'

if ARGV.empty?
  puts "  Please specify a data set directory,\n\
  and place in it a \"input_pars.yml\" file with\n\
  following details\n\
  (1) name of input sequences (fasta) and variant (vcf) files\n\
  (2) a name for the output folder,\n\
  (3) a threshold for ratio to discard the contigs\n\
  (4) a adjusting factor to calculate the ratio (1, 0.5, 0.1, 0.01...)\n\
  (5) kind of cross: back or out and any additional log folder details"
  exit
end

loc = ARGV[0]
FileUtils.cd("#{loc}")
#### Inputs
### sequences and variants from the shuffled genome
# filter parameter are to be read from a "input_pars.yml" file in the current folder
pars = YAML.load_file("./input_pars.yml")
fasta_shuffle = pars['fasta']
vcf_file = pars['vcf']
adjust = pars['ratio_adj'].to_f
threshold = pars['filter'].to_i
cross = pars['cross']
output_folder = "#{pars['outdir']}_#{threshold}_#{adjust}"
log_folder = "#{pars['logdir']}_#{threshold}_#{adjust}"

# Make Output directory
FileUtils.mkdir_p "#{output_folder}"
FileUtils.mkdir_p "#{log_folder}"

if threshold > 0
  puts "Filtering step on: #{threshold}% selected"
elsif threshold == 0
  puts 'Filtering step off. '
else
  puts 'Not valid filtering value, plese specify 0 to skip filtering or a positive integer to allow it'
  exit
end

puts "A factor of #{adjust} will be used to calculate the ratio"


# ###[1] Open VCF file
var_num, var_pos = Vcf.snps_in_vcf(vcf_file)
File.open("#{log_folder}/1_4_frag_pos.yml", "w") do |file|
  file.write var_pos.to_yaml
end

File.open("#{log_folder}/1_5_dic_hm.yml", "w") do |file|
  file.write var_num[:hom].to_yaml
end
File.open("#{log_folder}/1_6_dic_ht.yml", "w") do |file|
  file.write var_num[:het].to_yaml
end

# ###[2] Open FASTA files containing the unordered contigs
# #Create a hash with shuffled fragments seq ids - values are lengths and sequences
inseq, genome_length = FileRW.fasta_parse(fasta_shuffle)
ids = inseq[:len].keys
average_contig = genome_length / ids.length

input_frags = Vcf.varpos_aggregate(var_pos, inseq[:len], ids, adjust, 'no')
File.open("#{log_folder}/t_17_input_frags.yml", "w") do |file|
  file.write input_frags.to_yaml
end

# ###[3]
# #ratio of homozygous to heterozygous snps per each fragment is calculated (shuffled)
input_frags, ratios_hash = Ratio_filtering.selected_ratios(input_frags, threshold)
File.open("#{log_folder}/3_4_dic_ratios_inv_shuf.yml", "w") do |file|
  file.write ratios_hash.to_yaml
end


# ###[4] SDM
# #Iteration: look for the minimum value in the array of values, that will be 0 (fragments without SNPs) and put the fragments
# with this value in a list. Then, the list is cut by half and each half is added to a new array (right, that will be used
# to reconstruct the right side of the distribution, and left, for the left side)
sdm_frags, mut_frags, hyp_positions = SDM.arrange(ratios_hash, var_pos[:hom], cross, average_contig)
FileRW.write_txt("#{log_folder}/4_3_perm_ratio", sdm_frags)
FileRW.write_txt("#{log_folder}/4_4_mut", mut_frags)
FileRW.write_txt("#{log_folder}/4_5_hyp_positions", hyp_positions)

puts "Hypothetical positions carrying the causal mutation #{hyp_positions}"
FileRW.write_txt("#{output_folder}/hyp_positions", hyp_positions)

# ###[5]
# Calculate ratios in the contig permutation obtained from SDM
sdm_ratios = Ratio_filtering.get_ratios(input_frags, sdm_frags)
FileRW.write_txt("#{log_folder}/5_2_expected_ratios", sdm_ratios)

# ###[5] Outputs
# Create FASTA file for the contig permutation obtained from SDM
filename = "ordered_frags_thres#{threshold}.fasta"
FileRW.write_order(sdm_frags, inseq[:seq], filename)

region = average_contig * (sdm_frags.length)
puts "The length of the group of contigs that have a high Hom/het ratio is #{region.to_i} bp"
puts '______________________'

########## Test comparison inputs and analysis

# #Open FASTA files containing the ordered contigs
# #from the array take ids and lengths
fasta_file = "frags.fasta"
inseq_ok, genome_length = FileRW.fasta_parse(fasta_file)
ids_ok = inseq_ok[:len].keys

original = Vcf.varpos_aggregate(var_pos, inseq[:len], ids_ok, adjust)

# #Hashes with fragments ids and SNP positions for the correctly ordered genome
origin_pos = Vcf.frag_positions(original)
File.open("#{log_folder}/t_01_dic_pos_hm.yml", "w") do |file|
  file.write origin_pos[:hom].to_yaml
end
File.open("#{log_folder}/t_02_dic_pos_ht.yml", "w") do |file|
  file.write origin_pos[:het].to_yaml
end


# #ratio of homozygous to heterozygous snps per each fragment is calculated (ordered)
original, dic_ratios_inv  = Ratio_filtering.selected_ratios(original, threshold)
ratios = Ratio_filtering.get_ratios(original, original.keys)
FileRW.write_txt("#{log_folder}/t_08_ratios", ratios)
File.open("#{log_folder}/t_10_dic_ratios_inv.yml", "w") do |file|
  file.write dic_ratios_inv.to_yaml
end


outcome = Vcf.varpos_aggregate(var_pos, inseq[:len], sdm_frags, adjust)
File.open("#{log_folder}/t_13_original.yml", "w") do |file|
  file.write original.to_yaml
end
File.open("#{log_folder}/t_14_outcome.yml", "w") do |file|
  file.write outcome.to_yaml
end

out_original = File.open("#{log_folder}/t_15_original.txt", "w")
out_original.puts "Frag\thm\tht\tratio\tlen\thm_pos\tht_pos\n"
original.each_key { |key|
  hash = original[key]
  out_original.puts "#{key}\t#{hash[:hm]}\t#{hash[:ht]}\t#{hash[:ratio]}\t#{hash[:len]}\t#{hash[:hm_pos].join(" ")}\t#{hash[:ht_pos].join(" ")}\n"
}

out_outcome = File.open("#{log_folder}/t_16_outcome.txt", "w")
out_outcome.puts "Frag\thm\tht\tratio\tlen\thm_pos\tht_pos\n"
outcome.each_key { |key|
  hash = outcome[key]
  out_outcome.puts "#{key}\t#{hash[:hm]}\t#{hash[:ht]}\t#{hash[:ratio]}\t#{hash[:len]}\t#{hash[:hm_pos].join(" ")}\t#{hash[:ht_pos].join(" ")}\n"
}

# #Create arrays with the  SNP positions in the new ordered file.
hom_snps, het_snps = Vcf.varpositions(outcome)

FileRW.write_txt("#{output_folder}/perm_hm", hom_snps)
FileRW.write_txt("#{output_folder}/perm_ht", het_snps)

# ###[6] Plots

# #Plot expected vs SDM ratios, QQplots

Mutation.density_plots(average_contig.to_f, ratios, sdm_ratios, hom_snps, het_snps, region, genome_length, output_folder, mut_frags, var_pos[:hom], original, outcome)
