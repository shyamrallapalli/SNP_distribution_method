# encoding: utf-8

require_relative 'lib/fasta_handle'
require_relative 'lib/file_rw'
require_relative 'lib/mutation'
require_relative 'lib/ratio_filtering'
require_relative 'lib/SDM'
require_relative 'lib/stuff'
require_relative 'lib/vcf'

require 'pp'
require 'benchmark'
require 'csv'
require 'yaml'
require 'fileutils'

if ARGV.empty?
  puts "  Please specify a dataset directory,\n\
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
inseq, genome_length = FastaHandle.fasta_parse(fasta_shuffle)
ids = inseq[:len].keys
average_contig = genome_length / ids.length

input_frags = Vcf.varpos_aggregate(var_pos, inseq[:len], ids, adjust, "no")
File.open("#{log_folder}/t_17_input_frags.yml", "w") do |file|
  file.write input_frags.to_yaml
end

# ###[3]
# #ratio of homozygous to heterozygous snps per each fragment is calculated (shuffled)
ratios_shuf, input_frags, dic_ratios_inv_shuf = Ratio_filtering.selected_ratios(input_frags, threshold)
FileRW.write_txt("#{log_folder}/3_2_ratios_shuf", ratios_shuf)
File.open("#{log_folder}/3_4_dic_ratios_inv_shuf.yml", "w") do |file|
  file.write dic_ratios_inv_shuf.to_yaml
end

# #Calculate how many contigs were discarded
shuf_hm, shu_snps_hm = Vcf.define_snps(input_frags.keys, var_num[:hom])
File.open("#{log_folder}/3_6_shuf_hm.yml", "w") do |file|
  file.write shuf_hm.to_yaml
end
FileRW.write_txt("#{log_folder}/3_7_shu_snps_hm", shu_snps_hm)

# ###[4] SDM
# invert hash with ratios as keys and fragments as values
dic_shuf_hm_norm = Stuff.safe_invert(shuf_hm)
File.open("#{log_folder}/4_1_dic_shuf_hm_norm.yml", "w") do |file|
  file.write dic_shuf_hm_norm.to_yaml
end

# #Iteration: look for the minimum value in the array of values, that will be 0 (fragments without SNPs) and put the fragments
# with this value in a list. Then, the list is cut by half and each half is added to a new array (right, that will be used
# to reconstruct the right side of the distribution, and left, for the left side)
perm_hm, perm_ratio, mut, hyp_positions = SDM.calling_SDM(dic_shuf_hm_norm, dic_ratios_inv_shuf, var_pos[:hom], cross, average_contig)
FileRW.write_txt("#{log_folder}/4_2_perm_hm", perm_hm)
FileRW.write_txt("#{log_folder}/4_3_perm_ratio", perm_ratio)
FileRW.write_txt("#{log_folder}/4_4_mut", mut)
FileRW.write_txt("#{log_folder}/4_5_hyp_positions", hyp_positions)

puts "Hypothetical positions carrying the causal mutation #{hyp_positions}"
FileRW.write_txt("#{output_folder}/hyp_positions", hyp_positions)

# ###[5]
# thres = 0
# pp thres
# Calculate ratios in the contig permutation obtained from SDM
expected_ratios = []
perm_hm.each do | fragid |
  expected_ratios << input_frags[fragid][:ratio]
end
FileRW.write_txt("#{log_folder}/5_2_expected_ratios", expected_ratios)

# ###[5] Outputs
# Create FASTA file for the contig permutation obtained from SDM
fasta_perm = FastaHandle.create_perm_fasta(perm_hm, inseq[:seq])
File.open("ordered_frags_thres#{threshold}.fasta", 'w+') do |f|
  fasta_perm.each { |element| f.puts(element) }
end

region = average_contig * (perm_hm.length)
puts "The length of the group of contigs that have a high Hom/het ratio is #{region.to_i} bp"
puts '______________________'

########## Test comparison inputs and analysis

# #Open FASTA files containing the ordered contigs
# #from the array take ids and lengths
fasta_file = "frags.fasta"
inseq_ok, genome_length = FastaHandle.fasta_parse(fasta_file)
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

# #Assign the number of SNPs to each fragment in the ordered list (hash)
ok_hm, snps_hm = Vcf.define_snps(ids_ok, var_num[:hom])
ok_ht, snps_ht = Vcf.define_snps(ids_ok, var_num[:het])
File.open("#{log_folder}/t_03_ok_hm.yml", "w") do |file|
  file.write ok_hm.to_yaml
end
FileRW.write_txt("#{log_folder}/t_04_snps_hm", snps_hm)
File.open("#{log_folder}/t_05_ok_ht.yml", "w") do |file|
  file.write ok_ht.to_yaml
end
FileRW.write_txt("#{log_folder}/t_06_snps_ht", snps_ht)

# #ratio of homozygous to heterozygous snps per each fragment is calculated (ordered)
ratios, original, dic_ratios_inv  = Ratio_filtering.selected_ratios(original, threshold)
FileRW.write_txt("#{log_folder}/t_08_ratios", ratios)
File.open("#{log_folder}/t_10_dic_ratios_inv.yml", "w") do |file|
  file.write dic_ratios_inv.to_yaml
end


hm_sh = Ratio_filtering.important_pos(original.keys, origin_pos[:hom])
ht_sh = Ratio_filtering.important_pos(original.keys, origin_pos[:het])

FileRW.write_txt("#{output_folder}/hm_snps_short", hm_sh) # save the SNP distributions for the best permutation in the generation
FileRW.write_txt("#{output_folder}/ht_snps_short", ht_sh)


outcome = Vcf.varpos_aggregate(var_pos, inseq[:len], perm_hm, adjust)
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

candi_peak = Mutation.density_plots(average_contig.to_f, ratios, expected_ratios, hom_snps, het_snps, region, genome_length, output_folder, mut, var_pos[:hom], original, outcome)
puts "#{candi_peak}\n"
