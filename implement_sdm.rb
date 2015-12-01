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
snp_data, hm, ht, frag_pos = Vcf.snps_in_vcf(vcf_file)
File.open("#{log_folder}/1_1_snp_data.yml", "w") do |file|
  file.write snp_data.to_yaml
end
FileRW.write_txt("#{log_folder}/1_2_hm_positions", hm)
FileRW.write_txt("#{log_folder}/1_3_ht_positions", ht)
File.open("#{log_folder}/1_4_frag_pos.yml", "w") do |file|
  file.write frag_pos.to_yaml
end

# #Create hashes with the id of the fragment as the key and the NUMBER of SNPs as value
dic_hm = Stuff.create_hash_number(hm)
dic_ht = Stuff.create_hash_number(ht)
File.open("#{log_folder}/1_5_dic_hm.yml", "w") do |file|
  file.write dic_hm.to_yaml
end
File.open("#{log_folder}/1_6_dic_ht.yml", "w") do |file|
  file.write dic_ht.to_yaml
end

# ###[2] Open FASTA files containing the unordered contigs
# #Create a hash with shuffled fragments seq ids - values are lengths and sequences
inseq, genome_length = FastaHandle.fasta_parse(fasta_shuffle)
ids = inseq[:len].keys
average_contig = genome_length / ids.length

# #Assign the number of SNPs to each fragment in the shuffled list (hash)
# #If a fragment does not have SNPs, the value assigned will be 0.
shuf_hm, shuf_snps_hm = Vcf.define_snps(ids, dic_hm)
shuf_ht, shuf_snps_ht = Vcf.define_snps(ids, dic_ht)
File.open("#{log_folder}/2_1_shuf_hm.yml", "w") do |file|
  file.write shuf_hm.to_yaml
end
FileRW.write_txt("#{log_folder}/2_2_shuf_snps_hm", shuf_snps_hm)

File.open("#{log_folder}/2_3_shuf_ht.yml", "w") do |file|
  file.write shuf_ht.to_yaml
end
FileRW.write_txt("#{log_folder}/2_4_shuf_snps_ht", shuf_snps_ht)

# ###[3]
# #ratio of homozygous to heterozygous snps per each fragment is calculated (shuffled)
dic_ratios_shuf, ratios_shuf, ids_short_shuf, dic_ratios_inv_shuf = Ratio_filtering.important_ratios(shuf_snps_hm, shuf_snps_ht, ids, threshold, adjust)
File.open("#{log_folder}/3_1_dic_ratios_shuf.yml", "w") do |file|
  file.write dic_ratios_shuf.to_yaml
end
FileRW.write_txt("#{log_folder}/3_2_ratios_shuf", ratios_shuf)
FileRW.write_txt("#{log_folder}/3_3_ids_short_shuf", ids_short_shuf)
File.open("#{log_folder}/3_4_dic_ratios_inv_shuf.yml", "w") do |file|
  file.write dic_ratios_inv_shuf.to_yaml
end

# #Redefine the arrays of SNPs after discarding the contigs that fell below the threshold provided.
# #We refered to them as the "important contigs" and the SNPs on those are the "important positions"
shuf_short_ids = Ratio_filtering.important_ids(ids_short_shuf, ids)
FileRW.write_txt("#{log_folder}/3_5_shuf_short_ids", shuf_short_ids)

# #Calculate how many contigs were discarded
shuf_hm, shu_snps_hm = Vcf.define_snps(shuf_short_ids, dic_hm)
File.open("#{log_folder}/3_6_shuf_hm.yml", "w") do |file|
  file.write shuf_hm.to_yaml
end
FileRW.write_txt("#{log_folder}/3_7_shu_snps_hm", shu_snps_hm)

# ###[4] SDM
## Calculate scores (number of homozygous SNPs in each contig divided by fragment length)
dic_shuf_hm_norm = SDM.normalise_by_length(inseq[:len], shuf_hm)
File.open("#{log_folder}/4_1_dic_shuf_hm_norm.yml", "w") do |file|
  file.write dic_shuf_hm_norm.to_yaml
end

# #Iteration: look for the minimum value in the array of values, that will be 0 (fragments without SNPs) and put the fragments
# with this value in a list. Then, the list is cut by half and each half is added to a new array (right, that will be used
# to reconstruct the right side of the distribution, and left, for the left side)
perm_hm, perm_ratio, mut, hyp_positions = SDM.calling_SDM(dic_shuf_hm_norm, dic_ratios_inv_shuf, frag_pos[:hom], cross, average_contig)
FileRW.write_txt("#{log_folder}/4_2_perm_hm", perm_hm)
FileRW.write_txt("#{log_folder}/4_3_perm_ratio", perm_ratio)
FileRW.write_txt("#{log_folder}/4_4_mut", mut)
FileRW.write_txt("#{log_folder}/4_5_hyp_positions", hyp_positions)

puts "Hypothetical positions carrying the causal mutation #{hyp_positions}"
FileRW.write_txt("#{output_folder}/hyp_positions", hyp_positions)

# #Define SNPs in the r ordered array of fragments.
dic_or_hm, snps_hm_or = Vcf.define_snps(perm_hm, dic_hm)
dic_or_ht, snps_ht_or = Vcf.define_snps(perm_hm, dic_ht)
File.open("#{log_folder}/4_6_dic_or_hm.yml", "w") do |file|
  file.write dic_or_hm.to_yaml
end
FileRW.write_txt("#{log_folder}/4_7_snps_hm_or", snps_hm_or)

File.open("#{log_folder}/4_8_dic_or_ht.yml", "w") do |file|
  file.write dic_or_ht.to_yaml
end
FileRW.write_txt("#{log_folder}/4_9_snps_ht_or", snps_ht_or)

# ###[5]
# thres = 0
# pp thres
# Calculate ratios in the contig permutation obtained from SDM
dic_expected_ratios, expected_ratios, exp_ids_short, exp_inv_ratios = Ratio_filtering.important_ratios(snps_hm_or, snps_ht_or, perm_hm, threshold, adjust)
File.open("#{log_folder}/5_1_dic_expected_ratios.yml", "w") do |file|
  file.write dic_expected_ratios.to_yaml
end
FileRW.write_txt("#{log_folder}/5_2_expected_ratios", expected_ratios)
FileRW.write_txt("#{log_folder}/5_3_exp_ids_short", exp_ids_short)
File.open("#{log_folder}/5_4_exp_inv_ratios.yml", "w") do |file|
  file.write exp_inv_ratios.to_yaml
end

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

### Ordered genome and variants in ordered genome
fasta_file = "frags.fasta"
hm_list = FileRW.to_array_int("hm_snps.txt") # create arrays for SNP densities
ht_list = FileRW.to_array_int("ht_snps.txt")

# #Hashes with fragments ids and SNP positions for the correctly ordered genome
dic_pos_hm =  Vcf.dic_id_pos(hm, hm_list)
dic_pos_ht =  Vcf.dic_id_pos(ht, ht_list)
File.open("#{log_folder}/t_01_dic_pos_hm.yml", "w") do |file|
  file.write dic_pos_hm.to_yaml
end
File.open("#{log_folder}/t_02_dic_pos_ht.yml", "w") do |file|
  file.write dic_pos_ht.to_yaml
end

# #Open FASTA files containing the ordered contigs
# #from the array take ids and lengths
inseq_ok, genome_length = FastaHandle.fasta_parse(fasta_file)
ids_ok = inseq_ok[:len].keys

# #Assign the number of SNPs to each fragment in the ordered list (hash)
ok_hm, snps_hm = Vcf.define_snps(ids_ok, dic_hm)
ok_ht, snps_ht = Vcf.define_snps(ids_ok, dic_ht)
File.open("#{log_folder}/t_03_ok_hm.yml", "w") do |file|
  file.write ok_hm.to_yaml
end
FileRW.write_txt("#{log_folder}/t_04_snps_hm", snps_hm)
File.open("#{log_folder}/t_05_ok_ht.yml", "w") do |file|
  file.write ok_ht.to_yaml
end
FileRW.write_txt("#{log_folder}/t_06_snps_ht", snps_ht)

# #ratio of homozygous to heterozygous snps per each fragment is calculated (ordered)
dic_ratios, ratios, ids_short, dic_ratios_inv  = Ratio_filtering.important_ratios(snps_hm, snps_ht, ids_ok, threshold, adjust)
File.open("#{log_folder}/t_07_dic_ratios.yml", "w") do |file|
  file.write dic_ratios.to_yaml
end
FileRW.write_txt("#{log_folder}/t_08_ratios", ratios)
FileRW.write_txt("#{log_folder}/t_09_ids_short", ids_short)
File.open("#{log_folder}/t_10_dic_ratios_inv.yml", "w") do |file|
  file.write dic_ratios_inv.to_yaml
end


s_hm, s_snps_hm = Vcf.define_snps(ids_short, dic_hm)
File.open("#{log_folder}/t_11_s_hm.yml", "w") do |file|
  file.write s_hm.to_yaml
end
FileRW.write_txt("#{log_folder}/t_12_s_snps_hm", s_snps_hm)

hm_sh = Ratio_filtering.important_pos(ids_short, dic_pos_hm)
ht_sh = Ratio_filtering.important_pos(ids_short, dic_pos_ht)

FileRW.write_txt("#{output_folder}/hm_snps_short", hm_sh) # save the SNP distributions for the best permutation in the generation
FileRW.write_txt("#{output_folder}/ht_snps_short", ht_sh)



original = Vcf.varpos_aggregate(frag_pos, inseq[:len], ids_ok, adjust)
outcome = Vcf.varpos_aggregate(frag_pos, inseq[:len], perm_hm, adjust)
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

candi_peak = Mutation.density_plots(average_contig.to_f, ratios, expected_ratios, hom_snps, het_snps, region, genome_length, output_folder, mut, frag_pos[:hom], original, outcome)
puts "#{candi_peak}\n"
