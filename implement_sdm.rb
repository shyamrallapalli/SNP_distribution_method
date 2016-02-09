# encoding: utf-8

require_relative 'lib/file_rw'
require_relative 'lib/mutation'
require_relative 'lib/ratio_filter'
require_relative 'lib/fragments'
require_relative 'lib/vcf'
require_relative 'lib/pileup'

require 'pp'
require 'benchmark'
require 'csv'
require 'yaml'
require 'fileutils'
require 'bio-gngm'
require 'bio-samtools'

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
pars = YAML.load_file('./input_pars.yml')
fasta_shuffle = pars['fasta']
bamfile = pars['bam']
vcf_file = pars['vcf']
background = pars['background']
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
if background == ''
  var_pos = Vcf.get_vars(vcf_file)
else
  var_pos = Vcf.filtering(vcf_file, background)
end
File.open("#{log_folder}/1_4_frag_pos.yml", 'w') do |file|
  file.write var_pos.to_yaml
end

# ###[2] Open FASTA files containing the unordered contigs
# #Create a hash with shuffled fragments seq ids - values are lengths and sequences
# inseq = FileRW.fasta_parse(fasta_shuffle)
inseq = {}
temp  = Bio::DB::FastaLengthDB.new(:file => fasta_shuffle)
inseq[:len] = temp.instance_variable_get(:@seqs)
ids = inseq[:len].keys
genome_length = inseq[:len].values.inject { | sum, n | sum + n }
average_contig = genome_length / ids.length

input_frags = Vcf.varpos_aggregate(var_pos, inseq[:len], ids, adjust, 'no')
File.open("#{log_folder}/t_17_input_frags.yml", 'w') do |file|
  file.write input_frags.to_yaml
end

# ###[3]
# #ratio of homozygous to heterozygous snps per each fragment is calculated (shuffled)
ratios_hash = RatioFilter.selected_ratios(input_frags, threshold)
File.open("#{log_folder}/3_4_dic_ratios_inv_shuf.yml", 'w') do |file|
  file.write ratios_hash.to_yaml
end


# ###[4] SDM
# #Iteration: look for the minimum value in the array of values, that will be 0 (fragments without SNPs) and put the fragments
# with this value in a list. Then, the list is cut by half and each half is added to a new array (right, that will be used
# to reconstruct the right side of the distribution, and left, for the left side)
sdm_frags, mut_frags = Fragments.arrange(ratios_hash, cross, average_contig)
FileRW.write_txt("#{log_folder}/4_3_perm_ratio", sdm_frags)
FileRW.write_txt("#{log_folder}/4_4_mut", mut_frags)

sel_frags = Fragments.select_fragments(cross, ratios_hash, sdm_frags, adjust, threshold)
FileRW.write_txt("#{log_folder}/4_5_selected_frags", sel_frags)
mut_frags = sel_frags

sortfrags = Pileup.pick_frag_vars(bamfile,fasta_shuffle,sel_frags,input_frags)
File.open("#{log_folder}/4_6_sortfrags.yml", 'w') do |file|
  file.write sortfrags.to_yaml
end

# ###[5] Outputs
# Create FASTA file for the contig permutation obtained from SDM
filename = "ordered_frags_thres#{threshold}.fasta"
FileRW.write_order(sdm_frags, fasta_shuffle, filename)

region = average_contig * (sdm_frags.length)
puts "The length of the group of contigs that have a high Hom/het ratio is #{region.to_i} bp"
puts '______________________'

outcome = Vcf.varpos_aggregate(var_pos, inseq[:len], sdm_frags, adjust)


# ###[6] Plots

# #Plot expected vs SDM ratios, QQplots
# candidate_frag_vars = Mutation.get_candidates(mut_frags, var_pos[:hom])
File.open("#{output_folder}/mutation.txt", 'w+') do |f|
  # f.puts "The length of the group of contigs that form the peak of the distribution is #{region.to_i} bp"
  # f.puts "The mutation is likely to be found on the following contigs #{candidate_frag_vars}"
  f.puts "non_ref_ratio\tseq_id\tposition\tref_base\tcoverage\tbases\tbase_quals"
  sortfrags.keys.sort.reverse.each do | ratio_1 |
    sortfrags[ratio_1].each_key do | frag_1 |
      sortfrags[ratio_1][frag_1].each_key do | pos_1 |
        pileup = sortfrags[ratio_1][frag_1][pos_1].to_s
        f.puts "#{ratio_1}\t#{pileup}"
      end
    end
  end
end

Mutation.density_plot(outcome, average_contig.to_f, output_folder)


########## Test comparison inputs and analysis and comparison

if pars['test']
  # open file containing the ordered fragment ids and add them to an array
  frags_order = pars['frags_order']
  ids_ok = FileRW.to_array(frags_order)

  original = Vcf.varpos_aggregate(var_pos, inseq[:len], ids_ok, adjust)


  # #ratio of homozygous to heterozygous snps per each fragment is calculated (ordered)
  dic_ratios_inv  = RatioFilter.selected_ratios(original, threshold)
  File.open("#{log_folder}/t_10_dic_ratios_inv.yml", 'w') do |file|
    file.write dic_ratios_inv.to_yaml
  end

  # delete fragments which are discarded via ratio selection
  selected_frags = dic_ratios_inv.values.flatten
  original.each_key do | fragid |
    unless selected_frags.include?(fragid)
      original.delete(fragid)
    end
  end

  File.open("#{log_folder}/t_13_original.yml", 'w') do |file|
    file.write original.to_yaml
  end
  File.open("#{log_folder}/t_14_outcome.yml", 'w') do |file|
    file.write outcome.to_yaml
  end

  Mutation.compare_density(outcome, mut_frags, average_contig.to_f, genome_length, output_folder,  original)
end
