# encoding: utf-8

require_relative 'lib/file_rw'
require_relative 'lib/mutation'
require_relative 'lib/ratio_filter'
require_relative 'lib/fragments'
require_relative 'lib/vcf'
require_relative 'lib/pileup'
require_relative 'lib/poly_vcf'

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
pars = YAML.load_file("#{loc}/input_pars.yml")
fasta_shuffle = pars['fasta']
mut_bam = pars['mut_bam']
bg_bam = pars['bg_bam']

mut_vcf = pars['mut_vcf']
bg_vcf = pars['bg_vcf']

mut_pileup = pars['mut_pileup']
bg_pileup = pars['bg_pileup']

adjust = pars['ratio_adj'].to_f
threshold = pars['filter'].to_i
cross = pars['cross']

# Make Output directory
output_folder = "#{pars['outdir']}_#{threshold}_#{adjust}"
log_folder = "#{pars['logdir']}_#{threshold}_#{adjust}"
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
var_pos = ''
if bg_vcf == '' and bg_pileup == ''
  if mut_pileup != ''
    # do something with only mut pileup file
  elsif mut_vcf != ''
    var_pos = Vcf.get_vars(mut_vcf, :ht_low => 0.25, :ht_high => 0.75)
  else
    warn 'nothing to do here, provide me pileup or vcf file'
    exit
  end
else
  if mut_pileup != '' and bg_pileup != ''
    vars_bg = Polyploid.vars_in_pileup(bg_pileup)
    var_pos = Polyploid.filter_vars(mut_pileup, vars_bg)
  elsif mut_vcf != '' and bg_vcf != ''
    var_pos = Vcf.filtering(mut_vcf, bg_vcf)
  else
    warn 'nothing to do here, provide me with pileup or vcf file for both mutant and background'
    exit
  end
end

File.open("#{log_folder}/1_4_frag_pos.yml", 'w') do |file|
  file.write var_pos.to_yaml
end

# ###[2] Open FASTA files containing the unordered contigs
# Create a hash with shuffled fragments seq ids - values are lengths and sequences
# inseq = FileRW.fasta_parse(fasta_shuffle)
temp  = Bio::DB::FastaLengthDB.new(:file => fasta_shuffle)
inseq_len = temp.instance_variable_get(:@seqs)
ids = inseq_len.keys
genome_length = inseq_len.values.inject { | sum, n | sum + n }
average_contig = genome_length / ids.length

var_pos_orig = FileRW.deep_copy_hash(var_pos)
input_frags = ''
sdm_frags = ''
new_mut_frags = ''
repeat = 1
while repeat < 3 do
  input_frags = Vcf.varpos_aggregate(var_pos, inseq_len, ids, adjust, :cumulate => false)
  File.open("#{log_folder}/#{repeat}_t_17_input_frags.yml", 'w') do |file|
    file.write input_frags.to_yaml
  end

  # ###[3]
  # #ratio of homozygous to heterozygous snps per each fragment is calculated (shuffled)
  ratios_hash = RatioFilter.selected_ratios(input_frags, adjust, :only_frag_with_vars => true)
  File.open("#{log_folder}/#{repeat}_3_4_dic_ratios_inv_shuf.yml", 'w') do |file|
    file.write ratios_hash.to_yaml
  end


  # ###[4] SDM
  # #Iteration: look for the minimum value in the array of values, that will be 0 (fragments without SNPs) and put the fragments
  # with this value in a list. Then, the list is cut by half and each half is added to a new array (right, that will be used
  # to reconstruct the right side of the distribution, and left, for the left side)
  sdm_frags = Fragments.arrange(ratios_hash, input_frags)
  FileRW.write_txt("#{log_folder}/#{repeat}_4_3_perm_ratio", sdm_frags)

  sel_frags = Fragments.select_fragments(ratios_hash, sdm_frags, adjust, :cross => cross, :filter_out_low_hmes => true)
  FileRW.write_txt("#{log_folder}/#{repeat}_4_5_selected_frags", sel_frags)

  # remove below par snps and update var_pos hash for second iteration
  sortfrags, var_pos = Pileup.pick_frag_vars(mut_bam,fasta_shuffle,sel_frags,input_frags,var_pos, :bgbam => bg_bam,
                                                 :bq => 15, :mq => 20, :min_depth => 6, :min_non_ref_count => 3)
  File.open("#{log_folder}/#{repeat}_4_6_sortfrags.yml", 'w') do |file|
    file.write sortfrags.to_yaml
  end

  new_mut_frags = []
  mut_frags_pos = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
  File.open("#{output_folder}/#{repeat}_mutation.txt", 'w+') do |f|
    f.puts "HMEscore\tAlleleFreq\tseq_id\tposition\tref_base\tcoverage\tbases\tbase_quals"
    sortfrags.keys.sort.reverse.each do | ratio_1 |
      if ratio_1 >= 0.75
        sortfrags[ratio_1].each_key do | frag_1 |
          new_mut_frags << frag_1
          sortfrags[ratio_1][frag_1].each_key do | pos_1 |
            mut_frags_pos[frag_1][pos_1] = 1
            pileup = sortfrags[ratio_1][frag_1][pos_1].to_s
            f.puts "#{input_frags[frag_1][:ratio]}\t#{ratio_1}\t#{pileup}"
          end
        end
      end
    end
  end
  new_mut_frags.uniq!
  FileRW.write_txt("#{log_folder}/#{repeat}_6_7_final_selected_frags", new_mut_frags)

  # delete positions in the selected fragments that didn't pass the filtering
  mut_frags_pos.each_key do | fragment |
    selected_pos = mut_frags_pos[fragment].keys
    var_pos[:hom][fragment].each do | varpos |
      unless selected_pos.include?(varpos)
        warn "#{fragment}\t#{varpos}\n"
        var_pos[:hom][fragment].delete(varpos)
      end
    end
  end
  repeat += 1
end

# ###[5] Outputs
# Create FASTA file for the contig permutation obtained from SDM
filename = "ordered_frags_thres#{threshold}.fasta"
FileRW.write_order(sdm_frags, fasta_shuffle, filename)

region = average_contig * (sdm_frags.length)
puts "The length of the group of contigs that have a high Hom/het ratio is #{region.to_i} bp"
puts '______________________'



# ###[6] Plots

# do the pos aggregation after trimming filtered positions
outcome = Vcf.varpos_aggregate(var_pos, inseq_len, sdm_frags, adjust)

File.open("#{output_folder}/outcome_table.txt", 'w+') do |f|
  f.puts "Frag_id\tHMEscore\tnum hm\tnum ht\tlength\tcumulative length\thm positions"
  outcome.each_key do |key|
    f.puts "#{key}\t#{outcome[key][:ratio]}\t#{outcome[key][:hm]}\t#{outcome[key][:ht]}\t#{outcome[key][:len]}\t#{outcome[key][:cum_len]}\t#{outcome[key][:hm_pos]}"
  end
end

Mutation.density_plot(outcome, output_folder)

# Create FASTA file for the fragments selected to host mutation
filename = "selected_frags_thres#{threshold}.fasta"
FileRW.write_order(new_mut_frags, fasta_shuffle, filename)

########## Test comparison inputs and analysis and comparison

if pars['test']
  # open file containing the ordered fragment ids and add them to an array
  frags_order = pars['frags_order']
  ids_ok = FileRW.to_array(frags_order)

  original = Vcf.varpos_aggregate(var_pos_orig, inseq_len, ids_ok, adjust)


  # #ratio of homozygous to heterozygous snps per each fragment is calculated (ordered)
  dic_ratios_inv  = RatioFilter.selected_ratios(original, adjust, :only_frag_with_vars => false)
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

  Mutation.compare_density(outcome, new_mut_frags, average_contig.to_f, genome_length, output_folder,  original)
end
