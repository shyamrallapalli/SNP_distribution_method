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
pars = YAML.load_file('./input_pars.yml')
fasta_shuffle = pars['fasta']
mut_vcf = pars['mut_vcf']
mut_bam = pars['mut_bam']
bg_vcf = pars['bg_vcf']
bg_bam = pars['bg_bam']
adjust = pars['ratio_adj'].to_f
threshold = pars['filter'].to_i
log_folder = "#{pars['logdir']}_#{threshold}_#{adjust}"

# Make Output directory
FileUtils.mkdir_p "#{log_folder}"

# ###[1] Open VCF file
vars_mut = Polyploid.vars_in_file(mut_vcf, mut_bam, fasta_shuffle)
vars_bg = Polyploid.vars_in_file(bg_vcf, bg_bam, fasta_shuffle)
var_pos = Polyploid.filter_vars(vars_mut, vars_bg)

File.open("#{log_folder}/1_4_frag_pos.yml", 'w') do |file|
  file.write var_pos.to_yaml
end
