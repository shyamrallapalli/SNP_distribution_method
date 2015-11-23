#encoding: utf-8
require_relative 'locate_mutation'
require_relative 'snp_dist'
require_relative 'reform_ratio'
require_relative 'plot'
require_relative 'write_it'
require 'simple_stats'
require 'pp'

class Mutation
	def self.test_genomes_define(hm, ht, perm_hm, perm_ht, genome_length, ratios, expected_ratios)
		hm.pop
		# div = 100
		n = 2000
		hyp = SNPdist.hyp_snps(ratios, genome_length)
		peak =  LocateMutation.find_peak(hyp, n) # Find the peak in the approximated (hypothetical SNP) distribution
		causal = LocateMutation.closest_snp(peak, hm)
		perm_hyp = SNPdist.hyp_snps(expected_ratios, genome_length)

		perm_peak = LocateMutation.find_peak(perm_hyp, n)
		candidate = LocateMutation.closest_snp(perm_peak, perm_hm)
		normalised = (candidate - causal).abs
		percent = (normalised*100)/genome_length.to_f
		return causal, candidate, percent
	end

	def self.test_genomes_distribution_plot(genome_length, ratios, expected_ratios, dataset, perm)
		hm, ht, hyp, ylim_hm, ylim_ht, ylim_hyp = [],[],[],[],[],[]
		Dir.chdir(File.join(Dir.home, "SNP_distribution_method/Small_genomes/#{perm}")) do
			hom_snps = WriteIt.file_to_ints_array("hm_snps_short.txt")
			hm << hom_snps
			ylim_hm << SNPdist.get_ylim(hom_snps, genome_length, 'density')

			het_snps = WriteIt.file_to_ints_array("ht_snps_short.txt")
			ht << het_snps
			ylim_ht << SNPdist.get_ylim(het_snps, genome_length, 'density')

			hyp_snps = SNPdist.hyp_snps(expected_ratios, genome_length)
			hyp << hyp_snps
			ylim_hyp << SNPdist.get_ylim(hyp_snps, genome_length, 'density')
		end

		Dir.chdir(File.join(Dir.home, "SNP_distribution_method/Small_genomes/#{perm}")) do

			perm_hm = WriteIt.file_to_ints_array("perm_hm.txt")
			SNPdist.plot_snps(perm_hm, hm[0], "SNP_distribution_method/Small_genomes", "#{perm}", 1, genome_length, 'hm',
				'Homozygous SNP density', ylim_hm[0])

			perm_ht = WriteIt.file_to_ints_array("perm_ht.txt")
			SNPdist.plot_snps(perm_ht, ht[0], "SNP_distribution_method/Small_genomes", "#{perm}", 1, genome_length, 'ht',
				'Heterozygous SNP density', ylim_ht[0])

			perm_hyp = SNPdist.hyp_snps(ratios, genome_length)
			SNPdist.plot_snps(perm_hyp, hyp[0], "SNP_distribution_method/Small_genomes", "#{perm}", 1, genome_length, 'hyp',
				'Approximated ratio of homozygous to heterozygous SNP density', ylim_hyp[0])
		end
	end


	def self.candidate(mut, frag_pos_hm)
		candidate_mutations ={}
		mut.each do |frag|
			if frag_pos_hm.has_key?(frag)
				candidate_mutations.store(frag, frag_pos_hm[frag])
			end
		end
		return candidate_mutations
	end

	def self.density_plots(contig_size, ratios, expected_ratios, snps_hm, snps_ht, region, genome_len, file, mut, frag_pos_hm, original, outcome)
		n = 1048576*4
		average_positions = SNPdist.general_positions(contig_size, ratios)
		hyp_ratios = SNPdist.densities_pos(expected_ratios, average_positions)
		real_ratios = SNPdist.densities_pos(ratios, average_positions)
		WriteIt::write_txt("#{file}/hyp_ratios", hyp_ratios)
		WriteIt::write_txt("#{file}/ratios", real_ratios)
		peak =  LocateMutation.find_peak(hyp_ratios, n) # Find the peak in the approximated (hypothetical SNP) distribution
		ylim = Plot.get_ylim(hyp_ratios, region)
		candidate_peak = LocateMutation.closest_snp(peak, snps_hm)
		candidate_mutations = Mutation.candidate(mut, frag_pos_hm)
		File.open("#{file}/mutation.txt", "w+") do |f|
			f.puts "The length of the group of contigs that form the peak of the distribution is #{region.to_i} bp"
			f.puts "The mutation is likely to be found on the following contigs #{candidate_mutations}"
		end
		original_pos, outcome_pos = [], []
		candidate_mutations.each_key { |frag|
			original_pos << original[frag][:hm_pos] if original.key?(frag)
			outcome_pos << outcome[frag][:hm_pos] if outcome.key?(frag)
		}
		original_pos.flatten!
		outcome_pos.flatten!
		# calculate mean of orignal and outcome position to adj position on plots
		original_mean = original_pos.mean
		outcome_mean = outcome_pos.mean
		adj_mean = original_mean - outcome_mean
		outcome_pos.map! {|x| x + adj_mean.to_i }
		Plot::densities(snps_hm, snps_ht, hyp_ratios, region, file)
		Plot::comparison(real_ratios, hyp_ratios, genome_len, file, ylim, original_pos, outcome_pos)
		Plot::qqplot(snps_hm, file, "QQplot for hm density", "Theoretical normal distribution", "Hypothetical SNP density", "hm_snps")
		Plot::qqplot(hyp_ratios, file, "QQplot for the ratios", "Theoretical normal distribution", "Hypothetical ratios", "ratios")

		candidate_peak
	end

  # function to get cumulative variant positions from the order of fragments
  # input1: hash of frag ids with positions for homozygous and heterozygous variants
  # input2: hash of fragment lengths
  # input3: array of fragment order
  # input4: ratio adjustment factor
  # output: a hash of frag ids with all details and variant positions
  # are accumulated using length and order of fragments
  def self.varpos_aggregate(frag_info, frag_len, frag_order, ratio_adjust)
    details = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
    asmbly_len = 0
    frag_order.each { | frag |
      details[frag][:hm] = ratio_adjust
      details[frag][:ht] = ratio_adjust
      details[frag][:hm_pos] = []
      details[frag][:ht_pos] = []
      if frag_info[:hom].key?(frag)
        hm_pos = frag_info[:hom][frag]
        details[frag][:hm] += hm_pos.length
        details[frag][:hm_pos] = hm_pos.map { |position| position + asmbly_len }
      end
      if frag_info[:het].key?(frag)
        ht_pos = frag_info[:het][frag]
        details[frag][:ht] += ht_pos.length
        details[frag][:ht_pos] = ht_pos.map { |position| position + asmbly_len }
      end
      if details[frag][:hm] == ratio_adjust and details[frag][:ht] == ratio_adjust
        details[frag][:ratio] = 0.0
      else
        details[frag][:ratio] = details[frag][:hm]/details[frag][:ht]
      end
      details[frag][:len] = frag_len[frag].to_i
      asmbly_len += frag_len[frag].to_i
    }
    details
  end

end
