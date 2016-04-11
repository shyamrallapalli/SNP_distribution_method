# encoding: utf-8
require 'bio'
require 'bio-samtools'
require 'bio-gngm'
require_relative 'bfr'

class Pileup

  attr_accessor :defaults

  @defaults = {
      noise: 0.1,                  # noise level for read depth
      ht_low: 0.2,                 # min allele freq for heterozygosity
      ht_high: 0.9,                # max allele freq for heterozygosity
      min_depth: 6,                # minimum coverage for variant
      min_non_ref_count: 3,
      ignore_reference_n: true,
      min_indel_count_support: 3,
  }

  # check if the pileup has the parameters we are looking for
  def self.is_var?(pileup)
    ignore_reference_n = @defaults[:ignore_reference_n]
    min_depth  = @defaults[:min_depth]
    min_non_ref_count = @defaults[:min_non_ref_count]

    return false if pileup.ref_base == '*'
    return false if ignore_reference_n and pileup.ref_base =~ /^[nN]$/
    non_ref_count = get_nonref_count(pileup)
    return true if pileup.coverage >= min_depth and non_ref_count >= min_non_ref_count
    false
  end

  # get read bases from pileup object
  # removes mapping quality information
  def self.get_read_bases(pileup)
    read_bases = pileup.instance_variable_get(:@read_bases)
    # mapping quality after '^' symbol is substituted
    # to avoid splitting at non indel + or - characters
    # end of the read marking '$' symbol is substituted
    # insertion and deletion marking '*' symbol is substituted
    read_bases.gsub!(/\^./, '')
    read_bases.gsub!(/\$/, '')
    read_bases.gsub!(/\*/, '')
    # warn about reads with ambiguous codes
    # if read_bases.match(/[^atgcATGC,\.\+\-0-9]/)
    #   warn "Ambiguous nucleotide\t#{read_bases}"
    # end
    read_bases
  end

  def self.basehash_counts(read_bases)
    bases_hash = {}
    bases_hash[:ref] = read_bases.count('.,')
    bases_hash[:A] = read_bases.count('aA')
    bases_hash[:C] = read_bases.count('cC')
    bases_hash[:G] = read_bases.count('gG')
    bases_hash[:T] = read_bases.count('tT')
    bases_hash[:N] = read_bases.count('nN')
    bases_hash
  end

  # count number of indels and number non-indel base
  # and return a hash with bases and indel counts
  def self.indels_to_hash(read_bases, delimiter)
    indel_bases = 'acgtryswkmbdhvnACGTRYSWKMBDHVN'
    non_indel_bases = String.new
    array = read_bases.split(delimiter)
    non_indel_bases << array.shift
    array.each do |element|
      # get number of nucleotides inserted or deleted
      number = /^(\d+)[#{indel_bases}]/.match(element)[1].to_i
      # capture remaining nucleotides
      non_indel_bases << element.gsub(/^#{number}\w{#{number}}/, '')
    end
    bases_hash = basehash_counts(non_indel_bases)
    # check at least three reads are supporting indel
    indel_count = read_bases.count(delimiter)
    if indel_count >= @defaults[:min_indel_count_support]
      bases_hash[:indel] = indel_count
    end
    bases_hash
  end

  # count bases matching reference and non-reference
  # from snp variant and make a hash of bases with counts
  # for indels return the read bases information instead
  def self.read_bases_to_hash(pileup)
    read_bases = Pileup.get_read_bases(pileup)
    if read_bases =~ /\+/
      bases_hash = indels_to_hash(read_bases, '+')
    elsif read_bases =~ /\-/
      bases_hash = indels_to_hash(read_bases, '-')
    else
      bases_hash = basehash_counts(read_bases)
    end
    # some indels will have ref base in the read and using
    # sum of hash values is going to give wrond addtional coverage
    # from indels so including actual coverage from pileup
    # bases_hash keys are :A, :C, :G, :T, :N, :ref, :indel and :cov
    bases_hash[:cov] = pileup.coverage
    bases_hash
  end

  # count bases from indels
  # array of pileup bases is split at + / -
  # and number after each + / - is counted
  def self.count_indel_bases(read_bases, delimiter)
    array = read_bases.split(delimiter)
    number = 0
    array.shift
    array.each do |element|
      # deletions in reference could contain ambiguous codes,
      number += /^(\d+)[acgtryswkmbdhvnACGTRYSWKMBDHVN]/.match(element)[1].to_i
    end
    number
  end

  # count bases matching reference and non-reference
  # and calculate ratio of non_ref allele to total bases
  def self.get_nonref_count(pileup)
    read_bases = get_read_bases(pileup)
    if read_bases =~ /\+/
      non_ref_count = read_bases.count('atgcnATGCN')
      pluscounts = read_bases.count('+')
      indel_bases = count_indel_bases(read_bases, '+')
      non_ref_count += pluscounts - indel_bases
    elsif read_bases =~ /\-/
      non_ref_count = read_bases.count('acgtryswkmbdhvnACGTRYSWKMBDHVN')
      minuscounts = read_bases.count('-')
      indel_bases = count_indel_bases(read_bases, '-')
      non_ref_count += minuscounts - indel_bases
    else
      non_ref_count = read_bases.count('atgcATGC')
    end
    non_ref_count
  end

  # count bases matching reference and non-reference
  # and calculate ratio of non_ref allele to total bases
  def self.get_nonref_ratio(pileup)
    non_ref_count = get_nonref_count(pileup)
    non_ref_count.to_f / pileup.coverage.to_f
  end

end
