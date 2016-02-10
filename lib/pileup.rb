# encoding: utf-8
require 'bio'
require 'bio-samtools'
require 'bio-gngm'

class Pileup

  # count bases from indels
  # array of pileup bases is split at + / -
  # and number after each + / - is counted
  def self.count_indels(array)
    number = 0
    array.shift
    array.each do |element|
      num = /^(\d+)[atgcATGC]/.match(element)[1].to_i
      number += num
    end
    number
  end

  # count bases matching reference and non-reference
  # and calculate ratio of non_ref allele to total bases
  def self.read_base_hash(pileup)
    bases_hash = {}
    read_bases = pileup.instance_variable_get(:@read_bases)
    # ref_base = pileup.instance_variable_get(:@ref_base)
    # read_bases.gsub!(/[.,]/, ref_base)
    if read_bases =~ /\+/
      return bases_hash
    elsif read_bases =~ /\-/
      return bases_hash
    else
      bases_hash[:A] = read_bases.count('aA')
      bases_hash[:C] = read_bases.count('cC')
      bases_hash[:G] = read_bases.count('gG')
      bases_hash[:T] = read_bases.count('tT')
    end
    bases_hash
  end

  # count bases matching reference and non-reference
  # and calculate ratio of non_ref allele to total bases
  def self.get_nonref_ratio(pileup)
    read_bases = pileup.instance_variable_get(:@read_bases)
    ref_count = read_bases.count('.,')
    basecounts = read_bases.count('atgcATGC')
    non_ref_count = basecounts
    if read_bases =~ /\+/
      pluscounts = read_bases.count('+')
      indel_bases = count_indels(read_bases.split('+'))
      non_ref_count = basecounts + pluscounts - indel_bases
    elsif read_bases =~ /\-/
      minuscounts = read_bases.count('-')
      indel_bases = count_indels(read_bases.split('-'))
      non_ref_count = basecounts + minuscounts - indel_bases
    end
    ratio = non_ref_count.to_f / (ref_count.to_f + non_ref_count.to_f)
    ratio
  end

  def self.pick_frag_vars(bamfile,infasta,keyfrags,input_frags)
    bam = Bio::DB::Sam.new(:bam=>bamfile, :fasta=>infasta)
    bam.open
    sortfrags = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
    keyfrags.each do | selfrag |
      positions = input_frags[selfrag][:hm_pos]
      positions.each do | mutpos |
        pileups = get_pileup(bam,selfrag,mutpos)
        if pileups.empty?
          next
        end
        pileup = pileups[0]
        if pileup.is_snp?(:ignore_reference_n => true, :min_depth => 6, :min_non_ref_count => 3)
          base_hash = read_base_hash(pileup)
          ratio = get_nonref_ratio(pileup)
          sortfrags[ratio][selfrag][mutpos] = pileup
        end
      end
    end
    sortfrags
  end

  # a function to get pileup using a bam object and sequence information
  # if the pileup is giving zero in the information there is insertion or deletion upstream of it
  # there seems to be a location difference between vcf files and mpileup output
  # so moving upstream in position to get these var info
  # and not going further than 10bp upstream (10bp is arbitrary and this needs verification)
  def self.get_pileup(bamobject,id,pos1, bq=15, mq=20)
    # a variable to store initial position of var
    initial_pos = pos1
    pileuparray = []
    bamobject.mpileup(:r => "#{id}:#{pos1}-#{pos1}", :Q => bq, :q => mq) do |pileup|
      pileuparray << pileup
    end
    # calculate difference between var position and adjusted positions
    # and check if the pileup information has zeros and change position info
    # to get upstream position pileup information
    pos_diff = initial_pos - pos1
    if pileuparray[0].to_s =~ /^\t0/ and pos_diff <= 10 and pos1 > 1
      pos1 = pos1 - 1
      pileuparray = get_pileup(bamobject,id,pos1)
    # if the difference is larger than 10 bp then set up pileup information as null
    elsif pileuparray[0].to_s =~ /^\t0/ and pos_diff > 10
      pileuparray = []
    # might be best to increase the number to 25 instead of 1
    # to reduce possible errors in variant calls at ends of sequences
    elsif pos1 == 1
      pileuparray = []
    end
    pileuparray
  end

end