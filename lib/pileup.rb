# encoding: utf-8
require 'bio'
require 'bio-samtools'
require 'bio-gngm'

class Pileup

  # count bases from indels
  # array of pileup bases is split at + / -
  # and number after each + / - is counted
  def self.count_indels(read_bases, delimiter)
    array = read_bases.split(delimiter)
    number = 0
    array.shift
    array.each do |element|
      num = /^(\d+)[atgcATGC]/.match(element)[1].to_i
      number += num
    end
    number
  end

  # get read bases from pileup object
  # removes mapping quality information
  def self.get_read_bases(pileup)
    read_bases = pileup.instance_variable_get(:@read_bases)
    # mapping quality after '^' symbol is substituted
    # to avoid splitting at non indel + or - characters
    read_bases.gsub!(/\^./, '')
    read_bases
  end

  # count bases matching reference and non-reference
  # from snp variant and make a hash of bases with counts
  # for indels return the read bases information instead
  def self.read_bases_to_hash(read_bases)
    bases_hash = {}
    if read_bases =~ /\+/
      return read_bases
    elsif read_bases =~ /\-/
      return read_bases
    else
      bases_hash[:ref] = read_bases.count('.,')
      bases_hash[:A] = read_bases.count('aA')
      bases_hash[:C] = read_bases.count('cC')
      bases_hash[:G] = read_bases.count('gG')
      bases_hash[:T] = read_bases.count('tT')
    end
    bases_hash
  end

  # count bases matching reference and non-reference
  # and calculate ratio of non_ref allele to total bases
  def self.get_nonref_ratio(read_bases)
    ref_count = read_bases.count('.,')
    basecounts = read_bases.count('atgcATGC')
    non_ref_count = basecounts
    if read_bases =~ /\+/
      pluscounts = read_bases.count('+')
      indel_bases = count_indels(read_bases, '+')
      non_ref_count += pluscounts - indel_bases
    elsif read_bases =~ /\-/
      minuscounts = read_bases.count('-')
      indel_bases = count_indels(read_bases, '-')
      non_ref_count += minuscounts - indel_bases
    end
    ratio = non_ref_count.to_f / (ref_count.to_f + non_ref_count.to_f)
    ratio
  end

  def self.get_bg_ratio(bg_bam,selfrag,mutpos)
    bg_ratio = ''
    bg_pileups = get_pileup(bg_bam,selfrag,mutpos)
    if bg_pileups.empty?
      return bg_ratio
    end
    pileup = bg_pileups[0]
    read_bases = get_read_bases(pileup)
    bg_ratio = get_nonref_ratio(read_bases)
    bg_ratio
  end

  def self.pick_frag_vars(mutbam,infasta,keyfrags,input_frags, bgbam='')
    mut_bam = Bio::DB::Sam.new(:bam=>mutbam, :fasta=>infasta)
    mut_bam.open
    if bgbam != ''
      bg_bam = Bio::DB::Sam.new(:bam=>bgbam, :fasta=>infasta)
      bg_bam.open
    end
    sortfrags = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
    keyfrags.each do | selfrag |
      positions = input_frags[selfrag][:hm_pos]
      positions.each do | mutpos |
        pileups = get_pileup(mut_bam,selfrag,mutpos)
        if pileups.empty?
          next
        end
        pileup = pileups[0]
        if pileup.is_snp?(:ignore_reference_n => true, :min_depth => 6, :min_non_ref_count => 3)
          read_bases = get_read_bases(pileup)
          ratio = get_nonref_ratio(read_bases)
          if bgbam != ''
            bg_ratio = get_bg_ratio(bg_bam,selfrag,mutpos)
            if bg_ratio == ''
              sortfrags[ratio][selfrag][mutpos] = pileup
            elsif bg_ratio <= 0.35
              sortfrags[ratio][selfrag][mutpos] = pileup
            else
              warn "#{selfrag}\t#{mutpos}\t#{ratio}\t#{bg_ratio}\n"
              next
            end
          else
            sortfrags[ratio][selfrag][mutpos] = pileup
          end
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
