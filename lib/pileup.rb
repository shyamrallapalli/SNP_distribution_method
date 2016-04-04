# encoding: utf-8
require 'bio'
require 'bio-samtools'
require 'bio-gngm'

class Pileup

  attr_accessor :min_indel_count_support

  DEFAULT = {
      bgbam: '',
      bq: 15,
      mq: 20,
      ignore_reference_n: true,
      min_depth: 6,
      min_non_ref_count: 3,
      min_indel_count_support: 3,
  }

  # count bases from indels
  # array of pileup bases is split at + / -
  # and number after each + / - is counted
  def self.count_indels(read_bases, delimiter)
    array = read_bases.split(delimiter)
    number = 0
    array.shift
    array.each do |element|
      # deletions in reference could contain ambiguous codes,
      number += /^(\d+)[acgtryswkmbdhvnACGTRYSWKMBDHVN]/.match(element)[1].to_i
    end
    number
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
    indel_count = array.length - 1
    non_indel_bases << array.shift
    array.each do |element|
      # get number of nucleotides inserted or deleted
      number = /^(\d+)[#{indel_bases}]/.match(element)[1].to_i
      # capture remaining nucleotides
      non_indel_bases << element.gsub(/^#{number}\w{#{number}}/, '')
    end
    bases_hash = basehash_counts(non_indel_bases)
    # check at least three reads are supporting indel
    if indel_count >= @min_indel_count_support
      bases_hash[:indel] = indel_count
    end
    bases_hash
  end

  # count bases matching reference and non-reference
  # from snp variant and make a hash of bases with counts
  # for indels return the read bases information instead
  def self.read_bases_to_hash(read_bases)
    if read_bases =~ /\+/
      bases_hash = indels_to_hash(read_bases, '+')
    elsif read_bases =~ /\-/
      bases_hash = indels_to_hash(read_bases, '-')
    else
      bases_hash = basehash_counts(read_bases)
    end
    bases_hash
  end

  # count bases matching reference and non-reference
  # and calculate ratio of non_ref allele to total bases
  def self.get_nonref_ratio(read_bases)
    ref_count = read_bases.count('.,')
    if read_bases =~ /\+/
      non_ref_count = read_bases.count('atgcnATGCN')
      pluscounts = read_bases.count('+')
      indel_bases = count_indels(read_bases, '+')
      non_ref_count += pluscounts - indel_bases
    elsif read_bases =~ /\-/
      non_ref_count = read_bases.count('acgtryswkmbdhvnACGTRYSWKMBDHVN')
      minuscounts = read_bases.count('-')
      indel_bases = count_indels(read_bases, '-')
      non_ref_count += minuscounts - indel_bases
    else
      non_ref_count = read_bases.count('atgcATGC')
    end
    ratio = non_ref_count.to_f / (ref_count.to_f + non_ref_count.to_f)
    ratio
  end

  def self.get_bg_ratio(bg_bam,selfrag,mutpos, opts = {})
    opts = DEFAULT.merge(opts)
    ignore_reference_n = opts[:ignore_reference_n]
    min_depth  = opts[:min_depth]
    min_non_ref_count = opts[:min_non_ref_count]
    bq = opts[:bq]
    mq = opts[:mq]

    bg_ratio = ''
    bg_pileups = get_pileup(bg_bam,selfrag,mutpos, bq, mq)
    if bg_pileups.empty?
      return bg_ratio
    end
    pileup = bg_pileups[0]
    if pileup.is_snp?(:ignore_reference_n => ignore_reference_n, :min_depth => min_depth, :min_non_ref_count => min_non_ref_count)
      read_bases = get_read_bases(pileup)
      bg_ratio = get_nonref_ratio(read_bases)
    end
    bg_ratio
  end

  def self.pick_frag_vars(mutbam,infasta,keyfrags,input_frags, var_pos, opts = {})
    opts = DEFAULT.merge(opts)
    ignore_reference_n = opts[:ignore_reference_n]
    min_depth  = opts[:min_depth]
    min_non_ref_count = opts[:min_non_ref_count]
    bgbam = opts[:bgbam]
    bq = opts[:bq]
    mq = opts[:mq]
    @min_indel_count_support = opts[:min_indel_count_support]

    mut_bam = Bio::DB::Sam.new(:bam=>mutbam, :fasta=>infasta)
    mut_bam.open
    if bgbam == ''
      bg_bam = ''
    else
      bg_bam = Bio::DB::Sam.new(:bam=>bgbam, :fasta=>infasta)
      bg_bam.open
    end
    sortfrags = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
    keyfrags.each do | selfrag |
      positions = input_frags[selfrag][:hm_pos]
      positions.each do | mutpos |
        pileups = get_pileup(mut_bam,selfrag,mutpos, bq, mq)
        if pileups.empty?
          next
        end
        pileup = pileups[0]
        if pileup.is_snp?(:ignore_reference_n => ignore_reference_n, :min_depth => min_depth, :min_non_ref_count => min_non_ref_count)
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
              var_pos[:hom][selfrag].delete(mutpos)
              next
            end
          else
            sortfrags[ratio][selfrag][mutpos] = pileup
          end
        end
      end
    end
    [sortfrags, var_pos]
  end

  # a function to get pileup using a bam object and sequence information
  # if the pileup is giving zero in the information there is insertion or deletion upstream of it
  # there seems to be a location difference between vcf files and mpileup output
  # so moving upstream in position to get these var info
  # and not going further than 10bp upstream (10bp is arbitrary and this needs verification)
  def self.get_pileup(bamobject,id,pos1, bq, mq)
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
      pileuparray = get_pileup(bamobject,id,pos1, bq, mq)
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

  # form hash of base information, [ATGC] counts for snp
  # a hash of base proportion is calculated
  # base proportion hash below a selected depth is empty
  # base proportion below or equal to a noise factor are discarded
  def self.get_var_base_prop(hash)
    snp_hash = {}
    coverage = hash.values.inject { | sum, n | sum + n.to_f }
    return snp_hash if coverage < @min_depth
    # calculate proportion of each base in coverage
    hash.each_key do | base |
      next if base == :ref
      freq = hash[base].to_f/coverage
      next if freq <= @noise
      snp_hash[base] = freq
    end
    snp_hash
  end

  # calculate var zygosity for non-polyploid variants
  # increased range is used for heterozygosity for RNA-seq data
  def self.var_mode(ratio)
    mode = ''
    if ratio.between?(@ht_low, @ht_high)
      mode = :het
    elsif ratio > @ht_high
      mode = :hom
    end
    mode
  end

  # getting vars from pileup file
  # each var is checked from pileup information
  # added to a hash to return
  def self.vars_in_pileup(pileupfile, opts = {})
    opts = DEFAULT.merge(opts)
    ignore_reference_n = opts[:ignore_reference_n]
    min_depth  = opts[:min_depth]
    min_non_ref_count = opts[:min_non_ref_count]

    # hash of frag ids with respective variant positions and their base hash info
    # only snps have base hash info and indels base hash is read bases
    vars_hash = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }

    # read mpileup file and process each variant
    File.open(pileupfile, 'r').each do |line|
      pileup = Bio::DB::Pileup.new(line)
      if pileup.is_snp?(:ignore_reference_n => ignore_reference_n, :min_depth => min_depth,
                        :min_non_ref_count => min_non_ref_count) and
          pileup.consensus != pileup.ref_base
        read_bases = Pileup.get_read_bases(pileup)
        basehash = Pileup.read_bases_to_hash(read_bases)
        vars_hash[pileup.ref_name][pileup.pos] = basehash
        # puts "#{pileup.ref_name}\t#{pileup.pos}\t#{pileup.consensus}\t#{basehash}\n"
      end
    end
    vars_hash
  end

  def self.filter_vars(mut_pileup, vars_hash_bg, opts = {})
    opts = DEFAULT.merge(opts)
    ignore_reference_n = opts[:ignore_reference_n]
    min_non_ref_count = opts[:min_non_ref_count]
    @min_depth  = opts[:min_depth]
    @noise = opts[:noise]
    @polyploidy = opts[:polyploidy]
    @ht_low = opts[:ht_low]
    @ht_high  = opts[:ht_high]

    vars_hash = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
    # read mpileup file and process each variant
    File.open(mut_pileup, 'r').each do |line|
      pileup = Bio::DB::Pileup.new(line)
      if pileup.is_snp?(:ignore_reference_n => ignore_reference_n, :min_depth => @min_depth, :min_non_ref_count => min_non_ref_count) and
          pileup.consensus != pileup.ref_base
        read_bases = Pileup.get_read_bases(pileup)
        data1 = Pileup.read_bases_to_hash(read_bases)
        frag = pileup.ref_name
        pos = pileup.pos
        mut_bases = get_var_base_prop(data1)
        if vars_hash_bg[frag].key?(pos)
          bg_bases = get_var_base_prop(vars_hash_bg[frag][pos])
          vars_hash = push_base_hash(mut_bases, vars_hash, frag, pos, bg_bases)
        else
          vars_hash = push_base_hash(mut_bases, vars_hash, frag, pos)
        end
      end
    end
    vars_hash
  end

  def self.push_base_hash(base_hash, store_hash, frag, pos, background='')
    # we are only dealing with single element hashes
    # so discard hashes with more than one element and empty hashes
    # empty hash results from position below selected coverage or bases freq below noise
    return store_hash if base_hash.empty?
    if base_hash.length > 1
      vartype = multi_var_hash(base_hash, background)
      # warn "#{frag}\t#{pos}\t#{base_hash}\t#{background}\n"
      if vartype == ''
        return store_hash
      else
        store_hash = Vcf.push_to_hash(store_hash, frag, pos, vartype)
      end
    else
      base = base_hash.keys[0]
      mut_type = var_mode(base_hash[base])
      if background != '' and background.key?(base)
        bg_type = var_mode(background[base])
        # if both have the same base type then return original hash
        return store_hash if mut_type == bg_type
        store_hash = Vcf.push_to_hash(store_hash, frag, pos, mut_type)
      else
        store_hash = Vcf.push_to_hash(store_hash, frag, pos, mut_type)
      end
    end
    store_hash
  end

  # method to compare base hash between and background and mutant
  # returns the var_type if single base left after background subtraction
  # otherwise returns empty for multiple vars
  def self.multi_var_hash(base_hash, background='')
    var_type = ''
    mut_vars = base_hash.keys.sort
    # ignore complex variant locations
    # return var_type if mut_vars.length > 2
    if background != ''
      bg_vars = background.keys.sort
      return var_type if mut_vars == bg_vars
      mut_vars.delete_if { |base| bg_vars.include?(base) }
      if mut_vars.length == 1
        var_type = var_mode(base_hash[mut_vars[0]])
      else
        ratio = polybase_proportion(mut_vars, base_hash)
        var_type = var_mode(ratio)
      end
    else
      ratio = polybase_proportion(mut_vars, base_hash)
      var_type = var_mode(ratio)
    end
    var_type
  end

  # get total proportion of bases in hash
  def self.polybase_proportion(vars, hash)
    # if polyploidy set then take combination of proportions
    # if not then take maximum value
    if @polyploidy
      prop = 0.0
      hash.each_key { | key |
        if vars.include?(key)
          prop += hash[key]
        end
      }
      # warn "polybase\t#{hash}\t#{vars}\t#{prop}\n"
    else
      prop = hash.values.max
      # warn "non_polybase\t#{hash}\t#{vars}\t#{prop}\n"
    end
    prop
  end

end
