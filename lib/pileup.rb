# encoding: utf-8
require 'bio'
require 'bio-samtools'
require 'bio-gngm'

class Pileup

  attr_accessor :defaults

  @defaults = {
      bq: 15,                      # base quality
      mq: 20,                      # mapping quality
      noise: 0.1,                  # noise level for read depth
      ht_low: 0.2,                 # min allele freq for heterozygosity
      ht_high: 0.9,                # max allele freq for heterozygosity
      min_depth: 6,                # minimum coverage for variant
      min_non_ref_count: 3,
      ignore_reference_n: true,
      polyploidy: false,
      min_indel_count_support: 3,
      bgbam: '',                   # bam path for background bulk
      parent_hemi_hash: '',        # hash of hemi snps from parents
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

  def self.get_bg_ratio(bg_bam, selfrag, mutpos)
    bg_ratio = ''
    bg_pileups = get_pileup(bg_bam, selfrag, mutpos)
    return bg_ratio if bg_pileups.empty?
    if is_var?(bg_pileups[0])
      bg_ratio = get_nonref_ratio(bg_pileups[0])
    end
    bg_ratio
  end

  def self.pick_frag_vars(mutbam,infasta,keyfrags,input_frags, var_pos, opts = {})
    @defaults.merge!(opts)
    bgbam = @defaults[:bgbam]

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
        pileups = get_pileup(mut_bam, selfrag, mutpos)
        if pileups.empty?
          next
        end
        pileup = pileups[0]
        if is_var?(pileup)
          ratio = get_nonref_ratio(pileup)
          if bgbam != ''
            bg_ratio = get_bg_ratio(bg_bam, selfrag, mutpos)
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
  def self.get_pileup(bamobject, id, pos1)
    bq = @defaults[:bq]
    mq = @defaults[:mq]

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
      pileuparray = get_pileup(bamobject, id, pos1)
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
  def self.get_var_base_frac(hash)
    snp_hash = {}
    coverage = hash[:cov]
    return snp_hash if coverage < @defaults[:min_depth]
    # calculate proportion of each base in coverage
    hash.each_key do | base |
      next if base == :cov
      freq = hash[base].to_f/coverage.to_f
      next if freq <= @defaults[:noise]
      snp_hash[base] = freq
    end
    snp_hash
  end

  # calculate var zygosity for non-polyploid variants
  # increased range is used for heterozygosity for RNA-seq data
  def self.var_mode(ratio)
    ht_low = @defaults[:ht_low]
    ht_high = @defaults[:ht_high]
    mode = ''
    if ratio.between?(ht_low, ht_high)
      mode = :het
    elsif ratio > ht_high
      mode = :hom
    end
    mode
  end

  # getting vars from pileup file
  # each var is checked from pileup information
  # added to a hash to return
  def self.vars_in_pileup(pileupfile, opts = {})
    @defaults.merge!(opts)

    # hash of frag ids with respective variant positions and their base hash info
    # only snps have base hash info and indels base hash is read bases
    vars_hash = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }

    # read mpileup file and process each variant
    File.open(pileupfile, 'r').each do |line|
      pileup = Bio::DB::Pileup.new(line)
      if is_var?(pileup)
        basehash = Pileup.read_bases_to_hash(pileup)
        vars_hash[pileup.ref_name][pileup.pos] = basehash
        # puts "#{pileup.ref_name}\t#{pileup.pos}\t#{pileup.consensus}\t#{basehash}\n"
      end
    end
    vars_hash
  end

  def self.filter_vars(mut_pileup, bg_bulk_pileup_hash, opts = {})
    @defaults.merge!(opts)

    vars_hash = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
    # read mpileup file and process each variant
    File.open(mut_pileup, 'r').each do |line|
      pileup = Bio::DB::Pileup.new(line)
      if is_var?(pileup)
        vars_hash = compare_bulk_pileups(pileup, bg_bulk_pileup_hash, vars_hash)
      end
    end
    vars_hash
  end

  def self.compare_bulk_pileups(mut_pileup, bg_pileup_hash, out_hash)
    parent_hemi_hash = @defaults[:parent_hemi_hash]
    data1 = Pileup.read_bases_to_hash(mut_pileup)
    frag = mut_pileup.ref_name
    pos = mut_pileup.pos
    mut_bases = get_var_base_frac(data1)
    if @defaults[:polyploidy]
      if parent_hemi_hash != '' and parent_hemi_hash[frag].key?(pos)
        bg_bases = ''
        if bg_pileup_hash[frag].key?(pos)
          bg_bases = get_var_base_frac(bg_pileup_hash[frag][pos])
        end
        bfr = Bfr.get_bfr(mut_bases, :bg_hash => bg_bases)
        out_hash[:hemi][frag][pos] = bfr
      else
        out_hash = wrapper_to_push_base_hash(mut_bases, frag, pos, bg_pileup_hash, out_hash)
      end
    else
      out_hash = wrapper_to_push_base_hash(mut_bases, frag, pos, bg_pileup_hash, out_hash)
    end
    out_hash
  end

  def self.wrapper_to_push_base_hash(mut_bases, frag, pos, bg_pileup_hash, out_hash)
    if bg_pileup_hash[frag].key?(pos)
      bg_bases = get_var_base_frac(bg_pileup_hash[frag][pos])
      out_hash = push_base_hash(mut_bases, out_hash, frag, pos, bg_bases)
    else
      out_hash = push_base_hash(mut_bases, out_hash, frag, pos)
    end
    out_hash
  end

  def self.push_base_hash(base_hash, store_hash, frag, pos, background='')
    # we are only dealing with single element hashes
    # so discard hashes with more than one element and empty hashes
    # empty hash results from position below selected coverage or bases freq below noise
    base_hash.delete(:ref)
    return store_hash if base_hash.empty?
    # we could ignore complex loci or
    # take the variant type based on predominant base
    if base_hash.length > 1
      mut_type = var_mode(base_hash.values.max)
      if background != ''
        bg_type = var_mode(background.values.max)
        return store_hash if mut_type == bg_type
      end
      store_hash[mut_type][frag][pos] = base_hash.values.max
    else
      base = base_hash.keys[0]
      mut_type = var_mode(base_hash[base])
      if background != '' and background.key?(base)
        bg_type = var_mode(background[base])
        # if both have the same base type then return original hash
        return store_hash if mut_type == bg_type
        store_hash[mut_type][frag][pos] = base_hash[base]
      else
        store_hash[mut_type][frag][pos] = base_hash[base]
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
    polyploidy = @defaults[:polyploidy]
    # if polyploidy set then take combination of proportions
    # if not then take maximum value
    if polyploidy
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
