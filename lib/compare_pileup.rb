# encoding: utf-8
require 'bio'
require 'bio-samtools'
require 'bio-gngm'
require 'progressbar'
require_relative 'bfr'
require_relative 'pileup'

class BamCompare

  attr_accessor :defaults

  @defaults = {
      bq: 15,                      # base quality
      mq: 20,                      # mapping quality
      bgbam: '',                   # bam path for background bulk
  }

  def self.get_bg_ratio(bg_bam, selfrag, mutpos)
    bg_ratio = ''
    bg_pileups = get_pileup(bg_bam, selfrag, mutpos)
    return [bg_ratio, ''] if bg_pileups.empty?
    if is_var?(bg_pileups[0])
      bg_ratio = Pileup.get_nonref_ratio(bg_pileups[0])
    end
    [bg_ratio, bg_pileups[0]]
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
          ratio = Pileup.get_nonref_ratio(pileup)
          if bgbam != ''
            bg_ratio, bg_pileup = get_bg_ratio(bg_bam, selfrag, mutpos)
            if bg_ratio == ''
              sortfrags[ratio][selfrag][mutpos] = pileup
            elsif bg_ratio <= 0.35
              sortfrags[ratio][selfrag][mutpos] = pileup
            else
              warn "#{selfrag}\t#{mutpos}\t#{ratio}\t#{bg_ratio}\t#{pileup}\t#{bg_pileup}"
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

  def self.vars_in_bam(opts = {})
    @defaults.merge!(opts)
    bam = @defaults[:bam]
    fasta = @defaults[:fasta]
    bq = @defaults[:bq]
    mq = @defaults[:mq]

    inbam = Bio::DB::Sam.new(:bam=>bam, :fasta=>fasta)
    vars_hash = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
    inbam.each_region do | region |
      inbam.mpileup_cached(:r => region, :q => mq, :Q => bq) do | pileup |
        if is_var?(pileup)
          basehash = Pileup.read_bases_to_hash(pileup)
          vars_hash[pileup.ref_name][pileup.pos] = basehash
        end
      end
      inbam.mpileup_clear_cache (region)
    end
    vars_hash
  end

end # class BamCompare

class PileupCompare

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
    # totallines = %x[wc -l #{pileupfile}].to_i
    # onepercent = totallines/100
    # pbar = ProgressBar.new("pileupfile", 100)
    @defaults.merge!(opts)

    # hash of frag ids with respective variant positions and their base hash info
    # only snps have base hash info and indels base hash is read bases
    vars_hash = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }

    # read mpileup file and process each variant
    # linenum = 0
    # step = 1
    File.foreach(pileupfile) do |line|
      pileup = Bio::DB::Pileup.new(line)
      if is_var?(pileup)
        basehash = Pileup.read_bases_to_hash(pileup)
        vars_hash[pileup.ref_name][pileup.pos] = basehash
        # puts "#{pileup.ref_name}\t#{pileup.pos}\t#{pileup.consensus}\t#{basehash}\n"
      end
      # linenum += 1
      # if linenum == step * onepercent
      #   pbar.set(step)
      #   step += 1
      # end
    end
    # pbar.finish
    vars_hash
  end

  def self.filter_vars(mut_pileup, bg_bulk_pileup_hash, opts = {})
    # totallines = %x[wc -l #{mut_pileup}].to_i
    # onepercent = totallines/100
    # pbar = ProgressBar.new("mut_pileup", 100)
    @defaults.merge!(opts)

    vars_hash = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
    # read mpileup file and process each variant
    # linenum = 0
    # step = 1
    File.foreach(mut_pileup) do |line|
      pileup = Bio::DB::Pileup.new(line)
      if is_var?(pileup)
        vars_hash = compare_bulk_pileups(pileup, bg_bulk_pileup_hash, vars_hash)
      end
      # linenum += 1
      # if linenum == step * onepercent
      #   pbar.set(step)
      #   step += 1
      # end
    end
    # pbar.finish
    vars_hash
  end

  def self.filter_vars_hash(mut_pileup, bg_bulk_pileup_hash, opts = {})
    @defaults.merge!(opts)

    vars_hash = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
    # read mpileup file and process each variant
    mut_pileup.each_key do |contig|
      mut_pileup[contig].each_key do |pos|
        pileup = mut_pileup[contig][pos]
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

end # class PileupCompare
