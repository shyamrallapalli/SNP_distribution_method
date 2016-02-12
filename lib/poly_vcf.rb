# encoding: utf-8
require 'bio'
require 'bio-samtools'
require 'bio-gngm'

class Polyploid

  # getting vars from vcf file
  # each var is checked using pileup information from bam file
  # added to a hash to return
  def self.vars_in_file(vcf_file, bamfile, fastafile)
    # hash of frag ids with respective variant positions and their base hash info
    # only snps have base hash info and indels base hash is empty
    vars_hash = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }

    # open bam object to access pileup information
    bam = Bio::DB::Sam.new(:bam=>bamfile, :fasta=>fastafile)
    bam.open

    # read vcf file and process each variant
    File.open(vcf_file, 'r').each do |line|
      next if line =~ /^#/
      v = Bio::DB::Vcf.new(line)
      # some variant callers like Freebays are including some non variants as vars
      # so check ref and alt differ
      if v.variant? and v.ref != v.alt
        pileups = Pileup.get_pileup(bam,v.chrom,v.pos)
        next if pileups.empty?
        basehash = Pileup.read_base_hash(pileups[0])
        vars_hash[v.chrom][v.pos] = basehash
      end
    end
    vars_hash
  end

  # form hash of base information, [ATGC] counts for snp
  # a hash of base proportion is calculated
  # base proportion hash below a selected depth is empty
  # base proportion below or equal to a noise factor are discarded
  def self.get_base_freq(hash, depth, noise)
    snp_hash = {}
    coverage = hash.values.inject { | sum, n | sum + n.to_f }
    return snp_hash if coverage < depth
    # calculate proportion of each base in coverage
    hash.each_key do | base |
      next if base == :ref
      freq = data[base].to_f/coverage
      next if freq <= noise
      snp_hash[base] = freq
    end
    snp_hash
  end

  # calculate var zygosity for non-polyploid variants
  # increased range is used for heterozygosity for RNA-seq data
  def self.var_mode(ratio, ht_low = 0.10, ht_high = 0.90)
    mode = ''
    if ratio.between?(ht_low, ht_high)
      mode = :het
    elsif ratio > ht_high
      mode = :hom
    end
    mode
  end


  def self.push_base_hash(base_hash, store_hash, frag, pos, background='')
    # we are only dealing with single element hashes
    if base_hash.keys.length > 1
      warn "#{frag}\t#{pos}\t#{base_hash}\t#{background}\n"
      return store_hash
    end
    base = base_hash.keys[0]
    if background == ''
      mut_type = var_mode(base_hash[base])
      store_hash = Vcf.push_to_hash(store_hash, frag, pos, mut_type)
    else
      if background.key?(base)
        mut_type = var_mode(base_hash[base])
        bg_type = var_mode(background[base])
        # if both have the same base type then return original hash
        return store_hash if mut_type == bg_type
        store_hash = Vcf.push_to_hash(store_hash, frag, pos, mut_type)
      else
        mut_type = var_mode(base_hash[base])
        store_hash = Vcf.push_to_hash(store_hash, frag, pos, mut_type)
      end
    end
    store_hash
  end

  def self.filter_vars(vars_hash_mut, vars_hash_bg, depth=6, noise=0.1)
    vars_hash = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
    vars_hash_mut.each_key do | frag |
      positions = vars_hash_mut[frag].keys
      positions.each do | pos |
        data1 = vars_hash_mut[frag][pos]
        # background bulk has variant at same position
        if vars_hash_bg[frag].key?(pos)
          data2 = vars_hash_bg[frag][pos]
          if data1.instance_of? Hash
            mut_bases = get_base_freq(data1, depth, noise)
            if data2.instance_of? Hash
              bg_bases = get_base_freq(data2, depth, noise)
              vars_hash = push_base_hash(mut_bases, vars_hash, frag, pos, bg_bases)
            else
              vars_hash = push_base_hash(mut_bases, vars_hash, frag, pos)
            end
          elsif data1.instance_of? String
            mut_ratio = Pileup.get_nonref_ratio(data1)
            mut_type = var_mode(mut_ratio)
            if data2.instance_of? String
              bg_ratio = Pileup.get_nonref_ratio(data2)
              bg_type = var_mode(bg_ratio)
              next if mut_type == bg_type
              vars_hash = Vcf.push_to_hash(vars_hash, frag, pos, mut_type)
            else
              vars_hash = Vcf.push_to_hash(vars_hash, frag, pos, mut_type)
            end
          end
        # only mut bulk has variant at this position
        else
          # var is snp
          if data1.instance_of? Hash
            mut_bases = get_base_freq(data1, depth, noise)
            vars_hash = push_base_hash(mut_bases, vars_hash, frag, pos)
          # var is indel
          elsif data1.instance_of? String
            mut_ratio = Pileup.get_nonref_ratio(data1)
            mut_type = var_mode(mut_ratio)
            vars_hash = Vcf.push_to_hash(vars_hash, frag, pos, mut_type)
          end
        end
      end
    end
    vars_hash
  end

end
