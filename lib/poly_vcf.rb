# encoding: utf-8
require 'bio'
require 'bio-samtools'
require 'bio-gngm'

class Polyploid

  attr_accessor :polyploidy, :ht_low, :ht_high, :min_depth
  attr_accessor :noise, :bfr_adj, :parent_hemi_hash

  DEFAULT = {
      ignore_reference_n: true,
      min_depth: 6,
      min_non_ref_count: 3,
      noise: 0.1,
      polyploidy: false,
      ht_low: 0.1,
      ht_high: 0.9,
      bfr_adj: 0.01,
      parent_hemi_hash: '',
  }

  def self.mark_hemisnps_in_parent(mut_parent_hash, bg_parent_hash)
    out_hash = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
    # mark all the hemi snp based on both parents
    mut_parent_hash.each_key do |frag|
      mut_parent_hash[frag].each_key do |pos|
        mut_bases = get_var_base_frac(mut_parent_hash[frag][pos])
        if bg_parent_hash[frag].key?(pos)
          bg_bases = get_var_base_frac(bg_parent_hash[frag][pos])
          if mut_bases.length == 2 and mut_bases.key?(:ref)
            out_hash[frag][pos] = 'hemi'
            # calculate bfr
          elsif bg_bases.length == 2 and bg_bases.key?(:ref)
            out_hash[frag][pos] = 'hemi'
          end
          bg_parent_hash[frag].delete(pos)
        else
          if mut_bases.length == 2 and mut_bases.key?(:ref)
            out_hash[frag][pos] = 'hemi'
          end
        end
      end
      mut_parent_hash.delete(frag)
    end

    # now include all hemi snp unique background parent
    bg_parent_hash.each_key do |frag|
      bg_parent_hash[frag].each_key do |pos|
        bg_bases = get_var_base_frac(bg_parent_hash[frag][pos])
        if bg_bases.length == 2 and bg_bases.key?(:ref)
          unless out_hash[frag].key?(pos)
            out_hash[frag][pos] = 'hemi'
          end
        end
      end
      bg_parent_hash.delete(frag)
    end
    out_hash
  end

end
