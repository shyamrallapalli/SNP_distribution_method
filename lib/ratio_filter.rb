#encoding: utf-8
require 'csv'

class RatioFilter

  DEFAULT = {
      only_frag_with_vars: true,
      ratio_adjust: 0.5,
      ratio_type: 'ratio',
  }

  def self.prune_hash(inhash, adjust, ratio_type)
    first = inhash.length
    inhash.each_key do | frag |
      numhm = inhash[frag][:hm]
      numht = inhash[frag][:ht]
      numbfr = inhash[frag][:bfr]
      # selecting fragments which have a variant
      if ratio_type == :ratio and numht + numhm <= 2 * adjust
        inhash.delete(frag)
      # in polyploidy scenario selecting fragments with at least one bfr position
      elsif ratio_type == :brf_rat and numbfr <= 0
        inhash.delete(frag)
      end
    end
    last = inhash.length
    warn "Discarded #{first-last} out of #{first} fragments, which lack any variant\n"
    inhash
  end

  def self.selected_ratios(inhash, opts = {})
    opts = DEFAULT.merge(opts)
    only_frag_with_vars = opts[:only_frag_with_vars]
    adjust = opts[:ratio_adjust]
    ratio_type = opts[:ratio_type].to_sym
    ratios_hash = {}
    # select only fragments with variants
    # this is to discard fragments which may not contain
    # any information about causative mutations
    if only_frag_with_vars
      inhash = prune_hash(inhash, adjust, ratio_type)
    end
    inhash.keys.each do | frag |
      ratio = inhash[frag][ratio_type]
      if ratios_hash.key?(ratio)
        ratios_hash[ratio] << frag
      else
        ratios_hash[ratio] = []
        ratios_hash[ratio] << frag
      end
    end
    ratios_hash
  end

end
