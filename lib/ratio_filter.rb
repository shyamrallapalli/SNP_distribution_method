#encoding: utf-8
require 'csv'

class RatioFilter

  DEFAULT = {
      only_frag_with_vars: true,
  }

  def self.prune_hash(inhash, adjust)
    first = inhash.length
    inhash.each_key do | frag |
      numhm = inhash[frag][:hm]
      numht = inhash[frag][:ht]
      # selecting fragments which have a variant
      if numht + numhm <= 2 * adjust
        inhash.delete(frag)
      end
    end
    last = inhash.length
    warn "Discarded #{first-last} out of #{first} fragments, which lack any variant\n"
    inhash
  end

  def self.selected_ratios(inhash, adjust, opts = {})
    opts = DEFAULT.merge(opts)
    only_frag_with_vars = opts[:only_frag_with_vars]
    ratios_hash = {}
    # select only fragments with variants
    # this is to discard fragments which may not contain
    # any information about causative mutations
    if only_frag_with_vars
      inhash = prune_hash(inhash, adjust)
    end
    inhash.keys.each do | frag |
      ratio = inhash[frag][:ratio]
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
