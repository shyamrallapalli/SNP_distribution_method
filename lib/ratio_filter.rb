#encoding: utf-8
require 'csv'

class RatioFilter

  def self.get_ratios(inhash)
    ratios = []
    inhash.each_key do | frag |
      ratios << inhash[frag][:ratio]
    end
    ratios
  end

  def self.ratio_hash(inhash)
    hash = {}
    inhash.keys.each do | frag |
      ratio = inhash[frag][:ratio]
      if hash.key?(ratio)
        hash[ratio] << frag
      else
        hash[ratio] = []
        hash[ratio] << frag
      end
    end
    hash
  end

  def self.selected_ratios(inhash, threshold)
    initial = inhash.keys.length
    ratios = get_ratios(inhash)

    if threshold > 0
      thres = 100.0/threshold.to_f
      filter = (ratios.max.to_f)/thres
      # go through selection if the threshold of discarded fragments is low
      sel_ids = inhash.keys
      contigs_discarded = initial - sel_ids.length
      while sel_ids.length > 30 * contigs_discarded do
        puts "threshold #{threshold}%"
        # delete fragment information below selected threshold
        inhash.keys.each do | frag |
          if inhash[frag][:ratio] <= filter
            inhash.delete(frag)
          end
        end
        sel_ids = inhash.keys
        contigs_discarded = initial - sel_ids.length
        puts "#{contigs_discarded} contigs out of #{initial} discarded"
        # increase threshold in steps of 2% until while condition is reached
        threshold = threshold + 2
      end
      ratios_hash = ratio_hash(inhash)
    else
      ratios_hash = ratio_hash(inhash)
    end
    return inhash, ratios_hash
  end

end
