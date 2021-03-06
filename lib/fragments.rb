# encoding: utf-8

class Fragments

  DEFAULT = {
      filter_out_low_hmes: false,
      cross: 'back',
  }

  #shuffle the contigs with the minimum homozygous scores on the two halves of the expected normal distribution.
  #(1) If the number of contigs is even, half of the array goes to the left and half goes to the right part of the distribution
  #(2) If the number of contigs is odd, we define two situations:
  ##if there's only 1 contig, we randomly assign the right or left destination for it
  ##if there are more than 1 contigs, we take the first contig in the array and randomly assign the right or left destination for it.
  # The remaining array have a even number of elements, so we proceed as described in (1)
  def self.split(ratios_hash, left, right, ratio_keys, input_frags, right_or_left)
    minimum = ratio_keys.min # get minimum score
    contigs_at_min = ratios_hash.values_at(minimum).flatten # look for the contigs with the minimum score

    hmnum_hash = {}
    contigs_at_min.each do |frag|
      hm_num = input_frags[frag][:hm]
      if hmnum_hash.key?(hm_num)
        hmnum_hash[hm_num] << frag
      else
        hmnum_hash[hm_num] = []
        hmnum_hash[hm_num] << frag
      end
    end
    hmnum_hash.keys.sort.each do |number|
      num_frags = hmnum_hash[number].length
      half = num_frags/2
      ## odd number of contigs, pick odd one out
      if num_frags % 2 == 1
        one_element = hmnum_hash[number].shift
        ## odd contig, so choose randomly between left or right
        if right_or_left == 0
          left << one_element
        else
          right << one_element
        end
      end
      # if the array is empty proceed
      unless hmnum_hash[number].empty?
        left << hmnum_hash[number][0, half]
        right << hmnum_hash[number][half..-1]
      end
      # update right or left value for next iteration
      right_or_left = (right_or_left - 1).abs
    end

    #delete entries for contigs with the minimum homozygous score from the original hash
    ratio_keys.delete(minimum)
    [left, right, ratio_keys, right_or_left]
  end


  # Input for sorting: inverted hash containing the normalised homozygous scores as key and the fragments' ids as value.
  # Output from sorting: perm is the permutation array containing the candidate order of contigs
  # and mut is the array of candidate contigs taken from the central part of perm
  def self.arrange(ratios_hash, input_frags)
    right_or_left = 0
    left, right = [], []
    ratio_keys= ratios_hash.keys.to_a
    iterations = (ratio_keys.length.to_f/2.0).round
    iterations.times do #repeat the sorting process until the original hash is sorted.
      left, right, ratio_keys, right_or_left = Fragments.split(ratios_hash, left, right, ratio_keys, input_frags, right_or_left)
      if ratio_keys.length >= 1
        left, right, ratio_keys, right_or_left = Fragments.split(ratios_hash, left, right, ratio_keys, input_frags, right_or_left)
      end
    end

    # flatten the split arrays and reverse right side array
    # to bring the key fragments to the middle
    left.flatten!
    # right array is flattened first before reversing
    (right.flatten!).reverse!
    # instead of reversing first before flattening
    # (right.reverse!).flatten!
    # both should work, but would to lead change of order
    # for fragments with same ratio value

    # combine together both sides of the distribution
    perm = left + right
    perm
  end

  def self.select_fragments(ratios_hash, perm, adjust, opts = {})
    opts = DEFAULT.merge(opts)
    filter_out_low_hmes = opts[:filter_out_low_hmes]
    cross = opts[:cross]
    # set minimum cut off ratio to pick fragments with variants
    # calculate min hme score for back or out crossed data
    # if no filtering applied set cutoff to 1
    if filter_out_low_hmes
      if cross == 'back'
        cutoff = (1.0/adjust) + 1.0
      else
        cutoff = (2.0/adjust) + 1.0
      end
    else
      cutoff = 1.1
    end

    frags_to_keep = []
    ratios_hash.each_key { | ratio |
      # store fragment id which have more homozygosity
      if ratio >= cutoff
        frags_to_keep << ratios_hash[ratio]
      end
    }
    frags_to_keep.flatten!

    mutfrags = []
    perm.each do |id|
      if frags_to_keep.include?(id)
        mutfrags << id
      end
    end
    mutfrags
  end

end
