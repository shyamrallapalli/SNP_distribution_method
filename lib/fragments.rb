#encoding: utf-8

class Fragments

  #shuffle the contigs with the minimum homozygous scores on the two halves of the expected normal distribution.
  #(1) If the number of contigs is even, half of the array goes to the left and half goes to the right part of the distribution
  #(2) If the number of contigs is odd, we define two situations:
  ##if there's only 1 contig, we randomly assign the right or left destination for it
  ##if there are more than 1 contigs, we take the first contig in the array and randomly assign the right or left destination for it.
  # The remaining array have a even number of elements, so we proceed as described in (1)
  def self.split(ratios_hash, left, right, ratio_keys, dest)
    minimum = ratio_keys.min # get minimum score
    contigs_at_min = ratios_hash.values_at(minimum).flatten # look for the contigs with the minimum score
    num_frags = contigs_at_min.length
    half = num_frags/2
    ## even number of contigs
    if num_frags % 2 == 0
      left << contigs_at_min[0, half]
      right << contigs_at_min[half..-1]
    ## odd number of contigs
    else
      ## odd and more than 2 contigs
      if contigs_at_min.length.to_i > 2
        # add odd fragment to left and then add equal halves to both
        left << contigs_at_min.shift
        left << contigs_at_min[0, half]
        right << contigs_at_min[half..-1]
      else
        ## 1 contig, so choose randomly between left or right
        if dest == 0
          left << contigs_at_min
        elsif dest == 1
          right << contigs_at_min
        end
      end
    end
    #delete entries for contigs with the minimum homozygous score from the original hash
    ratio_keys.delete(minimum)
    [left, right, ratio_keys]
  end


  ##Input for sorting: inverted hash containing the normalised homozygous scores as key and the fragments' ids as value.
  ##Output from sorting: perm is the permutation array containing the candidate order of contigs
  # and mut is the array of candidate contigs taken from the central part of perm
  def self.arrange(ratios_hash, cross, average_contig)
    left, right = [], []
    ratio_keys= ratios_hash.keys.to_a
    iterations = (ratio_keys.length.to_f/2.0).round
    iterations.times do #repeat the sorting process until the original hash is sorted.
      left, right, ratio_keys = Fragments.split(ratios_hash, left, right, ratio_keys, 0)
      if ratio_keys.length >= 1
        left, right, ratio_keys = Fragments.split(ratios_hash, left, right, ratio_keys, 1)
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

    # check which of right or left are smaller array
    rlen = right.length
    llen = left.length
    smaller = rlen <= llen ? right : left
    num = smaller.length

    mut = []
    if average_contig.to_i < 10000
      # In case of a back-cross, 20 contigs in the middle part of the permutation taken
      if cross == 'back'
        if num > 10
          mut = left[-10, 10] + right[0, 10]
        else
          mut = left[-llen, llen] + right[0, rlen]
        end
      # In case of a out-cross, 40 contigs in the middle part of the permutation taken
      elsif cross == 'out'
        if num > 20
          mut = left[-20, 20] + right[0, 20]
        # If a strong filtering step reduces the total number of contigs to a number lower than 20,
        # perm.length/2 contigs on the right and perm.length/2 on the left side of the middle point are taken.
        else
          mut = left[-llen, llen] + right[0, rlen]
        end
      end
    else
      if num > 6
          num = 6
      end
      mut = left[-num, num] + right[0, num]
    end
    [perm, mut]
  end

  def self.select_fragments(cross, ratios_hash, perm, adjust, threshold)

    # set minimum cut off ratio to pick fragments with variants
    # calculate min ratio for back or out crossed data
    # and apply any threshold use applied
    if cross == 'back'
      min = (1.0/adjust) + 1.0
    else
      min = 2 * ((1/adjust) + 1)
    end
    maximum = ratios_hash.keys.max.to_f
    range = maximum - min
    thres = 100.0/threshold.to_f
    cutoff = min + (range/thres)

    # delete fragments which are below selected cutoff ratio
    frags_to_keep = []
    ratios_hash.each_key { | ratio |
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
