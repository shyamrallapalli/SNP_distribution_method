#encoding: utf-8
require_relative 'stuff'

class SDM
  ##Inputs: hashes with IDs as keys and the SNP density as value
  ##Divide absolute number of SNPs by the length of the given fragment.
  ##Output: hashes with IDs as keys and the normalised SNP density as value
  def self.normalise_by_length(lengths, shuf)
    shuf_norm = {}
    shuf.each_key { | frag |
      # each fragment snp frequency is normalized to the length of the fragment
      perbp = shuf[frag].to_f / lengths[frag].to_f
      shuf_norm.store(frag, perbp)
    }
    shuf_norm = Stuff.safe_invert(shuf_norm)
    return shuf_norm
  end

  def self.divide_array(dic_hm_inv, right, left, keys_hm, dest, cross)
    contigs_at_min = []
    minimum = keys_hm.min #define minimum score
    contigs_at_min << dic_hm_inv.values_at(minimum) #look for the contigs with the minimum score
    contigs_at_min.flatten!
    keys_hm.delete(minimum) #delete entries for contigs with the minimum homozygous score from the original hash

    #shuffle the contigs with the minimum homozygous scores on the two halfs of the expected normal distribution.
    #(1) If the number of contigs is even, half of the array goes to the left and half goes to the right part of the distribution
    #(2) If the number of contigs is odd, we define two situations:
    ##if there's only 1 contig, we randomly assign the right or left destination for it
    ##if there are more than 1 contigs, we take the first contig in the array and randomly assign the right or left destination for it.
    # The remaining array have a even number of elements, so we proceed as described in (1)
    ## even number of contigs
    if contigs_at_min.length.to_i % 2 == 0
      half = contigs_at_min.each_slice(contigs_at_min.length/2.to_i).to_a
      right << half[0]
      left << half[1]
    ## odd number of contigs
    else
      ## odd and more than 2 contigs
      if contigs_at_min.length.to_i > 2
        object = contigs_at_min.shift
        other_half = contigs_at_min.each_slice(contigs_at_min.length/2.to_i).to_a
        right << other_half[0]
        left << other_half[1]
        right << object
      else
        ## 1 contig, so choose randomly between left or right
        if dest == 0
          right << contigs_at_min
        elsif dest == 1
          left << contigs_at_min
        end
      end
    end
    return right, left, keys_hm
  end


  ##Input for sorting: inverted hash containing the normalised homozygous scores as key and the fragments' ids as value.
  ##Output from sorting: perm is the permutation array containing the candidate order of contigs
  # and mut is the array of candidate contigs taken from the central part of perm
  def self.sorting(dic_hm_inv, cross, average_contig)
    left, right = [], []
    keys= dic_hm_inv.keys.to_a
    (keys.length/2).times do #repeat the sorting process until the original hash is sorted.
      right, left, keys = SDM.divide_array(dic_hm_inv, right, left, keys, 0, cross)
      right, left, keys = SDM.divide_array(dic_hm_inv, right, left, keys, 1, cross)
    end
    right.compact.flatten!
    left.compact.flatten!

    # check which of right or left are smaller array
    smaller = []
    if right.length <= left.length
      smaller = right
    else
      smaller = left
    end

    perm = right + left.reverse # combine together both sides of the distribution
    perm.flatten!
    mut = []
    # pl = perm.length.to_i
    # r = pl/50
    if average_contig.to_i < 10000
      # In case of a back-cross, 20 contigs in the middle part of the permutation taken
      if cross == "back"
        if smaller.length > 10
          mut << right[-10, 10]
          mut << left[-10, 10].reverse
          mut.flatten!
        else
          mut << right[-right.length, right.length]
          mut << left[-left.length, left.length].reverse
          mut.flatten!
        end
      # In case of a out-cross, 50 contigs in the middle part of the permutation taken
      elsif cross == "out"
        if smaller.length > 20
          mut << right[-20, 20]
          mut << left[-20, 20].reverse
          mut.flatten!
        else #If a strong filtering step reduces the total number of contigs to a number lower than 20, perm.length/2 contigs on the right and perm.length/2 on the left side of the middle point are taken.
          mut << right[-right.length, right.length]
          mut << left[-left.length, left.length].reverse
          mut.flatten!
        end
      end
    else
      num = smaller.length
      if num > 6
          num = 6
      end
      mut << right[-num, num]
      mut << left[-num, num].reverse
      mut.flatten!
    end
    return perm, mut
  end

  ###Input for calling: inverted hashes with homozygous densities and ratios as key and ids as value, the type of cross and the hash containing the ids and the hm SNP positions.
  ###Output from calling: contig permutation based on homozygous SNP, contig permutation based on ratios, array of candidate contigs taken from the central part of both permutations (mut) and
  #candidate positions for causal mutations  (hyp) - this is used to prove the efficiency of SDM as we know the correct order and the position of the mutation,
  #in a real experiment this will not be known -
  def self.calling_SDM (dic_hm_inv, dic_ratios_inv, dic_pos_hm, cross, average_contig)
    mut, number_of_snps = [], []
    or_pos = {}
    perm_hm, mut  = SDM.sorting(dic_hm_inv, cross, average_contig) #sorting step based on homozygous SNP density -score-
    perm_ratio, mut_ratio = SDM.sorting(dic_ratios_inv, cross, average_contig) #repeat the sorting step based on the ratios
    mut << mut_ratio #merge together both candidate contig arrays into one, and remove duplications.
    mut.flatten!
    mut = mut.uniq
    #***Testing step***
    mut.each { |frag| #identify the SNP positions in the candidate contigs contained in mut
        if dic_pos_hm.has_key?(frag)
          or_pos.store(frag, dic_pos_hm[frag])
        end
      }
    or_pos.each do |frag, array|
        number_of_snps << array.length
    end
    #or_pos.delete_if { |id, array|  array.length < (number_of_snps.max - 1) }
    hyp = or_pos.values.sort.flatten!
    return perm_hm, perm_ratio, mut, hyp
  end
end
