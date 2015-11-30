#encoding: utf-8
require 'bio'
require 'bio-samtools'

##Open the vcf file and create lists of heterozygous and homozygous SNPs

class Stuff

	def self.safe_invert(hash)
    	hash.each_with_object( {} ) { |(key, value), out| ( out[value] ||= [] ) << key }
 	end

	##Input: Lists of hm and ht SNPs
	##Output: dictionaries with the id of the fragment as key and the absolute number of SNPs as value
	def self.create_hash_number(array)
    array.each_with_object(Hash.new(0)){|string, hash| hash[string] += 1}
	end

end
