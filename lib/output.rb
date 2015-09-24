#encoding: utf-8

##Input1: permutation array after SDM
##Input2: fasta file converted to array
##Input3: list of fragment ids from the shuffled fasta file
##Output: permutation of fragments after SDM with the data (lengths, etc) obtained from the original fasta file
class Output
	def self.create_perm_fasta(perm, frags)
		fasta_perm = []
    fragments = {}
		frags.each do |i|
      fragments[i.entry_id] = i.entry
    end
    perm.each { |frag|
      if fragments.has_key?(frag)
        fasta_perm << Bio::FastaFormat.new(fragments[frag].to_s)
      end
    }
		return fasta_perm
	end
end
