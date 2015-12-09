# encoding: utf-8
require 'bio'

class FileRW

  # Input 0: File to create/add a line to
  # Input 1: Something to add to a new line of the file
  def self.append(filename, line)
    File.open(filename, 'a') do |file|
      file.puts line
    end
  end

  # Input 0: Filename by which to save an array with filetype extension, one value per line
  # Input 1: Array to save
  def self.write(filename, array)
    File.open("#{filename}", 'w+') do |f|
      array.each { |i| f.puts(i) }
    end
  end

  # Input 0: Filename by which to save an array to a .txt file, one value per line
  # Input 1: Array to save
  def self.write_txt(filename, array)
    write("#{filename}.txt", array)
  end

  # Input: file location
  # Output: array with each line of the file as an entry (strings)
  def self.to_array(file)
    IO.foreach(file).collect(&:chomp)
  end

  # Input: file location (file must have integers on each line)
  # Output: array with each line of the file as an entry, converted to integer
  def self.to_array_int(file)
    to_array(file).collect { |l| Integer(Float(l)) }
  rescue ArgumentError => e
    $stderr.puts "Not all lines in file can be converted to ints: #{e}"
    exit
  end

  # Input: file location (file must have floats on each line)
  # Output: array with each line of the file as an entry, converted to float
  def self.to_array_float(file)
    to_array(file).collect { |l| Float(l) }
  rescue ArgumentError => e
    $stderr.puts "Not all lines in file can be converted to floats: #{e}"
    exit
  end

  # def self.safe_invert(hash)
  #   hash.each_with_object( {} ) { |(key, value), out| ( out[value] ||= [] ) << key }
  # end

  ##Input: Lists of hm and ht SNPs
  ##Output: dictionaries with the id of the fragment as key and the absolute number of SNPs as value
  def self.create_hash_number(array)
    array.each_with_object(Hash.new(0)){|string, hash| hash[string] += 1}
  end

  # Input: FASTA file
  # Output: hash of sequence ids with lengths and sequences and total assembly length
  def self.fasta_parse(fasta_file)
    sequences = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
    assembly_len = 0
    Bio::FastaFormat.open(fasta_file).each do |inseq|
      sequences[:seq][inseq.entry_id] = inseq.entry
      sequences[:len][inseq.entry_id] = inseq.length
      assembly_len += inseq.length
    end
    return sequences, assembly_len
  end

  # Input1: permutation array of frag ids after SDM
  # Input2: sequences hash of frag ids and a filename to write
  # Output will be written to file and if no filename is given
  # output will be written to "ordered_frags.fasta"
  def self.write_order(perm, frags, filename='ordered_frags.fasta')
    File.open(filename, 'w+') do |f|
      perm.each do |frag|
        element = Bio::FastaFormat.new(frags[frag].to_s)
        seqout = Bio::Sequence::NA.new(element.seq).upcase
        f.puts seqout.to_fasta(element.definition, 80)
      end
    end
  end

end
