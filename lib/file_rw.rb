# encoding: utf-8
require 'bio'
require 'bio-samtools'
require 'bio-gngm'
require 'yaml'

class FileRW


  # Input: file location
  # Output: array with each line of the file as an entry (strings)
  def self.to_array(file)
    IO.foreach(file).collect {|l| l.chomp }
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

  # deep copy hash
  def self.deep_copy_hash(in_hash)
    tempname = Time.now.to_f.to_s + '.yml'
    File.open("#{tempname}", 'w') do |file|
      file.write in_hash.to_yaml
    end
    out_hash = YAML.load_file(tempname)
    %x[rm #{tempname}]
    out_hash
  end

  # Input: FASTA file
  # Output: hash of sequence ids with lengths and sequences
  def self.fasta_parse(fasta_file)
    sequences = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
    Bio::FastaFormat.open(fasta_file).each do |inseq|
      sequences[:seq][inseq.entry_id] = inseq.entry
      sequences[:len][inseq.entry_id] = inseq.length
    end
    sequences
  end

  # Input1: permutation array of frag ids after SDM
  # Input2: sequences hash of frag ids and a filename to write
  # Output will be written to file and if no filename is given
  # output will be written to "ordered_frags.fasta"
  def self.write_order(array, fasta_file, filename='ordered_frags.fasta')
    sequences = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
    Bio::FastaFormat.open(fasta_file).each do |inseq|
      sequences[inseq.entry_id] = inseq.entry
    end
    File.open(filename, 'w+') do |f|
      array.each do |frag|
        element = Bio::FastaFormat.new(sequences[frag].to_s)
        seqout = Bio::Sequence::NA.new(element.seq).upcase
        f.puts seqout.to_fasta(element.definition, 80)
        #command = "#{samtools} faidx #{fasta_file} \"#{frag}\""
        #f.puts `#{command}`
      end
    end
  end

end
