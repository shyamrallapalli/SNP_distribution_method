# encoding: utf-8

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

end