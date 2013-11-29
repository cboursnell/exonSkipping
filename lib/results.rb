#!/usr/bin/env ruby

#require 'rubygems'
require 'trollop'

opts = Trollop::options do
  version "v0.0.1a"
  opt :input, "Output from skipping2.rb", :required => true, :type => String
  opt :verbose, "Be verbose"
end

Trollop::die :input, "must exist" if !File.exist?(opts[:input]) if opts[:input]


hash = Hash.new
File.open("#{opts.input}", "r").each_line do |line|
  line.chomp!
  cols = line.split(/\t/)
  transcript = cols[0]
  genename = transcript.split(".").first
  cell = cols[1]
  if hash.has_key?(genename)
    if hash[genename] == "35S" and cell=="GDC"
      hash[genename] = "both"
    elsif hash[genename] == "GDC" and cell=="35S"
      hash[genename] = "both"
    elsif hash[genename] == "" and (cell=="35S" or cell=="GDC")
      hash[genename] = cell
    end
  else
    # if cell!=""
    hash[genename] = cell
    # end
  end
end

hash.each_pair do |key, value|
  puts "#{key}\t#{value}"
end