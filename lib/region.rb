

class Region
  attr_accessor :start, :stop, :chromosome, :splicename, :cell, :genestart, :genestop, :strand

  def initialize(start, stop, chromosome, splicename, cell, genestart, genestop, strand)
    @start = start.to_i
    @stop = stop.to_i
    @chromosome = chromosome
    @splicename = splicename
    @cell = cell
    @genestart = genestart.to_i
    @genestop = genestop.to_i
    @strand = strand
  end

  def to_s
    "#{@chromosome} #{@splicename} #{@start}..#{@stop} #{@cell}"
  end
end
