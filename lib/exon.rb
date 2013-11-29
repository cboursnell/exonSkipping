
class Exon
  attr_accessor :start, :stop, :type, :cell

  def initialize(start, stop, type)
    @start = start.to_i
    @stop = stop.to_i
    @type = type
    @cell = ""
  end

  def contains?(position)
    if position >= @start and position < @stop
      return true
    end
  end

  # a.overlap(b)
  def overlap(exon)
    if @start == exon.start and @stop == exon.stop
      return 0
    elsif @start < exon.start
      if @stop > exon.start
        if @stop>exon.stop
          return 4
        else
          return 2
        end
      else
        return 5
      end
    else
      if exon.stop > @start
        if exon.stop >= @stop
          return 3
        else
          return 1
        end
      else
        return 6
      end
    end
  end

  def length
    return stop - start + 1
  end

  def to_s
    "#{@start}..#{@stop} #{@type}"
  end
end