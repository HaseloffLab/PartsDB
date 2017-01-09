from CoordinateMapper import CoordinateMapper
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.SeqRecord import SeqRecord

class RangeCoordinateMapper(CoordinateMapper):
	def __init__(self, selist, length, startOffset, endOffset):
		selist.location.parts.sort(key = lambda x: x.start)

		if selist.strand == -1:
			selist = selist._flip(length)
			self.strand = -1
		else:
			self.strand = 1
		super(RangeCoordinateMapper, self).__init__(selist)
		self.length = length
		self.startOffset = startOffset
		self.endOffset = endOffset
		
	# start - end in GenBank notation
	def rc2g(self, start, end, strand):

		locations = []

		
		start = self.c2g(start - self.startOffset, dialect = 'GenBank').to_genbank()
		end   = self.c2g(end - self.startOffset, dialect = 'GenBank').to_genbank()
		
		for exon in sorted(self._exons.parts, key = lambda x: x.start):
			eStart = int(exon.start) + 1
			eEnd   = int(exon.end)

			if start <= eStart:
				if end >= eStart and end <= exon.end:
					locations.append( FeatureLocation(eStart-1, end,strand) )
					break
				elif end > eEnd:
					locations.append( FeatureLocation(eStart-1, exon.end, strand) )
			
			elif start <= eEnd:
				if end <= eEnd:
					locations.append( FeatureLocation(start-1, end, strand) )
					break
				else:
					locations.append( FeatureLocation(start-1, eEnd, strand) )
		
		if strand == 1:
			locations.sort(key = lambda x: x.start)
		if strand == -1:
			locations.sort(key = lambda x: x.start, reverse = True)

		if len(locations) == 1:
			location = locations[0]
		else:
			location = CompoundLocation(locations)

		if self.strand == -1:
			return location._flip(self.length)

		return location
