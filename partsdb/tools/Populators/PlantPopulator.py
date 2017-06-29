from . import Populator

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from BCBio.GFF import GFFExaminer
from BCBio import GFF
from pprint import pprint
from Bio.SeqFeature import SeqFeature, CompoundLocation, FeatureLocation
import sys

from ..CoordinateMapper import RangeCoordinateMapper

class PlantPopulator(Populator):

	def _locationToCoordinates(self, location):
		if isinstance(location, CompoundLocation):
			return ";".join( [ "{0},{1},1".format( part.start, part.end ) for part in sorted(location.parts, key = lambda part: part.start) ] )

		elif isinstance(location, FeatureLocation):
			return  "{0},{1},1".format( location.start, location.end )

	def populate(self, *args, **kwargs):
		if isinstance(args[1],list):
			self.populateFromList(args[0])
		elif isinstance(args[1], str):
			self.populateFromFile(*args)

	def populateFromFile(self, transMapFileName, transFileName, proteinFileName, genomeFileName):
		transMapFile 	= open(transMapFileName)
		transFile 		= open(transFileName)
		proteinFile   	= open(proteinFileName)
		genomeFile 		= open(genomeFileName)
		
		scaffoldDict = SeqIO.to_dict( SeqIO.parse(genomeFile, "fasta") )
		genomeFile.close()

		transMapRecords = GFF.parse(transMapFile, base_dict = scaffoldDict)

		transcripts = {}

		# Extracting transcript locations

		for record in transMapRecords:
			for feature in record.features:
				target = feature.qualifiers["Target"][0]
				transcriptName  = target.split()[0]

				feature.qualifiers 	= { }
				feature.type 		= "mRNA"
				feature.id 			= transcriptName 
				if not transcriptName in transcripts:
					transcript = SeqRecord( record.seq, annotations = {"Locus" : record.id} )
					transcript.features = [ feature ]
					transcripts[transcriptName] = transcript
					transcripts[transcriptName].annotations["minTarget"] = int(target.split()[1])
					transcripts[transcriptName].annotations["maxTarget"] = int(target.split()[2])
				else:
					transcripts[transcriptName].features[0].location += feature.location
					transcripts[transcriptName].annotations["minTarget"] = min([ transcripts[transcriptName].annotations["minTarget"], int(target.split()[1]) ])
					transcripts[transcriptName].annotations["maxTarget"] = max([ transcripts[transcriptName].annotations["maxTarget"], int(target.split()[2]) ])

		transMapFile.close()


		transcripRecords = SeqIO.to_dict(SeqIO.parse(transFileName, "fasta"))

		# Getting transcript lengths and setting offset

		for transcriptName in transcripRecords:
			if transcriptName in transcripts:
				length = len(transcripRecords[transcriptName])
				transcripts[transcriptName].annotations["startOffset"] = transcripts[transcriptName].annotations["minTarget"] - 1
				transcripts[transcriptName].annotations["endOffset"]   = length - transcripts[transcriptName].annotations["maxTarget"]

		genes = {}
 		bad = 0

		for line in proteinFile:
			if line.startswith('>'):
				substring = line.split()[-1]
				transcriptName = substring.split(':')[0]
				pepName 	   = line.split()[0][1:]
				pepStart = int(substring.split(':')[1].split('-')[0])
				pepEnd   = int(substring.split(':')[1].split('-')[1].split('(')[0] )
				pepStrand = substring.split('(')[1][0]

				if pepStrand == '-':
					pepStrand = -1
				else:
					pepStrand = 1

				if transcriptName in transcripts:
					transcript = transcripts[transcriptName]

					gene = transcript

					#Cutting the gene
					start = gene.features[0].location.start
					end   = gene.features[0].location.end

					cutStart = max( 0, start - 3000 )
					cutEnd  = min( len(gene), end + 3000 )

					gene = gene[cutStart:cutEnd]
					gene.annotations = transcript.annotations
					gene.annotations["LocusCoordinates"] = "{0}:{1}-{2}".format( gene.annotations["Locus"], cutStart, cutEnd )
					

					# Annotating CDS
					exons = gene.features[0]

					rcm = RangeCoordinateMapper(exons, len(gene), gene.annotations["startOffset"], gene.annotations["endOffset"] )
					
					if rcm.startOffset + 1 > pepStart or len(exons) < pepEnd:
						continue   

					location = rcm.rc2g(pepStart, pepEnd, pepStrand)	

					cdsFeature = SeqFeature(type = 'cds', location = location )

					protSeq = cdsFeature.extract(gene.seq).translate()
					
					if transcriptName in genes:
						if len( genes[transcriptName].features[0] ) >= len(cdsFeature):
							continue

					if len(protSeq) == 0:
						continue

					if not (protSeq[0] == 'M' and protSeq.find('*') == len(protSeq)-1):
						bad = bad + 1
						continue

					gene.features.append(cdsFeature)

					if pepStart > 1:
						location = rcm.rc2g( rcm.startOffset + 1, pepStart-1, pepStrand )
						utrType = 'utr5' if pepStrand == 1 else 'utr3'
						utrLFeature = SeqFeature( type = utrType, location = location )
					else:
						utrLFeature = None


					if pepEnd < len(exons):
						location = rcm.rc2g( pepEnd+1, len(exons), pepStrand )
						utrType = 'utr3' if pepStrand == 1 else 'utr5'
						utrRFeature = SeqFeature( type = utrType, location = location )
					else:
						utrRFeature = None


					gene.features.append(utrRFeature)
					gene.features.append(utrLFeature)
					
					# Annotating Promoter / Terminator
					location  = FeatureLocation( 0, exons.location.start, cdsFeature.location.strand ) 
					partType  = 'promoter' if cdsFeature.location.strand == 1 else 'terminator'
					rightPart = SeqFeature(type = partType, location = location)

					location  = FeatureLocation( exons.location.end, len(gene), cdsFeature.location.strand )
					partType  = 'terminator' if cdsFeature.location.strand == 1 else 'promoter'
					leftPart  = SeqFeature(type = partType, location = location)

					gene.features.append(rightPart)
					gene.features.append(leftPart)

					genes[transcriptName] = gene

		proteinFile.close()
		print "Number of protein sequences that failed to import: ", bad
		self.populateFromList(genes)

	def populateFromList(self, genes):
		assert 'locus' 		in self.db.classes
		assert 'promoter' 	in self.db.classes
		assert 'terminator' in self.db.classes
		
		Locus 		= self.db.classes['locus']
		Promoter	= self.db.classes['promoter']
		Terminator	= self.db.classes['terminator']
		Gene 		= self.db.classes['gene']

		for geneName, gene in genes.iteritems():

			locusCoordinates = gene.annotations["LocusCoordinates"]

			locus =  self.db.session.query( Locus ).filter( Locus.coordinates == locusCoordinates ).first()
			if not locus:
				locus = self.db.addPart('locus', coordinates = locusCoordinates)

			mRNAStart 	= min(gene.features[0].location)
			mRNAStop 	= max(gene.features[0].location)

			parts = {}

			strand = 0

			for feature in gene.features:
				if feature:

					if feature.type == 'cds':
						strand = feature.location.strand

					if feature.location.strand == 1:
						location = feature.location._shift( -feature.location.start )
						seq = str(gene.seq[ feature.location.start : feature.location.end ])
					else:
						location = feature.location._flip( len(gene) )
						location = location._shift( -location.start )
						seq = str(gene.seq[ feature.location.start : feature.location.end ].reverse_complement())
					parts[ feature.type ] = { "seq" : seq.upper(), "coordinates" : self._locationToCoordinates(location) }


			gene = self.db.session.query(Gene).filter(Gene.locusID == locus.id, Gene.locusStrand == strand).first()

			if gene:
				promoter 	= gene.promoter
				terminator 	= gene.terminator
				if not promoter:
					print 'No promoter for gene ', gene.id
					sys.exit()
			else:
				promoter = self.db.addPart('promoter', seq = parts['promoter']['seq'])
				terminator = self.db.addPart('terminator', seq = parts['terminator']['seq'])


			cds  = self.db.addPart('cds', seq = parts['cds']['seq'], coordinates = parts['cds']['coordinates'] )
			
			utr5 = None
			utr3 = None

			if 'utr5' in parts:
				utr5 = self.db.addPart('utr5', seq = parts['utr5']['seq'], coordinates = parts['utr5']['coordinates'] )
			if 'utr3' in parts:
				utr3 = self.db.addPart('utr3', seq = parts['utr3']['seq'], coordinates = parts['utr3']['coordinates'] )

			newGene = self.db.addPart('gene', cds = cds, promoter = promoter, terminator = terminator, locus = locus, locusStrand = strand, transcriptName = geneName, utr5 = utr5, utr3 = utr3)			

		self.db.commit()	