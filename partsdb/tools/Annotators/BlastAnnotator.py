from . import Annotator

class BlastAnnotator(Annotator):

	def annotate(self, db, **kwargs):
		self._annotate(db, kwargs['fileName'])
		
	def _annotate(self, db, fileName):
		inFile = open(fileName)
		
		BUFF_SIZE = 1000000
		n = 0

		buff = inFile.readlines(BUFF_SIZE)
		
		while buff:
			n += 1
			print buff[0] 
			
			# Read hits
			for line in buff:
				tabs = line.rstrip().split('\t')
		
				names = ["cdsID", "uniID", "qlen", "slen", "qstart", "qend", "tstart", "tend", "qcovs", "pident", "evalue", "proteinName", "origin",  "geneName"]

				data = dict( zip(names, tabs) )
				data["uniID"] = data["uniID"].split('|')[2]

				
				hit = self.cls()
				hit.uniID					= data["uniID"]
				hit.coverage				= data["qcovs"]
				hit.qLen	 				= data["qlen"]
				hit.tLen	 				= data["slen"]
				hit.coordinates				= ""
				hit.eVal	 				= float(data["evalue"])
				hit.proteinName	 			= data["proteinName"]
				hit.geneName	 			= data["geneName"]
				hit.origin	 				= data["origin"]
				hit.eVal 					= float(data["evalue"])

				hit.coordinates = "{0}:{1},{2}:{3},{4}".format(data["qstart"], data["qend"], data["tstart"], data["tend"], data["pident"])
				
				cds = db.session.query(self.cls.__targetclass__).filter(self.cls.__targetclass__.dbid == data['cdsID']).first()

				if cds:
					hit.target = cds
					db.session.add(hit)
				else:
					print "Failed to locate {0}".format(hitName.split('_')[0])
			db.commit()
			buff = inFile.readlines(BUFF_SIZE)