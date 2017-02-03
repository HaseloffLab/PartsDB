from . import Annotator

class BlastAnnotator(Annotator):

	def annotate(self, db, **kwargs):
		self._annotate(db, kwargs['fileName'])
		
	def _annotate(self, db, fileName):
		
		hits = {}
		session = db.Session()

		with open(fileName) as inFile:
			for line in inFile:
				tabs = line.rstrip().split('\t')
				
				names = ["cdsID", "uniID", "qlen", "slen", "qstart", "qend", "tstart", "tend", "qcovs", "pident", "evalue", "proteinName", "origin",  "geneName"]

				data = dict( zip(names, tabs) )
				data["uniID"] = data["uniID"].split('|')[2]

				# name = "{0}_{1}".format(data["cdsID"], data["uniID"])

				cds = session.query(self.cls.__targetclass__).filter(self.cls.__targetclass__.dbid == data["cdsID"]).first()

				if not cds:
					print "Failed to locate {0}".format(cds.dbid)
					continue

				hit = session.query(self.cls).filter(self.cls.uniID == data["uniID"])\
											.filter(self.cls.targetID == cds.id).first()

				if not hit:
					hit = self.cls()
					hit.uniID				= data["uniID"]
					hit.coverage			= data["qcovs"]
					hit.qLen	 			= data["qlen"]
					hit.tLen	 			= data["slen"]
					hit.coordinates			= ""
					hit.eVal	 			= float(data["evalue"])
					hit.proteinName	 		= data["proteinName"]
					hit.geneName	 		= data["geneName"]
					hit.origin	 			= data["origin"]
					hit.target 				= cds
					session.add(hit)


				hit.eVal = min( hit.eVal, float(data["evalue"]) )

				coordinates = "{0}:{1},{2}:{3},{4};".format(data["qstart"], data["qend"], data["tstart"], data["tend"], data["pident"])
				hit.coordinates += coordinates
		session.commit()