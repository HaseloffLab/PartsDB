
hits = {}

def annotateBlastp(cls, fileName, session):

	inFile = open(fileName)

	for line in inFile:
		tabs = line.rstrip().split('\t')
		
		names = ["cdsID", "uniID", "qlen", "slen", "qstart", "qend", "tstart", "tend", "qcovs", "pident", "evalue", "proteinName", "origin",  "geneName"]

		data = dict( zip(names, tabs) )
		data["uniID"] = data["uniID"].split('|')[2]

		name = "{0}_{1}".format(data["cdsID"], data["uniID"])
		
		if not name in hits:
			hits[name] = cls()
			hits[name].uniID				= data["uniID"]
			hits[name].coverage				= data["qcovs"]
			hits[name].qLen	 				= data["qlen"]
			hits[name].tLen	 				= data["slen"]
			hits[name].coordinates			= ""
			hits[name].eVal	 				= float(data["evalue"])
			hits[name].proteinName	 		= data["proteinName"]
			hits[name].geneName	 			= data["geneName"]
			hits[name].origin	 			= data["origin"]
			
		
		hits[name].eVal = min( hits[name].eVal, float(data["evalue"]) )

		coordinates = "{0}:{1},{2}:{3},{4};".format(data["qstart"], data["qend"], data["tstart"], data["tend"], data["pident"])
		hits[name].coordinates += coordinates

	inFile.close()

	for hitName, hit in hits.iteritems():
		print hitName
		cdsID = hitName.split('_')[0]
		uniID = hitName.split('_')[1]

		cds = session.query(cls.__targetclass__).filter(cls.__targetclass__.dbid == cdsID).first()

		if cds:
			hit.cds = cds
			session.add(hit)
		else:
			print "Failed to locate {0}".format(hitName.split('_')[0])

	fileName.close()