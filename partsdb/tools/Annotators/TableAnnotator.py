from . import Annotator

class TableAnnotator(Annotator):
	def annotate(self, db, **kwargs):
		self._annotate(db, kwargs['fileName'])

	def _annotate(self, db, fileName):

		inFile = open(fileName)
		BUFF_SIZE = 100000

		header = inFile.readline().rstrip()
		columns = header.split('\t')

		buff = inFile.readlines(BUFF_SIZE)

		while buff:

			for line in buff:
				if not line.startswith("#"):
					tabs = line.rstrip().split('\t')
					
					if len(tabs) == len(columns):
						dbid = tabs[0]

						target = db.session.query(self.cls.__targetclass__).filter(self.cls.__targetclass__.dbid == dbid).first()
						
						if target:
							row = self.cls()
							row.target = target
							for key, val in zip(columns[1:], tabs[1:]):
								setattr(row, key, val)
							db.session.add(row)
			db.commit()
			buff = inFile.readlines(BUFF_SIZE)
		inFile.close()
