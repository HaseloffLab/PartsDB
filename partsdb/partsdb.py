from sqlalchemy.orm import sessionmaker, mapper
from sqlalchemy import Table, MetaData, create_engine
from sqlalchemy.ext.declarative import declarative_base, declared_attr

from system.Tables import BaseMixIn, Base, Sys
from system.IDGenerator import nextIDGenerator


class PartsDB(object):

	Session 	= sessionmaker()
	newParts	= []

	def __init__(self, address, Base, clean = False, idGenerator = nextIDGenerator):
		self.engine = create_engine(address, echo = False)
		self.Session.configure(bind=self.engine)
		self.Base = Base

		self.session = self.Session()

		self.classes = {}

		for cls in Base.__subclasses__():
			if issubclass(cls, BaseMixIn):
				self.classes[cls.__tablename__] = cls

		self.idGenerator = idGenerator

		if clean:
			self.Base.metadata.drop_all(self.engine)

		Base.metadata.create_all(self.engine, checkfirst=True)

	def setup(self, **kwargs):
		session = self.Session()
		
		for key, value in kwargs.iteritems():
			val = session.query(Sys).filter(Sys.variable == key).first()
			if val:
				val.value = value
			else:
				val = Sys( variable = key, value = value )
			session.add(val)
		session.commit()
		session.close()

	def _getSysVal(self, key):
		session = self.Session()
		val = session.query(Sys.value).filter(Sys.variable == key).first()
		session.close()
		return val[0]

	def annotate(self, tableName, **kwargs):
		annotator = self.classes[tableName].annotator
		annotator.annotate(self, **kwargs)

	def addPart(self, table, *args, **kwargs):
		if args:
			return self._addPartFromCls(args[0])
		else:
			newPart = self.classes[table](**kwargs)
			self.session.add(newPart)
			return newPart

	def _addPartFromCls(self, part):
		self.session.add(part)
		return part

	def commit(self):
		self.session.commit()
		try:
			prefix = self._getSysVal('prefix')
		except:
			prefix = 'partsdb'
		
		for clsName, cls in self.classes.iteritems():
			unNamed = self.session.query(cls).filter(cls.dbid == None).all()
			for obj in unNamed:
				obj.dbid = "{0}.{1}.{2}".format(prefix, clsName, obj.id)
				self.session.add(obj)
		self.session.commit()

