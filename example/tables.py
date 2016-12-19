from sqlalchemy.ext.declarative import declared_attr
from sqlalchemy import Column, Integer, String, ForeignKey, Text, Float
from sqlalchemy.orm import relationship
from partsdb.system.Tables import Base, BaseMixIn, PartMixIn


class Promoter(Base,BaseMixIn,PartMixIn):
	pass

class CDS(Base,BaseMixIn,PartMixIn):
	pass

class Terminator(Base,BaseMixIn,PartMixIn):
	pass

class Gene(Base,BaseMixIn):
	name 			= Column( String(100) )
	promoterID  	= Column( Integer, ForeignKey('promoter.id') )
	cdsID  			= Column( Integer, ForeignKey('cds.id') )
	terminatorID  	= Column( Integer, ForeignKey('terminator.id') )

	promoter 		= relationship(Promoter, 	enable_typechecks=False)
	cds 			= relationship(CDS, 		enable_typechecks=False)
	terminator 		= relationship(Terminator, 	enable_typechecks=False)