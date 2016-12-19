# PartsDB
SQLAlchemy extension for registries of biological parts

# Installation

The easiest way to install is by using pip:

```
pip install .
```

# Example

PartsDB assumes you already have some DB server installed. If you don't have one, try [PostgreSQL](http://postgresql.org).

Once installed, create a new database. In Linux or MacOS, this can be done by:

```
createdb [dbname]
```

Then you need to define some tables for your database. The system/Tables.py file has some useful [MixIns](http://docs.sqlalchemy.org/en/latest/orm/extensions/declarative/mixins.html) which can be used to define tables for your database. For example:

``` python

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

``` 

defines three tables for parts (uses PartMixIn): Promoter, CDS, Terminator; and a Gene table (only uses BaseMixIn), which is a particular combination of one of each.

Then, communication with your database is easy as:

1.	Importing PartsDB and your table definitions

``` python
from partsdb.partsdb import PartsDB
from tables import *
```

2.	Creating a connection to your database

``` python
partsdb = PartsDB('postgresql:///dbname', clean = True, Base = Base)
```
Where `dbname` is a name of the database you have created. The `clean = True` flag erases all rows in all entries.

3.	Adding some parts

``` python
promoter 	= partsdb.addPart('promoter', seq = "TTACGTA")
cds		 	= partsdb.addPart('cds', seq = "ATGTAA")
terminator 	= partsdb.addPart('terminator', seq = "GCGCGTA")

gene = partsdb.addPart('gene', promoter = promoter, cds = cds, terminator = terminator)

partsdb.commit()
```

Note, that the `seq` column is defined in `PartsMixIn` in system/Tables.py. The `commit()` is essential for new records to be added to the databse.

4.	Querying the databse

``` python
session = partsdb.Session()
query = session.query(Gene).first()
session.close()
```

5.	Accessing record information

``` python
print query.dbid
print query.promoter.seq + query.cds.seq + query.terminator.seq
```

Have a try with files in the example folder.
