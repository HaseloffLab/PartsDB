from partsdb.partsdb import PartsDB
from tables import *

partsdb = PartsDB('postgresql:///partsdb', clean = True, Base = Base)

promoter 	= partsdb.addPart('promoter', seq = "TTACGTA")
cds		 	= partsdb.addPart('cds', seq = "ATGTAA")
terminator 	= partsdb.addPart('terminator', seq = "GCGCGTA")

gene = partsdb.addPart('gene', promoter = promoter, cds = cds, terminator = terminator)

partsdb.commit()

session = partsdb.Session()

query = session.query(Gene).first()

session.close()

print query.dbid
print query.promoter.seq + query.cds.seq + query.terminator.seq
