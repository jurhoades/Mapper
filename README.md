# Mapper
A utility for mapping a genomic coordinate to a specific position in a corresponding CDS or AA sequence.

## Requirements
- refFlat.txt file - UCSC annotation file specific to each species
- RefSeq ID - mRNA RefSeq Identifier (ex NM_146145)
- Genomic Coordinate - nucloetide location on chromosome (ex 101153495)

## Uses
Command Line Tool:
```
./mapper.py NM_146145 101153495 data/refFlatMm10.txt
```

Python Object:
```
from mapper import Mapper

mapper = Mapper("data/refFlatMm10.txt")
cds_pos, aa_pos = mapper.map(101153495, "NM_146145")
```

Tested using Python 2.7 and Python 3.5
