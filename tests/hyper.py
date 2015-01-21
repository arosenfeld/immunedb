import sys
from sldb.identification.v_genes import VGermlines

germs = VGermlines(sys.argv[1])
for v in germs:
    print v, '|'.join(sorted(germs.get_ties([v], 300, .05)))
