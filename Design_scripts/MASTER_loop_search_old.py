## code to iterate the MASTER loop searching automatically
# Example syntax: "python MASTER_loop_search.py #####.pds 10 21 1.5" will search for gapLen's starting from 10 and will run the wGap until length 21 with RMSD 1.5

import numpy as np
import sys
import math
import os

query=(sys.argv[1])
start=int(sys.argv[2]) # grabs the initial starting gapLen
extend=int(sys.argv[3]) # grabs the length to query
rmsdCut=(sys.argv[4]) # grabs the rmsdCut per normal MASTER search
gapLen_list=list(range(start,(extend+1)))
for i in gapLen_list:
	gapLen=i
	os.system("~/Nick/master --query %s --targetList ~/Nick/biounits/20180719/pds/list --rmsdCut %s --gapLen %d --matchOut %d/match.txt --structOut %d --outType wgap --topN 200" % (query, rmsdCut, gapLen, gapLen, gapLen))
