####
## for parents, change GTid range
####

import sys
markerfile = sys.argv[1]
gtfile = sys.argv[2]
hapsize = sys.argv[3]

from itertools import izip
marker = open(str(markerfile))
gt = open(str(gtfile))

for lineM, lineGT in izip(marker,gt):
  bitsM = lineM.split()
  if bitsM[15] == 'nochange':
    print lineGT,
    continue
  elif bitsM[15] == 'change':
    bitsGT = lineGT.split()
    for GTid in range(9,9+int(hapsize)):
      if bitsGT[GTid] == '0':
        bitsGT[GTid] = '1'
      elif bitsGT[GTid] == '1':
        bitsGT[GTid] = '0'
      else:
        print 'check gt format'
        break
    print '\t'.join(bitsGT)
  else:
    print 'check change/nochange'
    break

    
