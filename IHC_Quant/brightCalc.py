from __future__ import division
import os, sys
import cv2
import numpy as np

print "IHC start"

CELL_CUTOFF = 209
STAIN_CUTOFF = 147

# total cell range count
def count(img, cutoff):
	img = img[:,:,2]
	sumd = len(img[img<cutoff])
	return sumd

# loops thru each dir
fd = open("percs.txt", 'w')
results = ""
ind = 0
di = os.walk('.').next()[1]
s = [st.split(" ")[1] for st in di]
for d in [y for (x, y) in sorted(zip(s,di))]:

	results += "\n" + d + "\n"
	# sort on scan number
	files = os.walk("./" + d).next()[2]
	end = [e.split(" ")[1] for e in files]

	# loops thru each file
	for im in [f for (e,f) in sorted(zip(end,files))]:
		img = cv2.imread(d + "/" + im)
		img = cv2.cvtColor(img, cv2.COLOR_BGR2HSV)

		st = count(img, STAIN_CUTOFF)
		ce = count(img, CELL_CUTOFF)
		per = st / ce

		results += im + "  " + str(per) + "\n"

	ind += 1
	print("done with " + str(ind) + " of " + str(len(di)))

fd.write(results)
fd.close()
