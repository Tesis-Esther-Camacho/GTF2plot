#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from __future__ import division
import sys
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
from collections import defaultdict


syntax = """------------------------------------------------------------------------------------------
Usage: python gtf2plot.py file.gtf*

*File containing gene (transcripts) annotations in gtf format (ensembl-like)
------------------------------------------------------------------------------------------
"""

#--Parameters control
if len(sys.argv) != 2:
	print syntax
	sys.exit()
	
#--Read the file and store transcripts information
gtf = open (sys.argv[1], 'r')
coord = [] #--To fix graph limits according to the gene length
trans = defaultdict(list) 
names = defaultdict(str)

new_id = ""
id = ""
exons = []
a = 1

for line in gtf:
	line = line.rstrip('\n')
	col = line.split('\t')[:-1]
	att = line.split('\t')[-1].split('"')[1::2]

	if '#' in line[0]:
		continue
	
	if col[2] == 'exon':
		new_id = att[1]

		if not id: #--Iniciador
			id = att[1]
			
		if new_id != id:
			#print a, id, exons
			trans[a] = exons #--List of tuples
			names[a] = id
			exons = []
			a += 1

		id = att[1]
		exons.append((int(col[3]), int(col[4])))
		coord.append(int(col[3]))
		coord.append(int(col[4]))
		
#--The last transcript
trans[a] = exons
names[a] = id

xo = sorted(coord)[0] #--gen start
xf = sorted(coord)[-1] #--gen end
gtf.close()

#print trans
#print xo, xf

#--text-box propierties
bbox_props = dict(boxstyle='rarrow', ec='c', fc='c')

#--Actions required to built a box
codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO,Path.CLOSEPOLY,]

#	Square coordinates order
#
#	x2,y2-->x3,y3	 
#	^			|
#	|			v
#	x1,y1<--x4,y4	
#

#--Figure definition
fig = plt.figure()
fig.set_size_inches(12,0.5*(len(trans)+1))
#fig.set_figheight(100)
#fig.set_figwidth(15)


ax = fig.add_subplot(111)

level = 1 #--Defines the graph height to draw each transcript
level = len(trans) #--Defines the graph height to draw each transcript
for t in sorted(trans.keys()):#, reverse = True):
	exons = sorted(trans[t]) #--exons sorted by coordinates
	#print t, names[t], exons

	#--Drawing lines between exones
	for i in range(0, len(exons)-1):
		
		#--We substrat gene start coordinate (xo) to all coordinates to set the origin of the graph
 		x = [exons[i][1]-xo, exons[i+1][0]-xo]
 		y = [level, level]
 		ax.plot(x, y, 'r-',linewidth=1)
		
		
	#--Drawing exons as boxes
	for i in range(0, len(exons)):
		
		#--We substrat gene start coordinate (xo) to all coordinates to set the origin of the graph
		x = [exons[i][0]-xo,exons[i][0]-xo,exons[i][1]-xo,exons[i][1]-xo, 0] #--exons[i][0], exons[i][1] = exon start and end coordinates
		y = [level-0.25, level+0.25, level+0.25,level-0.25,0] 
		v = zip(x,y) 
		
		path = Path(v, codes)
		patch = patches.PathPatch(path, facecolor='g', lw=0) 
		ax.add_patch(patch) 
		
	#--Adding the transcript name
	ax.text(xf-xo + 1000, level, names[t], size=10, ha="left", va="center", rotation=0)#, bbox = bbox_props)
		
	


	level = level - 1 #--Defines the graph height to draw each transcript

	
ax.set_ylim(0,len(trans)+1) #--Global height depends on the number of transcripts to draw
ax.set_xlim(-1000, xf-xo+5000) #--limits of the graph according to gene size, plus-minus a margin
ax.set_yticks([])
ax.set_xticks([])

plt.savefig(sys.argv[1].split('.')[0] + '.ps')
#plt.show()
