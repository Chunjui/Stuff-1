#!/usr/bin/env python

# From: http://www.michael-noll.com/tutorials/writing-an-hadoop-mapreduce-program-in-python/

import sys
import math
def trunc(num, digits):
	sp = str(num).split('.')
	return (sp[0]+'.'+sp[1][:digits])

# input comes from STDIN (standard input)
for line in sys.stdin:
	# remove leading and trailing whitespace
	line = line.strip()
	# split the line into words
	words = line.split('\t')
	#if float(words[0]) > 0:
	#	xtmp=trunc(float(words[0]),1)
	#	#xtmp1=str(float(xtmp)+0.1)
	#else :
	#	xtmp=trunc(float(words[0])-0.1,1)
	#	#xtmp1=str(float(xtmp)+0.1)

	#if float(words[1]) > 0:
	#	ytmp=trunc(float(words[1]),1)
	#	#ytmp1=str(float(ytmp)+0.1)
	#else :
	#	ytmp=trunc(float(words[1])-0.1,1)
	#	#ytmp1=str(float(ytmp)+0.1)
	xtmp=math.floor(float(words[0])*10)/10
	ytmp=math.floor(float(words[1])*10)/10
	print '%.1f_%.1f\t%s' % (xtmp,ytmp,1)
	# increase counters
	#for word in words:
		# write the results to STDOUT (standard output);
		# what we output here will be the input for the
		# Reduce step, i.e. the input for reducer.py
		#
		# tab-delimited; the trivial word count is 1
		#xtmp=trunc(float(word))
		#print '%s\t%s' % (xtmp, 1)

