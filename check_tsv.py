# -*- coding: utf-8 -*-

import sys,os

if __name__ == '__main__':
	
	filename = sys.argv[1]
	
	with open(filename,'r') as f:
		line = f.readline()
		sys.stdout.write(line)
		L = len(line.split('\t'))
		for line in f:
			l = len(line.split('\t'))
			if l != L:
				sys.stderr.write("ERROR:col=" + str(l) + ":" + line)
			else:
				sys.stdout.write(line)
	exit(0)
