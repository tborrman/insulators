#!/usr/bin/env python

import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Compute the Dixon et al. (2012) directionality index')
parser.add_argument('-i', help='input matrix', type=str, required=True, dest='i')
parser.add_argument('-w', help='size of upstream/downstream window', type=int, default=2000000, dest='w')
parser.add_argument('-b', help='bin size (resolution) of input matrix', type=int, default=40000, dest='b')
args = parser.parse_args()

def enough_space(i, w, l):
	'''
	Check if enough space upstream and downstream to calculate DI
	'''
	if i- w < 0 or i + w >= l:
		return False
	else:
		return True

def calc_DI(A, B, E):
	if B-A == 0 or E == 0:
		return 'nan'
	else:
		DI = ( (B-A)/abs(B-A) ) * ( (((A-E)**2)/E) + (((B-E)**2)/E) )
		return DI




def main():
	x = np.array([[np.nan,1,4,2,1],
			  [1,3,1,1,2],
			  [4,1,4,0,3],
			  [2,1,0,2,4],
			  [1,2,3,4,5]])

	for i, row in enumerate(x):

		if enough_space(i, args.w, len(row)):
				if not np.isnan(row).any():
					# Calculate upstream (A) and downstream (B) interactions
					A = np.sum(row[i-args.w:i])
					B = np.sum(row[i+1 :i+args.w + 1])
					# Expected number of reads (E)
					E = (A + B) / 2		
					# Calculate directionality index (DI)
					DI = calc_DI(A, B, E)
					print DI
				else:
					print 'nan'

		else:
			print 'nan'








if __name__ == '__main__':
	main()