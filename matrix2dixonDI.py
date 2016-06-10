#!/usr/bin/env python

import argparse
import numpy as np
import re

parser = argparse.ArgumentParser(description='Compute the Dixon et al. (2012) directionality index')
parser.add_argument('-i', help='input matrix', type=str, required=True, dest='i')
parser.add_argument('-w', help='size of upstream/downstream window', type=int, default=2000000, dest='w')
parser.add_argument('-b', help='bin size (resolution) of input matrix', type=int, default=40000, dest='b')
parser.add_argument('-n', help='should NAs be zero', type=bool, default=False, dest='n')
args = parser.parse_args()

def load_matrix(matrix):
	'''
	Parse input Dekker format matrix
	'''
	IN = open(matrix, 'r')
	raw_x= []
	is_header = True
	for line in IN:
		if '#' not in line:
			if is_header:
				header = line.split('\t')[1:]
				is_header = False
			else:
				raw_x.append(map(lambda x: np.nan if 'nan' in x or 'NA' in x else float(x), line.split('\t')[1:]))
	X = np.array(raw_x)

	return X, header

def get_chr_start_end(colname):
	'''
	Get chromosome and bp start and end
	'''
	searchObj = re.search('(chr[\dXY]+):(\d+)-(\d+)', colname)
	return searchObj.group(1), searchObj.group(2), searchObj.group(3)

	

def enough_space(i, w, l):
	'''
	Check if enough space upstream and downstream to calculate DI
	'''
	if i- w < 0 or i + w >= l:
		return False
	else:
		return True

def calc_DI(A, B, E):
	'''
	Directionality Index
	'''
	if B-A == 0 or E == 0:
		return 'nan'
	else:
		DI = ( (B-A)/abs(B-A) ) * ( (((A-E)**2)/E) + (((B-E)**2)/E) )
		return DI

def nan_in_window(i, bin_window, row):
	'''
	Check for missing values in upstream and downstream windows
	'''
	if np.isnan(row[i-bin_window:i]).any() or np.isnan(row[i+1:i+bin_window+1]).any():
		return True
	else:
		return False

def main():
	OUT = open('_'.join([args.i[:-7], 'dixonDI', 'w-' + str(args.w)]), 'w')
	OUT.write('\t'.join(['chrom', 'start', 'end', 'DI']) + '\n')
	bin_window =  args.w / args.b
	X, header = load_matrix(args.i)
	if args.n:
		X = np.nan_to_num(X)
	for i, row in enumerate(X):
		chrom, start, end = get_chr_start_end(header[i])
		if enough_space(i, bin_window, len(row)):
				if not nan_in_window(i, bin_window, row):
					# Calculate upstream (A) and downstream (B) interactions
					A = np.sum(row[i-bin_window:i])
					print A
					B = np.sum(row[i+1:i+bin_window+1])
					print B
					# Expected number of reads (E)
					E = (A + B) / 2		
					# Calculate directionality index (DI)
					DI = calc_DI(A, B, E)
					if DI != 'nan':
						OUT.write('\t'.join([chrom, start, end, str(round(DI,4))]) + '\n')
					else:
						OUT.write('\t'.join([chrom, start, end, 'nan']) + '\n')			
				else:
					OUT.write('\t'.join([chrom, start, end, 'nan']) + '\n')
		else:
			OUT.write('\t'.join([chrom, start, end, 'nan']) + '\n')


if __name__ == '__main__':
	main()