#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Profile the Matcher Module
"""
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
import os
import re
import pysam
import array
import csv

import numpy as np

DIR_PATH = os.path.dirname(os.path.realpath(__file__))
TEST_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),"..", "tests"))

# - Create directories if not exists
SEQUENCE_FOLDER = os.path.join(DIR_PATH, 'seq')
SUFFIX_ARRAY_FOLDER = os.path.join(DIR_PATH, 'sarray')
if not os.path.exists(SEQUENCE_FOLDER):
	os.makedirs(SEQUENCE_FOLDER)
if not os.path.exists(SUFFIX_ARRAY_FOLDER):
	os.makedirs(SUFFIX_ARRAY_FOLDER)

""" 
	Split every chromosome sequence into chunks before processing.
	This is the maximum length of bases for each chunk
"""
TRUNK_SIZE = 3000000

GENOME_FILE  = '../tests/data/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz'

def create_suffix_array():
	""" read genome sequence """
	seq_files = os.listdir(SEQUENCE_FOLDER)
	seq_files.sort()
	for seq_file in seq_files:
		seq_file_full_path = os.path.join(SEQUENCE_FOLDER, seq_file)
		f = open(seq_file_full_path, 'r')
		sequence = f.read()
		f.close()
		seq_len = len(sequence)
		name = re.search(r'chromosome_(.+)_', seq_file).group(1)	
		print "sequence length : %s" % (seq_len)
		print "processing chromosome %s" % (name)	
	
		"""
			METHOD 1
			
			library ksa by Mark Mazumder
		"""
		from sa import ksa
		sequence_len = len(sequence)

		# to cover missing checks for tails in every trunks (except the final one)
		
		trunk_transit_base_size = 150
		counter = 1

		while TRUNK_SIZE * (counter-1) < sequence_len:
			end = TRUNK_SIZE * counter + trunk_transit_base_size
			start = TRUNK_SIZE * (counter-1)
			if end <= sequence_len:
				ranks = ksa(sequence[start:end])
			else:
				ranks = ksa(sequence[start:])
			""" store the suffix array in hard disk in compressed txt format """ 
			
			filename = 'chromosome_' 
			filename += str(int(name))
			filename += '_partition_' 
			filename += str(counter)
			filename += '_suffix_array.txt.gz'

			full_file_name = os.path.join(SUFFIX_ARRAY_FOLDER, filename)
			np.savetxt(full_file_name, ranks)

			print "Suffix array for chunk %s in chromosome %s saved" % (str(counter), 
																		str(name), )
			
			counter += 1
		print "Chromosome %s suffix array written" % (name, )

		"""
			METHOD 2
			
			SA-IS algorithm
			
			source: https://github.com/AlexeyG/PySAIS
		"""
		# import pysais
		# sa = pysais.sais(sequence)
		# filename = 'chromosome_'+str(name)+'_suffix_array.txt.gz'
		# full_file_name = os.path.join(DIR_PATH, 
		# 								'sarray', 
		# 								filename)
		# np.savetxt(full_file_name, sa)
		# print "Chromosome %s suffix array written" % (name, )

def binary_search(pattern, SA, text):
	"""
		This is a binary search for only one letter (base)
	"""
	l = 0; r = len(SA) - 1
	while l < r:
		mid = int((l+r) / 2)
		if pattern > text[SA[mid]]:
			l = mid + 1
		else:
			r = mid
	s = l; r = len(SA) - 1
	while l < r:
		mid = int((l+r) / 2)
		if pattern < text[SA[mid]]:
			r = mid
		else:
			l = mid + 1
	return (s, r)

def separate_chromosome_sequence():
	"""
		This is to separate human genome sequence 
		from different chromosomes into different files
	"""
	with pysam.FastxFile(GENOME_FILE) as genfile:
		for gen in genfile:
			name = gen.name
			sequence = gen.sequence
			filename = 'chromosome_'+str(name)+'_sequence.txt.gz'
			full_file_name = os.path.join(SEQUENCE_FOLDER, filename)
			text_file = open(full_file_name, "w")
			text_file.write(sequence)
			text_file.close()
			print "Chromosome %s get" % (name, )

def profile_genome_search(pattern, maxmismatch):
	# - Preset variables
	pat_len = len(pattern)
	chrom_pos_tuple = []
	# - read genome sequence
	seq_files = os.listdir(SEQUENCE_FOLDER)
	seq_files.sort()
	for seq_file in seq_files:
		# - search in every chromosome
		seq_file_full_path = os.path.join(SEQUENCE_FOLDER, seq_file)
		f = open(seq_file_full_path, 'r')	
		sequence = f.read()
		f.close()
		name = re.search(r'chromosome_(.+)_', seq_file).group(1)
		must_have_words_in_file = 'chromosome_'+str(name)
		files = []
		# - find all suffix array files for this chromosome
		for file in os.listdir(SUFFIX_ARRAY_FOLDER):
			if must_have_words_in_file in file:
				files.append(file)
		# - for each truck of chromosome suffix array, match patterns
		for file in files:
			partition = int(re.search(r'partition_(\d+)_', file).group(1))
			# - read suffix array
			file = str(os.path.join(SUFFIX_ARRAY_FOLDER, file))
			SA = np.loadtxt(file,dtype=np.int32)
			# - start and end position
			start_pos = TRUNK_SIZE * (partition-1)
			end_pos = start_pos + len(SA)
			# - get sequence
			sequence_part = sequence[start_pos:end_pos]
			"""
				# - special treatment for last trunk
				# - empty space calculated in SA 
				# - just remove it
			"""
			if len(sequence_part) < len(SA):
				SA = SA[1:]
			""" 
				*** Initializaiton ***
				Mismatch counter array
		 	"""
			mismatches = np.array([0] * (len(sequence_part)))
			counter = 0
			
			# - Look for candidate interval character by character in pattern
			for base in pattern:
				""" Binary search """
				(left, right) = binary_search(base, SA, sequence_part)
				print "%s th base %s search" % (counter+1, base)
				print "range (%s, %s)" % (left, right)
				""" 
					Update mismatch counter 
					# - binary search for current base in SA that is transformed
					# - before every comparison SA corresponding to 1:counter position
						in original sequence are removed because they represent no
						real SA representing this base
					# - when result returned positions have to be shifted back (minus counter)
						to keep synchronised with original SA
				"""
				mismatches[SA[:left] - counter] = mismatches[SA[:left] - counter] + 1
				mismatches[SA[right:] - counter] = mismatches[SA[right:] - counter] + 1
				""" get suffix array corresponding to next position of base """
				SA = np.delete(SA, np.where(SA == counter))
				counter += 1
			""" 
				Remove last positions with length of pattern 
				because they are too short to match.
			"""
			mismatches = mismatches[:len(mismatches)-pat_len]
			# - remove candidates that exceed maxmismatch
			result = np.where(mismatches <= maxmismatch)[0]
			"""
				Output to a list of tuples
				(chromosome name, mismatch number, start position, sequence)
				Remember to add one to output position as python count from 0.
			"""
			
			for pos in result:
				chrom_pos_tuple += [
										(
											name, 
											mismatches[pos], 
											pos+start_pos+1, 
											sequence[pos+start_pos:pos+start_pos+pat_len]
										)
									]

			print "Chromosome %s chunk %s have been searched" % (str(name),
																partition)
		print "Chromosome %s have all been searched" % (str(name))

	filename = 'Pattern_Matches_For_%s.txt' % (pattern, )
	full_file_path = os.path.join(TEST_PATH, filename)
	with open(full_file_path, 'w') as csvfile:
		writer = csv.writer(csvfile)
		for match_pair in chrom_pos_tuple:
			writer.writerow(match_pair)



def main():
	separate_chromosome_sequence()
	create_suffix_array()
	profile_genome_search(pattern=b'TGGATGTGAAATGAGTCAAG', maxmismatch=3)

if __name__ == '__main__':
	main()