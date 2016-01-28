"""
Class built for easy access to VCF file info, including any VEP annotations.

Similar to PyVCF interface.
"""

__author__ = 'bernie'

import re
import gzip
import argparse
import sys
import os
import subprocess as sp
import errno
from collections import defaultdict
from minimal_representation import get_minimal_representation 

class vcf_reader:

	# read through VCF header, stopping after first variant line 
	def __init__(self, vcf, output=None, version=None):
		self.file_name = vcf
		self.vcf = gzip.open(vcf) if vcf.endswith('.gz') else open(vcf)

		if not output:
			self.output = output 
		elif output == sys.stdout:
			self.output = output
		elif output.endswith('.gz'):
			self.output = gzip.open(output, 'w')
		else:
			self.output = open(output, 'w')

		if version and not re.match('\d\.\d$', version):
			raise Exception, "Version %s invalid." % version
		
		self.metainfo = {}
		self.metaformat = {}
		# read up to, but not including first variant line of vcf (i.e. process header and stop at beginning)
		for line in self.vcf:
			if version and 'fileformat=VCF' in line:
				line = '##fileformat=VCFv%s\n' % version
			if output:
				self.output.write(line)
			
			line = line.strip()
			self.line = line
			# Reading header lines to get VEP and individual arrays
			if line.startswith('#'):
				line = line.lstrip('#')
				if line.startswith('FORMAT') or line.startswith('INFO'):
					info_tag = re.search('ID=(.*?)[,>]', line).group(1)
					num = re.search('Number=(.*?)[,>]', line).group(1)
					infotype = re.search('Type=(.*?)[,>]', line).group(1)
					descr = re.search('Description=(.*?)[,>]', line).group(1)
					if line.startswith('INFO'):
						self.metainfo[info_tag] = {'number' : num, 'type' : infotype, 'description' : descr}
					else:
						self.metaformat[info_tag] = {'number' : num, 'type' : infotype, 'description' : descr}
				if line.find('ID=CSQ') > -1:
					self.vep_field_names = line.split('Format: ')[-1].strip('">').split('|')
				if line.startswith('CHROM'):
					header = re.split("\t", line)
					self.sample_names = header[9:]
					self.header = dict(zip(header, range(len(header))))
					break								
				continue

		# read first variant line
		self.read_line()
		self.reading = True

	# return info from primary VCF columns (i.e. CHROM, POS, REF, ALT)
	def __getitem__(self, key):
		return self.fields[self.header[key]]

	# change info in primary VCF columns (i.e. CHROM, POS, REF, ALT)
	def __setitem__(self, key, val):
		self.fields[self.header[key]] = val

	def write_line(self, line=None):
		if line:
			self.output.write(line + '\n')
		else:
			self.output.write('\t'.join(self.fields) + '\n')


	def read_line(self):
		try:
			self.line = self.vcf.next().strip()
		except StopIteration:
			self.vcf.close()
			if self.output:
				self.output.close()
			self.reading =  False

		# safe splitting, happy with either tab character or white space
		# self.fields = filter(lambda f: f != "", re.split("\t", self.line))

		self.fields = re.split("\t", self.line)
		if len(self.header.keys()) != len(self.fields):
			raise Exception, "VCF header does not correspond with fields. Is your VCF file tab-seperated?"

		if 'FORMAT' in self.header:
			format_string = self.fields[self.header['FORMAT']].split(':')
			self.format = dict(zip(format_string, range(len(format_string))))
		else:
			self.format = None

		# store info field and VEP annotation info in their own seperate variables, if present in vcf
		if 'INFO' in self.header:
			info = self.header['INFO']
			if info == '.':
				self.info_field = None
				self.annotations = None
			else:
				self.info_kv = [(x.split('=', 1)) if '=' in x else (x, x) for x in re.split(';(?=\w)', self.fields[info])]
				self.info_field = dict(self.info_kv)
				if 'CSQ' not in self.info_field: 
					self.annotations = None
				else:
					# self.annotations = [dict(zip(self.vep_field_names, x.split('|'))) for x in self.info_field['CSQ'].split(',')]	
					self.annotations = [dict(zip(self.vep_field_names, x.split('|'))) for x in self.info_field['CSQ'].split(',') if len(self.vep_field_names) == len(x.split('|'))]
		else:
			self.info_field = None
			self.annotations = None

	
	def read(self, reset=False):
		while self.reading:
			# skip commented lines
			if self.line.startswith('#'):
				try:
					self.line = self.vcf.next().strip()
				except StopIteration:
					self.reading = False
					if self.output:
						self.output.close()
				continue
			yield self
			self.read_line()
		if reset:
			self.__init__(self.file_name)


	def get_variants(self):
		variants = []
		for site in self.read():
			for alt in site['ALT'].split(','):
				v = (site['CHROM'], site['POS'], site['REF'], alt)
				variants.append(v)
		return variants

	"""
	Requires that file be tabix-ed.
	"""
	def go_to(self, region):
		p = sp.Popen(['tabix', self.file_name, region], stdout=sp.PIPE, stderr=sp.PIPE)
		stdout, stderr = p.communicate()
		if p.returncode:
			raise SystemExit, "VCF file must be tabix-ed."
		if stdout == "":
			return
		self.vcf = (line for line in stdout.strip().split("\n"))

	"""
	Store individual gentoypes for current site in a dictionary
	"""
	def genotype(self):
		self.gt = {}
		for sample in self.sample_names:
			sample_info = self.fields[self.header[sample]]
			gt = sample_info.split(':')[self.format['GT']]
			self.gt[sample] = gt


	"""
	Helper function to locate a given variant within the vcf. Assumes nothing about representation of variant.
	Returns allele number, and sets current state of reader object to variant of interest.
	"""
	def seek_out_variant(self, chrom, pos, ref, alt):
		pos = int(pos)
		leeway = max([len(ref), len(alt)])
		target = get_minimal_representation(pos, ref, alt)
		found = False
		reg = tuple(map(str,(chrom, pos-leeway, pos+leeway)))
		region = "%s:%s-%s" % reg
		self.go_to(region)
		for line in self.read():
			pos = line['POS']
			ref = line['REF']
			alts = line['ALT'].split(",")
			# get min rep of alleles in full vcf so they correspond to alleles in sites table
			candidates = [get_minimal_representation(pos, ref, alt) for alt in alts]
			if target in candidates:
				found = True
				allele_number = candidates.index(target) + 1
				break
		if found:
			return allele_number
		else:
			return False

	"""
	Get hom_refs, hets, and hom_alts for given variant.
	"""
	def get_genotype_for_variant(self, chrom, pos, ref, alt, gq_threshold=float("-inf"), dp_threshold=float("-inf")):
		hom_refs = set()
		hom_alts = set()
		hets = set()
		
		allele_number = self.seek_out_variant(chrom, pos, ref, alt)
		if not allele_number:
			return False
		if not self.sample_names or not self.format:
			raise Exception, "Unable to parse individual genotypes. Fields missing."

		# get carriers of variant
		num_missing = 0
		for sample in self.sample_names:
			# filter non-calls
			sample_info = self.fields[self.header[sample]]
			gt = sample_info.split(':')[self.format['GT']]
			if gt == './.':
				continue
			
			gq = sample_info.split(':')[self.format['GQ']]
			dp = sample_info.split(':')[self.format['DP']]			
			if gq == '.' or dp == '.':
				num_missing += 1
				continue
			if int(gq) < gq_threshold or int(dp) < dp_threshold:
				continue

			lofs = sum(1 for a in map(int, re.split('/|\|', gt)) if a == allele_number)
			if lofs == 0:
				hom_refs.add(sample)
			elif lofs == 1:
				hets.add(sample)
			elif lofs == 2:
				hom_alts.add(sample)

		return hom_refs, hets, hom_alts
	

	def get_allele_frequency(self, allele_num, mode='info_field', gq_threshold=float("-inf"), dp_threshold=float("-inf")):
		# if nof indivual genotypes are present, look for allele counts in the info field
		if mode == 'query_info':
			allele_idx = int(allele_num) - 1
			if not self.info_field:
				raise Exception, "Unable to parse INFO field because it doesn't exist."
			
			if "AC_Adj" in self.info_field:
				# print allele_num, allele_idx, self.info_field["AC_Adj"]
				ac = self.info_field["AC_Adj"].split(',')[allele_idx]
				an = self.info_field["AN_Adj"]
			elif "AC" in self.info_field:
				ac = self.info_field["AC"].split(',')[allele_idx]
				an = self.info_field["AN"]
			else:
				raise Exception
			
			if float(an) == 0:
				return None
			else:
				return float(ac)/float(an)
		
		# if individual genotypes are present, use them
		elif mode == 'count_alleles':
			if not (self.format or self.sample_names):
				raise Exception, "Unable to parse individual genotypes. Check formatting."			
			an = 0
			ac = 0
			# get carriers of variant
			for sample in self.sample_names:
				# filter non-calls
				sample_info = self.fields[self.header[sample]]
				gt = sample_info.split(':')[self.format['GT']]
				if gt == './.':
					continue

				gq = sample_info.split(':')[self.format['GQ']]
				dp = sample_info.split(':')[self.format['DP']]
				if gq < gq_threshold or dp < dp_threshold:
					continue
				
				for allele in re.split('/|\|', gt):
					if allele == str(allele_num):
						ac += 1
					if allele != ".":
						an += 1
			if float(an) == 0:
				return
			return float(ac)/float(an)
		else:
			raise Exception, "Improper mode %s given to get_allele_frequency function."


	



