"""
Subset VCF based on VEP annotations.
"""

__author__ = 'bernie'

import argparse
import os
import re
import sys
from parse_vep_vcf import vcf_reader

def subset(args):
	tree = StemNode('root')
	build_tree(args.expr, tree)
	ensg = read_in_genes(args.ensg)
	hgnc = read_in_genes(args.symbol)
	
	vcf = vcf_reader(args.ivcf, output=args.ovcf)
	for site in vcf.read():
		if not site.annotations:
			print >> sys.stderr, "Skipping chr%s:%s %s -> %s because no VEP annotations present." % (site['CHROM'], site['POS'], site['REF'], site['ALT'])
			continue
		for a in site.annotations:
			if tree.eval(a) and check_gene(a, hgnc, ensg):
				site.write_line()
				break
			else:
				continue


def check_gene(a, hgnc, ensg):
	if not (hgnc or ensg):
		return True
	elif hgnc:
		return a['SYMBOL'] in hgnc
	elif ensg:
		return a['Gene'] in ensg


def read_in_genes(arg):
	if not arg:
		return []
	elif os.path.isfile(arg):
		g = open(arg, 'r')
		genes = [line.strip() for line in g]
	else:
		genes = arg.split(',')
	return genes


def parse_expression(e):
	i = a = 0
	conditions = []
	nesting = 0
	for char in e:
		if char == '(':
			nesting += 1
		elif char == ')':
			nesting -= 1
		elif char in ['&', '|'] and nesting == 0:
			conditions.append(e[a:i])
			a = i
		i += 1
	
	conditions.append(e[a:i])	
	conjunctions = []
	newconds = []
	for cond in conditions:
		if re.match('[|&]', cond):
			conjunctions.append(cond[0])
			newconds.append(cond[1:])
		else:
			newconds.append(cond)
	return newconds, conjunctions


class StemNode:
	def __init__(self, operator):
		self.oper = operator
		self.children = []
	
	def add_child(self, node):
		self.children.append(node)

	def eval(self, d):
		if self.oper == 'root':
			return self.children[0].eval(d)
		elif len(self.children) < 2:
			raise Exception, "Invalid operation in -e. Fewer than 2 arguments given to %s operator." % self.oper
		if self.oper == '|':
			return any(child.eval(d) for child in self.children)
		elif self.oper == '&':
			return all(child.eval(d) for child in self.children)
		else:
			raise Exception, "Invalid logical operator: %s" % self.oper


class LeafNode:
	def __init__(self, expr):
		if not expr:
			self.fn = lambda x: True
			self.checked = True
		else:
			self.fn = self.define_operator(expr)
			self.checked = False

	def eval(self, d):
		if not self.checked:
			if self.key not in d:
				raise Exception, "\"%s\" not in annotation field." % self.key
			self.checked = True 
		return self.fn(d)

	def define_operator(self, e):
		match = re.search('(!)?(=|CONTAINS|IN|~)', e)
		if not match:
			raise Exception, "Invalid logical expression: \"%s\". Valid operators are =, !, ~, IN, CONTAINS." % e
		
		neg, oper = match.groups()
		key, val = map(lambda s: s.strip(), e.split(match.group(0)))
		self.key = key

		if '=' == oper:
			fn = lambda d: d[key] == val
		elif 'CONTAINS' in e:
			fn = lambda d: val in d[key].split(',')
		elif 'IN' in e:
			fn = lambda d: d[key] in val.split(',')
		elif '~' in e:
			fn = lambda d: re.match(val, d[key])

		if neg:
			return lambda x: not fn(x)
		else: 
			return fn


def build_tree(expr, parent):
	if not expr:
		parent.add_child(LeafNode(expr))
		return
	clauses, conjuncts = parse_expression(expr)
	if '&' in conjuncts and '|' in conjuncts:
		raise Exception, "Invalid logical expression: \"%s\". Consider using parentheses." % expr
	# base case
	if len(clauses) == 1:
		newNode = LeafNode(clauses[0])
		parent.add_child(newNode)
	else:
		oper = next(iter(set(conjuncts)))
		newParent = StemNode(oper)
		parent.add_child(newParent)
		for clause in clauses:
			build_tree(clause.strip('( )'), parent=newParent)


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--ivcf', '-i', dest='ivcf', help='Input VCF.', required=True)
	parser.add_argument('--ovcf', '-o', dest='ovcf', help='Subsetted VCF.', default=sys.stdout)
	parser.add_argument('-e', '--expr', dest='expr', help="Logical expression using VEP annotation fields as inputs.\n" \
						"Valid operators: =, !, ~, IN, CONTAINS, |, &. Parentheses are also valid\n" \
						"e.g. \'Gene=ABCA1 & (Consequence CONTAINS missense_variant | Consequence ~ splice_.*_variant)\'\n"\
						"At least one annotation must be satisfied for a given site for it to be included in the output.")
	parser.add_argument('--ensg', help='Text file or comma-delimited list of ENSG gene IDs.')
	parser.add_argument('--symbol', help='Text file or comma-delimited list of HGNC-style gene symbols.')
	args = parser.parse_args()

	if args.ensg and args.symbol:
		parser.error('Pick either --ensg or --symbol. Not sure how to deal with both.') 
	if args.expr and ('&&' in args.expr or '||' in args.expr):
		parser.error('Invalid logical expression. Use & and | for logical AND/OR operations.')

	subset(args)

	# Testing code
	# e = 'Gene=DMD | Gene IN ADSG,AA1'
	# root = StemNode('root')
	# build_tree(e, root)
	# ann = {'Gene' : 'AA1', 'Consequence' : 'stop_gained_variant,missense_variant', 'AC' : '5'}
	# print root.eval(ann)



