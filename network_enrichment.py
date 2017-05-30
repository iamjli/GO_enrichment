
import sys, os
import cPickle as pkl

import pandas as pd
import numpy as np

from enrichment import *



def load_interactome(pkl_file):

	if os.path.exists(pkl_file):
		print "Loading interactome gene list from pickled file:", pkl_file
		return pkl.load(open(pkl_file, 'rb'))

	else:
		print "Loading interactome gene list fresh. This will take a while."
		symbols = [f.split("\t")[0:2] for f in open("../../ALS/data/interactome/iRefIndex_v13_MIScore_interactome.txt", 'r').readlines()]
		symbols = list(set([item for sublist in symbols for item in sublist]))

		entrez = hgnc2entrez(symbols)
		pkl.dump(entrez, open(pkl_file, 'wb'))

		return entrez


def load_gene_list(node_attr_file, cutoff):

	node_attributes = pd.read_csv(node_attr_file, sep="\t")
	robust_nodes = node_attributes[node_attributes["robustness"] >= cutoff]["protein"].tolist()

	print "---%d genes found---" %len(robust_nodes)

	return hgnc2entrez(robust_nodes)


def main():

	node_attributes_file = sys.argv[1]
	cutoff 				 = float(sys.argv[2])
	base_dir 			 = os.path.dirname(node_attributes_file)

	# Load background file from interactome
	background = load_interactome("iRefIndex13_entrez.pkl")

	# Load gene list
	gene_list = load_gene_list(node_attributes_file, cutoff)

	# Perform GO enrichment analysis
	enr = Enrichment(background)
	enr.run_analysis(gene_list, "%s/GO_enrichment_%s.tsv" %(base_dir, str(cutoff)))



main()