import sys, os

import pandas as pd
import numpy as np

import mygene

from goatools.base import download_go_basic_obo, download_ncbi_associations
from goatools.obo_parser import GODag
from goatools.associations import read_ncbi_gene2go
from goatools.go_enrichment import GOEnrichmentStudy


def hgnc2entrez(l):
	"""
	Converts HGNC gene symbols to NCBI's Entrez 
	"""
	mg = mygene.MyGeneInfo()
	out = mg.querymany(l, scopes='symbol', fields='entrezgene', species='human')
	l_entrez = [entry["entrezgene"] for entry in out if "entrezgene" in entry]
	return l_entrez


class Enrichment():
	def __init__(self, background, sig_cutoff=0.05):
		self.background = background
		self.alpha = sig_cutoff

		self.obodag, self.geneid2gos_human = self.load_ontologies_and_associations()

		self.goeaobj = self.create_enrichment_study()

	def load_ontologies_and_associations(self):
		print "---LOADING ONTOLOGIES AND ASSOCIATIONS---"
		# Check if files exist and download if not
		obo_fname = download_go_basic_obo()
		gene2go = download_ncbi_associations()

		# Load ontologies and associations
		obodag = GODag(obo_fname)
		geneid2gos_human = read_ncbi_gene2go("gene2go", taxids=[9606])
		print "{N:,} annotated human genes".format(N=len(geneid2gos_human))

		return obodag, geneid2gos_human

	def create_enrichment_study(self):
		obj = GOEnrichmentStudy(
			self.background, 			# List of human protein-coding genes (Entrez IDs)
			self.geneid2gos_human, 		# Gene ID / GO associtations
			self.obodag, 				# Ontologies
			propagate_counts = False,
			alpha = self.alpha, 		# Significance cutoff
			methods = ['fdr_bh']		# Multiple hypothesis correction
		)

		return obj

	def run_analysis(self, gene_list, path):
		"""
		Returns enrichments below significance cutoff and writes to file
		"""
		goea_all = self.goeaobj.run_study(gene_list)
		goea_sig = [r for r in goea_all if r.p_fdr_bh < self.alpha]

		if path.endswith(".xlsx"): 	
			self.goeaobj.wr_xlsx(path, goea_sig)
		elif path.endswith(".tsv"):
			self.goeaobj.wr_tsv(path, goea_sig)
		else:		
			self.goeaobj.wr_txt(path, goea_sig)			
		
		return goea_sig
