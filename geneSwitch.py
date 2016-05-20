#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from __future__ import division


def load_csv(filename):
	dictGeneSwitch = {}
	with open(filename, 'r') as csvfile:
		lines = csv.reader(cgcsvFile, delimiter = '\t')
		next(lines, None)
		for line in lines:
			tissue = line[2].strip()
			if tissue in dictGeneSwitch:
				tdata = dictGeneSwitch[tissue]
			else:
				tdata = {}

			gene = line[0]
			tdata[gene] += [(int(line[1]), int(line[3]))]

			dictGeneSwitch[tissue] = tdata
	csvfile.close()

	for tissue in dictGeneSwitch:
		tdata = dictGeneSwitch[tissue]
		dictGeneSwitch[tissue] = sorted(tdata.iteritems(), key=lambda d:d[0])

	return(dictGeneSwitch)

def _has_open_switch(gene):
	for promoter in gene:
		if(promoter[1] == 0):
			return True
	return False

def _build_link(preGene, preGeneData, gene, promoter):
	link = '{"name":"flare.' + gene + '.P' + str(promoter[0]) + '",';
	imports = []
	for prePromoter in preGeneData:
		if prePromoter[1] == 0:
			imports += ['"flare.'+ preGene +'.P'+ str(prePromoter[0]) + '"']
	link += '"imports":[' + ','.join(imports) + ']}'
	return link

def generate_d3js_json(dictGeneSwitch):
	json = ''
	for tissue in dictGeneSwitch:
		tdata = dictGeneSwitch[tissue]
		genes = tdata.keys()
		genes.sort()
		for i, gene in enumerate(genes):
			
			# where this gene has open switches

			if !_has_open_switch(gene):
				continue

			# fine the previous gene

			preGene = None
			j = i - 1
			if(j < 0):
				j += len(genes)
			while j != i:

				if _has_open_switch(genes[j]):
					preGene = genes[j]
					break;

				j -= 1
				if(j < 0):
					j += len(genes)

			# build links

			links = []
			if(preGene):
				for promoter in tdata[gene]:
					if(promoter[1] == 0):
						links += [_build_link(gene, tdata[preGene], promoter)]

			return ('[' + ','.join(links) + ']')						
															
def write_csv(links, filename):
	try:
		csvFile = open(filename, 'w')
	except IOError:
		log.info('error: write to csv file "' + filename + '" failed!')
		sys.exit(-1)
	
	csvFile.write('\n'.join([format('%f\t%f' % (cg[0], cg[1])) for cg in cgdata]))		
	csvFile.close()

dictGeneSwitch = load_csv("/home/fsch/Project/CpGDenLowess/test/bonemarrow/multipromoter/geneSwitch.csv")
links = generate_d3js_json(dictGeneSwitch)

jsonFile = open("/home/fsch/Project/CpGDenLowess/test/bonemarrow/multipromoter/geneSwitch.d3js.json", 'w')
jsonFile.write(links)
jsonFile.close()