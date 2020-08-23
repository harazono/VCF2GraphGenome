#! /usr/bin/env python3
import sys
from optparse import OptionParser
import vcf
from Bio import SeqIO
from dataclasses import dataclass, field
from typing import List, Tuple, Dict

@dataclass
class node:
	pass

@dataclass
class node:
	identifier: int
	seq: str
	ref_chrom: str
	ref_chrom_pos: (int, int)
	edges: List[node]
	isRef: bool
	def node_print(self):
		print(f"identifier: {self.identifier}")
		print(f"ref_chrom: {self.ref_chrom}")
		print(f"ref_chrom_pos: {self.ref_chrom_pos}")
		print(f"isRef: {self.isRef}")
		print(f"edges: {[x.identifier for x in self.edges]}")
		print(f"seq: {self.seq}")
		print(f"\n")

def findAlleleinRef(node_array, pos):
	for each_node in node_array:
		if each_node.ref_chrom_pos[0] <= pos and each_node.ref_chrom_pos[1] >= pos:
			return each_node
	return None

def findRefNodeId(node_array):
	ret_array = []
	for each_node in node_array:
		if each_node.isRef:
			ret_array.append(each_node.identifier)
	ret_array = sorted(ret_array)
	return ret_array


def construct_graph(node_array, ref_dict, snp_dict):
	for each_chrom_id, refseq in ref_dict.items():
		snp_in_chrom = sorted(snp_dict[each_chrom_id], key = lambda x: x[0]) if each_chrom_id in snp_dict.keys() else None
		if snp_in_chrom is None:
			continue

		unique_snp_pos_list = sorted(set([x[0] for x in snp_in_chrom]))

		previous_snp_index = 0
		id_counter = 0
		for each_snp_pos in unique_snp_pos_list:
			pos = each_snp_pos - 1

			if previous_snp_index != pos:
				tmp_node1 = node(
						identifier = id_counter,
						seq = str(refseq[previous_snp_index:pos]),
						ref_chrom = each_chrom_id,
						ref_chrom_pos = (previous_snp_index, pos - 1),
						edges = [],
						isRef = True)
			else:
				tmp_node1 = node(
						identifier = id_counter,
						seq = str(refseq[previous_snp_index:pos + 1]),
						ref_chrom = each_chrom_id,
						ref_chrom_pos = (previous_snp_index, pos),
						edges = [], 
						isRef = True)

			node_array.append(tmp_node1)
			id_counter = id_counter + 1

			if previous_snp_index != pos:# In this case, this allele is not divided
				tmp_node2 = node(
						identifier = id_counter,
						seq = str(refseq[pos]),
						ref_chrom = each_chrom_id,
						ref_chrom_pos = (pos, pos), 
						edges = [], 
						isRef = True)
				id_counter = id_counter + 1
				node_array.append(tmp_node2)
			previous_snp_index = each_snp_pos


		last_node = node(identifier = id_counter, seq = refseq[previous_snp_index:], ref_chrom = each_chrom_id, ref_chrom_pos = (previous_snp_index, len(refseq)), edges = [], isRef = True)
		id_counter = id_counter + 1
		node_array.append(last_node)
		ref_node_length = len(node_array)

		for i in range(len(node_array) - 1):
			pass
			node_array[i].edges.append(node_array[i + 1])


		previous_snp_pos = unique_snp_pos_list[0]
		for each_snp in snp_in_chrom:
			pos = each_snp[0] - 2
			alt_seq = each_snp[2]
			previous_ref_node = findAlleleinRef(node_array, pos)
			next_ref_node = findAlleleinRef(node_array, pos + 2)
			tmp_node = node(identifier = id_counter, seq = alt_seq, ref_chrom = pos, ref_chrom_pos = (pos + 1, pos + 1), edges = [next_ref_node], isRef = False)
			previous_ref_node.edges.append(tmp_node)
			node_array.append(tmp_node)
			id_counter += 1
			previous_snp_pos = pos
		for i in range(ref_node_length, len(node_array)):
			pos = node_array[i].ref_chrom_pos[0]
			corresponding_ref_node = findAlleleinRef(node_array, pos)
			next_ref_node = findAlleleinRef(node_array, pos + 1)
			ref_node_len = len(corresponding_ref_node.edges)
			if ref_node_len > 1:
				node_array[i].edges.append(node_array[i + 1])


def convert2GFA1(node_array, filename = "otameshiOUTPUT.gfa"):
	try:
		with open(filename, "w") as f:
			f.write("H\tVN:Z:1.0\n")
			ref_node_id = ",".join([str(x) for x in findRefNodeId(node_array)])
			tmpstr = f"P\tREF\t{ref_node_id}\n"
			f.write(tmpstr)
			for each_node in node_array:
				if not isinstance(each_node.seq, list):
					seq = str(each_node.seq)
				else:
					seq = "".join([x.sequence for x in each_node.seq])
				tmpstr = f"S\t{each_node.identifier}\t{seq}\n"
				f.write(tmpstr)

				for each_edge in each_node.edges:
					tmpstr = f"L\t{each_node.identifier}\t+\t{each_edge.identifier}\t+\t0M\n"
					f.write(tmpstr)
	except FileNotFoundError as e:
		print(e)



def main():
	usage = f"usage: {__file__} なんか書く。"
	parser = OptionParser(usage = usage)
	parser.add_option("--vcf", dest="vcf_path", help = "Path to vcf file(required)")
	parser.add_option("--ref", dest="ref_path", help = "Path to reference genome fasta file(required)")
	(options, args) = parser.parse_args()
	vcf_path = options.vcf_path
	ref_path = options.ref_path
	if vcf_path is None:
		print(f"vcf_path is required")
		sys.exit(1)
	if ref_path is None:
		print(f"ref_path is required")
		sys.exit(1)

	ref_dict = {}
	for each_chrom in SeqIO.parse(ref_path, 'fasta'):
		ref_dict[each_chrom.id] = each_chrom.seq

	snp_dict = {}
	vcf_fh = vcf.Reader(open(vcf_path, 'r'))
	for each_snp in vcf_fh:
		if each_snp.CHROM in snp_dict.keys():
			snp_dict[each_snp.CHROM].append((each_snp.POS, each_snp.REF, each_snp.ALT))
			#print(type(each_snp.ALT[0]))
		else:
			snp_dict[each_snp.CHROM] = [(each_snp.POS, each_snp.REF, each_snp.ALT)]

	node_array = []
	construct_graph(node_array, ref_dict, snp_dict)
	convert2GFA1(node_array)



if __name__ == "__main__":
	main()