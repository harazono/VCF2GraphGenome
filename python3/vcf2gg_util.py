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
