#!/bin/env python

import sys,os,itertools
import numpy as np

class Protein:
	""" A class for reading and manipulating protein structures"""
	def __init__(self,filename=None):
		if filename:
			self.name = filename
			self.read_coords()
			self.reslist = self.resid.keys()
			self.nres = len(self.reslist)

	def read_coords(self):
		raw = filter(lambda x: x.split()[0] == "ATOM", 
				open(self.name,"r").readlines())
		resid = {}
		res = -1
		resseq = []
		for r in raw:
			rs = r.split()
			if rs[3:5] != resseq:
				resseq = rs[3:5]
				res +=1
				resid[res] = {}
				rtype = rs[3]
				resid[res]['type'] = rtype
				resid[res]['resseq'] = int(resseq[1])
				resid[res]['coords'] = []
			atype = rs[2]
			for atom in ['C','O','N','S']:
				if atom in atype:
					xyz = (float(rs[5]),float(rs[6]),float(rs[7]))
					resid[res]['coords'].append((atype,xyz))
					break
		self.resid = resid

	def calc_rco(self,cutoff):
		""" calculate contact order [Plaxco, Simons & Baker, JMB 1998] """
		rc2 = cutoff**2
		aco = 0
		ncont = 0
		for i in self.reslist:
			xyz_i = self.resid[i]['coords']
			jreslist = filter(lambda x: i - x >= 1 ,self.reslist)
			for j in jreslist:
				xyz_j = self.resid[j]['coords']
				for (x,y) in itertools.product(xyz_i,xyz_j):
					d2 = dist2(x[1],y[1])
					if d2 <= rc2:
						aco += i - j
						ncont +=1
		aco = float(aco)/ncont
		self.aco = aco
		self.rco = aco/self.nres
		return self.aco,self.rco

	def calc_lro(self,cutoff):
		""" calculate long range order [Gromiha & Selvaraj, JMB 1999] """
		rc2 = cutoff**2
		lro = 0
		ncont = 0
		for i in self.reslist:
			xyz_i = filter(lambda x: 'CA' in x,self.resid[i]['coords'])
			jlist = filter(lambda x: i - x >= 12 ,range(self.nres))
			for j in jlist:
				xyz_j = filter(lambda x: 'CA' in x,self.resid[j]['coords'])
				for (x,y) in itertools.product(xyz_i,xyz_j):
					d2 = dist2(x[1],y[1])
					if d2 <= rc2:
						#print ncont,i,j,np.sqrt(d2),xyz_i,xyz_j
						ncont +=1
		lro = float(ncont)/self.nres
		self.lro = lro
		return self.lro

	def contactmap(self,out,cutoff):
		""" print contact map in gnuplot format"""
		rc2 = cutoff**2
		fout = open(out,"w")
		self.cmap = np.zeros((self.nres,self.nres))
		for i in self.reslist:
			xyz_i = self.resid[i]['coords']
			for j in self.reslist:
				xyz_j = self.resid[j]['coords']
				contact = False
				for (x,y) in itertools.product(xyz_i,xyz_j):
					d2 = dist2(x[1],y[1])
					if d2 <= rc2:
						contact = True
				if contact:
					fout.write( "%i  %i  %i\n"%(i+1,j+1,1))
				else:
					fout.write ("%i  %i  %i\n"%(i+1,j+1,0))
			fout.write ("\n")

	def permutate(self,newNt,name):
		"""" do the circular permutation """
		permutant = self.__class__()
		permutant.name = name
		permutant.resid = {}
		# define residues after new start
		iNt = int(filter(lambda x: self.resid[x]['resseq']==newNt,self.reslist)[0])
		rlist = filter(lambda x: x >= iNt, self.reslist)
		res = 0
		for i in rlist:
			permutant.resid[res] = self.resid[i]
			res += 1
		# connect with WT N-terminus
		rlist = filter(lambda x: x < iNt ,self.reslist)
		for i in rlist:
			permutant.resid[res] = self.resid[i]
			res += 1
		permutant.reslist = permutant.resid.keys()
		permutant.nres = len(permutant.reslist)
		return permutant 

def dist2(v1,v2):
	d2 = (v1[0] - v2[0])**2 + (v1[1] - v2[1])**2 + (v1[2] - v2[2])**2 
	return d2
