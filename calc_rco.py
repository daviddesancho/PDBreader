#!/usr/bin/env python

import sys,os
import numpy as np
sys.path.append("/home/dd363/Projects/Proteins/Death/code")
#sys.path.append("/Users/daviddesancho/Research/Projects/Proteins/Death/code")
import PDBreader

# read PDB file
print " Permutants calculated from experimental structure"
name = "1E41_model5.pdb"
wt = PDBreader.Protein(name)
# contact order of WT experimental
print " Structure  ACO     RCO    LRO "
print " WT (Exp)   %6.4f  %6.4f" %wt.calc_rco(6),"%6.4f"%wt.calc_lro(8)
out = name[:name.rfind(".pdb")] + "_map.dat"
wt.contactmap(out,6)

# define series of permutants
cp1 = wt.permutate(110,"cp1_calc")
print " CP1 (Calc) %6.4f  %6.4f"%cp1.calc_rco(6)," %6.4f"%cp1.calc_lro(8)
out = cp1.name + "_map.dat"
cp1.contactmap(out,6)

cp2 = wt.permutate(154,"cp2_calc")
print " CP2 (Calc) %6.4f  %6.4f"%cp2.calc_rco(6)," %6.4f"%cp2.calc_lro(8)
out = cp2.name + "_map.dat"
cp2.contactmap(out,6)

cp3 = wt.permutate(121,"cp3_calc")
print " CP3 (Calc) %6.4f  %6.4f"%cp3.calc_rco(6)," %6.4f"%cp3.calc_lro(8)
out = cp3.name + "_map.dat"
cp3.contactmap(out,6)

cp4 = wt.permutate(135,"cp4_calc")
print " CP4 (Calc) %6.4f  %6.4f"%cp4.calc_rco(6)," %6.4f"%cp4.calc_lro(8)
out = cp4.name + "_map.dat"
cp4.contactmap(out,6)

cp5 = wt.permutate(171,"cp5_calc")
print " CP5 (Calc) %6.4f  %6.4f"%cp5.calc_rco(6)," %6.4f"%cp5.calc_lro(8)
out = cp5.name + "_map.dat"
cp5.contactmap(out,6)

# read filenames for simulated permutants
print "\n Permutants from MD simulations"
filenames = map(lambda x: x.split()[0],
		open("names.txt").readlines())
print " Structure  ACO     RCO    LRO "
for name in filenames:
	p = PDBreader.Protein(name)
	p.calc_lro(8)

# contact order of WT simulated 
	print " %s "%name, "  %6.4f  %6.4f" %p.calc_rco(6),"%6.4f"%p.calc_lro(8)
	out = name[:name.rfind(".pdb")] + "_map.dat"
	p.contactmap(out,6)

