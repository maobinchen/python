#!/usr/bin/env python

import sys
import re
from os import system
ROSETTA_PATH = '/work/binchen/rosetta/'
idealization_app = ROSETTA_PATH + "rosetta_source/bin/idealize_jd2.linuxgccrelease"
rosetta_database = ROSETTA_PATH + 'rosetta_database'
pdb = sys.argv[1]
m = re.match(r"(.*)\.pdb",pdb)
pdb = m.group(1)
cmd = idealization_app + " -in:file:s " + pdb + ".pdb -out:file:output_torsions -chemical:exclude_patches LowerDNA  UpperDNA Cterm_amidation SpecialRotamer VirtualBB ShoveBB VirtualDNAPhosphate VirtualNTerm CTermConnect sc_orbitals pro_hydroxylated_case1 pro_hydroxylated_case2 ser_phosphorylated thr_phosphorylated  tyr_phosphorylated tyr_sulfated lys_dimethylated lys_monomethylated  lys_trimethylated lys_acetylated glu_carboxylated cys_acetylated tyr_diiodinated N_acetylated C_methylamidated MethylatedProteinCterm -chemical:exclude_patches ActeylatedProteinNterm -chainbreaks -overwrite -database %s" % rosetta_database
system(cmd)
relax_app = ROSETTA_PATH + "rosetta_source/bin/relax.linuxgccrelease"
cmd = relax_app + " -in:file:s " + pdb + "_0001.pdb  -relax:fast -in:file:fullatom -ignore_unrecognized_res -out:file:output_torsions -chemical:exclude_patches LowerDNA  UpperDNA Cterm_amidation SpecialRotamer VirtualBB ShoveBB VirtualDNAPhosphate VirtualNTerm CTermConnect sc_orbitals pro_hydroxylated_case1 pro_hydroxylated_case2 ser_phosphorylated thr_phosphorylated  tyr_phosphorylated tyr_sulfated lys_dimethylated lys_monomethylated  lys_trimethylated lys_acetylated glu_carboxylated cys_acetylated tyr_diiodinated N_acetylated C_methylamidated MethylatedProteinCterm -chemical:exclude_patches ActeylatedProteinNterm -overwrite -database %s" % rosetta_database
system(cmd)
