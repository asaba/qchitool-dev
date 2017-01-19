'''
Created on 01/mar/2012

@author: asaba
'''
from .. import performancetest as pt

import geometry_tools as gt
debug = True
#gt.test("cml_out_file.cml", "cml_out_file3.cml", debug)
#gt.test("cml_out_file.cml", "cml_out_file3_wrong.cml", debug)
#gt.test("circumcoro_c.cml", "circumcoro_c_rotated1.cml", debug)
#gt.test("circumcoro_c.cml", "circumcoro_c_rotated1_wrong.cml", debug)
#gt.test("circumcoro_n.cml", "circumcoro_n_rotated2.cml", debug)
#gt.test("circumcoro_n.cml", "circumcoro_n_rotated2_wrong.cml", debug)
#gt.test("fullerene_n.cml", "fullerene_n_rotated1.cml", debug)
#gt.test("fullerene_n.cml", "fullerene_n_rotated1_wrong.cml", debug)
gt.test("antracene_631gs_2c.cml", "antracene_631gs_2c_nosym.cml", debug)
#gt.test("antracene_631gs_2c.cml", "antracene_631gs_2c_nvert.cml", debug)
#gt.test("anthracene_631gs_a_nosym.cml", "anthracene_631gs_a.cml", debug)
#gt.test("circumanthracene_631gs_c.cml", "circumanthracene_631gs_a.cml", debug)
#gt.test("bisanthene_631gs_a.cml", "bisanthene_631gs_n.cml", debug)
#gt.test("tetracene_631gs_a_nosym.cml", "tetracene_631gs_a.cml", debug)
#gt.rotate("circumovalene_631gs_a.cml", "circumovalene_631gs_a_rotated.xyz", None, traspose = False)
#gt.rotate("circumovalene_631gs_a_rotated.cml", "circumovalene_631gs_a_rotated_2.xyz", None, traspose = True)
#gt.rotate("circumovalene_631gs_a_rotated.cml", "circumovalene_631gs_a_antirotated_1.xyz", None, traspose = True)
gt.test("circumovalene_631gs_a_rotated.cml", "circumovalene_631gs_a_antirotated.cml", debug)
gt.test("tetracene_631gs_a_nosym.cml", "tetracene_631gs_a.cml", debug)
