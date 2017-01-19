'''
Created on 26/ago/2011

@author: asaba
'''
from ... import performancetest as pt


nwchem_task_theory = {"SCF": "Hartree-Fock",
                    "DFT": "Density functional theory for molecules",
                    "SODFT": "Spin-Orbit Density functional theory",
                    "MP2" : "MP2 using a semi-direct algorithm",
                    "DIRECT_MP2" : "MP2 using a full-direct algorithm",
                    "RIMP2" : "MP2 using the RI approximation",
                    "CCSD" : "Coupled-cluster single and double excitations",
                    "CCSD(T)" : "Coupled-cluster linearized triples approximation",
                    "CCSD+T(CCSD)#" : "Fourth order triples contribution", 
                    "MCSCF" : "Multiconfiguration SCF",
                    "SELCI" : "Selected configuration interaction with perturbation correction",
                    "MD"  : "Classical molecular dynamics simulation",
                    "PSPW" : "Pseudopotential plane-wave density functional theory for molecules and insulating solids using NWPW",
                    "BAND" : "Pseudopotential plane-wave density functional theory for solids using NWPW",
                    "TCE" : "Tensor Contraction Engine",
                    "TDDFT" : "time-dependent density functional theory",
                    "ERROR" : "ERROR: Undefined Theory" }
nwchem_task_operation = {"ENERGY" : "Evaluate the single point energy",
                         "GRADIENT" : "Evaluate the derivative of the energy with respect to nuclear coordinates",
                         "OPTIMIZE" : "Minimize the energy by varying the molecular structure. By default, this geometry optimization is presently driven by the Driver module, but the Stepper module may also be used",
                         "SADDLE" : "Conduct a search for a transition state (or saddle point) using either Driver module (the default) or Stepper",
                         "HESSIAN" : "Compute second derivatives. See hessian section for analytic hessians",
                         "FREQUENCIES" : "Compute second derivatives and print out an analysis of molecular vibrations. See vibration section for controls for vibration calculations",
                         "FREQ" : "Compute second derivatives and print out an analysis of molecular vibrations. See vibration section for controls for vibration calculations",
                         "PROPERTY" : "Calculate the properties for the wave function",
                         "DYNAMICS" : "Perform classical molecular dynamics",
                         "THERMODYNAMICS" : "Perform multi-configuration thermo-dynamic integration using classical MD",
                         "ERROR" : "ERROR: Undefined Operation" } 
nwchem_task_basis_sets = {"4-31G": "Pople basis set 4-31g",
                          "6-31G": "Pople basis set 6-31g",
                          "6-31G*" : "Pople basis set 6-31g*",
                          "6-31G**" : "Pople basis set 6-31g**",
                          "6-31+G*": "Pople basis set 6-31+g*",
                          "6-311G**" : "Pople basis set 6-311g**",
                          "6-311++G**": "Pople basis set 6-311++g**",
                          "CC-PVDZ": "Double-zeta",
                          "CC-PVTZ": "Triple-zeta",
                          "CC-PVQZ": "Quadruple-zeta",
                          "CC-PV5Z": "Quintuple-zeta",
                          "TZVP_(DFT_ORBITAL)" : "Pople basis set TZVP",
                          "ERROR" : "ERROR: Undefined basis set",
                          }
nwchem_electronicstates_description = {1:"S",
                                       2:"D",
                                       3:"T",
                                       4:"Q",
                                       5:"P",
                                       6:"E",}
fields_repr_by_list = ["data", "internucleardistances", "internuclearangles", "basis", "lines", "value_Z"]
fields_repr_by_dict = ["values"]
cutter_extremes = [["\n-+\nEAF file.*aoints", "Full wait time used for read and write\.\n\n"],
                   [" Parallel integral file used", "large values\n"],
                   ["\n-+\nEAF file.*gridpts", " Parallel grid_pts file used\s*[0-9]+ records\n"],
                   ["\n  itol2e modified to match", "convergence criterion\.\n"]] #couple of regular expressions used to remove unwanted output parts in nwc output file

nwchem_bin_path = "/usr/local/NWChem-6.0/bin/nwchem-serial-acml"

nwchem_units = {"energyunit": "eV",
                "rotationalunit": "cm-1",
                "dipolemomentunit": "Debye",
                "frequencyunit": "cm-1",
                "intensityunit": "km/mol",
                "muunit": "Debye" }