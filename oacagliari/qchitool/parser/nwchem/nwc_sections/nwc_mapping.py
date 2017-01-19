'''
Created on 30/ago/2011

@author: asaba
'''
import re
import nwc_level_0_sections, nwc_level_1_sections, nwc_level_2_sections
from .... import performancetest as pt

class MappingObject:
    def __init__(self, name, re_ex, relatedclass):
        self.name = name #description name of the section
        self.re_ex = re_ex #regular expression for the first line of the section
        self.relatedclass = relatedclass #Name of the class that manage this section
        self.re_ex_compiled = re.compile(self.re_ex) #compiled regular expression for faster parsing

class MappingSections:
        def __init__(self, filetype):
            self.levels = []
            if filetype == "NwChemOutput":
                self.levels.append([MappingObject("arguments", r"argument  1 =.*", nwc_level_0_sections.Arguments),
                               MappingObject("echo", r"[=]+ echo of input deck [=]+",  nwc_level_0_sections.InputDeck),
                               MappingObject("NWChem package version", r"Northwest Computational Chemistry Package \(NWChem\) [0-9]\.[0-9]",  nwc_level_0_sections.NwchemManifest),
                               MappingObject("NWChem Input Module", r"NWChem Input Module",  nwc_level_0_sections.NwchemInputModule), #more than one
                               MappingObject("NWChem Geometry Optimization", r"NWChem Geometry Optimization",  nwc_level_0_sections.NwchemGeometryOptimization),
                               MappingObject("NWChem Nuclear Hessian and Frequency Analysis", r"NWChem Nuclear Hessian and Frequency Analysis", nwc_level_2_sections.NwchemNuclearHessianAndFrequencyAnalysis),
                                                                                   
                               MappingObject("NWChem DFT Module", r"NWChem DFT Module",  nwc_level_1_sections.NwchemDftModule), #dft - Density functional theory for molecules
                               MappingObject("NWChem Property Module", r"NWChem Property Module", nwc_level_1_sections.NwchemPropertyModule),
                               ]) #(section name, regular expression to match the row, [start row, end row], related object name)
                self.levels.append([
                               MappingObject("hessian step", r"atom:.+ xyz: .+\(.+\) wall time:.+date:.+", nwc_level_0_sections.HessianStep), #more than one
                               MappingObject("Step", r"[\s]*Step[\s]*[0-9]+", nwc_level_0_sections.Step), #n (step index)     
                               MappingObject("NWChem SCF Module", r"NWChem SCF Module", nwc_level_1_sections.NwchemScfModule), #scf - Hartree-Fock
                               MappingObject("NWChem DFT Module", r"NWChem DFT Module",  nwc_level_1_sections.NwchemDftModule), #dft - Density functional theory for molecules
                               #("sodft - Spin-Orbit Density functional theory", r"", [0,0], None),
                               MappingObject("NWChem MP2 Semi-direct Energy/Gradient Module", r"NWChem MP2 Semi-direct Energy\/Gradient Module", None), #mp2 - MP2 using a semi-direct algorithm
                               #("direct_mp2 - MP2 using a full-direct algorithm", r"", None),
                               #("rimp2 - MP2 using the RI approximation", r"",  None),
                               #("ccsd - Coupled-cluster single and double excitations", r"",  None),
                               #("ccsd(t) - Coupled-cluster linearized triples approximation", r"",  None),
                               #("mcscf - Multiconfiguration SCF", r"", None),
                               #("selci - Selected configuration interaction with perturbation correction", r"",  None),
                               #("md - Classical molecular dynamics simulation", r"",  None),
                               #("pspw - Pseudopotential plane-wave density functional theory for molecules and insulating solids using NWPW", r"", None),
                               #("band - Pseudopotential plane-wave density functional theory for solids using NWPW", r"", None),
                               #("tce - Tensor Contraction Engine", r"", None),
                               MappingObject("NWChem CPHF Module", r"NWChem CPHF Module",   nwc_level_1_sections.NwchemCphfModule),
                               MappingObject("NWChem TDDFT Module", r"NWChem TDDFT Module", None),
                               MappingObject("NWChem DFT Gradient Module", r"NWChem DFT Gradient Module", nwc_level_1_sections.DftGradientModule),
                               MappingObject("Vibrational analysis via the FX method", r".*Vibrational analysis via the FX method.*", nwc_level_1_sections.VibrationalAnalysisViaTheFXMethod),
                                                                                         
                               ##("NWChem Geometry Optimization", r"NWChem Geometry Optimization",  nwc_level_0_sections.NwchemGeometryOptimization)#,
                               #("hessian pass", r"atom:.+ xyz: .+\(.+\) wall time:.+date:.+", None) #more than one
                               ]) #(section name, regular expression to match the row, [start row, end row], related object name)
                self.levels.append([MappingObject("Memory information", r"Memory information",  nwc_level_2_sections.MemoryInformation),
                               MappingObject("Directory information", r"Directory information",  nwc_level_2_sections.DirectoryInformation),
                               MappingObject("Previous task information", r"Previous task information", nwc_level_2_sections.PreviousTaskInformation),
                               MappingObject("Geometries in the database", r"Geometries in the database",  nwc_level_2_sections.GeometriesInTheDatabase),
                               MappingObject("Basis sets in the database", r"Basis sets in the database",  nwc_level_2_sections.BasisSetsInTheDatabase),    
                               MappingObject("DFT ENERGY GRADIENTS", r"DFT ENERGY GRADIENTS", nwc_level_2_sections.DftEnergyGradients),                                                             
                               MappingObject("Finite difference hessian delta", r"finite difference hessian delta.*", nwc_level_2_sections.FiniteDifferenceHessianDelta),
                               MappingObject("Finite difference derivative dipole", r"finite difference derivative dipole.*", nwc_level_2_sections.FiniteDifferenceDerivativeDipole),
                               MappingObject("Summary of", r"Summary of \".*\" -> \".*\" \(.+\)", nwc_level_2_sections.SummaryOf), # "" -> ""
                               MappingObject("NWChem Finite-difference Hessian", r"NWChem Finite\-difference Hessian", nwc_level_2_sections.NwchemFiniteDifferenceHessian),
                               MappingObject("Geometry", r"Geometry \".*\" -> \".*\"", nwc_level_2_sections.Geometry), # "" -> ""
                               #("Energy Minimization", r"Energy Minimization",  sub_sections.energy_minimization_class),
                               MappingObject("Atomic Mass", r"Atomic Mass",  nwc_level_2_sections.AtomicMass),
                               MappingObject("Nuclear Dipole moment (a.u.)", r"Nuclear Dipole moment \(.*\)", nwc_level_2_sections.NuclearDipoleMoment),
                               MappingObject("Symmetry information", r"Symmetry information",  nwc_level_2_sections.SymmetryInformation),
                               MappingObject("Z-matrix (autoz)", r"Z-matrix \(autoz\)",   nwc_level_2_sections.ZMatrix),
                               MappingObject("Basis", r"Basis \".*\" -> \".*\" \(.+\)", nwc_level_2_sections.BasisFrom), # "" -> ""
                               MappingObject("General Information", r"General Information",  nwc_level_2_sections.GeneralInformation),
                               MappingObject("XC Information", r"XC Information",  nwc_level_2_sections.XcInformation),
                               MappingObject("Grid Information", r"Grid Information",  nwc_level_2_sections.GridInformation),
                               MappingObject("Convergence Information", r"Convergence Information",  nwc_level_2_sections.ConvergenceInformation),
                               MappingObject("Screening Tollerance Information", r"Screening Tolerance Information",  nwc_level_2_sections.ScreeningToleranceInformation),
                               MappingObject("Superposition of Atomic Density Guess", r"Superposition of Atomic Density Guess",  nwc_level_2_sections.SuperpositionOfAtomicDensityGuess),
                               MappingObject("Non-variational initial energy", r"Non-variational initial energy",  nwc_level_2_sections.NonVariationalInitialEnergy),
                               MappingObject("Swapping alpha orbitals", r"Swapping alpha orbitals.*", nwc_level_2_sections.SwappingOrbitals),
                               MappingObject("Swapping beta orbitals", r"Swapping beta orbitals.*", nwc_level_2_sections.SwappingOrbitals),
                               MappingObject("Symmetry analysis of molecular orbitals - initial alpha", r"Symmetry analysis of molecular orbitals - initial alpha",  nwc_level_2_sections.SymmetryAnalysisOfMolecularOrbitals),
                               MappingObject("Symmetry analysis of molecular orbitals - initial beta", r"Symmetry analysis of molecular orbitals - initial beta",  nwc_level_2_sections.SymmetryAnalysisOfMolecularOrbitals),
                               MappingObject("Symmetry analysis of molecular orbitals - initial", r"Symmetry analysis of molecular orbitals - initial",  nwc_level_2_sections.SymmetryAnalysisOfMolecularOrbitals),
                               MappingObject("DTF Final Alpha Molecular Orbital Analysis", r"DFT Final Alpha Molecular Orbital Analysis",  nwc_level_2_sections.FinalMolecularOrbitalAnalysis),
                               MappingObject("DTF Final Beta Molecular Orbital Analysis", r"DFT Final Beta Molecular Orbital Analysis",  nwc_level_2_sections.FinalMolecularOrbitalAnalysis),
                               MappingObject("DFT Final Molecular Orbital Analysis", r"DFT Final Molecular Orbital Analysis",  nwc_level_2_sections.FinalMolecularOrbitalAnalysis),
                               MappingObject("alpha - beta orbital overlaps", r"alpha - beta orbital overlaps",  nwc_level_2_sections.OrbitalOverlaps),
                               MappingObject("Expectation value of S2:", r"Expectation value of S2:", nwc_level_2_sections.ExpectationValueOfS2),
                               MappingObject("center of mass", r"center of mass", nwc_level_2_sections.CenterOfMass),
                               MappingObject("moment of inertia", r"moments of inertia \(.*\)", nwc_level_2_sections.MomentOfInertia),
                               MappingObject("Multipole analysis of the density", r"Multipole analysis of the density wrt the origin", nwc_level_2_sections.MultipoleAnalysisOfTheDensityWRTOrigini),
                               MappingObject("Multipole analysis of the density", r"Multipole analysis of the density", nwc_level_2_sections.MultipoleAnalysisOfTheDensity),
                               MappingObject("Summary of allocated global arrays", r"Summary of allocated global arrays", nwc_level_2_sections.SummaryOfAllocatedGlobalArray),
                               MappingObject("GA Statistics for process", r"GA Statistics for process[\s]*[0-9]+", nwc_level_2_sections.GaStatisticsForProcess),
                               MappingObject("Total electron density", r"Total electron density", nwc_level_2_sections.TotalElectronDensity),
                               MappingObject("Total spin density", r"Total spin density", nwc_level_2_sections.TotalSpinDensity),
                               MappingObject("Job Information", r"Job information", nwc_level_2_sections.JobInformation),
                               MappingObject("Effective nuclear repulsion energy", r"Effective nuclear repulsion energy .*", nwc_level_2_sections.EffectiveNuclearRepulsionEnergy),
                               MappingObject("Final UHF  results", r"Final UHF  results", nwc_level_2_sections.FinalUhfResult),
                               MappingObject("Energy Minimization", r"Energy Minimization", nwc_level_2_sections.EnergyMinimization),
                               MappingObject("Mass-weighted nuclear Hessian", r"MASS[-]WEIGHTED NUCLEAR HESSIAN.*", nwc_level_2_sections.MassWeightedNuclearHessian),
                               MappingObject("Step Info", r"[\s]*Step[\s]*Energy[\s]*Delta E.*", nwc_level_2_sections.StepInfo),
                               MappingObject("Optimization converged", r".*Optimization converged.*", nwc_level_2_sections.OptimizationConverged),
                               MappingObject("Atom information", r"[-]+ Atom information [-]+",  nwc_level_2_sections.AtomInformation),
                               MappingObject("Normal Mode Eigenvectors", r".*NORMAL MODE EIGENVECTORS.*",  nwc_level_2_sections.NormalModeEigenvectors),
                               MappingObject("Normal Eigenvalue Projected Dipole Moment", r".*Normal Eigenvalue.*Projected Derivative Dipole Moments.*", nwc_level_2_sections.NormalEigenvalue_DipoleMoment),
                               MappingObject("Normal Eigenvalue Projected Infrared Intensities", r".*Normal Eigenvalue.*Projected Infra Red Intensities.*", nwc_level_2_sections.NormalEigenvalue_InfraRedIntensities),
                               MappingObject("Rotational Constants", r".*Rotational Constants.*", nwc_level_2_sections.Rot_Const_And_Energy),                    
                               ])
            elif filetype == "NwChemDatabase":
                self.levels.append([MappingObject("arguments", r"argument  1 =.*", nwc_level_0_sections.Arguments),
                               MappingObject("echo", r"[=]+ echo of input deck [=]+",  nwc_level_0_sections.InputDeck),
                               MappingObject("NWChem package version", r"Northwest Computational Chemistry Package \(NWChem\) [0-9]\.[0-9]",  nwc_level_0_sections.NwchemManifest),
                               MappingObject("NWChem Input Module", r"NWChem Input Module",  nwc_level_0_sections.NwchemInputModule), #more than one
                               MappingObject("NWChem DFT Module", r"NWChem DFT Module",  nwc_level_1_sections.NwchemDftModule), #dft - Density functional theory for molecules
                               MappingObject("NWChem Geometry Optimization", r"NWChem Geometry Optimization",  nwc_level_0_sections.NwchemGeometryOptimization),
                               MappingObject("RTDB Database", r".*Contents of RTDB .*", nwc_level_0_sections.RTDBDatabase)
                               ]) #(section name, regular expression to match the row, [start row, end row], related object name)
                self.levels.append([MappingObject("hessian step", r"atom:.+ xyz: .+\(.+\) wall time:.+date:.+", nwc_level_0_sections.HessianStep), #more than one
                                    MappingObject("Step", r"[\s]*Step[\s]*[0-9]+", nwc_level_0_sections.Step), #n (step index)
                                    MappingObject("NWChem SCF Module", r"NWChem SCF Module", None), #scf - Hartree-Fock
                               MappingObject("NWChem DFT Module", r"NWChem DFT Module",  nwc_level_1_sections.NwchemDftModule), #dft - Density functional theory for molecules
                               #("sodft - Spin-Orbit Density functional theory", r"", [0,0], None),
                               MappingObject("NWChem MP2 Semi-direct Energy/Gradient Module", r"NWChem MP2 Semi-direct Energy\/Gradient Module", None), #mp2 - MP2 using a semi-direct algorithm
                               #("direct_mp2 - MP2 using a full-direct algorithm", r"", None),
                               #("rimp2 - MP2 using the RI approximation", r"",  None),
                               #("ccsd - Coupled-cluster single and double excitations", r"",  None),
                               #("ccsd(t) - Coupled-cluster linearized triples approximation", r"",  None),
                               #("mcscf - Multiconfiguration SCF", r"", None),
                               #("selci - Selected configuration interaction with perturbation correction", r"",  None),
                               #("md - Classical molecular dynamics simulation", r"",  None),
                               #("pspw - Pseudopotential plane-wave density functional theory for molecules and insulating solids using NWPW", r"", None),
                               #("band - Pseudopotential plane-wave density functional theory for solids using NWPW", r"", None),
                               #("tce - Tensor Contraction Engine", r"", None),
                               MappingObject("NWChem CPHF Module", r"NWChem CPHF Module",   nwc_level_1_sections.NwchemCphfModule),
                               MappingObject("NWChem TDDFT Module", r"NWChem TDDFT Module", None),
                               MappingObject("NWChem Property Module", r"NWChem Property Module", nwc_level_1_sections.NwchemPropertyModule),
                               MappingObject("NWChem DFT Gradient Module", r"NWChem DFT Gradient Module", nwc_level_1_sections.DftGradientModule), 
                               #("DFT ENERGY GRADIENTS", r"DFT ENERGY GRADIENTS", super_sections.DftEnergyGradient),                                                             
                               MappingObject("NWChem Geometry Optimization", r"NWChem Geometry Optimization",  nwc_level_0_sections.NwchemGeometryOptimization)#,
                               #("hessian pass", r"atom:.+ xyz: .+\(.+\) wall time:.+date:.+", None) #more than one
                               ]) #(section name, regular expression to match the row, [start row, end row], related object name)
                self.levels.append([MappingObject("Memory information", r"Memory information",  nwc_level_2_sections.MemoryInformation),
                               MappingObject("Directory information", r"Directory information",  nwc_level_2_sections.DirectoryInformation),
                               MappingObject("Previous task information", r"Previous task information", nwc_level_2_sections.PreviousTaskInformation),
                               MappingObject("Geometries in the database", r"Geometries in the database",  nwc_level_2_sections.GeometriesInTheDatabase),
                               MappingObject("Basis sets in the database", r"Basis sets in the database",  nwc_level_2_sections.BasisSetsInTheDatabase),    
                               MappingObject("DFT ENERGY GRADIENTS", r"DFT ENERGY GRADIENTS", nwc_level_2_sections.DftEnergyGradients),                                                             
                               MappingObject("Summary of", r"Summary of \".*\" -> \".*\" \(.+\)", nwc_level_2_sections.SummaryOf), # "" -> ""
                               MappingObject("NWChem Nuclear Hessian and Frequency Analysis", r"NWChem Nuclear Hessian and Frequency Analysis", nwc_level_2_sections.NwchemNuclearHessianAndFrequencyAnalysis),
                               MappingObject("NWChem Finite-difference Hessian", r"NWChem Finite\-difference Hessian", nwc_level_2_sections.NwchemFiniteDifferenceHessian),
                               MappingObject("Geometry", r"Geometry \".*\" -> \".*\"", nwc_level_2_sections.Geometry), # "" -> ""
                               #("Energy Minimization", r"Energy Minimization",  sub_sections.energy_minimization_class),
                               MappingObject("Atomic Mass", r"Atomic Mass",  nwc_level_2_sections.AtomicMass),
                               MappingObject("Nuclear Dipole moment (a.u.)", r"Nuclear Dipole moment \(.*\)", nwc_level_2_sections.NuclearDipoleMoment),
                               MappingObject("Symmetry information", r"Symmetry information",  nwc_level_2_sections.SymmetryInformation),
                               MappingObject("Z-matrix (autoz)", r"Z-matrix \(autoz\)",   nwc_level_2_sections.ZMatrix),
                               MappingObject("Basis", r"Basis \".*\" -> \".*\" \(.+\)", nwc_level_2_sections.BasisFrom), # "" -> ""
                               MappingObject("General Information", r"General Information",  nwc_level_2_sections.GeneralInformation),
                               MappingObject("XC Information", r"XC Information",  nwc_level_2_sections.XcInformation),
                               MappingObject("Grid Information", r"Grid Information",  nwc_level_2_sections.GridInformation),
                               MappingObject("Convergence Information", r"Convergence Information",  nwc_level_2_sections.ConvergenceInformation),
                               MappingObject("Screening Tollerance Information", r"Screening Tolerance Information",  nwc_level_2_sections.ScreeningToleranceInformation),
                               MappingObject("Superposition of Atomic Density Guess", r"Superposition of Atomic Density Guess",  nwc_level_2_sections.SuperpositionOfAtomicDensityGuess),
                               MappingObject("Non-variational initial energy", r"Non-variational initial energy",  nwc_level_2_sections.NonVariationalInitialEnergy),
                               MappingObject("Symmetry analysis of molecular orbitals - initial alpha", r"Symmetry analysis of molecular orbitals - initial alpha",  nwc_level_2_sections.SymmetryAnalysisOfMolecularOrbitals),
                               MappingObject("Symmetry analysis of molecular orbitals - initial beta", r"Symmetry analysis of molecular orbitals - initial beta",  nwc_level_2_sections.SymmetryAnalysisOfMolecularOrbitals),
                               MappingObject("Symmetry analysis of molecular orbitals - initial", r"Symmetry analysis of molecular orbitals - initial",  nwc_level_2_sections.SymmetryAnalysisOfMolecularOrbitals),
                               MappingObject("DTF Final Alpha Molecular Orbital Analysis", r"DFT Final Alpha Molecular Orbital Analysis",  nwc_level_2_sections.FinalMolecularOrbitalAnalysis),
                               MappingObject("DTF Final Beta Molecular Orbital Analysis", r"DFT Final Beta Molecular Orbital Analysis",  nwc_level_2_sections.FinalMolecularOrbitalAnalysis),
                               MappingObject("DFT Final Molecular Orbital Analysis", r"DFT Final Molecular Orbital Analysis",  nwc_level_2_sections.FinalMolecularOrbitalAnalysis),
                               MappingObject("alpha - beta orbital overlaps", r"alpha - beta orbital overlaps",  nwc_level_2_sections.OrbitalOverlaps),
                               MappingObject("Expectation value of S2:", r"Expectation value of S2:", nwc_level_2_sections.ExpectationValueOfS2),
                               MappingObject("center of mass", r"center of mass", nwc_level_2_sections.CenterOfMass),
                               MappingObject("moment of inertia", r"moments of inertia \(.*\)", nwc_level_2_sections.MomentOfInertia),
                               MappingObject("Multipole analysis of the density", r"Multipole analysis of the density", nwc_level_2_sections.MultipoleAnalysisOfTheDensity),
                               MappingObject("Summary of allocated global arrays", r"Summary of allocated global arrays", nwc_level_2_sections.SummaryOfAllocatedGlobalArray),
                               MappingObject("GA Statistics for process", r"GA Statistics for process[\s]*[0-9]+", nwc_level_2_sections.GaStatisticsForProcess),
                               MappingObject("Total electron density", r"Total electron density", nwc_level_2_sections.TotalElectronDensity),
                               MappingObject("Total spin density", r"Total spin density", nwc_level_2_sections.TotalSpinDensity),
                               MappingObject("Job Information", r"Job information", nwc_level_2_sections.JobInformation),
                               MappingObject("Effective nuclear repulsion energy", r"Effective nuclear repulsion energy .*", nwc_level_2_sections.EffectiveNuclearRepulsionEnergy),
                               MappingObject("Final UHF  results", r"Final UHF  results", nwc_level_2_sections.FinalUhfResult),
                               MappingObject("Energy Minimization", r"Energy Minimization", nwc_level_2_sections.EnergyMinimization),
                               MappingObject("Mass-weighted nuclear Hessian", r"MASS[-]WEIGHTED NUCLEAR HESSIAN.*", nwc_level_2_sections.MassWeightedNuclearHessian),
                               MappingObject("Step Info", r"[\s]*Step[\s]*Energy[\s]*Delta E.*", nwc_level_2_sections.StepInfo),
                               MappingObject("Optimization converged", r".*Optimization converged.*", nwc_level_2_sections.OptimizationConverged),
                               MappingObject("Database Entry", r".*[a-z]+[\[][0-9]+[\]][\s]*[a-zA-Z]{3}[\s]*[a-zA-Z]{3}[\s]*[0-9]{1,2}[\s]*[0-9]{2}:[0-9]{2}:[0-9]{2}[\s]*[0-9]{4}.*", nwc_level_2_sections.DatabaseEntry),
                               MappingObject("Atom information", r"[-]+ Atom information [-]+",  nwc_level_2_sections.AtomInformation),
                               MappingObject("Normal Mode Eigenvectors", r".*NORMAL MODE EIGENVECTORS.*",  nwc_level_2_sections.NormalModeEigenvectors),
                               ])

#class Sections:
#    def __init__(self):
#        self.types = {"NwChemOutput": MappingSections("NwChemOutput"),
#                      "NwChemDatabase": MappingSections("NwChemDatabase")}
#section = Sections()