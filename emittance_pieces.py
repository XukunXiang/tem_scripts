import sys
import os
import argparse
import numpy as np
from scipy import constants
from phase_volume_6d import Phase6DVolume

parser = argparse.ArgumentParser(description='Calculates and reports the variance, emittance and correllation for input distributions.')
parser.add_argument('phase_volume_files', nargs='+', type=str, help='The path(s) to the simultion output files containing the phase space.')
parser.add_argument('-m','--number_of_electron_per_macroparticle', dest="number_of_electrons_per_macroparticle", type=int, help='The number of electrons per macroparticle for the simulation.  This defaults to 100 unless specified.', default=100)
parser.add_argument('--header', dest="header", action="store_true", help='Prints out the names of the columns in a separate line.  Default is off.', default=False)

args = parser.parse_args()

mass_of_electron = constants.physical_constants["electron mass energy equivalent in MeV"][0]
mass_of_macroparticle = args.number_of_electrons_per_macroparticle*mass_of_electron

phase_volume = Phase6DVolume(mass_of_electron)
for filepath in args.phase_volume_files:
  filename = os.path.basename(filepath)
  phase_volume.injectFile(filepath,mass=mass_of_electron,momentum_weight=args.number_of_electrons_per_macroparticle)

cov_matrix = phase_volume.getCovarianceMatrix()

if args.header:
  header = ["ex","varx","varpx","one_minus_cor_x","ey","vary","varpy","one_minus_cor_y","ez","varz","varpz","one_minus_cor_z"]
  print ",".join([str(h) for h in header])
  

output = []
output.append(np.sqrt(cov_matrix.getSubDeterminant(["x","px"]))/mass_of_electron)
var_x = cov_matrix.getCovarianceElement("x","x")
output.append(var_x)
var_px = cov_matrix.getCovarianceElement("px","px")
output.append(var_px)
output.append(1 - cov_matrix.getCovarianceElement("x","px")**2/(var_x*var_px))
output.append(np.sqrt(cov_matrix.getSubDeterminant(["y","py"]))/mass_of_electron)
var_y = cov_matrix.getCovarianceElement("y","y")
output.append(var_y)
var_py = cov_matrix.getCovarianceElement("py","py")
output.append(var_py)
output.append(1 - cov_matrix.getCovarianceElement("y","py")**2/(var_y*var_py))
output.append(np.sqrt(cov_matrix.getSubDeterminant(["z","pz"]))/mass_of_electron)
var_z = cov_matrix.getCovarianceElement("z","z")
output.append(var_z)
var_pz = cov_matrix.getCovarianceElement("pz","pz")
output.append(var_pz)
output.append(1 - cov_matrix.getCovarianceElement("z","pz")**2/(var_z*var_pz))
output.append(len(phase_volume)*args.number_of_electrons_per_macroparticle)
print ",".join([str(o) for o in output])
