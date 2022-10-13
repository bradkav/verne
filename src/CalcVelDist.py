import numpy as np
from scipy.interpolate import interp1d

import MaxwellBoltzmann as MB
import argparse

import verne as verne
from CalcVelDist_function import calcVelDist_full

try:
    from tqdm import tqdm
except:
    def tqdm(x):
        return x
    

#Parse the arguments!
parser = argparse.ArgumentParser(description='...')
parser.add_argument('-m_x','--m_x', help='DM mass in GeV', type=float,default = 1e5)
parser.add_argument('-sigma_p','--sigma_p', help='DM-nucleon cross section, sigma_p in cm^2', type=float, default = 1e-36)
parser.add_argument('-loc','--location', help='Detector location to consider.', type=str, default="surface")
parser.add_argument('-int', '--interaction', help='Interaction type: `SI`, `SD`, `hDP` or `Millicharge`.', type=str, default="SI")
parser.add_argument('-d', '--depth', help='Underground depth of detector [m]', type=float, default=0.0)
args = parser.parse_args()
m_x = args.m_x
sigma_p = args.sigma_p
loc = args.location
interaction = args.interaction
depth = args.depth

calcVelDist_full(m_x, sigma_p, loc, interaction, depth)


#if __name__ == '__main__':
#    calcVelDist_full(0.5, 0.521e-36, "full", interaction="SI", depth_in=1400)
