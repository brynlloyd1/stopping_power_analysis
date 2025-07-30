from EhrenfestAnalysis import EhrenfestAnalysis

import matplotlib.pyplot as plt
from logger_config import setup_logging
import logging


setup_logging()
logger = logging.getLogger(__name__)
logger.info("Starting Ehrenfest Analysis")

"""
SUPERCELL SIZE
322
622
333
422

TRAJECTORY
hyperchannelling
presampled
volumecapture

PROJECTILE SETUP
hydrogen -> just a hydrogen
proton -> hydrogen with external potential around the starting location
proton1 -> hydrogen with external potential that follows the proton
proton2 -> custom setup of hydrogen with valence state removed
proton3 -> custom setup of hydrogen with valence state removed + external potential around starting location

_ -> projectile is initialised inside the crystal
surface -> projectile is initialised in vacuum and then enters the crystal

TARGET SETUP
_ -> Al has 10 frozen core states and 3 valence states
7valence -> only did this one bc gpaw couldnt make me an 11valence one the first time i tried
11valence -> only n=1 states are frozen
"""


# parameters for the charge state analysis
parameters = {
    "size" : 5,
    "offset" : 0,
    "crop_low" : 4,
    "crop_high" : 64,
}

analysis = EhrenfestAnalysis()

directories = {
    # "322_hyperchannelling_hydrogen" : "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/322_hyperchannelling_hydrogen/",
    # "322_presampled1_hydrogen" : "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/322_presampled1_hydrogen/",
    # "333_presampled1_proton" : "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/333_presampled1_proton/",
    "333_presampled1_proton2" : "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/333_presampled1_proton2/",
    # "333_presampled1_proton2_7valence" : "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/333_presampled1_proton2_7valence/",
    # "333_presampled2_proton2_7valence_200timesteps": "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/333_presampled1_proton2_7valence_200timesteps/",
    # "333_presampled2_proton2_11valence_halftimestep": "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/333_presampled1_proton2_11valence_halftimestep/",
    "333_presampled2_proton2_11valence_quartertimestep": "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/333_presampled1_proton2_11valence_quartertimestep/",
    # "333_presampled2_proton2": "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/333_presampled2_proton2/",

    # "422_hyperchannelling_proton2_surface" : "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/422_hyperchannelling_proton2_surface/",
    # "422_hyperchannelling_proton3_surface" : "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/422_hyperchannelling_proton3_surface/",
    # "422_hyperchannelling_proton2_7valence" : "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/422_hyperchannelling_proton2_7valence/",

    # "622_hyperchannelling_hydrogen" : "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/622_hyperchannelling_hydrogen/",
    # "622_hyperchannelling_proton" : "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/622_hyperchannelling_proton/",
    "622_hyperchannelling_proton2" : "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/622_hyperchannelling_proton2/",

    # "ehh": "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/test_halftimestep/",
    # "ehhhh": "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/temp/",


    # "scarf_333_presampled3" : "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/scarf_data/333_presampled3/",
    "scarf_333_presampled4": "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/scarf_data/333_presampled4/",

}

names = {key: analysis.initialise_analysis(value) for key, value in directories.items()}
# this is where you would set_which_energies
# _ = [analysis.set_which_energies(name, which_energies = ["400 keV"]) for name in names.values()]
_ = [analysis.load_data(name) for name in names.values()]

"""
see fit to kinetic energy data and instantaneous stopping powers 
"""
# _ = [analysis.calculate_stopping_curve(name) for name in names.values()]
# _ = [analysis.view_fits(name) for name in names.values()]


"""
compare stopping data
"""
energy = "400 keV"
(analysis.compare_stopping_curves(list(names.values()), energy))

"""
compare with presampling
"""
# energy = "400 keV"
# _ = [analysis.compare_to_presampling(name, energy) for name in names.values()]

"""
get stopping curve for each and plot against SRIM
"""
# _ = [analysis.calculate_stopping_curve(name) for name in names.values()]
# analysis.compare_to_montecarlo(list(names.values()))


"""
projectile charge state
"""
# _ = [analysis.set_which_energies(name, which_energies=["400 keV"]) for name in names.values()]
# _ = [analysis.load_densities(name) for name in names.values()]
# _ = [analysis.analyse_proton_charge_state(name, parameters) for name in names.values()]



"""
visualise electron density
"""
# energy = "400 keV"
# _ = [analysis.visualise_electron_density(name, energy) for name in names.values()]

plt.show()