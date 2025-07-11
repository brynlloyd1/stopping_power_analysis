from EhrenfestAnalysis import EhrenfestAnalysis

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
ae -> all-electron (Al has only 4 core states + 7 valence) 
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
    # "333_presampled1_proton2_ae" : "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/333_presampled1_proton2_ae/",
    # "333_presampled2_proton2": "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/333_presampled2_proton2/",

    # "422_hyperchannelling_proton2_surface" : "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/422_hyperchannelling_proton2_surface/",
    # "422_hyperchannelling_proton3_surface" : "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/422_hyperchannelling_proton3_surface/",
    # "422_hyperchannelling_proton2_ae" : "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/422_hyperchannelling_proton2_ae/",

    # "622_hyperchannelling_hydrogen" : "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/622_hyperchannelling_hydrogen/",
    # "622_hyperchannelling_proton" : "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/622_hyperchannelling_proton/",
    # "622_hyperchannelling_proton2" : "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/622_hyperchannelling_proton2/",
}

names = {key: analysis.initialise_analysis(value) for key, value in directories.items()}
# this is where you would set_which_energies
analysis.set_which_energies(*list(names.values()), which_energies = ["400 keV"])
_ = [analysis.load_data(name) for name in names.values()]


"""
see fit to kinetic energy data and instantaneous stopping powers 
"""
_ = [analysis.calculate_stopping_curve(name) for name in names.values()]
_ = [analysis.view_fits(name) for name in names.values()]


"""
compare stopping data
"""
# energy = "400 keV"
# _ = [analysis.calculate_stopping_curve(name) for name in names.values()]
# analysis.compare_fits(list(names.values()), energy)




"""
get stopping curve for each and plot against SRIM
"""
# _ = [analysis.calculate_stopping_curve(name) for name in names.values()]
# analysis.compare_to_montecarlo(list(names.values()))


"""
visualise electron density
"""
# energy = "40 keV"
# _ = [analysis.visualise_electron_density(name, energy) for name in names.values()]

