from EhrenfestAnalysis import EhrenfestAnalysis

import matplotlib.pyplot as plt
from logger_config import setup_logging
import logging


setup_logging()
logger = logging.getLogger(__name__)
logger.info("Starting Ehrenfest Analysis")

analysis = EhrenfestAnalysis()

directories = {
    # "322_hyperchannelling_hydrogen" : "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/322_hyperchannelling_hydrogen/",
    # "322_presampled1_hydrogen" : "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/322_presampled1_hydrogen/",

    # "333_presampled1_proton" : "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/333_presampled1_proton/",
    # "333_presampled1_proton2" : "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/333_presampled1_proton2/",
    # "333_presampled1_proton2_7valence" : "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/333_presampled1_proton2_7valence/",
    # "333_presampled2_proton2_7valence_200timesteps": "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/333_presampled1_proton2_7valence_200timesteps/",
    "333_presampled2_proton2_11valence_halftimestep": "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/333_presampled1_proton2_11valence_halftimestep/",
    # "333_presampled2_proton2": "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/333_presampled2_proton2/",

    # "422_hyperchannelling_proton2_surface" : "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/422_hyperchannelling_proton2_surface/",
    # "422_hyperchannelling_proton3_surface" : "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/422_hyperchannelling_proton3_surface/",
    # "422_hyperchannelling_proton2_7valence" : "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/422_hyperchannelling_proton2_7valence/",

    # "622_hyperchannelling_hydrogen" : "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/622_hyperchannelling_hydrogen/",
    # "622_hyperchannelling_proton" : "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/622_hyperchannelling_proton/",
    # "622_hyperchannelling_proton2" : "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/622_hyperchannelling_proton2/",

    # "ehh": "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/test_halftimestep/",
    # "ehhhh": "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/temp/",
}

names = {key: analysis.initialise_analysis(value) for key, value in directories.items()}
# this is where you would set_which_energies
energy = "400 keV"
_ = [analysis.set_which_energies(name, which_energies = [energy]) for name in names.values()]
_ = [analysis.load_data(name) for name in names.values()]

name = list(names.values())[0]
processor = analysis.data_handlers[name].data_processor
data = analysis.data_handlers[name].all_data


rolling_stopping_powers = []
for i in range(len(data[energy].projectile_positions), 0, -1):

    crop = [0, i]

    try:
        fit = processor.calculate_stopping_powers(data, crop=crop)
    except:
        continue

    rolling_stopping_powers.append(-1e3 * fit[energy].fit[0])

plt.plot(rolling_stopping_powers)
plt.show()