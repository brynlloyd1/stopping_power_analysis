from EhrenfestAnalysis import EhrenfestAnalysis

import json

analysis = EhrenfestAnalysis()
directory = "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/322_hyperchannelling_hydrogen/"
name = analysis.initialise_analysis(directory)
# # analysis.set_timesteps(name, "::10")
analysis.load_data(name)
analysis.calculate_stopping_curve(name, crop=[11, None])
analysis.view_fits(name)
analysis.compare_to_geant4([name])

# data_loader = DataLoader(directory)

# data_loader.load_data()
# all_gpw_files = data_loader.get_files("gpw")
# print(json.dumps(all_gpw_files, indent=4))

