from EhrenfestAnalysis import EhrenfestAnalysis

analysis = EhrenfestAnalysis()

directory_622_hyperchannelling_hydrogen = "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/622_hyperchannelling_hydrogen/"
directory_622_hyperchannelling_proton = "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/622_hyperchannelling_proton/"

name_622_hyperchannelling_hydrogen = analysis.initialise_analysis(directory_622_hyperchannelling_hydrogen)
name_622_hyperchannelling_proton = analysis.initialise_analysis(directory_622_hyperchannelling_proton)

analysis.load_data(name_622_hyperchannelling_hydrogen)
# analysis.load_data(name_622_hyperchannelling_proton)

analysis.calculate_stopping_curve(name_622_hyperchannelling_hydrogen)
# analysis.calculate_stopping_curve(name_622_hyperchannelling_proton)

analysis.view_fits(name_622_hyperchannelling_hydrogen)
# analysis.view_fits(name_622_hyperchannelling_proton)

# analysis.compare_to_montecarlo([name_622_hyperchannelling_hydrogen, name_622_hyperchannelling_proton])


