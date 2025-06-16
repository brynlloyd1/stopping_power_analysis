from EhrenfestAnalysis import EhrenfestAnalysis

analysis = EhrenfestAnalysis()
directory = "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/322_hyperchannelling_hydrogen/"
name = analysis.initialise_analysis(directory)
analysis.set_timesteps(name, "::10")
analysis.load_data(name)

