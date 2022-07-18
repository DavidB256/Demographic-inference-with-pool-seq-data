import msprime
import moments

dem_params = {"ancestral_population_size": 2e5,
              "split_population_size": 1e5,
              "t_split": 10,
              "mig_rate": 1e-2,
              "mu": 2.9e-9,
              "rho": 1.25e-7,
              "seq_len": 1e6,
              "ploidy": 2}
random_seed = 1

two_pop_split_demography = msprime.Demography()
two_pop_split_demography.add_population(name="ancestral",
                                        initial_size=dem_params["ancestral_population_size"])
two_pop_split_demography.add_population(name="pop0",
                                        initial_size=dem_params["split_population_size"])
two_pop_split_demography.add_population(name="pop1",
                                        initial_size=dem_params["split_population_size"])
two_pop_split_demography.add_population_split(time=dem_params["t_split"],
                                              derived=["pop0", "pop1"], ancestral="ancestral")
two_pop_split_demography.set_symmetric_migration_rate(["pop0", "pop1"], dem_params["mig_rate"])
