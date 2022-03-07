import gillespy2 as gillespy
import numpy as np
import math

import base

# test with number of species = 5, 10, 20

# BD stands for Birth Death
class BD20(base.BaseCRNModel):
    """Class for BD network from the DeepCME paper with n = 20."""

    params = {
        'k': 10.,
        'gamma': 1.
    }

    def __init__(self, endtime, timestep):
        """
        Initialize the model.

        Parameters
        ----------
        endtime : endtime of simulations
        timestep : time-step of simulations
        """
        super().__init__(
            endtime=endtime,
            timestep=timestep,
            model_name="BD20",
        )

        k = gillespy.Parameter(name='k', expression=self.params['k'])
        gamma = gillespy.Parameter(name='gamma', expression=self.params['gamma'])
        self.add_parameter([k, gamma])

        # Species, n = 5, 10, 20 in the DeepCME paper
        S1 = gillespy.Species(name='S1', initial_value=0)
        S2 = gillespy.Species(name='S2', initial_value=0)
        S3 = gillespy.Species(name='S3', initial_value=0)
        S4 = gillespy.Species(name='S4', initial_value=0)
        S5 = gillespy.Species(name='S5', initial_value=0)
        S6 = gillespy.Species(name='S6', initial_value=0)
        S7 = gillespy.Species(name='S7', initial_value=0)
        S8 = gillespy.Species(name='S8', initial_value=0)
        S9 = gillespy.Species(name='S9', initial_value=0)
        S10 = gillespy.Species(name='S10', initial_value=0)
        S11 = gillespy.Species(name='S11', initial_value=0)
        S12 = gillespy.Species(name='S12', initial_value=0)
        S13 = gillespy.Species(name='S13', initial_value=0)
        S14 = gillespy.Species(name='S14', initial_value=0)
        S15 = gillespy.Species(name='S15', initial_value=0)
        S16 = gillespy.Species(name='S16', initial_value=0)
        S17 = gillespy.Species(name='S17', initial_value=0)
        S18 = gillespy.Species(name='S18', initial_value=0)
        S19 = gillespy.Species(name='S19', initial_value=0)
        S20 = gillespy.Species(name='S20', initial_value=0)

        self.add_species([
                S1, S2, S3, S4, S5, S6, S7, S8, S9, S10,
                S11, S12, S13, S14, S15, S16, S17, S18, S19, S20,
            ])

        # Reactions
        # births
        S1_birth = gillespy.Reaction(
            name='S1_birth',
            reactants={},
            products={S1: 1},
            rate=k,
            massaction=True
        )
        S2_birth = gillespy.Reaction(
            name='S2_birth',
            reactants={},
            products={S2: 1},
            rate=k,
            massaction=True
        )
        S3_birth = gillespy.Reaction(
            name='S3_birth',
            reactants={},
            products={S3: 1},
            rate=k,
            massaction=True
        )
        S4_birth = gillespy.Reaction(
            name='S4_birth',
            reactants={},
            products={S4: 1},
            rate=k,
            massaction=True
        )
        S5_birth = gillespy.Reaction(
            name='S5_birth',
            reactants={},
            products={S5: 1},
            rate=k,
            massaction=True
        )
        S6_birth = gillespy.Reaction(
            name='S6_birth',
            reactants={},
            products={S6: 1},
            rate=k,
            massaction=True
        )
        S7_birth = gillespy.Reaction(
            name='S7_birth',
            reactants={},
            products={S7: 1},
            rate=k,
            massaction=True
        )
        S8_birth = gillespy.Reaction(
            name='S8_birth',
            reactants={},
            products={S8: 1},
            rate=k,
            massaction=True
        )
        S9_birth = gillespy.Reaction(
            name='S9_birth',
            reactants={},
            products={S9: 1},
            rate=k,
            massaction=True
        )
        S10_birth = gillespy.Reaction(
            name='S10_birth',
            reactants={},
            products={S10: 1},
            rate=k,
            massaction=True
        )
        S11_birth = gillespy.Reaction(
            name='S11_birth',
            reactants={},
            products={S11: 1},
            rate=k,
            massaction=True
        )
        S12_birth = gillespy.Reaction(
            name='S12_birth',
            reactants={},
            products={S12: 1},
            rate=k,
            massaction=True
        )
        S13_birth = gillespy.Reaction(
            name='S13_birth',
            reactants={},
            products={S13: 1},
            rate=k,
            massaction=True
        )
        S14_birth = gillespy.Reaction(
            name='S14_birth',
            reactants={},
            products={S14: 1},
            rate=k,
            massaction=True
        )
        S15_birth = gillespy.Reaction(
            name='S15_birth',
            reactants={},
            products={S15: 1},
            rate=k,
            massaction=True
        )
        S16_birth = gillespy.Reaction(
            name='S16_birth',
            reactants={},
            products={S16: 1},
            rate=k,
            massaction=True
        )
        S17_birth = gillespy.Reaction(
            name='S17_birth',
            reactants={},
            products={S17: 1},
            rate=k,
            massaction=True
        )
        S18_birth = gillespy.Reaction(
            name='S18_birth',
            reactants={},
            products={S18: 1},
            rate=k,
            massaction=True
        )
        S19_birth = gillespy.Reaction(
            name='S19_birth',
            reactants={},
            products={S19: 1},
            rate=k,
            massaction=True
        )
        S20_birth = gillespy.Reaction(
            name='S20_birth',
            reactants={},
            products={S20: 1},
            rate=k,
            massaction=True
        )

        # deaths
        S1_death = gillespy.Reaction(
            name='S1_death',
            reactants={S1: 1},
            products={},
            rate=gamma,
            massaction=True
        )
        S2_death = gillespy.Reaction(
            name='S2_death',
            reactants={S2: 1},
            products={},
            rate=gamma,
            massaction=True
        )
        S3_death = gillespy.Reaction(
            name='S3_death',
            reactants={S3: 1},
            products={},
            rate=gamma,
            massaction=True
        )
        S4_death = gillespy.Reaction(
            name='S4_death',
            reactants={S4: 1},
            products={},
            rate=gamma,
            massaction=True
        )
        S5_death = gillespy.Reaction(
            name='S5_death',
            reactants={S5: 1},
            products={},
            rate=gamma,
            massaction=True
        )
        S6_death = gillespy.Reaction(
            name='S6_death',
            reactants={S6: 1},
            products={},
            rate=gamma,
            massaction=True
        )
        S7_death = gillespy.Reaction(
            name='S7_death',
            reactants={S7: 1},
            products={},
            rate=gamma,
            massaction=True
        )
        S8_death = gillespy.Reaction(
            name='S8_death',
            reactants={S8: 1},
            products={},
            rate=gamma,
            massaction=True
        )
        S9_death = gillespy.Reaction(
            name='S9_death',
            reactants={S9: 1},
            products={},
            rate=gamma,
            massaction=True
        )
        S10_death = gillespy.Reaction(
            name='S10_death',
            reactants={S10: 1},
            products={},
            rate=gamma,
            massaction=True
        )
        S11_death = gillespy.Reaction(
            name='S11_death',
            reactants={S11: 1},
            products={},
            rate=gamma,
            massaction=True
        )
        S12_death = gillespy.Reaction(
            name='S12_death',
            reactants={S12: 1},
            products={},
            rate=gamma,
            massaction=True
        )
        S13_death = gillespy.Reaction(
            name='S13_death',
            reactants={S13: 1},
            products={},
            rate=gamma,
            massaction=True
        )
        S14_death = gillespy.Reaction(
            name='S14_death',
            reactants={S14: 1},
            products={},
            rate=gamma,
            massaction=True
        )
        S15_death = gillespy.Reaction(
            name='S15_death',
            reactants={S15: 1},
            products={},
            rate=gamma,
            massaction=True
        )
        S16_death = gillespy.Reaction(
            name='S16_death',
            reactants={S16: 1},
            products={},
            rate=gamma,
            massaction=True
        )
        S17_death = gillespy.Reaction(
            name='S17_death',
            reactants={S17: 1},
            products={},
            rate=gamma,
            massaction=True
        )
        S18_death = gillespy.Reaction(
            name='S18_death',
            reactants={S18: 1},
            products={},
            rate=gamma,
            massaction=True
        )
        S19_death = gillespy.Reaction(
            name='S19_death',
            reactants={S19: 1},
            products={},
            rate=gamma,
            massaction=True
        )
        S20_death = gillespy.Reaction(
            name='S20_death',
            reactants={S20: 1},
            products={},
            rate=gamma,
            massaction=True
        )
        
        self.add_reaction([
            S1_birth, S2_birth, S3_birth, S4_birth, S5_birth,
            S6_birth, S7_birth, S8_birth, S9_birth, S10_birth,
            S11_birth, S12_birth, S13_birth, S14_birth, S15_birth,
            S16_birth, S17_birth, S18_birth, S19_birth, S20_birth,
            S1_death, S2_death, S3_death, S4_death, S5_death,
            S6_death, S7_death, S8_death, S9_death, S10_death,
            S11_death, S12_death, S13_death, S14_death, S15_death,
            S16_death, S17_death, S18_death, S19_death, S20_death,
        ])
        nb_of_steps = int(math.ceil((endtime / timestep))) + 1
        self.timespan(np.linspace(0, endtime, nb_of_steps))

    def set_species_initial_value(self, species_initial_value):
        """
        Set initial values to species.

        Parameters
        ----------
        species_initial_value : list or 1d array of values, size should be equal
            to the number of species, the order should be coherent with theo order
            of species returned by get_initial_state method.

        Returns
        -------
        None

        """
        # YIKES
        for i, name in enumerate(self.get_species_names()):
            # self.listOfSpecies[name].initial_value = species_initial_value[i]
            self.listOfSpecies[name].initial_value = 0.0
    @staticmethod
    def get_species_names():
        """
        Returns list of species names.

        Returns
        -------
        list of all species names. The order of names should be coherent
        with the list returned by get_initial_state method.

        """
        return [
            'S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8', 'S9', 'S10',
            'S11', 'S12', 'S13', 'S14', 'S15', 'S16', 'S17', 'S18', 'S19', 'S20'
        ]

    @staticmethod
    def get_initial_state():
        pass

    @classmethod
    def get_n_species(cls):
        """Total number of species."""
        return len(cls.get_species_names())

    @classmethod
    def get_initial_settings(cls, n_settings, sigm=None):
        """
        Generate a set of (random) initial states.

        Parameters
        ----------
        n_settings : number of settings to produce.
        sigm : not used

        Returns
        -------
        array of n_settings initial states

        """
        n_species = cls.get_n_species()
        settings = np.zeros((n_settings, n_species))
        return settings

    @classmethod
    def get_histogram_bounds(cls, species_names_list=None):
        """
        Returns bounds for species histograms.

        Parameters
        ----------
        species_names_list: not used

        Returns
        -------
        histogram_bounds: list of [min, max] values for species

        """
        n_species_for_histogram = len(cls.get_species_for_histogram())
        histogram_bounds = [[0.5, 200.5]] * n_species_for_histogram
        return histogram_bounds

    @staticmethod
    def get_species_for_histogram():
        """
        Returns subset of species of interest.

        Returns
        -------
        list of species names.

        """
        return [
            'S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8', 'S9', 'S10',
            'S11', 'S12', 'S13', 'S14', 'S15', 'S16', 'S17', 'S18', 'S19', 'S20'
        ]

    @classmethod
    def get_randomized_parameters(cls, param_names, n_settings, sigm=0.5):
        randomized = {}
        for name in param_names:
            if name not in cls.params:
                raise KeyError(f"Could not find param {name} in {cls.__name__} class `params` dict.")
            val = float(cls.params[name])
    
            randomized[name] = np.random.uniform(val * (1. - sigm), val * (1. + sigm), n_settings)
        return randomized
    #
    # def set_parameters(self, params_dict):
    #     for name, val in params_dict.items():
    #         if name not in self.listOfParameters:
    #             raise KeyError(
    #                 f"Could not find {name} parameter in {self.__class__.__name__} model listOfParameters.")
    #         self.set_parameter(name, str(val))
