import gillespy2 as gillespy
import numpy as np
import math

import base

# test with number of species = 5, 10, 20

# BD stands for Birth Death
class BD5(base.BaseCRNModel):
    """Class for BD network from the DeepCME paper with n = 5."""

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
            model_name="BD5",
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

        self.add_species([
                S1, S2, S3, S4, S5
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
        
        self.add_reaction([
            S1_birth, S2_birth, S3_birth, S4_birth, S5_birth,
            S1_death, S2_death, S3_death, S4_death, S5_death
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
        self.listOfSpecies['S1'].initial_value = species_initial_value[0]
        self.listOfSpecies['S2'].initial_value = species_initial_value[1]
        self.listOfSpecies['S3'].initial_value = species_initial_value[2]
        self.listOfSpecies['S4'].initial_value = species_initial_value[3]
        self.listOfSpecies['S5'].initial_value = species_initial_value[4]

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
            'S1', 'S2', 'S3', 'S4', 'S5'
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
            'S1', 'S2', 'S3', 'S4', 'S5'
        ]

    # @classmethod
    # def get_randomized_parameters(cls, param_names, n_settings, sigm=0.5):
    #     randomized = {}
    #     for name in param_names:
    #         if name not in cls.params:
    #             raise KeyError(f"Could not find param {name} in {cls.__name__} class `params` dict.")
    #         val = float(cls.params[name])
    #
    #         randomized[name] = np.random.uniform(val * (1. - sigm), val * (1. + sigm), n_settings)
    #     return randomized
    #
    # def set_parameters(self, params_dict):
    #     for name, val in params_dict.items():
    #         if name not in self.listOfParameters:
    #             raise KeyError(
    #                 f"Could not find {name} parameter in {self.__class__.__name__} model listOfParameters.")
    #         self.set_parameter(name, str(val))
