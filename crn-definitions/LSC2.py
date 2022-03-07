import gillespy2 as gillespy
import numpy as np
import math

import base

# test with number of species = 2, 5, 10

# BD stands for Linear Signalling Cascade
class LSC2(base.BaseCRNModel):
    """Class for LSC2 network from the DeepCME paper with n = 2."""

    params = {
        'beta_0': 10.,
        'k': 5.,
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
            model_name="LSC2",
        )

        beta_0 = gillespy.Parameter(name='beta_0', expression=self.params['beta_0'])
        k = gillespy.Parameter(name='k', expression=self.params['k'])
        gamma = gillespy.Parameter(name='gamma', expression=self.params['gamma'])
        self.add_parameter([beta_0, k, gamma])

        # Species
        S1 = gillespy.Species(name='S1', initial_value=0)
        S2 = gillespy.Species(name='S2', initial_value=0)

        self.add_species([
                S1, S2
            ])

        # Reactions
        # births
        S1_birth = gillespy.Reaction(
            name='S1_birth',
            reactants={},
            products={S1: 1},
            rate=beta_0,
            massaction=True
        )

        # transitions
        S1_to_S2 = gillespy.Reaction(
            name='S1_to_S2',
            reactants={S1: 1},
            products={S1:1, S2: 1},
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

        
        self.add_reaction([
            S1_birth, S1_to_S2, S1_death, S2_death
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
            'S1', 'S2'
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
            'S1', 'S2'
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
