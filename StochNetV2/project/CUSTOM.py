import gillespy2 as gp
import numpy as np
import math

import base


class CUSTOM(base.BaseCRNModel):
    """Class for SIR network."""

    params = {
        'alpha': '2.', # recovery rate
        'beta': '1.', # infection rate
        'gamma': '0.1' # death rate
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
            model_name="CUSTOM",
        )

        alpha = gp.Parameter(name='alpha', expression=self.params['alpha'])
        beta = gp.Parameter(name='beta', expression=self.params['beta'])
        gamma = gp.Parameter(name='gamma', expression=self.params['gamma'])
        self.add_parameter([alpha, beta, gamma])

        # Species
        S = gp.Species(name='S', initial_value=50)
        I = gp.Species(name='I', initial_value=100)
        D = gp.Species(name='D', initial_value=0)
        self.add_species([S, I, D])

        # Reactions
        infect = gp.Reaction(
            name='infect',
            reactants={S: 1, I: 1},
            products={I: 2},
            propensity_function='beta*S*I/(S+I+D)',
        )
        recover = gp.Reaction(
            name='recover',
            reactants={I: 1},
            products={S: 1},
            rate=alpha,
        )
        die = gp.Reaction(
            name='die',
            reactants={I: 1},
            products={D: 1},
            rate=gamma,
        )
        self.add_reaction([infect, recover, die])
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
        self.listOfSpecies['S'].initial_value = species_initial_value[0]
        self.listOfSpecies['I'].initial_value = species_initial_value[1]
        self.listOfSpecies['D'].initial_value = species_initial_value[2]

    @staticmethod
    def get_species_names():
        """
        Returns list of species names.

        Returns
        -------
        list of all species names. The order of names should be coherent
        with the list returned by get_initial_state method.

        """
        return ['S', 'I', 'D']

    @staticmethod
    def get_initial_state():
        return [50, 100, 0]

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
        settings = np.random.randint(low=10, high=300, size=(n_settings, n_species-1))
        settings = np.hstack([settings, np.zeros((n_settings, 1))])
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
        return ['S', 'I', 'D']

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
