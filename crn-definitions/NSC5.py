import gillespy2 as gillespy
import numpy as np
import math

import base

# test with number of species = 2, 5, 10

# BD stands for Nonlinear Signalling Cascade
class NSC5(base.BaseCRNModel):
    """Class for NSC5 network from the DeepCME paper with n = 2."""

    params = {
        'b': 1.,
        'k_m': 100.,
        'k_0': 10.,
        'H': 1.,
        'beta_0': 10.,
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
            model_name="NSC5",
        )

        b = gillespy.Parameter(name='b', expression=self.params['b'])
        k_m = gillespy.Parameter(name='k_m', expression=self.params['k_m'])
        k_0 = gillespy.Parameter(name='k_0', expression=self.params['k_0'])
        H = gillespy.Parameter(name='H', expression=self.params['H'])
        beta_0 = gillespy.Parameter(name='beta_0', expression=self.params['beta_0'])
        gamma = gillespy.Parameter(name='gamma', expression=self.params['gamma'])
        self.add_parameter([b, k_m, k_0, H, beta_0, gamma])

        # Species
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
            rate=beta_0,
            massaction=True
        )

        # transitions
        S1_to_S2 = gillespy.Reaction(
            name='S1_to_S2',
            reactants={S1: 1},
            products={S1:1, S2: 1},
            propensity_function='b + (k_m * pow(S1,H))/(k_0 + pow(S1,H))'
        )
        S2_to_S3 = gillespy.Reaction(
            name='S2_to_S3',
            reactants={S2: 1},
            products={S2: 1, S3: 1},
            propensity_function='b + (k_m * pow(S2,H))/(k_0 + pow(S2,H))'
        )
        S3_to_S4 = gillespy.Reaction(
            name='S3_to_S4',
            reactants={S3: 1},
            products={S3: 1, S4: 1},
            propensity_function='b + (k_m * pow(S3,H))/(k_0 + pow(S3,H))'
        )
        S4_to_S5 = gillespy.Reaction(
            name='S4_to_S5',
            reactants={S4: 1},
            products={S4: 1, S5: 1},
            propensity_function='b + (k_m * pow(S4,H))/(k_0 + pow(S4,H))'
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
            S1_birth, S1_to_S2, S2_to_S3, S3_to_S4, S4_to_S5,
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
        # self.listOfSpecies['S1'].initial_value = species_initial_value[0]
        # self.listOfSpecies['S2'].initial_value = species_initial_value[1]
        # self.listOfSpecies['S3'].initial_value = species_initial_value[2]
        # self.listOfSpecies['S4'].initial_value = species_initial_value[3]
        # self.listOfSpecies['S5'].initial_value = species_initial_value[4]
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

    @classmethod
    def get_randomized_parameters(cls, param_names, n_settings, sigm=0.1):
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
