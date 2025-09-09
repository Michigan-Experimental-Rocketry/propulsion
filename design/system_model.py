import matplotlib.pyplot as plt
import numpy as np
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary
import os
from pint import UnitRegistry, Quantity

os.environ['RPPREFIX'] = r'/home/jasonyc/REFPROP-cmake/build'
ureg = UnitRegistry()
Q_ = ureg.Quantity
SATURATION_DELTA = Q_(1, 'Pa')  # Small pressure difference to move between vapor or liquid phase


class State:
    """
    Static class variables for running REFPROP across all instances.
    NOTE: All internal member variables of all classes will be in base SI units.
    """
    _RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])
    _RP.SETPATHdll(os.environ['RPPREFIX'])
    _MASS_BASE_SI = _RP.GETENUMdll(0,"MASS BASE SI").iEnum
    @property
    def RP(self):
        return type(self)._RP
    @property
    def MASS_BASE_SI(self):
        return type(self)._MASS_BASE_SI


class GasNode(State):
    def __init__(self, P: Quantity, T: Quantity, V: Quantity, species: str = 'nitrogen'):
        self.species = species
        self.P = P
        self.T = T
        self.V = V
        self.rho = Q_(self.RP.REFPROPdll(self.species,"PT","D",self.MASS_BASE_SI,0,0,
                                         self.P.to('Pa').magnitude,self.T.to('K').magnitude,[1.0]).Output[0], 'kg/m^3')
        self.e = Q_(self.RP.REFPROPdll(self.species,"PT","E",self.MASS_BASE_SI,0,0,
                                       self.P.to('Pa').magnitude,self.T.to('K').magnitude,[1.0]).Output[0], 'J/kg')
        self.n = self.rho*self.V
        self.s = Q_(self.RP.REFPROPdll(self.species,"PT","S",self.MASS_BASE_SI,0,0,
                                       self.P.to('Pa').magnitude,self.T.to('K').magnitude,[1.0]).Output[0], 'J/kg/K')

    def update_state(self, delta_mass: Quantity):
        """Isentropic expansion (mass removal)"""
        self.n += delta_mass
        self.rho = self.n/self.V
        self.s = self.s

        self.P = Q_(self.RP.REFPROPdll(self.species,"DS","P",self.MASS_BASE_SI,0,0,
                                       self.rho.to('kg/m^3').magnitude,self.s.to('J/kg/K').magnitude,[1.0]).Output[0], 'Pa')
        self.T = Q_(self.RP.REFPROPdll(self.species,"DS","T",self.MASS_BASE_SI,0,0,
                                       self.rho.to('kg/m^3').magnitude,self.s.to('J/kg/K').magnitude,[1.0]).Output[0], 'K')
        self.e = Q_(self.RP.REFPROPdll(self.species,"DS","E",self.MASS_BASE_SI,0,0,
                                       self.rho.to('kg/m^3').magnitude,self.s.to('J/kg/K').magnitude,[1.0]).Output[0], 'J/kg')

    def get_state_string(self) -> str:
        return (f"[STATE] P: {self.P.to('psi').magnitude:.2f} psi | T: {self.T.to('K').magnitude:.2f} K | D: {self.rho.to('kg/m^3').magnitude:.2f} kg/m^3")


class LiquidNode(State):
    """Models an incompressible liquid node"""
    def __init__(self, P: Quantity, T: Quantity, V: Quantity, species: str = 'water'):
        self.species = species
        self.P = P
        self.T = T
        self.V = V
        self.rho = Q_(self.RP.REFPROPdll(self.species,"PT","D",self.MASS_BASE_SI,0,0,
                                         self.P.to('Pa').magnitude,self.T.to('K').magnitude,[1.0]).Output[0], 'kg/m^3')
        self.n = self.rho*self.V

    @classmethod
    def from_quality(cls, T: Quantity, V: Quantity, quality: float = 0.0, species: str = 'water'):
        """Constructor for a liquid, assuming it is at saturation temperature"""
        state = State()
        sat_pressure = Q_(state.RP.REFPROPdll(species,"TQ","P",state.MASS_BASE_SI,0,0,
                                              T.to('K').magnitude,quality,[1.0]).Output[0], 'Pa')
        return cls(sat_pressure + SATURATION_DELTA, T, V, species)

    def remove_mass(self, delta_mass: Quantity):
        """Remove mass at constant temperature, and therefore density"""
        self.n -= delta_mass
        delta_V = delta_mass/self.rho
        self.V -= delta_V
        return delta_V
    
    def get_latent_heat_of_vaporization(self):
        """Assuming we're at saturation conditions, how much enthalpy it take to phase change"""
        h_vap = Q_(self.RP.REFPROPdll(self.species,"TQ","H",self.MASS_BASE_SI,0,0,
                                      1.0,self.T.to('K').magnitude,[1.0]).Output[0], 'J/kg')
        h_liq = Q_(self.RP.REFPROPdll(self.species,"TQ","H",self.MASS_BASE_SI,0,0,
                                      0.0,self.T.to('K').magnitude,[1.0]).Output[0], 'J/kg')
        return h_vap - h_liq


if __name__ == "__main__":
    # Target saturation: 257.33 K, 2.0313 MPa
    liquid = LiquidNode.from_quality(Q_(258, 'K'), Q_(3, 'L'), quality=0.0, species='nitrous oxide')
    ullage = GasNode((liquid.P - 2*SATURATION_DELTA), Q_(258, 'K'), Q_(1, 'L'), species='nitrous oxide')
