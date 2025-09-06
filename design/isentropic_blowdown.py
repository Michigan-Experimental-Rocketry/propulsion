import os
import numpy as np
import plotly.graph_objects as go
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary

# Set REFPROP path
os.environ['RPPREFIX'] = r'/home/jasonyc/REFPROP-cmake/build'
RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])
RP.SETPATHdll(os.environ['RPPREFIX'])
MOLAR_BASE_SI = RP.GETENUMdll(0,"MOLAR BASE SI").iEnum

# Tank parameters
V_tank = 0.002  # m^3
T0 = 298        # K
P0 = 1e7        # Pa (10 MPa initial pressure)
P_atm = 101325  # Pa

fluid = "Nitrogen"
composition = [1.0]

# Initial state
r0 = RP.REFPROPdll(fluid, "TP", "D;S", MOLAR_BASE_SI, 0, 0, T0, P0, composition)
rho0 = r0.Output[0]  # mol/L
S0 = r0.Output[1]    # J/mol-K

# Convert density to mol/m^3
rho0_m3 = rho0 * 1000
n0 = rho0_m3 * V_tank  # initial moles

# Time stepping
dt = 0.1  # seconds
t = [0]
P = [P0]
n = [n0]

while P[-1] > P_atm:
    # Assume small amount of mass leaves (e.g., 0.1% per step)
    n_new = n[-1] * 0.999
    # Find new pressure at constant entropy and tank volume
    # Use REFPROP: Given S, V, and n, find P
    # Guess pressure, iterate to match S0
    P_guess = P[-1]
    for _ in range(20):
        # Calculate density
        rho_guess = n_new / V_tank / 1000  # mol/L
        r = RP.REFPROPdll(fluid, "DS", "P", MOLAR_BASE_SI, 0, 0, rho_guess, S0, composition)
        P_guess = r.Output[0]
    t.append(t[-1] + dt)
    P.append(P_guess)
    n.append(n_new)

# Plot using Plotly
fig = go.Figure()
fig.add_trace(go.Scatter(x=t, y=np.array(P)/1e5, mode='lines', name='Pressure'))
fig.update_layout(
    title='Isentropic Blowdown of GN2 Tank',
    xaxis_title='Time (s)',
    yaxis_title='Pressure (bar)',
    template='plotly_white'
)
fig.show()
