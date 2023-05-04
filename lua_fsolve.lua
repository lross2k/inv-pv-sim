newlib = require 'newlib'

-- Constant values
q   = 1.6e-19 -- elementary charge
K   = 1.4e-23 -- Boltzmann constant
T_n = 298 -- Kelvin

-- Datasheet values
n     = 1.3 -- Diode factor
N_s   = 60 -- Amount of PV cells in series
N_p   = 1 -- Amount of PV cells in parallel
i_sc  = 8.42 -- short circuit current (A)
k_i   = 0.0032 -- short circuit thermal gain
v_oc  = 36.6 -- open circuit voltage of PV module (V) 
R_s   = 0.221 -- series resistance (Ohm)
R_sh  = 415.405 -- shunt resistance (Ohm)
E_go  = 1.1 -- bandgap semiconductor (eV)

-- input values
T     = 25 -- degrees C
G     = 1000 -- solar irradiance

-- Calculated values
T = T + 273.15
-- photocurrent
i_ph = (i_sc + k_i * (T - 298)) * G / 1000
-- reverse saturation current
i_rs = i_sc / (math.exp(q * v_oc / (n * N_s * K * T)) - 1)
-- saturation current
i_0 = i_rs * (T / T_n) ^ 3 * math.exp(q * E_go * (1 / T_n - 1 / T) / (n * K))

v_pv = 3

--print(v_pv,R_s,R_sh)
--print(i_ph,i_0,q,v_pv,R_s,n,K,N_s,T)

a = v_pv/R_sh
b = R_s/R_sh
c = -i_0 * math.exp(q * v_pv / (n * K * N_s * T))
d = q * R_s / (n * K * N_s * T)
e = i_0 + i_ph

print(newlib.c_pv_calc(a,b,c,d,e,0,0))
