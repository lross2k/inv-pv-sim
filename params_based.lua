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

-- Array declaration
i = {}; p = {}; v = {}

-- NUMERICAL METHOD NEEDED HERE PLEASE

-- shunt current
i_sh = (v + i_pv * R_s) / R_sh
-- photovoltaic current
i_pv = i_ph - i_0 * (math.exp(q * (v_pv + i_pv * R_s) / (n * K * N_s * T)) - 1) - i_sh

-- Print generated values of current and power
--print(table.concat(i,", "))
--print(table.concat(p,", "))
--print(table.concat(v,", "))

--[[
require('plplotluac')
-- Parse and process command line arguments
--pl.parseopts( arg, pl.PL_PARSE_FULL )

pl.init() -- Initialize plplot
pl.env(0, 40, 0, 260, 0, 0) -- Create a labelled box to hold the plot.
pl.lab("Tension (V)", "Potencia (W)", "")
pl.line(v, p) -- Plot the data that was prepared above.
pl.plend() -- Close PLplot library

pl.init() -- Initialize plplot
pl.env(0, 40, 0, 9, 0, 0) -- Create a labelled box to hold the plot.
pl.lab("Tension (V)", "Corriente (A)", "")
pl.line(v, i) -- Plot the data that was prepared above.
pl.plend() -- Close PLplot library
]]
