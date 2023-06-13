--params = require 'params'

-- Constant values
q   = 1.6e-19 -- elementary charge
K   = 1.4e-23 -- Boltzmann constant

-- Datasheet values
G_n   = 1000
T_n = 298 -- Kelvin
a_1   = 1 -- ideality factor diode (assumed)
a_2   = 2 -- ideality factor diode 2 (assumed)
N_s   = 36 -- Amount of PV cells in series
i_scn = 3.8 -- short circuit current (A)
k_i   = 0.003 -- short circuit thermal gain (A/K)
k_v   = -80e-3  -- (V/K)
v_oc  = 21.1 -- open circuit voltage of PV module (V)
i_mp  = 3.5 -- maximum power point current (A)
v_mp  = 17.1 -- maximum power point voltage (V)
p_max = 60 -- maximum power (W)

-- input values
T     = 25 -- degrees C
G     = 1000 -- solar irradiance

-- Calculated values
v_oc_mp = v_oc
v_oc = v_oc / N_s
T = T + 273.15 -- Kelvin
V_T = N_s * K * T / q
R_s = 0
I_pvn = 1 / i_mp -- Starting value equals to 1 / i_mp
DeltaT = T - T_n
I_0 = (i_mp + k_i * DeltaT) / math.exp((v_oc_mp + k_v * DeltaT) / (a_1 * V_T) - 1)
I_01 = I_0 * (i_mp + k_i * DeltaT) / math.exp(v_mp + k_v * DeltaT / ((a_1 + a_2) * V_T) - 1)
P_maxe = i_mp * v_mp
I_pv = (I_pvn + k_i * DeltaT) * (G/G_n)
I_sc = (i_scn + k_i * DeltaT) * (G/G_n)  
-- Solve the system with newton rhapson

I = I_pv - I_01 * (math.exp((V + I * R_s)/(a1 * V_T)) - 1) - I_02 * (math.exp((V + I * R_s)/(a2 * V_T)) - 1) - (V + I * R_s)/R_p







-- photocurrent
i_ph = (i_sc + k_i * (T - 298)) * G / 1000
-- reverse saturation current
i_rs = i_sc / (math.exp(q * v_oc / (n * N_s * K * T)) - 1)
-- saturation current
i_0 = i_rs * (T / T_n) ^ 3 * math.exp(q * E_go * (1 / T_n - 1 / T) / (n * K))

-- Array declaration
i = {}; p = {}; v = {}

v_oc = v_oc / N_s
b = R_s/R_sh
d = q * R_s / (n * K * N_s * T)
e = i_0 + i_ph
-- Generate values for a range of the equation
index = 1
for v_pv=0,52,v_oc do
  a = v_pv/R_sh
  c = -i_0 * math.exp(q * v_pv / (n * K * N_s * T))
  i_pv,i_sh = params.solve_currents(a,b,c,d,e,0,0)
	if i_pv >= 0 then
    	i[index] = i_pv
	    p[index] = i_pv * v_pv
	    v[index] = v_pv
	    index = index + 1
	end
end

-- Print generated values of current and power
print(table.concat(i,", "))
print(table.concat(p,", "))
print(table.concat(v,", "))

-- Using file I/O
file = io.open("params.txt", "w")
for iter=1,index-1 do
  file:write(v[iter])
  file:write("\t")
  file:write(i[iter])
  file:write("\t")
  file:write(p[iter])
  file:write("\n")
end

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
