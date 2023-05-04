-- Constant values
q   = 1.6e-19 -- elementary charge
k   = 1.38e-23 -- Boltzmann constant
T_r = 25 + 273.15 -- degrees C

-- Datasheet values
A     = 1.2 -- Ideality factor
n_p   = 60 -- Amount of PV cells in PV module
i_sc  = 8.49 -- short circuit current (A)
alpha = 6 -- mA / C
T     = 25 + 273.15 -- degrees C
S     = 1000 -- solar irradiance
v_oc  = 37.5 -- open circuit voltage of PV module (V) 

-- Calculated values
v_oc = v_oc / n_p -- open circuit voltage of PV Cell (V) 
i_rr = math.exp(math.log(i_sc) - v_oc * q / (A * k * T))
E_g  = 1.16 - 7.02e-4 * (T ^ 2) / (T - 1108)
i_rs = i_rr * (T / T_r) ^ 3 * math.exp(q * E_g / (k * A) * (1 / T_r - 1 / T))
n_s  = n_p
i_ph = (i_sc + alpha * (T - T_r)) * S / 1e3 -- Photovoltaic current source

-- Array declaration
i = {}; p = {}; v = {}

-- Generate values for a range of the equation
index = 0
for v_p=0,52,v_oc do
    i_calc = i_ph - i_rs * (math.exp(q * v_p / (A * k * T * n_s)) - 1)
	if i_calc >= 0 then
    	i[index] = i_calc
	    p[index] = i_calc * v_p
	    v[index] = v_p
	    index = index + 1
	else
	    break
	end
end

-- Print generated values of current and power
--print(table.concat(i,", "))
--print(table.concat(p,", "))
--print(table.concat(v,", "))

-- Using file I/O
file = io.open("single_diode.txt", "w")
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
--]]
