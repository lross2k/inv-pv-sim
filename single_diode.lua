-- Constantes fisicas
q    = 1.6e-19 -- elementary charge
k    = 1.38e-23 -- Boltzmann constant

-- Valores correspondientes al panel
n_p  = 72
n_s  = 72 

A    = 1 -- Ideality factor
i_sc = 14.04 -- short circuit current (A)
T_r  = 45 + 273.15 -- degrees C
alpha = i_sc / 1000 / T_r -- mA / C

T    = 45 + 273.15 -- degrees C
S    = 1000 -- solar irradiance
v_oc = 46.73 -- open circuit voltage (V)
--i_rr = math.exp(math.log(i_sc) - v_oc * q / (A * k * T))
i_rr = 5.39e-5

-- current source due to irradiance
i_ph = S * (i_sc + alpha / 1e3 * (T - T_r))
E_g = 1.16 - 7.02 * 10 ^ (-4) * (T ^ 2) / (T - 1108)
i_rs = i_rr * (T / T_r) ^ 3 * math.exp(q * E_g / (k * A) * (1 / T_r - 1 / T))

for v_p=0,50,10 do
	-- power output
	p = v_p * n_p * i_ph - v_p * n_s * i_rs * (math.exp(q * v_p / (A * k * T * n_s)) - 1) 
	print(p)
end

-- https://electricalacademia.com/renewable-energy/photovoltaic-pv-cell-working-characteristics/#:~:text=Iscn%20is%20the%20nominal%20short-circuit%20current%20Ki%20is,ratings%20listed%20for%20commercial%20PV%20cells%20or%20panels.
