using DifferentialEquations

p = 10 # Pressure in Pa
T_g = 300 # Gas temperature in K
const k_B = 1.38e-23 # Boltzmann constant
const eps0 = 8.85e-12 # Permittivity of vacuum
const c = 3e8 # speed of light
const h = 6.62e-34 # Planck constant

const me = 9.11e-31 # Electron mass in kg
const mAr = 1.67e-27 * 40 # Ar mass in kg
const e = 1.6e-19 # electron charge
P_d = 0.1 # power deposition per unit volume Wm-3

n_g = p/k_B/T_g # ideal gas law for gas density


#k_laser = a * exp(-(t-b)^2/2/gaussc^2) * y[3]

# The function "chemistry" below defines the system of ordinary differential equations to be solved
# Here, the densities of electrons, Ar, Ar+, Ar*, Ar2* and Ar2+ are solved as a function of time
# The electron temperature, Te, is also solved for
# The values of each are given by y[1, 2,..]
# The main requirement for the solution are source terms for each species density and the electron temperature
# these are defined as S_e, S_ee, S_Ar etc
# These source terms constitute dy[1, 2, ...]

function chemistry(dy,y,p,t)

 # Below, the rate coefficients for electron-heavy particle and heavy particle - heavy particle
 # reactions are defined
 # electron impact impact reactions are a function of electron temperature  
 # These rate coefficients define the rates of production and loss of particles and energy in the system i.e. the source terms
 # R1 e + Ar -> e + Ar (elastic collisions)    
 keg = 2.336e-14*(y[2]/11604)^1.609

 # R2 e + Ar -> e + Ar*
 kx = 2.48e-14*(y[2]/11604)^0.33*exp(-12.78/(y[2]/11604))

 # R3 e + Ar -> Ar+ + e + e
 kiz = 2.34e-14 * (y[2]/11604)^0.59 * exp(-17.44/(y[2]/11604))

 # R4 e + e + Ar+ -> Ar + e
 krec = 5e-39 * (y[2]/11604)^4.5

 # R5 Ar* + Ar* -> Ar+ + Ar + e
 kmm_iz = 1.2e-15*(T_g/300)^0.5

 # R6 Ar* + Ar + Ar -> Ar2* + Ar
 km_assoc3 = 1.14e-44*(T_g/300)^-1

 # R7 Ar* + Ar* -> Ar2+ + e
 km_assoc2 = 5.7e-16*(T_g/300)^0.5

 # R8 Ar+ + Ar + Ar -> Ar2+ + Ar
 kp_assoc = 2.5e-43*(T_g/300)^-1

 # R9 Ar2* + Ar2* -> Ar2+ + Ar + Ar + e
 km2_assoc = 5e-16*(T_g/300)^0.5

 # R10 e + Ar2+ -> Ar + Ar
 krec_2 = 2.69e-14*(y[2]/11604)^-0.67
 
 # R11 e + Ar* -> Ar+ + e + e
 kem = 2.71e-13*(y[2]/11604)^0.26*exp(-4.49/(y[2]/11604))
 
 # hypothetical rate of production of species by laser excitation (not used here)
 k_laser = 0.0

 # prefactor needed for electron energy balance equation
 C_e = 3 / 2 * k_B * y[1]
 
 # loss of electron energy in elastic collisions (R1)
 Q_eg = 3 / 2 * k_B * y[1] * keg * y[3] * (y[4] - T_g) * 2 * me / mAr

 # loss of electron energy due to electron impact ionization (R3)
 Q_iz = e * 15.76 * kiz * y[3] * y[1]

 # loss of electron energy during electron impact excitation of Ar (R2)
 Q_x = e * 12.14 * kx * y[3] * y[1]

 # hypothetical gain or loss of energy of electrons due to laser (not used here)
 Q_laser = e * 1 * k_laser

 # Rate of production and loss of electrons due to all elementary reactions
 S_e = y[1]*y[3]*kiz - y[1]* y[1]*y[4]*krec + k_laser + y[5]*y[5]*kmm_iz +
       y[1]*y[5]*kem + y[5]*y[5]*km_assoc2 + y[6]*y[6]*km2_assoc -  y[1]*y[7]*krec_2

 # Rate of production and loss of electron energy due to all elementary reactions
 # Notice that the power per unit volume, P_d, acts as the source term for electron energy
 # We divide the source terms by the pre-factor C_e, defined above
 S_ee = (P_d - Q_iz - Q_x - Q_eg + Q_laser) / C_e

 # Rate of production and loss of Ar due to all elementary reactions
 S_Ar =  y[1]* y[1]*y[4]*krec -  y[1]*y[3]*kiz - k_laser + y[5]*y[5]*kmm_iz +
       y[5]*y[3]*y[3]*km_assoc3 - y[4]*y[3]*y[3]*kp_assoc + 2* y[1]*y[6]*krec_2 +
       2*y[6]*y[6]*km2_assoc

 # Rate of production and loss of Ar+ due to all elementary reactions      
 S_Arp =  y[1]*y[3]*kiz -  y[1]* y[1]*y[4]*krec + k_laser + y[5]*y[5]*kmm_iz +
        y[1]*y[5]*kem - y[4]*y[3]*y[3]*kp_assoc

 # Rate of production and loss of Ar* due to all elementary reactions       
 S_Arm =  y[1]*y[3]*kx - y[5]*y[5]*kmm_iz -  y[1]*y[5]*kem - y[5]*y[3]*y[3]*km_assoc3

 # Rate of production and loss of Ar2* due to all elementary reactions
 S_Ar2m = y[5]*y[3]*y[3]*km_assoc3 - y[6]*y[6]*km2_assoc

 # Rate of production and loss of Ar2+ due to all elementary reactions
 S_Ar2p = y[5]*y[5]*km_assoc2 + y[4]*y[3]*y[3]*kp_assoc -  y[1]*y[7]*krec_2


 # Set derivatives equal to source functions defined above
 dy[1]  = S_e
 dy[2]  = S_ee
 dy[3]  = S_Ar
 dy[4]  = S_Arp
 dy[5]  = S_Arm
 dy[6]  = S_Ar2m
 dy[7]  = S_Ar2p
end

# define initial conditions of all quantities
# For Ar, this should be initial gas density, if we assume plasma is sustained in pure Ar
# For species densities, negligible initial densities are to be assumed at the beginning
# In this way, when we apply power, at the start of the simulation, the plasma can "ignite"
# For electron temperature, it also makes sense to assume a lower value than would be expected
# when the plasma turns on
u0 = zeros(7)
u0[1] = 1e6
u0[2] = 1000
u0[3] = n_g
u0[4] = 1e6
u0[5] = 1e6
u0[6] = 1e6
u0[7] = 1e6

# Here, the function and initial conditions are converted into an ODEProblem using the DifferentialEquations package
prob_ode_chemistry = ODEProblem(chemistry,
                                u0,(0.0,1.0))

# Here, the ODEProblem is solved using the DifferentialEquations package
sol = solve(prob_ode_chemistry)
