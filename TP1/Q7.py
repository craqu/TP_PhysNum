import sympy as sp

# Define symbols
T = sp.Symbol('T')  # Temperature (or any other dependent variable)
re = sp.Symbol('r_e', real=True, positive=True)  # Classical electron radius
me = sp.Symbol('m_e', real=True, positive=True)  # Electron mass
c = sp.Symbol('c', real=True, positive=True)  # Speed of light
ne = sp.Symbol('n_e', real=True, positive=True)  # Electron number density
mp = sp.Symbol('m_p', real=True, positive=True)  # Proton mass
gamma = sp.Symbol('\gamma', real=True, positive=True)  # Lorentz factor
a = sp.Symbol('a', real=True, positive=True)
I = sp.Symbol('I', real=True, positive=True)
b = sp.Symbol('b', real=True, positive=True)
delta = sp.Symbol('\delta', real=True, positive=True)

# Define components of the equation
ln_term = sp.log((a**2 * (gamma**2 - 1)**2) / (I**2 * (b + delta * gamma)))
numerator = ((-2 * gamma**3 / (gamma**2 - 1)**2 + 2 * gamma / (gamma**2 - 1)) * ln_term
             + (1 / (gamma**2 - 1)) * (4 * gamma - delta / (b + delta * gamma)))
denominator = (ln_term * (gamma**2 / (gamma**2 - 1)) - 2)**2

# Full expression
dS_col_inv_dT = (-1 / (2 * sp.pi * re**2 * me * c**4 * ne * mp)) * (numerator / denominator)

# Display the expression
sp.pprint(dS_col_inv_dT, use_unicode=True)

