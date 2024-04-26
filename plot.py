import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def snells_law(beta1, n2):
    beta2 = np.arcsin((1/n2) * np.sin(beta1))
    gamma1 = np.arcsin(n2 * np.sin(np.pi/3 - beta2))
    return beta1 + gamma1 - np.pi/3

def n(L, A1, A2):
    return A1 + A2 / L**2

csv_file_path = '/Users/edvin/git/waves_and_optics/optik/data.csv'

with open(csv_file_path, mode='r') as file:
    beta_deg_lst = []
    x_lst = []

    csv_reader = csv.reader(file, delimiter='\t')
    next(csv_reader)
    
    for row in csv_reader:
        beta, x = row  # Unpack the values from the row
        beta_deg_lst.append(float(beta))
        x_lst.append(float(x))

beta = [np.deg2rad(k) for k in beta_deg_lst]
delta = [np.arctan(i/14.5) for i in x_lst] # Distance between prism and image ~ 14.5 cm

fit_vals = curve_fit(snells_law, beta, delta, p0=[1.6])[0]
x_new = np.linspace(min(beta), max(beta), 100)

plt.plot(x_new, snells_law(x_new, *fit_vals), color='red', label='Fitted Data')
plt.scatter(beta, delta, marker='o', color='blue')
plt.xlabel('Incoming Angle (beta) [rad]', size=12)
plt.ylabel('Deviation (delta) [rad]', size=12)
plt.show()

delta_min = min(snells_law(x_new, *fit_vals))
n_prism = np.sin((delta_min + np.pi/3)/2) / np.sin(np.pi/6)

# Compare the calculated n-value with glass properties of some known glass, for the specific wavelength
glass_dict = {
    "Fused silica": (1.4580, 0.00354),
    "Borosilicate glass BK7": (1.5046, 0.00420),
    "Hard crown glass K5": (1.5220, 0.00459),
    "Barium crown glass BaK4": (1.5690, 0.00531),
    "Barium flint glass BaF10": (1.6700, 0.00743),
    "Dense flint glass SF10": (1.7280, 0.01342)
}

wavelength = 641 # nm (nanometers)
n_values = {}

for glass, A_value in glass_dict.items():
    n_values[glass] = n(wavelength, A_value[0], A_value[1])

closest_glass = min(n_values, key=lambda glass: abs(n_values[glass] - n_prism))

print(f'The prism with a calculated n-value of approximately {round(n_prism, 3)} is most similar to {closest_glass} with an n-value of {round(n_values[closest_glass], 3)}')
print(f'The minimum deviation was {round(np.rad2deg(delta_min), 3)} degrees')