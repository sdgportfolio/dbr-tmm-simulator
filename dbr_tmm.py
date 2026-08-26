"""
DBR-TMM: Distributed Bragg Reflector Simulator
================================================
Transfer Matrix Method (TMM) simulation of reflectance and transmittance spectra for a DBR stack.
The dependence of reflectance on the number of pairs can be visualized. 
The stopband edges and bandwidth are also calculated.

Author: Shroyon Dasgupta
GitHub: github.com/sdgportfolio/dbr-tmm-simulator
"""


#Importing the libraries
import numpy as np
import matplotlib.pyplot as plt


#--- INPUTS ---
#Input the target wavelength
lam0 = 1550e-9 #Target wavelength

#Input the refractive indices
n_high = 2.5 #Refractive index of the higher index layer
n_low = 1.5 #Refractive index of the lower index layer
ni = 1 #Refractive index of the incident medium
ns = 1.5 #Refractive index of the substrate

#Input the number of pairs of layers
N_lp = 15 #Number of pairs in the DBR

#Input the angle of incidence and polarization
theta = 0 #Angle of incidence in degrees
polarization = 's' #Polarization ('s'/'p' or 'TE'/'TM')


#--- CALCULATIONS ---
t1 = lam0/(4 * n_high) #Quarter-wave condition of the higher index layer
t2 = lam0/(4 * n_low) #Quarter-wave condition of the lower index layer

bandwidth = lam0 * (4/np.pi) * np.arcsin((n_high - n_low)/(n_high + n_low)) #Bandwidth of an infinite DBR

#Transfer matrix of a single layer
def single_layer_matrix(n, t, lam, theta_i, n_i, pol):
    #n of the layer, thickness of the layer, wavelength, angle of incidence, refractive index of incident medium, polarization
    sin_theta_n = (n_i/n) * np.sin(theta_i) #Calculating the sine of the angle in the layer
    cos_theta_n = np.cos(np.arcsin(sin_theta_n + 0j)) #Calculating the cosine of the angle in the layer

    delta = 2 * np.pi * n * t * cos_theta_n/lam #Calculating the phase accumulated in the layer

    if pol == 's' or pol == 'TE': #Transfer matrix for s-polarization or TE polarization
        M_sl = np.array([
            [np.cos(delta), 1j * np.sin(delta)/(n * cos_theta_n)],
            [1j * n * cos_theta_n * np.sin(delta), np.cos(delta)]
        ], dtype=complex)
    elif pol == 'p' or pol == 'TM': #Transfer matrix for p-polarization or TM polarization
        M_sl = np.array([
            [np.cos(delta), 1j * np.sin(delta) * cos_theta_n/n],
            [1j * n * np.sin(delta)/cos_theta_n, np.cos(delta)]
        ], dtype=complex) 

    return M_sl #returning the transfer matrix of a single layer

def dbr_matrix(nh, nl, th, tl, lam, N, theta_i, n_i, n_s, pol):
    #n of high index layer, n of low index layer, thickness of high index layer, thickness of low index layer, wavelength, number of pairs, angle of incidence, refractive index of incident medium, refractive index of substrate, polarization
    wavelengths = np.linspace(0.5 * lam, 1.5 * lam, 1000) #Creating an array of wavelengths from 0.5*lam to 1.5*lam with 1000 points
    R = np.zeros_like(wavelengths, dtype=float) #Creating the array to store reflectance values
    T = np.zeros_like(wavelengths, dtype=float) #Creating the array to store transmittance values

    sin_i = np.sin(theta_i) #sin of the angle of incidence
    cos_i = np.cos(theta_i) #cos of the angle of incidence
    sin_s = (n_i/n_s) * sin_i #sin of the angle in the substrate
    cos_s = np.cos(np.arcsin(sin_s + 0j)) #cos of the angle in the substrate

    #Calculating polarization-based optical admittance
    if pol == 's' or pol == 'TE':
        eta_i = n_i * cos_i
        eta_s = n_s * cos_s
    elif pol == 'p' or pol == 'TM':
        eta_i = n_i/cos_i
        eta_s = n_s/cos_s

    for i, wav in enumerate(wavelengths): #Looping through the wavelengths
        M_total = np.eye(2, dtype=complex) #Initializing the total transfer matrix as identity matrix
        M_high = single_layer_matrix(nh, th, wav, theta_i, n_i, pol) #Calculating the transfer matrix of the high index layer
        M_low = single_layer_matrix(nl, tl, wav, theta_i, n_i, pol) #Calculating the transfer matrix of the low index layer
        M_pair = M_high @ M_low #Calculating the transfer matrix of a pair of layers

        for j in range(N): #Looping through the number of pairs
            M_total = M_total @ M_pair #Multiplying the total transfer matrix with the transfer matrix of the pair

        #Assigning variables to matrix elements for ease of calculation
        A = M_total[0, 0]
        B = M_total[0, 1]
        C = M_total[1, 0]
        D = M_total[1, 1]

        R[i] = np.abs((eta_i * A + eta_i * eta_s * B - C - eta_s * D) / (eta_i * A + eta_i * eta_s * B + C + eta_s * D))**2 #Calculating the reflectance at a particular wavelength
        T[i] = (eta_s.real/eta_i.real) * np.abs((2 * eta_i)/(eta_i * A + eta_i * eta_s * B + C + eta_s * D))**2 #Calculating the transmittance at a particular wavelength

    return wavelengths, R, T #returning the reflectance and transmittance arrays corresponding to the wavelengths


#--- STOPBAND CALCULATION ---
def stopband_edges(wavelengths, R, threshold):
    if np.max(R) < threshold: #Checking if the maximum reflectance is less than the threshold
        print("Maximum reflectance is less than the threshold. ") #Printing a message if no stopband edges are found
        return None, None #Returning None for both stopband edges
    else:
        idx = np.where(R >= threshold)[0] #Creating an array of indices where the reflectance is greater than or equal to the threshold
        return wavelengths[idx[0]], wavelengths[idx[-1]] #Returning the wavelengths corresponding to the extreme values

#--- PLOTS ---
def plot_dbr_wavelength_spectrum(wavelengths, R, T, lambda0):
    fig, ax = plt.subplots(figsize=(8, 5)) #Plotting the reflectance and transmittance spectra in a single figure
    ax.plot(wavelengths * 1e9, R, label="Reflectance (R)", color="red") #Reflectance plot vs wavelength
    ax.plot(wavelengths * 1e9, T, label="Transmittance (T)", color="blue", ls="--") #Transmittance plot vs wavelength
    ax.axvline(lambda0 * 1e9, color="black", ls=":", label=r"$\lambda_0$ = {:.0f} nm".format(lambda0*1e9)) #Adding a vertical line at the target wavelength
    ax.set_xlabel("Wavelength (nm)", fontsize=12) #Setting the x-axis label
    ax.set_ylabel("Reflectance / Transmittance", fontsize=12) #Setting the y-axis label
    ax.set_title("DBR reflectance and transmittance wavelength dependence", fontsize=14) #Setting the title of the plot
    ax.legend(loc="upper right") #Adding plot legend
    plt.show() #Displaying the plot

def plot_dbr_frequency_spectrum(wavelengths, R, T, lambda0):
    freq = 3e8/wavelengths #Converting wavelengths to frequencies
    f0 = 3e8/lambda0 #Calculating central frequency

    #Reversing the frequency, reflectance, and transmittance arrays to sort it correctly in ascending order of frequency
    freq = freq[::-1]
    R= R[::-1]
    T = T[::-1]

    #Restrict to the range 0.5*f0 to 1.5*f0 for a cleaner plot
    mask = (freq >= 0.5 * f0) & (freq <= 1.5 * f0)
    freq = freq[mask]
    R = R[mask]
    T = T[mask]

    fig, ax = plt.subplots(figsize=(8, 5)) #Plotting the reflectance and transmittance spectra in a single figure
    ax.plot(freq / 1e12, R, label="Reflectance (R)", color="red") #Reflectance plot vs frequency
    ax.plot(freq / 1e12, T, label="Transmittance (T)", color="blue", ls="--") #Transmittance plot vs frequency
    ax.axvline(f0 / 1e12, color="black", ls=":", label=r"$\nu_0$ = {:.1f} THz".format(f0/1e12)) #Adding a vertical line at the central frequency
    ax.set_xlabel("Frequency (THz)", fontsize=12) #Setting the x-axis label
    ax.set_ylabel("Reflectance / Transmittance", fontsize=12) #Setting the y-axis label
    ax.set_title("DBR reflectance and transmittance frequency dependence", fontsize=14) #Setting the title of the plot
    ax.legend(loc="upper right") #Adding plot legend
    plt.show() #Displaying the plot

def plot_reflectance_and_N_pairs(nh, nl, th, tl, lam, N, theta_i, n_i, n_s, pol):
    fig, ax = plt.subplots(figsize=(8, 5)) #Plotting the reflectance and DBR pair dependence in a single figure

    N_val = [N//5, N//4, N//2, N] #Creating an array of number of pairs
    N_val = np.unique(N_val) #Removing duplicates from the array of number of pairs
    colors = ['blue', 'orange', 'green', 'red'] #Creating an array of colors

    for i, n in enumerate(N_val): #Looping through the number of pairs
        wavelengths, R, T = dbr_matrix(nh, nl, th, tl, lam, n, theta_i, n_i, n_s, pol) #Getting the reflectance spectra for the pair
        ax.plot(wavelengths * 1e9, R, label=f"N = {n}", color=colors[i]) #Plotting the reflectance spectra for the pair

    ax.axvline(lam0 * 1e9, color="black", ls=":", label=r"$\lambda_0$") #Adding a vertical line at the target wavelength
    ax.set_xlabel("Wavelength (nm)", fontsize=12) #Setting the x-axis label
    ax.set_ylabel("Reflectance (R)", fontsize=12) #Setting the y-axis label
    ax.set_title("Reflectance Spectrum vs. N (Wavelength Domain)", fontsize=14) #Setting the title of the plot
    ax.legend(loc="upper right", fontsize=8, ncol=2) #Adding plot legend
    plt.show() #Displaying the plot


#--- EXPORT ---
def export_to_csv(wavelengths, R, T, filename="dbr_tmm_results.csv"):
    freq = 3e8 / wavelengths #Converting wavelengths to frequency for reference
    data = np.column_stack([wavelengths * 1e9, freq / 1e12, R, T]) #Stacking columns: wavelength(nm), freq(THz), R, T
    header = "wavelength_nm,frequency_THz,reflectance,transmittance" #Creating the header for the CSV file
    np.savetxt(filename, data, delimiter=",", header=header, comments="") #Saving the data to a CSV file with the specified filename
    print(f"Saved results to {filename}") #Printing a message indicating that the results have been saved to the specified CSV file


#--- OUTPUTS ---
wavelengths, R, T = dbr_matrix(n_high, n_low, t1, t2, lam0, N_lp, np.radians(theta), ni, ns, polarization) #Getting the reflectance and transmittance spectra
stopband_min, stopband_max = stopband_edges(wavelengths, R, threshold=0.95) #Calculating the stopband edges where reflectance is greater than or equal to 95%
print(f"Thickness of high index layer: {t1*1e9:.2f} nm") #Printing the thickness of the high index layer in nanometers
print(f"Thickness of low index layer: {t2*1e9:.2f} nm") #Printing the thickness of the low index layer in nanometers
if stopband_min is not None: #Checking if stopband exists
    print(f"Stopband edges: {stopband_min*1e9:.2f} nm to {stopband_max*1e9:.2f} nm") #Printing the stopband edges in nanometers
    print(f"Bandwidth: {(stopband_max - stopband_min)*1e9:.2f} nm") #Printing the bandwidth of the stopband in nanometers
else:
    print("No stopband edges found.") #Printing a message if no stopband edges are found
print(f"Bandwidth: {bandwidth*1e9:.2f} nm") #Printing the theoretical bandwidth of an infinite DBR in nanometers
plot_dbr_wavelength_spectrum(wavelengths, R, T, lam0) #Plotting the reflectance and transmittance spectra in wavelength domain
plot_dbr_frequency_spectrum(wavelengths, R, T, lam0) #Plotting the reflectance and transmittance spectra in frequency domain
plot_reflectance_and_N_pairs(n_high, n_low, t1, t2, lam0, N_lp, np.radians(theta), ni, ns, polarization) #Plotting the reflectance spectra for different numbers of pairs
#export_to_csv(wavelengths, R, T) #Exporting the results to a CSV file
