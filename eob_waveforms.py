import numpy as np
from pyseobnr.generate_waveform import generate_modes_opt
import subprocess
import os
import matplotlib.pyplot as plt
#import grav_waldo.evaltoolkit as waldo_mm
from scipy.signal import windows
from scipy.interpolate import interp1d

def fourier(t, y):

    # interpolation:
    interp = interp1d(t, y, kind='cubic')
    t = np.linspace(t[0], t[-1], t.size)
    y = interp(t)

    Nt = t.size
    dt = t[1] - t[0]

    # windowing:
    y *= windows.tukey(Nt, 0.065)

    # Fourier Transform:
    yTil = dt*np.fft.fft(y)
    f = np.fft.fftfreq(Nt, dt)

    yTil = np.fft.fftshift(yTil)
    f = np.fft.fftshift(f)
    
    return f, yTil


def inner(f, Sn, y, g=None):

    if type(g) == type(None): num = y*np.conj(y)
    else: num = y*np.conj(g)

    integ = 4*(num/Sn).real

    return np.trapz(integ, x=f, dx=(f[1]-f[0]))


def Mismatch(t, h1, h2):
    
    f, h1 = fourier(t, h1)
    f, h2 = fourier(t, h2)
    
    Sn = np.ones(f.size)

    N1 = inner(f, Sn, h1)
    N2 = inner(f, Sn, h2)
    
    I = [inner(f, Sn, h1, h2*np.exp(-1j*phi)) 
         for phi in np.linspace(0, 2*np.pi, 1000)]
    
    match = max(I)/np.sqrt(N1*N2)

    return 1 - match
def get_nrpy(q,chi1,chi2,M,omegaref):
    parfile_text = f"""
#### seobnrv5_aligned_spin_inspiral BH@H parameter file. NOTE: only commondata CodeParameters appear here ###
###########################
###########################
### Module: nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_codegen_library
chi1 = {chi1}                   # (REAL)
chi2 = {chi2}                  # (REAL)
dt = 2.4627455127717882e-05  # (REAL)
initial_omega = {omegaref}      # (REAL)
mass_ratio = {q}               # (REAL)
total_mass = {M}              # (REAL)
"""
    executable = "./project/seobnrv5_aligned_spin_inspiral/seobnrv5_aligned_spin_inspiral"
    parfile = open("nrpy_parfile.par","w")
    parfile.write(parfile_text)
    parfile.close()
    args = ["nrpy_parfile.par"]
    with open("output.txt","w") as outfile:
        subprocess.run([executable] + args, check=True, stdout=outfile, stderr=subprocess.PIPE)
    data = np.loadtxt("output.txt")
    times = data[:,0]
    strains = data[:,1] + 1j*data[:,2]
    return times, strains

def get_pyseobnr(q,chi1,chi2,M,omegaref):
    times , modes = generate_modes_opt(q,chi1,chi2,omegaref,settings = {"M":M},debug = False)
    strains = modes['2,2']
    return times, strains


def individual_mismatch(input_params):
    q , chi1 , chi2 , M , omegaref = input_params
    times_pyseobnr , strains_pyseobnr = get_pyseobnr(q , chi1 , chi2 , M , omegaref)
    times_nrpy , strains_nrpy = get_nrpy(q , chi1 , chi2 , M , omegaref)
    # interpolate nrpy strains for waldo
    amp_nrpy = np.abs(strains_nrpy)
    phase_nrpy = np.unwrap(np.angle(strains_nrpy))
    amp_interped = np.interp(times_pyseobnr,times_nrpy,amp_nrpy)
    phase_interped = np.interp(times_pyseobnr,times_nrpy,phase_nrpy)
    strains_nrpy_interped = amp_interped*np.exp(1j*phase_interped)
    mismatch_total = Mismatch(times_pyseobnr,strains_pyseobnr,strains_nrpy_interped)
    return mismatch_total

def mismatches(n_samples):
    Qs = np.linspace(1.,8.,n_samples)
    np.random.shuffle(Qs)
    chi1s = np.linspace(-0.9,0.9,n_samples)
    np.random.shuffle(chi1s)
    chi2s = np.linspace(-0.9,0.9,n_samples)
    np.random.shuffle(chi2s)
    M_sol = 4.925491025543575903411922162094833998e-6
    #fiducial frequency given by 1 solar mass binary at 10Hz
    fiducial_omega = M_sol*np.pi*10
    Ms = np.linspace(20,50,n_samples)
    np.random.shuffle(Ms)
    print(Ms)
    omegas = Ms*fiducial_omega
    mm_list = []
    for i in range(n_samples):
        input_params = [Qs[i],chi1s[i],chi2s[i],Ms[i],omegas[i]]
        mm_list.append(individual_mismatch(input_params))
    return mm_list

n_samples = 500
mm_list = mismatches(n_samples)

logbins = np.arange(-5.,0.,0.2)
bins = 10**logbins
plt.hist(mm_list,bins, label = "v5HM vs v5HM_BOB")
plt.xscale('log')
plt.ylabel('# of simulations')
plt.xlabel('Mismatch')
plt.legend()
plt.savefig('histogram_mismatches_EOBOB.png',dpi = 300)