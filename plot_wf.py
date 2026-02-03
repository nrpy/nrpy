import matplotlib.pyplot as plt
import numpy as np

filename = "project/sebobv2/waveform.txt"

try:
    # Load data: Time, Real(h), Imag(h)
    data = np.loadtxt(filename)
except ValueError:
    print(f"Error loading {filename}. Open it and check if there are text lines mixed with numbers.")
    exit(1)

t = data[:, 0]
h_plus = data[:, 1]  # Real part
h_cross = data[:, 2] # Imaginary part

# Calculate Amplitude: |h| = sqrt(h+^2 + hx^2)
amp = np.sqrt(h_plus**2 + h_cross**2)

# Calculate Frequency: omega = d(phase)/dt
phase = np.unwrap(np.angle(h_plus - 1j*h_cross)) # -1j convention usually
# Simple finite difference for frequency
dt = np.diff(t)
omega = np.diff(phase) / dt

plt.figure(figsize=(12, 12))

# 1. Waveform (h+)
plt.subplot(4, 1, 1)
plt.plot(t, h_plus, label=r'$h_+$', linewidth=1)
plt.ylabel("Strain")
plt.title("SEBOBv2 Waveform")
plt.legend(loc='upper left')
plt.grid(alpha=0.3)

# 2. Amplitude (Zoomed in on Merger)
plt.subplot(4, 1, 2)
plt.plot(t, amp, color='orange', label=r'$|h|$')
plt.ylabel("Amplitude")
plt.legend(loc='upper left')
plt.grid(alpha=0.3)

# Auto-zoom to the merger part (peak of amplitude)
peak_idx = np.argmax(amp)
t_peak = t[peak_idx]
plt.xlim(t_peak - 100, t_peak + 100)
plt.axvline(t_peak, color='k', linestyle='--', alpha=0.5, label='Peak')

# 3. Frequency
plt.subplot(4, 1, 3)
plt.scatter(t[:-1], np.abs(omega), color='green', label=r'$\omega$',s=5.)
plt.plot(t[:-1], np.abs(omega), color='green', label=r'$\omega$')
plt.ylabel("Frequency")
plt.xlim(t_peak - 100, t_peak + 100)
plt.legend(loc='upper left')
plt.grid(alpha=0.3)

# 4. Phase
plt.subplot(4, 1, 4)
plt.plot(t, phase, color='purple', label=r'$\phi$')
plt.ylabel("Phase")
plt.xlabel("Time (M)")
plt.xlim(t_peak - 100, t_peak + 100)
plt.legend(loc='upper left')
plt.grid(alpha=0.3)



plt.tight_layout()
plt.show()