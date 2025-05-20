# Signal Processing Projects using MATLAB 

This repository contains a collection of MATLAB-based simulations and experiments on analog modulation techniques and signal filtering. It includes practical implementations of DSB-LC, DSB-SC, SSB-SC, and FM, along with noise analysis and the threshold effect.

---

## Project Contents:

### Bandlimited Audio Filtering
- Load and interpolate real speech signals: F, B, D, S, M, N
- Apply low-pass filters with cutoff frequencies ranging from 300 Hz to 3400 Hz
- Observe how intelligibility is affected by the bandwidth
- Plot signals in time and frequency domains
- Playback results and analyze the clarity

###  Amplitude Modulation
####  DSB-LC (Double Sideband with Large Carrier)
- Modulate with carrier frequency = 48 kHz, modulation index μ = 0.8
- Envelope and square-law detection
- Time/frequency analysis and energy evaluation

#### DSB-SC (Suppressed Carrier)
- Multiply the message by the carrier for modulation
- Coherent demodulation using the same carrier
- Analyze sensitivity to frequency mismatch

####  (Single Sideband Suppressed Carrier)
- Use the Hilbert transform to generate an analytic signal
- Suppress one sideband by modulation with a complex exponential
- Recover message using coherent detection
- Explore the frequency offset effect on intelligibility

###  Frequency Modulation (FM)
- FM modulation of real audio signals and 3 kHz tone
- Vary modulation index β = [0.1, 0.5, 1, 3, 5, 10]
- Observe bandwidth expansion and tone clarity
- Use direct demodulation method: differentiation + envelope detection

### Noise and Threshold Effects
- Add Gaussian noise at powers 0.0001, 0.01, 0.1
- Measure SNR before and after demodulation
- Explore FM threshold effect by increasing β
- Compare noise resilience of different modulation schemes

---

## Technical Notes & Insights

### Sampling and Reconstruction:
- **Sampling rate**: Number of samples per second. Must satisfy **Nyquist theorem**: `f_s > 2 × f_max`.
- **Aliasing** occurs if `f_s` is too low, causing different signals to appear identical at sampled points.
- **44.1 kHz** is standard for audio since it’s > 2 × 20 kHz (human hearing range).
- Sampling in time = multiplication by a delta train → convolution in frequency, creating copies of the spectrum.
- If the sampling frequency is not high enough, these copies **overlap**, causing **irreversible loss of information**.

### Recording Notes:
- **16-bit audio**: Bit depth defining amplitude resolution. MATLAB stores audio in `double`, normalized from `-1` to `+1`.
- **Clipping**: Happens when amplitude exceeds representable limits, leading to harsh distortion
