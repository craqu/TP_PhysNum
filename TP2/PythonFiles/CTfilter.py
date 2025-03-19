#!/usr/bin/env python
# -*- coding: utf-8 -*-
# TP reconstruction TDM (CT)
# Prof: Philippe Després
# programme: Dmitri Matenine (dmitri.matenine.1@ulaval.ca)


# libs
import numpy as np
from scipy.fft import fft, ifft, fftfreq
## filtrer le sinogramme
## ligne par ligne

def filterSinogram(sinogram):
    for i in range(sinogram.shape[0]):
        sinogram[i] = filterLine(sinogram[i])

## filter une ligne (projection) via FFT
def filterLine(projection):
    n = len(projection)
    # Transformée de Fourier de la ligne
    projection_fft = fft(projection)

    # Créer le filtre rampe dans le domaine fréquentiel
    freq = fftfreq(n).astype(np.float32)
    ramp_filter = np.abs(freq)

    # Appliquer le filtre rampe
    filtered_fft = projection_fft * ramp_filter

    # Revenir dans le domaine spatial avec la transformée inverse
    filtered_projection = np.real(ifft(filtered_fft))

    return filtered_projection
    # votre code ici
    # un filtre rampe est suffisant
