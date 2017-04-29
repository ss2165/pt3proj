"""
Based on plots.ipynb
"""

import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LogNorm, Normalize

import h5py

def discrete_mass(jet_image):
    '''
    Calculates the jet mass from a pixelated jet image
    Args:
    -----
        jet_image: numpy ndarray of dim (1, 25, 25)
    Returns:
    --------
        M: float, jet mass
    '''
    Px = np.sum(jet_image * np.cos(phi), axis=(1, 2))
    Py = np.sum(jet_image * np.sin(phi), axis=(1, 2))

    Pz = np.sum(jet_image * np.sinh(eta), axis=(1, 2))
    E = np.sum(jet_image * np.cosh(eta), axis=(1, 2))

    PT2 = np.square(Px) + np.square(Py)
    M2 = np.square(E) - (PT2 + np.square(Pz))
    M = np.sqrt(M2)
    return M

def discrete_pt(jet_image):
    '''
    Calculates the jet transverse momentum from a pixelated jet image
    Args:
    -----
        jet_image: numpy ndarray of dim (1, 25, 25)
    Returns:
    --------
        float, jet transverse momentum
    '''
    Px = np.sum(jet_image * np.cos(phi), axis=(1, 2))
    Py = np.sum(jet_image * np.sin(phi), axis=(1, 2))
    return np.sqrt(np.square(Px) + np.square(Py))

def dphi(phi1, phi2):
    '''
    Calculates the difference between two angles avoiding |phi1 - phi2| > 180 degrees
    '''
    import math
    return math.acos(math.cos(abs(phi1 - phi2)))

def _tau1(jet_image):
    '''
    Calculates the normalized tau1 from a pixelated jet image
    Args:
    -----
        jet_image: numpy ndarray of dim (1, 25, 25)
    Returns:
    --------
        float, normalized jet tau1
    '''
    # find coordinate of most energetic pixel, then use formula to compute tau1
    tau1_axis_eta = eta.ravel()[np.argmax(jet_image)]
    tau1_axis_phi = phi.ravel()[np.argmax(jet_image)]
    tau1 = np.sum(jet_image *
            np.sqrt(np.square(tau1_axis_eta - eta) + np.square([dphi(tau1_axis_phi, p) for p in phi.ravel()]).reshape(25, 25))
                 )
    return tau1 / np.sum(jet_image) # normalize by the total intensity


def _tau2(jet_image):
    '''
    Calculates the normalized tau2 from a pixelated jet image
    Args:
    -----
        jet_image: numpy ndarray of dim (1, 25, 25)
    Returns:
    --------
        float, normalized jet tau2
    Notes:
    ------
        slow implementation
    '''
    proto = np.array(zip(jet_image[jet_image != 0],
                         eta[jet_image != 0],
                         phi[jet_image != 0]))

    while len(proto) > 2:
        candidates = [
            (
                (i, j),
                (min(pt1, pt2) ** 2) * ((eta1 - eta2) ** 2 + (phi1 - phi2) ** 2)
            )
            for i, (pt1, eta1, phi1) in enumerate(proto)
            for j, (pt2, eta2, phi2) in enumerate(proto)
            if j > i
        ]

        index, value = zip(*candidates)
        pix1, pix2 = index[np.argmin(value)]
        if pix1 > pix2:
            # swap
            pix1, pix2 = pix2, pix1

        (pt1, eta1, phi1) = proto[pix1]
        (pt2, eta2, phi2) = proto[pix2]

        e1 = pt1 / np.cosh(eta1)
        e2 = pt2 / np.cosh(eta2)
        choice = e1 > e2

        eta_add = (eta1 if choice else eta2)
        phi_add = (phi1 if choice else phi2)
        pt_add = (e1 + e2) * np.cosh(eta_add)

        proto[pix1] = (pt_add, eta_add, phi_add)

        proto = np.delete(proto, pix2, axis=0).tolist()

    (_, eta1, phi1), (_, eta2, phi2) = proto
    np.sqrt(np.square(eta - eta1) + np.square(phi - phi1))

    grid = np.array([
        np.sqrt(np.square(eta - eta1) + np.square(phi - phi1)),
        np.sqrt(np.square(eta - eta2) + np.square(phi - phi2))
    ]).min(axis=0)

    return np.sum(jet_image * grid) / np.sum(jet_image) # normalize by the total intensity

def tau21(jet_image):
    '''
    Calculates the tau21 from a pixelated jet image using the functions above
    Args:
    -----
        jet_image: numpy ndarray of dim (1, 25, 25)
    Returns:
    --------
        float, jet tau21
    Notes:
    ------
        slow implementation
    '''
    tau1 = _tau1(jet_image)
    if tau1 <= 0:
        return 0
    else:
        tau2 = _tau2(jet_image)
        return tau2 / tau1

training_file = 'data/prepared.hdf'
with h5py.File(training_file, 'r') as f:
    real_images = f['image'][:]
    real_labels = f['signal'][:]
    real_images[real_images < 1e-3] = 0.0  # everything below 10^-3 is unphysical and due to instabilities in the rotation

n_jets = real_images.shape[0]
latent_space = 200 # size of the vector z

gen_weights = 'models/params_generator_epoch_049.hdf5'
disc_weights = 'models/params_discriminator_epoch_049.hdf5'

from models.networks.lagan import generator as build_generator
g = build_generator(latent_space, return_intermediate=False)
g.load_weights(gen_weights)

noise = np.random.normal(0, 1, (n_jets, latent_space))
sampled_labels = np.random.randint(0, 2, n_jets)
generated_images = g.predict(
    [noise, sampled_labels.reshape(-1, 1)], verbose=False, batch_size=64)
from jetimage.analysis import average_image, plot_jet

generated_images *= 100
generated_images = np.squeeze(generated_images)
# generated_images = generated_images.reshape(generated_images.shape[:-1])
signal_images = generated_images[sampled_labels==1]
noise_images = generated_images[sampled_labels==0]
av_gen = average_image(signal_images)
av_noise = average_image(noise_images)

# for i in range(10):
#     plot_jet(noise_images[i])
# plot_jet(av_gen)
# plot_jet(av_noise)

outdir = 'plots'
grid = 0.5 * (np.linspace(-1.25, 1.25, 26)[:-1] + np.linspace(-1.25, 1.25, 26)[1:])
eta = np.tile(grid, (25, 1))
phi = np.tile(grid[::-1].reshape(-1, 1), (1, 25))
##PIXEL INTENSITY
# fig, ax = plt.subplots(figsize=(6, 6))
#
# _, bins, _ = plt.hist(real_images.ravel(),
#            bins=np.linspace(0, 300,50), histtype='step', label='Pythia', color='purple')
# _ = plt.hist(generated_images.ravel(),
#              bins=bins, histtype='step', label='GAN', color='green')
#
# plt.xlabel('Pixel Intensity')
# plt.ylabel('Number of Pixels')
# plt.yscale('log')
# plt.legend(loc='upper right')
#
# plt.savefig(os.path.join(outdir, 'pixel_intensity.pdf'))

##MASS
fig, ax = plt.subplots(figsize=(6, 6))
bins = np.linspace(40, 120, 50)
_ = plt.hist(discrete_mass(generated_images[sampled_labels == 1]),
             bins=bins, histtype='step', label=r"generated ($W' \rightarrow WZ$)", normed=True, color='red')
_ = plt.hist(discrete_mass(real_images[real_labels == 1]),
             bins=bins, histtype='step', label=r"Pythia ($W' \rightarrow WZ$)", normed=True, color='red', linestyle='dashed')

_ = plt.hist(discrete_mass(generated_images[sampled_labels == 0]),
             bins=bins, histtype='step', label=r'generated (QCD dijets)', normed=True, color='blue')
_ = plt.hist(discrete_mass(real_images[real_labels == 0]),
             bins=bins, histtype='step', label=r'Pythia (QCD dijets)', normed=True, color='blue', linestyle='dashed')

plt.xlabel(r'Discretized $m$ of Jet Image')
plt.ylabel(r'Units normalized to unit area')
plt.legend()
plt.ylim(0, 0.11)
plt.savefig(os.path.join(outdir, 'mass.pdf'))

##PT
fig, ax = plt.subplots(figsize=(6, 6))
bins = np.linspace(200, 340, 50)
_ = plt.hist(discrete_pt(generated_images[sampled_labels == 1]),
             bins=bins, histtype='step', label=r"generated ($W' \rightarrow WZ$)", normed=True, color='red')
_ = plt.hist(discrete_pt(real_images[real_labels == 1]),
             bins=bins, histtype='step', label=r"Pythia ($W' \rightarrow WZ$)", normed=True, color='red', linestyle='dashed')

_ = plt.hist(discrete_pt(generated_images[sampled_labels == 0]),
             bins=bins, histtype='step', label=r'generated (QCD dijets)', normed=True, color='blue')
_ = plt.hist(discrete_pt(real_images[real_labels == 0]),
             bins=bins, histtype='step', label=r'Pythia (QCD dijets)', normed=True, color='blue', linestyle='dashed')
plt.xlabel(r'Discretized $p_T$ of Jet Image')
plt.ylabel(r'Units normalized to unit area')
plt.legend()
plt.ylim(0, 0.045)
plt.savefig(os.path.join(outdir, 'pt.pdf'))