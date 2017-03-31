import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LogNorm

from readjet import delta_phi


def plot_diff(phiran, etaran, image_array, output_name=None, extr=None, title='', cmap=cm.PRGn_r):
    """
    Adapted from Luke de Oliveira, & Michela Paganini. (2017). lukedeo/adversarial-jets: Initial Release [Data set]. Zenodo. http://doi.org/10.5281/zenodo.400708

        Function to help you visualize the difference between two sets of jet images on a linear scale
        Args:
        -----
           content : numpy array of dimensions 25x25, first arg to imshow, content of the image
                     e.g.: generated_images.mean(axis=0) - real_images.mean(axis=0) --> difference between avg generated and avg Pythia image
                           etc...
           output_name : string, name of the output file where the plot will be saved. Note: it will be located in ../plots/
           extr : (default = None) float, magnitude of the upper and lower bounds of the pixel intensity scale before saturation (symmetric around 0)
           title : (default = '') string, title of the plot, to be displayed on top of the image
           cmap : (default = matplotlib.cm.PRGn_r) matplotlib colormap, ideally white in the middle
        Outputs:
        --------
           no function returns
           saves file in ../plots/output_name
        """
    fig, ax = plt.subplots(figsize=(6, 6))
    extent = phiran + etaran
    if extr is None:
        extr = max(abs(image_array.min()), abs(image_array.max()))
    im = ax.imshow(image_array,
                   interpolation='nearest', norm=Normalize(vmin=-extr, vmax=+extr), extent=extent,
                   cmap=cmap)
    plt.colorbar(im, fraction=0.05, pad=0.05)
    plt.xlabel(r'[Transformed] Pseudorapidity $(\eta)$')
    plt.ylabel(r'[Transformed] Azimuthal Angle $(\phi)$')
    plt.title(title)
    if output_name is None:
        plt.show()
    else:
        plt.savefig(output_name)


def discrete_pT(image_array, phivals):
    #calculate the image pT
    cos_term = image_array * np.cos(phivals)
    sin_term = image_array * np.sin(phivals)
    return np.sqrt(np.sum(cos_term)**2 + np.sum(sin_term)**2)


def discrete_m(image_array, etavals, phivals):
    #calculate jet mass from image
    sinh_term = image_array * np.sinh(etavals)
    return np.sqrt(np.sum(image_array) ** 2 - discrete_pT(image_array, phivals) - np.sum(sinh_term) ** 2)


def tau1(image_array, etavals, phivals, dim):
    """"
    Calculates the normalized tau1 from a pixelated jet image
    Args:
    -----
        jet_image: numpy ndarray of dim (1, 25, 25)
    Returns:
    --------
        float, normalized jet tau1
    """
    # find coordinate of most energetic pixel, then use formula to compute tau1
    eta = etavals
    phi = phivals
    tau1_axis_eta = eta.ravel()[np.argmax(image_array)]
    tau1_axis_phi = phi.ravel()[np.argmax(image_array)]
    t1 = np.sum(image_array*
                  np.sqrt(np.square(tau1_axis_eta - eta) + np.square(
                      [delta_phi(tau1_axis_phi, p) for p in phi.ravel()]).reshape(dim))
                  )
    return t1 / np.sum(image_array)  # normalize by the total intensity


def tau2(image_array, etavals, phivals):
    """
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
    """
    jet_image = image_array
    eta= etavals
    phi= phivals
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

    return np.sum(jet_image * grid) / np.sum(jet_image)  # normalize by the total intensity


def tau21(image_array, etavals, phivals, dim):
    """
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
    """
    t1 = tau1(image_array, etavals, phivals, dim)
    if t1 <= 0:
        return 0
    else:
        t2 = tau2(image_array, etavals, phivals)
        return t2 / t1


def plot_jet(image_array, etaran=(-1.25,1.25), phiran=(-1.25,1.25), output_name=None, vmin=1e-6, vmax=300, title=''):
    """
        Adapted from Luke de Oliveira, & Michela Paganini. (2017). lukedeo/adversarial-jets: Initial Release [Data set]. Zenodo. http://doi.org/10.5281/zenodo.400708
        Function to help you visualize a jet image on a log scale
        Args:
        -----
           content : numpy array of dimensions 25x25, first arg to imshow, content of the image
                     e.g.: generated_images.mean(axis=0) --> the average generated image
                           real_images.mean(axis=0) --> the average Pythia image
                           generated_images[aux_out == 1].mean(axis=0) --> the average generated image labeled as real by the discriminator
                           etc...
           output_name : string, name of the output file where the plot will be saved. Note: it will be located in ../plots/
           vmin : (default = 1e-6) float, lower bound of the pixel intensity scale before saturation
           vmax : (default = 300) float, upper bound of the pixel intensity scale before saturation
           title : (default = '') string, title of the plot, to be displayed on top of the image
        Outputs:
        --------
           no function returns
           saves file in ../plots/output_name
        """
    fig, ax = plt.subplots(figsize=(7, 6))
    extent = phiran + etaran
    im = ax.imshow(image_array, interpolation='nearest', norm=LogNorm(vmin=vmin, vmax=vmax), extent=extent)
    cbar = plt.colorbar(im, fraction=0.05, pad=0.05)
    cbar.set_label(r'Pixel $p_T$ (GeV)', y=0.85)
    plt.xlabel(r'[Transformed] Pseudorapidity $(\eta)$')
    plt.ylabel(r'[Transformed] Azimuthal Angle $(\phi)$')
    plt.title(title)
    if output_name is None:
        plt.show()
    else:
        plt.savefig(output_name)


def average_image(images_array):
    #avereage image from array of images
    return np.mean(images_array, axis=0)



def maxim(images):
    #return the JetImage object with highest discrete pT in list
    return max(images, key=lambda item: discrete_pT(item.image_array, item.phivals))


def image_set(images):
    #convert list of JetImages to array of image arrays shape (N_images, Nphi, Neta)
    return np.array([image.image_array for image in images])