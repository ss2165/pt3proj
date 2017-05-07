from plots.plot_tools import *
import h5py
import os


def total_dist(pythia, delphes, title='Title'):
    # fig, ax = plt.subplots(figsize=(6, 6))
    bins = np.linspace(200, 1000, 50)
    _ = plt.hist(delphes,
                 bins=bins, histtype='step', label=r"P+D", normed=True)
    _ = plt.hist(pythia,
                 bins=bins, histtype='step', label=r"P", normed=True)

    plt.xlabel(r'Discretized $p_T$ of Jet Image')
    plt.ylabel(r'Units normalized to unit area')
    plt.legend()
    plt.ylim(0, 0.05)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    plt.title(title)

d_images, d_labels = load_images(training_file)
# generated_images, sampled_labels = load_images(generated_file)

with h5py.File(os.path.abspath('data/pythia_24k.hdf5'), 'r') as f:
    p_images = f['image'][:]
    p_labels = f['signal'][:]
    p_images[p_images < 1e-3] = 0.0


av_d = average_image(d_images)
av_p = average_image(p_images)

# total_p = np.sum(p_images.reshape((p_images.shape[0], 625)), axis=1)
# total_d = np.sum(d_images.reshape((d_images.shape[0], 625)), axis=1)
# total_dist(total_p, total_d)

# f, ax  = plt.subplots(1, 3, figsize=(10, 4))

# plt.subplot(131)
# plot_jet(av_p, title='Pythia only')
# plt.subplot(132)
# plot_diff_jet_image(np.clip(av_p - av_d, -10, 10), title='Difference')

# plt.subplot(133)
# plot_jet(av_d, title='Pythia+Delphes')

plt.semilogy(phi[:, 0], av_p[:, 12])
plt.semilogy(phi[:, 0], av_d[:, 12])
plt.show()


# with h5py.File(os.path.abspath('data/pythia_24k.hdf5'), 'w') as f:
#     dset = f.create_dataset('image', data=p_images)
#     sigs = f.create_dataset('signal', data=p_labels)