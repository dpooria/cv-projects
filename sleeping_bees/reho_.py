import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from scipy.ndimage import convolve


def preproc(inp):
  inp = np.array(inp)

  weights = np.array([[0, 0, 1, 0, 0],
                      [0, 2, 4, 2, 0],
                      [1, 4, 8, 4, 1],
                      [0, 2, 4, 2, 0],
                      [0, 0, 1, 0, 0]],
                     dtype=np.float32)

  weights = weights / np.sum(weights[:])

  inpsm = np.zeros(inp.shape)

  for i in range(len(inp)):
    inpsm[i] = convolve(inp[i], weights, mode='constant')

  # Apply [ t - mean(t) / mean (t) ] to each pixel

  inpen = -(inpsm - np.mean(inpsm, axis=0)) / np.mean(inpsm, axis=0)

  # -0.5 < threshold < 0.5
  inpen[inpen >= 0.5] = 0
  inpen[inpen <= -0.5] = 0

  return inpen


def correlation_pearson(data: np.ndarray, neighs=5):

  IMSIDE = data.shape[1]
  ncorr = np.zeros((IMSIDE, IMSIDE), dtype=np.float32)

  for i in range(neighs, IMSIDE-neighs):
      for j in range(neighs, IMSIDE-neighs):
          ncorr[i, j] = pearsonr(data[:, i, j], np.mean(
              data[:, i-neighs:i+neighs, j-neighs:j+neighs], axis=(1, 2)))[0]

  return ncorr


def plot_corr(ncorr, file_name="correlation.jpg"):
  fig = plt.imshow(ncorr, cmap='RdBu_r', vmin=-1, vmax=1)
  cbar = plt.colorbar(fig)
  cbar.set_label('Correlation')
  plt.savefig(file_name)
  plt.show()
