"""Usage:
    extractloss.py <in_file> 
    extractloss.py -h | --help

Arguments:
    <in_file>   train.py output text file to extract loss information from

"""


import os
from docopt import docopt
import numpy as np
import matplotlib.pyplot as plt

# Set Latex font for figures
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def lines2num(line):
    n1 = float(line[25:29])
    n2 = float(line[32:36])
    n3 = float(line[50:54])

    return [n1, n2, n3]

def main(fname):
    with open(os.path.abspath(fname), "r") as f:
        lines = f.readlines()

    epochs = int(lines[4][11:])
    i1 = lines.index('Testing for epoch 1:\n')
    i2 = lines.index('Testing for epoch 2:\n')

    period = i2 - i1

    g_trai = np.zeros((epochs, 3))
    g_test = np.zeros((epochs, 3))
    d_trai = np.zeros((epochs, 3))
    d_test = np.zeros((epochs, 3))

    for epoch in range(epochs):
        i = i1+period*epoch
        for j, loss in enumerate((g_trai, g_test, d_trai, d_test)):
            loss[epoch] = lines2num(lines[i+j+3])

    xaxis = np.arange(1,epochs+1)
    plt.plot(xaxis, g_trai[:, 2], label=r'G (train)')
    plt.plot(xaxis, g_test[:, 2], label=r'G (test)')
    plt.plot(xaxis, d_trai[:, 2], label=r'D (train)')
    plt.plot(xaxis, d_test[:, 2], label=r'D (test)')
    plt.xlabel(r'Epoch')
    plt.ylabel(r'Auxillary Loss')
    plt.xlim([0, 50])
    plt.legend()
    plt.grid()
    plt.show()


if __name__ == '__main__':
    arguments = docopt(__doc__, help=True)
    main(arguments['<in_file>'])