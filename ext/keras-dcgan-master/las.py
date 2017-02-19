#uses code from https://github.com/MayankSingal/Generative-Adversarial-Networks-Lasagne

from __future__ import print_function
import sys
import numpy as np
import theano
import theano.tensor as T

from scipy.stats import norm

import lasagne
import matplotlib.pyplot as plt


seed = 42
np.random.seed(seed)

# def generator_model(input_var=None, hid_n=5):
#     l_in = lasagne.layers.InputLayer(shape=(None, 1), input_var=input_var)
#
#     l_hid1 = lasagne.layers.DenseLayer(
#         l_in, num_units=hid_n,
#         nonlinearity=lasagne.nonlinearities.softmax,
#         W = lasagne.init.GlorotUniform())
#     l_out = lasagne.layers.DenseLayer(l_hid1, num_units=1, nonlinearity=lasagne.nonlinearities.linear)
#
#     return l_out


def generator_model(input_var=None, hid_n=5):
    g_in = lasagne.layers.InputLayer(shape=(None, 1), input_var=input_var)
    g_fw_1 = lasagne.layers.DenseLayer(g_in, num_units=6, nonlinearity=lasagne.nonlinearities.tanh,
                                       W=lasagne.init.GlorotUniform(),
                                       b=lasagne.init.Constant(0.0))
    g_fw_2 = lasagne.layers.DenseLayer(g_fw_1, num_units=5, nonlinearity=lasagne.nonlinearities.tanh,
                                       W=lasagne.init.GlorotUniform(),
                                       b=lasagne.init.Constant(0.0))
    g_out = lasagne.layers.DenseLayer(g_fw_2, num_units=1, nonlinearity=lasagne.nonlinearities.tanh,
                                      W=lasagne.init.GlorotUniform(),
                                      b=lasagne.init.Constant(0.0))


    return g_out

def discriminator_model(input_var=None, hid_n=5):
    d_in = lasagne.layers.InputLayer(shape=(None, 1), input_var=input_var)
    d_fw_1 = lasagne.layers.DenseLayer(d_in, num_units=6, nonlinearity=lasagne.nonlinearities.tanh,
                                       W=lasagne.init.GlorotUniform(),
                                       b=lasagne.init.Constant(0.0))
    d_fw_2 = lasagne.layers.DenseLayer(d_fw_1, num_units=5, nonlinearity=lasagne.nonlinearities.tanh,
                                       W=lasagne.init.GlorotUniform(),
                                       b=lasagne.init.Constant(0.0))
    d_out = lasagne.layers.DenseLayer(d_fw_2, num_units=1, nonlinearity=lasagne.nonlinearities.tanh,
                                      W=lasagne.init.GlorotUniform(),
                                      b=lasagne.init.Constant(0.0))

    return d_out

def discriminator_g_model(D=None, input_var=None, hid_n=5):
    lays = lasagne.layers.get_all_layers(D)
    lays = [(lay.W, lay.b) for lay in lays[1:]]

    dg_fw_1 = lasagne.layers.DenseLayer(input_var, num_units=6, nonlinearity=lasagne.nonlinearities.tanh,
                                        W=lays[0][0],
                                        b=lays[0][1])
    dg_fw_2 = lasagne.layers.DenseLayer(dg_fw_1, num_units=5, nonlinearity=lasagne.nonlinearities.tanh,
                                        W=lays[1][0],
                                        b=lays[1][1])
    dg_out = lasagne.layers.DenseLayer(dg_fw_2, num_units=1, nonlinearity=lasagne.nonlinearities.tanh,
                                       W=lays[2][0],
                                       b=lays[2][1])

    return dg_out

# def discriminator_model(input_var=None, hid_n=5):
#     l_in = lasagne.layers.InputLayer(shape=(None, 1), input_var=input_var)
#     hid_n *= 2
#     l_hid1 = lasagne.layers.DenseLayer(
#         l_in, num_units=hid_n,
#         nonlinearity=lasagne.nonlinearities.tanh,
#         W = lasagne.init.GlorotUniform())
#     l_hid2 = lasagne.layers.DenseLayer(
#         l_hid1, num_units=hid_n,
#         nonlinearity=lasagne.nonlinearities.tanh)
#     l_hid3 = lasagne.layers.DenseLayer(
#         l_hid2, num_units=hid_n,
#         nonlinearity=lasagne.nonlinearities.tanh)
#     l_out = lasagne.layers.DenseLayer(l_hid3, num_units=1, nonlinearity=lasagne.nonlinearities.sigmoid)
#
#     return l_out

# def discriminator_g_model(D=None, input_var=None, hid_n=5):
#
#     lays = lasagne.layers.get_all_layers(D)
#     lays = [lay.W for lay in lays[1:]]
#
#     hid_n *= 2
#     l_hid1 = lasagne.layers.DenseLayer(
#         input_var, num_units=hid_n,
#         nonlinearity=lasagne.nonlinearities.tanh,
#         W = lays[0])
#     l_hid2 = lasagne.layers.DenseLayer(
#         l_hid1, num_units=hid_n,
#         nonlinearity=lasagne.nonlinearities.tanh,
#         W = lays[1])
#     l_hid3 = lasagne.layers.DenseLayer(
#         l_hid2, num_units=hid_n,
#         nonlinearity=lasagne.nonlinearities.tanh,
#         W = lays[2])
#     l_out = lasagne.layers.DenseLayer(l_hid3, num_units=1, nonlinearity=lasagne.nonlinearities.sigmoid,
#                                       W=lays[3])
#
#     return l_out

def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def nz(batch_size, rang):
    # return np.linspace(-rang, rang, batch_size) + np.random.random(batch_size) * 0.01
    a = np.random.normal(0.1, 0.2, batch_size)

    return np.sort(a)

def main(num_epochs=500, batch_size=12, batch_n=10):
    print("main starting")
    hid = 10
    mu, sigma = 0.1, 0.2

    input_var_d = T.col('input_d')
    target_var_d = T.col('output_d')
    input_var_g = T.col('input_g')
    target_var_g = T.col('output_g')

    print("Pre train")

    d_pre = discriminator_model(input_var=input_var_d, hid_n=hid)

    prediction = lasagne.layers.get_output(d_pre)
    loss = loss.mean()
    loss = lasagne.objectives.squared_error(prediction, target_var_d)

    params = lasagne.layers.get_all_params(d_pre, trainable=True)
    updates = lasagne.updates.momentum(loss, params, learning_rate=0.03)

    train = theano.function([input_var_d, target_var_d], loss, updates=updates, allow_input_downcast=True)

    output = theano.function([input_var_d], prediction, allow_input_downcast=True)

    lh = np.zeros(1000)
    for i in range(1000):
        d = (np.random.random(batch_size))
        d = np.reshape(d, (batch_size, 1))
        labels = norm.pdf(d, loc=mu, scale=sigma)
        labels = np.reshape(labels, (batch_size, 1))
        lh[i] = train(d, labels)
        if i % 100 == 0:
            print(lh[i])

    print("Building models")
    generator = generator_model(input_var=input_var_g, hid_n=hid)
    prediction_g = lasagne.layers.get_output(generator)

    discriminator = discriminator_model(input_var=input_var_d, hid_n=hid)
    prediction_d = lasagne.layers.get_output(discriminator)

    discriminator_g = discriminator_g_model(D=discriminator, input_var=generator, hid_n=hid)
    prediction_dg = lasagne.layers.get_output(discriminator_g)

    params_d_g = lasagne.layers.get_all_params(discriminator_g)

    params_d = lasagne.layers.get_all_params(discriminator)

    params_g = lasagne.layers.get_all_params(generator)

    params_pretrained_d_values = lasagne.layers.get_all_param_values(d_pre)

    lasagne.layers.set_all_param_values(discriminator, params_pretrained_d_values)



    obj_d = T.mean(T.log(prediction_d) + T.log(1 - prediction_dg))
    obj_g = T.mean(T.log(prediction_dg))

    updates_d = lasagne.updates.momentum(1 - obj_d, params_d, learning_rate=0.005)
    updates_g = lasagne.updates.momentum(1 - obj_g, params_d_g, learning_rate=0.005)

    print("Compiling")
    train_d = theano.function([input_var_g, input_var_d], obj_d, updates=updates_d, allow_input_downcast=True)

    train_g = theano.function([input_var_g], obj_g, updates=updates_g, allow_input_downcast=True)

    out_d = theano.function([input_var_d], prediction_d, allow_input_downcast=True)
    out_dg = theano.function([input_var_g], prediction_dg, allow_input_downcast=True)
    out_g = theano.function([input_var_g], prediction_g, allow_input_downcast=True)

    xs = np.linspace(-1, 3, 1000)
    y = gaussian(xs, mu, sigma)
    gos = np.random.normal(0, 1.0, 1000)

    k = 1

    print("Begin training")
    for epoch in range(num_epochs):
        for batch in range(batch_n):

            for j in range(k):
                x = np.random.normal(mu, sigma, batch_size)
                x.sort()
                z = nz(batch_size, 1.0)
                histd = train_d(np.reshape(z, (batch_size, 1)), np.reshape(x, (batch_size, 1)))
                # print(histd)

            z = nz(batch_size, 1.0)
            histg = train_g(np.reshape(z, (batch_size, 1)))
            # print(histg)

        plt.clf()
        r = 1000
        # plt.hist(gos, bins=50, histtype='step')
        plt.plot(xs, y)
        z = nz(r, 1.0).reshape((r, 1))
        histc, edges = np.histogram(out_g(z), bins=10)
        # plt.hist(out_g(z), bins=50, histtype='step') #only pass to g in batch size?
        print(sum(histc/float(r)))
        plt.plot((edges[-1]+edges[1:])/2, histc/float(r))
        # plt.show()
        plt.savefig('las_{}.png'.format(epoch))
        print(epoch)

if __name__ == '__main__':
    kwargs = {}
    if len(sys.argv) > 1:
        kwargs['num_epochs'] = int(sys.argv[1])
    if len(sys.argv) > 2:
        kwargs['batch_size'] = int(sys.argv[2])
    if len(sys.argv) > 3:
        kwargs['batch_n'] = int(sys.argv[3])
    main(**kwargs)
