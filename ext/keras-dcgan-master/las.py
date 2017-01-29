#uses code from https://github.com/MayankSingal/Generative-Adversarial-Networks-Lasagne

from __future__ import print_function
import sys
import numpy as np
import theano
import theano.tensor as T

from scipy.stats import norm

import lasagne
# import matplotlib.pyplot as plt


seed = 42
np.random.seed(seed)

def generator_model(input_var=None, hid_n=5):
    l_in = lasagne.layers.InputLayer(shape=(None, 1), input_var=input_var)

    l_hid1 = lasagne.layers.DenseLayer(
        l_in, num_units=hid_n,
        nonlinearity=lasagne.nonlinearities.softmax,
        W = lasagne.init.GlorotUniform())
    l_out = lasagne.layers.DenseLayer(l_hid1, num_units=1, nonlinearity=lasagne.nonlinearities.linear)

    return l_out

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

def discriminator_model(D=None, input_var=None, hid_n=5):
    if D is not None:
        lays = lasagne.layers.get_all_layers(D)
    else:
        lays  = [lasagne.init.GlorotUniform()]*5
    l_in = lasagne.layers.InputLayer(shape=(None, 1), input_var=input_var)
    hid_n *= 2
    l_hid1 = lasagne.layers.DenseLayer(
        l_in, num_units=hid_n,
        nonlinearity=lasagne.nonlinearities.tanh,
        W = lays[1].W)
    l_hid2 = lasagne.layers.DenseLayer(
        l_hid1, num_units=hid_n,
        nonlinearity=lasagne.nonlinearities.tanh,
        W = lays[2].W)
    l_hid3 = lasagne.layers.DenseLayer(
        l_hid2, num_units=hid_n,
        nonlinearity=lasagne.nonlinearities.tanh,
        W = lays[3].W)
    l_out = lasagne.layers.DenseLayer(l_hid3, num_units=1, nonlinearity=lasagne.nonlinearities.sigmoid,
                                      W=lays[4].W)

    return l_out

def main(num_epochs=500, batch_size=12, batch_n=10):
    hid = 10

    train_iters = 500
    M = 200

    input_var_d = T.col('input_d')
    target_var_d = T.col('output_d')
    input_var_g = T.col('input_g')
    target_var_g = T.col('output_g')

    generator = generator_model(input_var=input_var_g, hid_n=hid)
    prediction_g = lasagne.layers.get_output(generator)

    discriminator = discriminator_model(input_var=input_var_d, hid_n=hid)
    prediction_d = lasagne.layers.get_output(discriminator)

    discriminator_g = discriminator_model(D=discriminator, input_var=generator, hid_n=hid)
    prediction_dg = lasagne.layers.get_output(discriminator_g)

    params_d_g = lasagne.layers.get_all_params(discriminator_g)

    params_d = lasagne.layers.get_all_params(discriminator)

    params_g = lasagne.layers.get_all_params(generator)

    obj_d = T.mean(T.log(prediction_d) + T.log(1 - prediction_dg))
    obj_g = T.mean(T.log(prediction_dg))

    updates_d = lasagne.updates.momentum(1 - obj_d, params_d, learning_rate=0.01)
    updates_g = lasagne.updates.momentum(1 - obj_g, params_d_g, learning_rate=0.01)

    train_d = theano.function([input_var_g, input_var_d], obj_d, updates=updates_d, allow_input_downcast=True)

    train_g = theano.function([input_var_g], obj_g, updates=updates_g, allow_input_downcast=True)

    out_d = theano.function([input_var_d], prediction_d, allow_input_downcast=True)
    out_dg = theano.function([input_var_g], prediction_dg, allow_input_downcast=True)
    out_g = theano.function([input_var_g], prediction_g, allow_input_downcast=True)

    mu, sigma = -1, 0.2
    xs = np.linspace(-1, 1, batch_size)

    k = 1
    batches = (norm.pdf(xs, loc=mu, scale=sigma))*batch_n
    for epoch in range(num_epochs):
        for batch in batches:

            for j in range(k):
                x = np.random.normal(mu, sigma, batch_size)
                x.sort()
                z = np.linspace(-1.0, 1.0, batch_size) + np.random.random(batch_size) * 0.01
                histd = train_d(np.reshape(z, (batch_size, 1)), np.reshape(x, (batch_size, 1)))

            z = np.linspace(-1.0, 1.0, batch_size) + np.random.random(batch_size) * 0.01
            histg = train_g(np.reshape(z, (batch_size, 1)))

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
