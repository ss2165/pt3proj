from keras.models import Sequential
from keras.layers import Dense
from keras.layers.core import Activation
from keras.optimizers import SGD
import numpy as np
import matplotlib.pyplot as plt

seed = 42
np.random.seed(seed)
import theano.tensor as T
from theano import function

def discriminator_model(input, hidden_size):
    model = Sequential()
    model.add(Dense(hidden_size*2, input_dim=input))
    model.add(Activation('tanh'))
    model.add(Dense(hidden_size*2))
    model.add(Activation('tanh'))
    # model.add(Dense(hidden_size*2))
    # model.add(Activation('tanh'))
    model.add(Dense(1))
    model.add(Activation('linear'))
    return model

def generator_model(input, hidden_size):
    model = Sequential()
    model.add(Dense(hidden_size, input_dim=input))
    model.add(Activation('softplus'))
    model.add(Dense(1))
    model.add(Activation('linear'))
    return model

def generator_containing_discriminator(generator, discriminator):
    model = Sequential()
    model.add(generator)
    discriminator.trainable = False
    model.add(discriminator)
    return model

def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

# def make_log():
#     yt = T.matrix('yt')
#     yp = T.matrix('yp')
#     z = T.mean(T.log)
#
# def noise_batch(bsize, no):

def train():

    X_train = np.random.normal(0, 1, 10000)
    X_train = X_train.reshape(X_train.shape[0],-1)
    # Y_train = np.greater(X_train, 0.25).astype(int)
    # data = np.stack((X_train, Y_train), axis=1)
    # np.random.shuffle(data)
    x1 = np.linspace(-2, 2, 1000)
    yd = gaussian(x1, 0, 1)



    discriminator = discriminator_model(1, 5)
    generator = generator_model(1, 5)
    discriminator_on_generator = \
        generator_containing_discriminator(generator, discriminator)

    d_optim = SGD(lr=0.001, momentum=0.9, nesterov=True)
    g_optim = SGD(lr=0.001, momentum=0.9, nesterov=True)
    generator.compile(loss='binary_crossentropy', optimizer="SGD")
    discriminator_on_generator.compile(loss='binary_crossentropy', optimizer=g_optim)
    discriminator.trainable = True
    discriminator.compile(loss='binary_crossentropy', optimizer=d_optim)
    BATCH_SIZE = 10
    noise = np.zeros((BATCH_SIZE, 1))

    for epoch in range(100):
        print("Epoch is", epoch)
        batch_no = int(X_train.shape[0]/BATCH_SIZE)
        print("Number of batches", batch_no)
        for index in range(batch_no):
            # for i in range(BATCH_SIZE):
            noise = np.linspace(-1, 1, BATCH_SIZE) + np.random.random(BATCH_SIZE) * 0.01
            train_batch = X_train[index * BATCH_SIZE:(index + 1) * BATCH_SIZE]

            generated = generator.predict(noise, verbose=0)
            X = np.concatenate((train_batch, generated))
            y = [1] * BATCH_SIZE + [0] * BATCH_SIZE
            d_loss = discriminator.train_on_batch(X, y)
            print("batch %d d_loss : %f" % (index, d_loss))

            # for i in range(BATCH_SIZE):
            noise = np.linspace(-1, 1, BATCH_SIZE) + np.random.random(BATCH_SIZE) * 0.01
            discriminator.trainable = False
            g_loss = discriminator_on_generator.train_on_batch(
                noise, [1] * BATCH_SIZE)
            discriminator.trainable = True
            print("batch %d g_loss : %f" % (index, g_loss))

        nz = np.random.uniform(-1, 1, 100)
        yg = generator.predict(nz, verbose=0)
        plt.clf()
        plt.plot(x1, yd)
        # plt.plot(nz, yg, 'x')
        plt.hist(yg, bins=50, histtype='step')
        plt.savefig("class_{}.png".format(epoch))




    # x = np.random.uniform(-1, 1, 100)
    # y = discriminator.predict(x, verbose=0)
    # plt.clf()
    # plt.plot(x, y, 'x')
    # plt.plot([0.25, 0.25], [0, 1.0])
    # plt.savefig('class.png')

train()


