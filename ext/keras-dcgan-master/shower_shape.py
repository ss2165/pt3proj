from keras.models import Sequential
from keras.layers import Dense
from keras.layers.core import Activation
from keras.optimizers import SGD
import numpy as np
import matplotlib.pyplot as plt
import argparse

def real_rad(size, width):
    return np.random.exponential(width, size)

def real_angle(size):
    return np.random.uniform(0, 2*np.pi, size)

def real_logE(size, width, mean):
    return np.random.normal(mean, width, size)

def generator_model():
    model = Sequential()
    model.add(Dense(12, input_dim=3))
    model.add(Activation('relu'))
    model.add(Dense(1))
    model.add(Activation('tanh'))
    return model

def discriminator_model():
    model = Sequential()
    # model.add(Dense(12, input_dim=1))
    # model.add(Activation('tanh'))
    model.add(Dense(1, input_dim=1))
    model.add(Activation('sigmoid'))
    return model

def generator_containing_discriminator(generator, discriminator):
    model = Sequential()
    model.add(generator)
    discriminator.trainable = False
    model.add(discriminator)
    return model

def train(BATCH_SIZE):
    # (X_train, y_train), (X_test, y_test) = mnist.load_data()
    # X_train = (X_train.astype(np.float32) - 127.5) / 127.5
    t_size = 10000
    r_width = 100
    E_width = 3
    E_mean = -15
    # X_train = np.stack((real_rad(t_size, r_width), real_angle(t_size), real_logE(t_size, E_width, E_mean)), axis=-1)
    X_train = real_logE(t_size, 1, 0)
    X_train = X_train.reshape(X_train.shape[0],-1)
    discriminator = discriminator_model()
    generator = generator_model()
    discriminator_on_generator = \
        generator_containing_discriminator(generator, discriminator)

    d_optim = SGD(lr=0.0005, momentum=0.9, nesterov=True)
    g_optim = SGD(lr=0.0005, momentum=0.9, nesterov=True)

    generator.compile(loss='binary_crossentropy', optimizer="SGD")
    discriminator_on_generator.compile(
        loss='binary_crossentropy', optimizer=g_optim)
    discriminator.trainable = True
    discriminator.compile(loss='binary_crossentropy', optimizer=d_optim)
    noise = np.zeros((BATCH_SIZE, 3))
    for epoch in range(200):
        print("Epoch is", epoch)
        batch_no = int(X_train.shape[0]/BATCH_SIZE)
        print("Number of batches", batch_no)
        for index in range(batch_no):
            for i in range(BATCH_SIZE):
                noise[i, :] = np.random.uniform(-1, 1, 3)
            train_batch = X_train[index * BATCH_SIZE:(index + 1) * BATCH_SIZE]

            generated = generator.predict(noise, verbose=0)

            if index % int(batch_no*BATCH_SIZE) == 0:
                save_hist(train_batch, generated, epoch, index)

            X = np.concatenate((train_batch, generated))
            y = [1] * BATCH_SIZE + [0] * BATCH_SIZE
            d_loss = discriminator.train_on_batch(X, y)
            print("batch %d d_loss : %f" % (index, d_loss))
            for i in range(BATCH_SIZE):
                noise[i, :] = np.random.uniform(-1, 1, 3)
            discriminator.trainable = False
            g_loss = discriminator_on_generator.train_on_batch(
                noise, [1] * BATCH_SIZE)
            discriminator.trainable = True
            print("batch %d g_loss : %f" % (index, g_loss))
            if index % 10 == 9:
                generator.save_weights('generator', True)
                discriminator.save_weights('discriminator', True)


def save_hist(dist1, dist2, epoch, index):
    plt.clf()
    fig = plt.figure()
    plt.hist(dist1, 50, histtype="step")
    plt.hist(dist2, 50, histtype="step")
    plt.savefig('r_dist_{}_{}.png'.format(epoch, index))

def combine_images(generated_images):
    generated_images=generated_images.reshape((generated_images.shape[0], 28, 28))
    num = generated_images.shape[0]
    width = int(np.sqrt(num))
    height = int(np.ceil(float(num)/width))
    shape = generated_images.shape[1:]
    image = np.zeros((height*shape[0], width*shape[1]),
                     dtype=generated_images.dtype)
    for index, img in enumerate(generated_images):
        i = int(index/width)
        j = index % width
        image[i*shape[0]:(i+1)*shape[0], j*shape[1]:(j+1)*shape[1]] = \
            img
    return image

def generate(BATCH_SIZE):
    generator = generator_model()
    generator.compile(loss='binary_crossentropy', optimizer="SGD")
    generator.load_weights('generator')
    noise = np.zeros((BATCH_SIZE, 3))
    for i in range(BATCH_SIZE):
        noise[i, :] = np.random.uniform(-1, 1, 3)
    generated = generator.predict(noise, verbose=1)
    save_hist(generated, generated, "G", "G")

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", type=str)
    parser.add_argument("--batch_size", type=int, default=128)
    parser.add_argument("--nice", dest="nice", action="store_true")
    parser.set_defaults(nice=False)
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = get_args()
    if args.mode == "train":
        train(BATCH_SIZE=args.batch_size)
    elif args.mode == "generate":
        generate(BATCH_SIZE=args.batch_size)
