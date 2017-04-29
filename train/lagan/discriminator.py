#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
file: discriminators.py
description: discrimination submodel for [arXiv/1701.05927]
author: Luke de Oliveira (lukedeoliveira@lbl.gov)
"""

import keras.backend as K
from keras.layers import (Input, Dense, Reshape, Flatten, Lambda, merge,
                          Dropout, BatchNormalization)
from keras.layers.advanced_activations import LeakyReLU
from keras.layers.convolutional import (UpSampling2D, Conv2D, ZeroPadding2D,
                                        AveragePooling2D)
from keras.layers.local import LocallyConnected2D
from keras.models import Model

from .ops import minibatch_discriminator, minibatch_output_shape, Dense3D

K.set_image_dim_ordering('tf')


def discriminator():

    image = Input(shape=(25, 25, 1))

    # block 1: normal 5x5 conv,
    # *NO* batchnorm (recommendation from [arXiv/1511.06434])
    x = Conv2D(32, 5, 5, border_mode='same')(image)
    x = LeakyReLU()(x)
    x = Dropout(0.2)(x)

    # block 2: 'same' bordered 5x5 locally connected block with batchnorm and
    # 2x2 subsampling
    x = ZeroPadding2D((2, 2))(x)
    x = LocallyConnected2D(8, 5, 5, border_mode='valid', subsample=(2, 2))(x)
    x = LeakyReLU()(x)
    x = BatchNormalization()(x)
    x = Dropout(0.2)(x)

    # block 2: 'same' bordered 5x5 locally connected block with batchnorm
    x = ZeroPadding2D((2, 2))(x)
    x = LocallyConnected2D(8, 5, 5, border_mode='valid')(x)
    x = LeakyReLU()(x)
    x = BatchNormalization()(x)
    x = Dropout(0.2)(x)

    # block 3: 'same' bordered 3x3 locally connected block with batchnorm and
    # 2x2 subsampling
    x = ZeroPadding2D((1, 1))(x)
    x = LocallyConnected2D(8, 3, 3, border_mode='valid', subsample=(2, 2))(x)
    x = LeakyReLU()(x)
    x = BatchNormalization()(x)
    x = Dropout(0.2)(x)

    x = AveragePooling2D((2, 2))(x)
    h = Flatten()(x)

    dnn = Model(image, h)

    image = Input(shape=(25, 25, 1))

    dnn_out = dnn(image)

    # nb of features to obtain
    nb_features = 20

    # dim of kernel space
    vspace_dim = 10

    # creates the kernel space for the minibatch discrimination
    K_x = Dense3D(nb_features, vspace_dim)(dnn_out)

    minibatch_featurizer = Lambda(minibatch_discriminator,
                                  output_shape=minibatch_output_shape)

    # concat the minibatch features with the normal ones
    features = merge([
        minibatch_featurizer(K_x),
        dnn_out
    ], mode='concat')

    # fake output tracks binary fake / not-fake, and the auxiliary requires
    # reconstruction of latent features, in this case, labels
    fake = Dense(1, activation='sigmoid', name='generation')(features)
    aux = Dense(1, activation='sigmoid', name='auxiliary')(features)

    return Model(input=image, output=[fake, aux])
