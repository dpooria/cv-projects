#! /usr/bin/python

"""

solving an ode with neural networks

f'(x) = exp(-x/5) * cos(x) - f(x) / 5
f(0) = 0

"""

import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt


def label_function(x, psi):
    return tf.exp(-x/5) * tf.cos(x) - psi / 5


def build_model():
    inputs = tf.keras.Input((1, ))
    x = tf.keras.layers.Dense(12, activation="relu", name="dense_relu_1")(inputs)
    x = tf.keras.layers.Dense(12, activation="relu", name="dense_relu_2")(x)
    x = tf.keras.layers.Dense(1, activation="relu", name="dense_relu_3")(x)
    outputs = tf.keras.layers.Dense(
        1, activation="linear", name="dense_linear_1")(x)
    model = tf.keras.Model(inputs=inputs, outputs=outputs, name="ode_model")
    return model


Model = build_model()
print(Model.summary())
#tf.keras.utils.plot_model(Model, "ode_model.png", show_shapes=True)
# hyper parameters
SAMPLE_SIZE = 1000000
LR = 0.01
BATCH_NUM = 16
BATCH_SIZE = SAMPLE_SIZE // BATCH_NUM
EPOCHS = 40

# preparing the dataset
X = tf.constant(np.linspace(0, 2, SAMPLE_SIZE), dtype=tf.float32, shape=(SAMPLE_SIZE, 1))
train_dataset = tf.data.Dataset.from_tensor_slices(X)
train_dataset = train_dataset.shuffle(buffer_size=1024).batch(BATCH_SIZE)

# optimizer and loss
opt = tf.keras.optimizers.SGD(learning_rate=LR)
#opt = tf.keras.optimizers.Adam(learning_rate=LR)
loss_fn = tf.keras.losses.MeanSquaredError(reduction="auto", name="mean_squared_error")

@tf.function
def train_step(X):
    with tf.GradientTape() as g:
        g.watch(X)
        N = Model(X)
    dN = g.gradient(N, X)
    with tf.GradientTape() as tape:
        N = Model(X, training=True)
        f = label_function(X, X * N)
        loss = loss_fn(f, N + X * dN)

    grads = tape.gradient(loss, Model.trainable_variables)
    opt.apply_gradients(zip(grads, Model.trainable_variables))
    return loss


# training loop
for epoch in range(EPOCHS):
    print(f"------Epoch number {epoch}-----")
    for step, x_batch in enumerate(train_dataset):
        loss_value = train_step(x_batch)
        # Log every 10 batches.
        if step % 10 == 0:
            print(
                "Training loss (for one batch) at step %d: %f" % (
                    step, loss_value.numpy())
            )
            print("Seen so far: %s samples" % ((step + 1) * BATCH_SIZE))


print("ok")
psi = X.numpy() * Model(X).numpy()
plt.plot(X, psi)
plt.savefig("test_ode.png")

