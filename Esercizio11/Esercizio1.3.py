import numpy as np
import matplotlib.pyplot as plt
import tensorflow as tf
from tensorflow import keras
from keras import Sequential
from keras.layers import Dense, Activation
from keras import backend as K
from keras.utils import get_custom_objects
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error

# Function to generate data for the trigonometric 2D function
def generate_data(num_samples, noise_level):
    x = np.linspace(-3/2, 3/2, num_samples)
    y = np.linspace(-3/2, 3/2, num_samples)
    X, Y = np.meshgrid(x, y)
    Z = np.sin(X**2 + Y**2) + np.random.normal(0, noise_level, (num_samples, num_samples))
    return X, Y, Z

# Function to build the model with fixed configurations
def build_model():
    model = Sequential()
    model.add(Dense(20, input_dim=2, activation='relu'))
    model.add(Dense(20, activation='relu'))
    model.add(Dense(20, activation='relu'))
    model.add(Dense(1, activation='linear'))
    model.compile(loss='mean_squared_error', optimizer='nadam')
    return model

# Function to train the model
def train_model(model, X_train, y_train, num_epochs, x_valid, y_valid):
    # Custom callback to plot the current estimate during training
    class PlotCurrentEstimate(tf.keras.callbacks.Callback):
        def __init__(self, x_valid, y_valid):
            self.x_valid = x_valid
            self.y_valid = y_valid
            self.iter = 0

        def on_epoch_end(self, epoch, logs={}):
            temp = self.model.predict(self.x_valid, batch_size=None, verbose=False, steps=None)
            self.y_curr = temp.flatten()

            self.iter += 1
            if self.iter % 10 == 0:
                plt.clf()
                plt.scatter(self.x_valid[:, 0], self.y_curr, color="blue", s=4, marker="o", label="estimate")
                plt.scatter(self.x_valid[:, 0], self.y_valid, color="red", s=4, marker="x", label="valid")
                plt.legend()
                plt.pause(0.05)

    # Custom callback for early stopping based on validation loss
    early_stopping = keras.callbacks.EarlyStopping(monitor='val_loss', patience=20, restore_best_weights=True)

    history = model.fit(X_train, y_train, epochs=num_epochs, validation_data=(x_valid, y_valid),
                        callbacks=[PlotCurrentEstimate(x_valid, y_valid), early_stopping], verbose=0)
    return history

# Function to evaluate the model
def evaluate_model(model, X_test, y_test):
    y_pred = model.predict(X_test)
    mse = mean_squared_error(y_test, y_pred)
    return mse

# Parameters
num_epochs = 200
num_train_samples = 50
noise_level = 0.1

# Generate noisy data for the trigonometric 2D function
X, Y, Z = generate_data(num_train_samples, noise_level)

# Flatten the 2D arrays to use as input for the model
X_train = np.column_stack((X.ravel(), Y.ravel()))
y_train = Z.ravel()

# Split the data into training and validation sets
X_train, X_valid, y_train, y_valid = train_test_split(X_train, y_train, test_size=0.2)

# Build the model with fixed configurations
model = build_model()

# Train the model
train_model(model, X_train, y_train, num_epochs, X_valid, y_valid)

# Plot the results
X_plot, Y_plot = np.meshgrid(np.linspace(-3/2, 3/2, 100), np.linspace(-3/2, 3/2, 100))
Z_true = np.sin(X_plot**2 + Y_plot**2)

Z_pred_flat = model.predict(np.column_stack((X_plot.ravel(), Y_plot.ravel())))
Z_pred = Z_pred_flat.reshape(X_plot.shape)

fig = plt.figure(figsize=(12, 6))

# True Function and Predicted Function Surface Plot
ax = fig.add_subplot(1, 1, 1, projection='3d')
true_plot = ax.plot_surface(X_plot, Y_plot, Z_true, cmap='viridis', alpha=0.7, label='True Function')
pred_plot = ax.plot_surface(X_plot, Y_plot, Z_pred, cmap='plasma', alpha=0.7, label='Predicted Function')

# Add a colorbar to indicate the mapping of the true and predicted functions
fig.colorbar(true_plot, ax=ax, shrink=0.6)

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('f(X, Y)')
ax.set_title('True Function vs. Predicted Function')

plt.show()