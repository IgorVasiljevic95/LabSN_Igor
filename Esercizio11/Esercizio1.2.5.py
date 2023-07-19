
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

def generate_data(num_samples, noise_level):
    X = np.linspace(-1, 1, num_samples)
    y = 4 - 3*X - 2*X**2 + 3*X**3 + np.random.normal(0, noise_level, num_samples)
    return X, y

def build_model():
    model = Sequential()
    model.add(Dense(30, input_dim=1, activation='tanh'))
    model.add(Dense(30, activation='tanh'))
    model.add(Dense(30, activation='tanh'))
    model.add(Dense(1, activation='linear'))
    model.compile(loss='mean_squared_error', optimizer='rmsprop')
    return model

def train_model(model, X_train, y_train, num_epochs, x_valid, y_valid):
    # Custom callback to plot the current model estimate during training
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
                plt.scatter(self.x_valid, self.y_curr, color="blue", s=4, marker="o", label="estimate")
                plt.scatter(self.x_valid, self.y_valid, color="red", s=4, marker="x", label="valid")
                plt.legend()
                plt.pause(0.05)

    # Custom callback for early stopping based on validation loss
    early_stopping = keras.callbacks.EarlyStopping(monitor='val_loss', patience=20, restore_best_weights=True)

    history = model.fit(X_train, y_train, epochs=num_epochs, validation_data=(x_valid, y_valid),
                        callbacks=[PlotCurrentEstimate(x_valid, y_valid), early_stopping], verbose=0)
    return history

def evaluate_model(model, X_test, y_test):
    y_pred = model.predict(X_test)
    mse = mean_squared_error(y_test, y_pred)
    return mse

# Parameters
num_epochs = 200
num_train_samples = 500
noise_level = 0.3

# Generate noisy data
X, y = generate_data(num_train_samples, noise_level)

# Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)

# Build the model
model = build_model()

# Train the model
train_model(model, X_train, y_train, num_epochs, X_test, y_test)

# Evaluate the model on the test data
mse_test = evaluate_model(model, X_test, y_test)

# Print the MSE on the test data
print("Mean Squared Error on Test Data:", mse_test)

# Plot the results for the best model on the training data and test data
plt.figure(figsize=(10, 6))
plt.scatter(X, y, color='blue', label='True Data')
plt.scatter(X, model.predict(X), color='red', label='Predicted Data (Training Range)')
plt.xlabel('X')
plt.ylabel('y')
plt.title('Regression using Neural Network')
plt.legend()
plt.show()

