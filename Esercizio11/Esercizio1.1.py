
import numpy as np
import matplotlib.pyplot as plt

# compose the NN model
import tensorflow as tf
from tensorflow import keras

from keras import Sequential
from keras.layers import Dense, Activation
from keras import backend as K
from keras.utils import get_custom_objects

from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from keras.models import Sequential
from keras.layers import Dense

def generate_data(num_samples, noise_level):
    X = np.linspace(0, 1, num_samples)
    y = 5*X + 2 + np.random.normal(0, noise_level, num_samples)
    return X, y

def build_model():
    model = Sequential()
    model.add(Dense(10, input_dim=1, activation='relu'))
    model.add(Dense(1, activation='linear'))
    model.compile(loss='mean_squared_error', optimizer='adam')
    return model

def train_model(model, X_train, y_train, num_epochs):
    history = model.fit(X_train, y_train, epochs=num_epochs, verbose=0)
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
train_model(model, X_train, y_train, num_epochs)

# Evaluate the model
mse = evaluate_model(model, X_test, y_test)
print("Mean Squared Error:", mse)

# Plot the results
plt.scatter(X_test, y_test, color='blue', label='True Data')
plt.scatter(X_test, model.predict(X_test), color='red', label='Predicted Data')
plt.xlabel('X')
plt.ylabel('y')
plt.title('Regression using Neural Network')
plt.legend()
plt.show()