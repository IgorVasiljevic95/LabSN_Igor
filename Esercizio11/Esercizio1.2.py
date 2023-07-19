
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

def build_model(num_layers, num_neurons, activation):
    model = Sequential()
    model.add(Dense(num_neurons, input_dim=1, activation=activation))
    for _ in range(num_layers - 1):
        model.add(Dense(num_neurons, activation=activation))
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

# Model configurations to explore
num_layers_list = [1, 2, 3]
num_neurons_list = [10, 20, 30]
activation_list = ['relu', 'tanh', 'sigmoid']
optimizer = 'adam'

# Train and evaluate models for different configurations
results = {}
for num_layers in num_layers_list:
    for num_neurons in num_neurons_list:
        for activation in activation_list:
            model = build_model(num_layers, num_neurons, activation)
            train_model(model, X_train, y_train, num_epochs)
            mse = evaluate_model(model, X_test, y_test)
            results[(num_layers, num_neurons, activation)] = mse

# Find the best model based on MSE
best_config = min(results, key=results.get)
best_mse = results[best_config]
best_model = build_model(*best_config)

print("Best Configuration (Layers, Neurons, Activation):", best_config)
print("Best Mean Squared Error:", best_mse)

# Train the best model on the entire dataset for the final fit
best_model.fit(X, y, epochs=num_epochs, verbose=0)

# Plot the results for the best model
plt.figure(figsize=(10, 6))
plt.scatter(X, y, color='blue', label='True Data')
plt.scatter(X, best_model.predict(X), color='red', label='Predicted Data')
plt.xlabel('X')
plt.ylabel('y')
plt.title('Regression using Neural Network (Best Model)')
plt.legend()
plt.show()

# Generate test data points within and outside the training range
X_test_within_range, y_test_within_range = generate_data(1000, noise_level)
X_test_outside_range, y_test_outside_range = generate_data(500, noise_level)

# Evaluate the best model on the test data within the training range
mse_within_range = evaluate_model(best_model, X_test_within_range, y_test_within_range)

# Evaluate the best model on the test data outside the training range
mse_outside_range = evaluate_model(best_model, X_test_outside_range, y_test_outside_range)

print("MSE on Test Data Within Training Range:", mse_within_range)
print("MSE on Test Data Outside Training Range:", mse_outside_range)

# Plot the results for the best model on the training data and test data
plt.figure(figsize=(10, 6))
plt.scatter(X, y, color='blue', label='True Data')
plt.scatter(X, best_model.predict(X), color='red', label='Predicted Data (Training Range)')
plt.scatter(X_test_outside_range, best_model.predict(X_test_outside_range), color='green', label='Predicted Data (Outside Training Range)')
plt.xlabel('X')
plt.ylabel('y')
plt.title('Regression using Neural Network (Best Model)')
plt.legend()
plt.show()
