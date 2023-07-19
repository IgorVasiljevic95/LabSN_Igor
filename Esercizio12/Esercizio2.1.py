from keras.models import Sequential
from keras.layers import Dense, Dropout
import matplotlib.pyplot as plt
import numpy as np
from keras.datasets import mnist
from keras.utils import to_categorical
from keras.optimizers import SGD, Adam, RMSprop
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'

def create_DNN():
    model = Sequential()
    model.add(Dense(400, input_shape=(img_rows * img_cols,), activation='relu'))
    model.add(Dense(100, activation='relu'))
    model.add(Dropout(0.5))
    model.add(Dense(num_classes, activation='softmax'))
    return model

def compile_model(optimizer):
    model = create_DNN()
    model.compile(loss='categorical_crossentropy',
                  optimizer=optimizer,
                  metrics=['accuracy'])
    return model


# Load MNIST dataset
(X_train, y_train), (X_test, y_test) = mnist.load_data()

# Reshape and normalize the input data
img_rows, img_cols = 28, 28
num_classes = 10

X_train = X_train.reshape(X_train.shape[0], img_rows * img_cols)
X_test = X_test.reshape(X_test.shape[0], img_rows * img_cols)
X_train = X_train.astype('float32') / 255
X_test = X_test.astype('float32') / 255

# Convert the class labels to categorical format
Y_train = to_categorical(y_train, num_classes=num_classes)
Y_test = to_categorical(y_test, num_classes=num_classes)

# Training parameters
batch_size = 32
epochs = 20  # Increase the number of epochs

# Create and compile the models with different optimizers
model_DNN_sgd = compile_model(SGD())
model_DNN_adam = compile_model(Adam())
model_DNN_rmsprop = compile_model(RMSprop())

# Train the DNNs and store training info in history
history_sgd = model_DNN_sgd.fit(X_train, Y_train,
                                batch_size=batch_size,
                                epochs=epochs,
                                verbose=1,
                                validation_data=(X_test, Y_test))

history_adam = model_DNN_adam.fit(X_train, Y_train,
                                  batch_size=batch_size,
                                  epochs=epochs,
                                  verbose=1,
                                  validation_data=(X_test, Y_test))

history_rmsprop = model_DNN_rmsprop.fit(X_train, Y_train,
                                        batch_size=batch_size,
                                        epochs=epochs,
                                        verbose=1,
                                        validation_data=(X_test, Y_test))

# Plot the training and validation accuracy for each optimizer
plt.plot(history_sgd.history['accuracy'])
plt.plot(history_sgd.history['val_accuracy'])
plt.plot(history_adam.history['accuracy'])
plt.plot(history_adam.history['val_accuracy'])
plt.plot(history_rmsprop.history['accuracy'])
plt.plot(history_rmsprop.history['val_accuracy'])
plt.title('Model Accuracy')
plt.ylabel('Accuracy')
plt.xlabel('Epoch')
plt.legend(['SGD_train', 'SGD_val', 'Adam_train', 'Adam_val', 'RMSprop_train', 'RMSprop_val'], loc='best')
plt.show()

# Plot the training and validation loss for each optimizer
plt.plot(history_sgd.history['loss'])
plt.plot(history_sgd.history['val_loss'])
plt.plot(history_adam.history['loss'])
plt.plot(history_adam.history['val_loss'])
plt.plot(history_rmsprop.history['loss'])
plt.plot(history_rmsprop.history['val_loss'])
plt.title('Model Loss')
plt.ylabel('Loss')
plt.xlabel('Epoch')
plt.legend(['SGD_train', 'SGD_val', 'Adam_train', 'Adam_val', 'RMSprop_train', 'RMSprop_val'], loc='best')
plt.show()
