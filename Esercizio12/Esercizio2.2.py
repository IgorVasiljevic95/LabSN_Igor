from keras.models import Sequential
from keras.layers import Conv2D, MaxPooling2D, Dropout, Flatten, Dense
import matplotlib.pyplot as plt
import numpy as np
from keras.datasets import mnist
from keras.utils import to_categorical
from keras.optimizers import SGD, Adam, RMSprop
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'

def create_CNN():
    model = Sequential()
    model.add(Conv2D(32, kernel_size=(3, 3), activation='relu', input_shape=(img_rows, img_cols, 1)))
    model.add(MaxPooling2D(pool_size=(2, 2)))
    model.add(Dropout(0.25))
    model.add(Conv2D(64, kernel_size=(3, 3), activation='relu'))
    model.add(MaxPooling2D(pool_size=(2, 2)))
    model.add(Dropout(0.25))
    model.add(Flatten())
    model.add(Dense(128, activation='relu'))
    model.add(Dropout(0.5))
    model.add(Dense(num_classes, activation='softmax'))
    return model

def compile_model(optimizer):
    model = create_CNN()
    model.compile(loss='categorical_crossentropy',
                  optimizer=optimizer,
                  metrics=['accuracy'])
    return model


# Load MNIST dataset
(X_train, y_train), (X_test, y_test) = mnist.load_data()

# Reshape and normalize the input data
img_rows, img_cols = 28, 28
num_classes = 10

X_train = X_train.reshape(X_train.shape[0], img_rows, img_cols, 1)
X_test = X_test.reshape(X_test.shape[0], img_rows, img_cols, 1)
X_train = X_train.astype('float32') / 255
X_test = X_test.astype('float32') / 255

# Convert the class labels to categorical format
Y_train = to_categorical(y_train, num_classes=num_classes)
Y_test = to_categorical(y_test, num_classes=num_classes)

# Training parameters
batch_size = 32
epochs = 20  # Increase the number of epochs

# Create and compile the models with different optimizers
model_CNN_sgd = compile_model(SGD())
model_CNN_adam = compile_model(Adam())
model_CNN_rmsprop = compile_model(RMSprop())

# Train the CNNs and store training info in history
history_cnn_sgd = model_CNN_sgd.fit(X_train, Y_train,
                                    batch_size=batch_size,
                                    epochs=epochs,
                                    verbose=1,
                                    validation_data=(X_test, Y_test))

history_cnn_adam = model_CNN_adam.fit(X_train, Y_train,
                                      batch_size=batch_size,
                                      epochs=epochs,
                                      verbose=1,
                                      validation_data=(X_test, Y_test))

history_cnn_rmsprop = model_CNN_rmsprop.fit(X_train, Y_train,
                                            batch_size=batch_size,
                                            epochs=epochs,
                                            verbose=1,
                                            validation_data=(X_test, Y_test))

# Plot the training and validation accuracy for each optimizer
plt.plot(history_cnn_sgd.history['accuracy'])
plt.plot(history_cnn_sgd.history['val_accuracy'])
plt.plot(history_cnn_adam.history['accuracy'])
plt.plot(history_cnn_adam.history['val_accuracy'])
plt.plot(history_cnn_rmsprop.history['accuracy'])
plt.plot(history_cnn_rmsprop.history['val_accuracy'])
plt.title('CNN Model Accuracy')
plt.ylabel('Accuracy')
plt.xlabel('Epoch')
plt.legend(['SGD_train', 'SGD_val', 'Adam_train', 'Adam_val', 'RMSprop_train', 'RMSprop_val'], loc='best')
plt.show()

# Plot the training and validation loss for each optimizer
plt.plot(history_cnn_sgd.history['loss'])
plt.plot(history_cnn_sgd.history['val_loss'])
plt.plot(history_cnn_adam.history['loss'])
plt.plot(history_cnn_adam.history['val_loss'])
plt.plot(history_cnn_rmsprop.history['loss'])
plt.plot(history_cnn_rmsprop.history['val_loss'])
plt.title('CNN Model Loss')
plt.ylabel('Loss')
plt.xlabel('Epoch')
plt.legend(['SGD_train', 'SGD_val', 'Adam_train', 'Adam_val', 'RMSprop_train', 'RMSprop_val'], loc='best')
plt.show()
