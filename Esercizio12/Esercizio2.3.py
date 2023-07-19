from tensorflow import keras
from keras.models import Sequential
from keras.layers import Conv2D, MaxPooling2D, Flatten, Dense, Dropout
from keras.datasets import mnist
from keras.utils import to_categorical
from keras.models import Model
import matplotlib.pyplot as plt
from PIL import Image
import numpy as np
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'

def create_CNN():
    model = Sequential()
    model.add(Conv2D(10, kernel_size=(5, 5),
                     activation='relu',
                     input_shape=input_shape))
    model.add(MaxPooling2D(pool_size=(2, 2)))
    model.add(Dropout(0.25))
    model.add(Conv2D(20, kernel_size=(3, 3), activation='relu'))
    model.add(Flatten())
    model.add(Dense(128, activation='relu'))
    model.add(Dropout(0.5))
    model.add(Dense(10, activation='softmax'))
    model.compile(loss=keras.losses.categorical_crossentropy,
                  optimizer='adam',
                  metrics=['accuracy'])
    return model

# Function to display the activations
def display_activation(activations, col_size, row_size, layer_index):
    activation = activations[layer_index]
    activation_index=0
    fig, ax = plt.subplots(row_size, col_size, figsize=(row_size*3,col_size*3))
    for row in range(0,row_size):
        for col in range(0,col_size):
            ax[row][col].imshow(activation[0, :, :, activation_index], cmap='gray')
            activation_index += 1

# Load MNIST dataset
(X_train, y_train), (X_test, y_test) = mnist.load_data()

# Reshape the input data
if keras.backend.image_data_format() == 'channels_first':
    X_train = X_train.reshape(X_train.shape[0], 1, 28, 28)
    X_test = X_test.reshape(X_test.shape[0], 1, 28, 28)
    input_shape = (1, 28, 28)
else:
    X_train = X_train.reshape(X_train.shape[0], 28, 28, 1)
    X_test = X_test.reshape(X_test.shape[0], 28, 28, 1)
    input_shape = (28, 28, 1)

# Convert the data type to float32 and normalize the pixel values
X_train = X_train.astype('float32') / 255
X_test = X_test.astype('float32') / 255

# Convert the class labels to categorical format
Y_train = to_categorical(y_train, num_classes=10)
Y_test = to_categorical(y_test, num_classes=10)

# Training parameters
batch_size = 32
epochs = 30

# Create the deep conv net
model_CNN = create_CNN()

# Train CNN
model_CNN.fit(X_train, Y_train,
              batch_size=batch_size,
              epochs=epochs,
              verbose=1,
              validation_data=(X_test, Y_test))

# Load the images
image_paths = ['digit0.png', 'digit1.png', 'digit2.png', 'digit3.png', 'digit4.png',
               'digit5.png', 'digit6.png', 'digit7.png', 'digit8.png', 'digit9.png']

images = []
for path in image_paths:
    image = Image.open(path).convert('L')  # Convert to grayscale
    image = image.resize((28, 28))  # Resize to 28x28
    image = np.array(image)  # Convert to numpy array
    image = image / 255.0  # Normalize to have pixel values between 0 and 1
    images.append(image)

# Convert the images to a numpy array
images = np.array(images)

# Reshape the images to add a channel dimension
images = np.expand_dims(images, axis=-1)

# Use the trained CNN model to predict the handwritten digits
predictions = model_CNN.predict(images)

# Get the predicted labels
predicted_labels = np.argmax(predictions, axis=1)

# Print the predicted labels
for i, label in enumerate(predicted_labels):
    print(f"Image {i}: Predicted Label = {label}")

# Index of the test image to visualize activations
test_index = 0

# Create an activation model
layer_outputs = [layer.output for layer in model_CNN.layers]
activation_model = Model(inputs=model_CNN.input, outputs=layer_outputs)

# Get the activations for the test image
activations = activation_model.predict(X_test[test_index].reshape(1,28,28,1))

plt.imshow(X_test[test_index][:,:,0], cmap='gray')
plt.show()
display_activation(activations, 4, 2, 0)


