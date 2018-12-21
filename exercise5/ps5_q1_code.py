# CM226: PS5, Q1

from keras.models import Sequential
from keras.layers import Dense, Dropout
from keras.optimizers import SGD

from keras.utils import to_categorical

import pandas as pd
import numpy as np

np.random.seed(314159)

# number of SNP features
M = 20

# read in the data
genos = pd.read_csv('ps5_q1.genos', sep='\t')
phenos = pd.read_csv('ps5_q1.phenos', sep='\t', header=None)

# split data into training and testing
training_indices = range(0,800)
testing_indices = range(801,1000)

X_train = genos.iloc[training_indices].values
X_test  = genos.iloc[testing_indices].values

y_train = phenos.iloc[training_indices].values.ravel()
y_test = phenos.iloc[testing_indices].values.ravel()

y_train = to_categorical(y_train)
y_test = to_categorical(y_test)

# declare the NN architecture
model = Sequential()
model.add(Dense(64, activation='relu', input_dim=M))
model.add(Dense(16, activation='relu'))
model.add(Dense(32, activation='relu'))
model.add(Dropout(0.1))
model.add(Dense(2, activation='softmax'))

# use stochastic gradient descent for optimization
# lr: learning rate
sgd = SGD(lr=0.01, decay=1e-6, momentum=0.9, nesterov=True)
model.compile(loss='categorical_crossentropy',
              optimizer=sgd,
              metrics=['accuracy'])

model.fit(X_train, y_train, batch_size=8, epochs=10)
testing_loss, accuracy = model.evaluate(X_test, y_test, batch_size=8)

print('Testing accuracy: {}'.format(accuracy))
