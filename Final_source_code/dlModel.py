#Import the required packages

import os
from pathlib import Path
import pandas as pd
import numpy as np

import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout, Activation, Flatten, Conv2D, MaxPooling2D
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras.models import load_model
from sklearn.preprocessing import LabelEncoder
from sklearn import preprocessing

import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt

import Init as init

import warnings
warnings.filterwarnings("ignore")

meta = init.meta
manifest = init.manifest

current_dir = os.getcwd() #current directory
data_dir = os.path.join(current_dir, "data")

fileName = os.path.join(data_dir, "merged_meta.tsv")
meta=pd.read_csv(os.path.join(data_dir,manifest['merged_meta']),sep="\t",dtype={'SMUBRID':object,'SEX':object,'DTHHRDY':object})

meta=meta[~(meta['AGE'].isnull())] # removes all samples without age

#Extract only the tissue types with count > 200
counts=pd.DataFrame(meta['SMTS'].value_counts())
df=meta[meta['SMTS'].isin(counts[counts['SMTS'] > 200].index)]

#Identify the unique tissue types
tissue_types = df['SMTS'].unique()

#print("Unique tissue types in the GTEx data: ", [t for t in tissue_types])

from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout, Activation, Flatten, Conv2D, MaxPooling2D
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras.models import load_model

class keras_model:
    DATA_DIR = "keras_models"
    
    def __init__(self):
        self.early_stopping_monitor = EarlyStopping(patience=3) #Hyperparameter tuning
        
    def construct_model(self, x_train, y_train):
        self.model = tf.keras.models.Sequential() #Sequential model
        self.model.add(tf.keras.layers.Flatten())
        self.model.add(tf.keras.layers.Dense(1024, input_dim=x_train.shape[1], activation=tf.nn.relu))
        self.model.add(tf.keras.layers.Dense(512, activation=tf.nn.relu))
        #self.model.add(tf.keras.layers.Dense(256, activation=tf.nn.relu))
        self.model.add(tf.keras.layers.Dense(128, activation=tf.nn.relu))
        self.model.add(tf.keras.layers.Dense(y_train.shape[1], activation=tf.nn.softmax))
        #Note: Output layer is designed to hold the number of neurons equivalent to the number of classes of age groups

        self.model.compile(optimizer='adam', loss='categorical_crossentropy', metrics=['accuracy'])
        
    def model_train(self, x_train, y_train):
        if x_train.shape[0] == y_train.shape[0]: #Checking shape since there are tissues with missing gene expressions
            try:
                
                self.construct_model(x_train, y_train) #Constructing the model
                self.model.fit(x_train, 
                               y_train, 
                               batch_size=32, #Setting batch size for ease of processing in local machines
                               epochs=30, #Maximum of 30 epochs
                               validation_split=0.1,
                               verbose = 0,
                               callbacks=[self.early_stopping_monitor]) #Early stopping hyperparameter
                
                #print("Evaluating training accuracy...")
                loss, accuracy = self.model.evaluate(x_train, y_train, verbose=0)
                
                try:
                    #Persisting the model trained for the corresponding tissue type
                    fileName = TISSUE + "_keras_model.h5"
                    filePath = os.path.join(self.DATA_DIR, fileName)
                    self.model.save(filePath)
                    return accuracy, fileName
                except Exception as e:
                    print("Exception while saving the model.")
                    print(e)
                    return accuracy, None
            except:
                print("Exception while processing!")
                return -1, None
        else:
            #print("Shape mismatch encountered!")
            return -1, None
        
tissue_specific_path = "tissue-specific"

tissue_type = []
tissue_model_persist = []
tissue_model_accuracy = []

for tissue in tissue_types:
    k_model = keras_model()
    TISSUE=tissue
    infiles=os.listdir(tissue_specific_path)
    TISSUE_files=[f for f in infiles if  TISSUE in f]
    for entry in TISSUE_files:
        if "_cpm" in entry: #Identify the file with _cpm suffix; cpm stands for Counts Per Million
            pdd = pd.read_csv(os.path.join(tissue_specific_path,entry), sep='\t')
            
            #print("Tissue type: ", TISSUE)

            #Dropping the gene id colunm since it plays no role in classification
            pdd = pdd.drop(pdd.columns[0], axis='columns')

            #Min_max normalization
            min_max_scaler = preprocessing.MinMaxScaler()
            np_scaled = min_max_scaler.fit_transform(pdd)
            pdd = pd.DataFrame(np_scaled)
            
            numpy_matrix = pdd.as_matrix()
            #Categorizing the target column and performing one-hot encoding
            tissue_meta=meta[meta['SMTS']==TISSUE]
            encoder = LabelEncoder()
            age_y = tissue_meta['AGE']
            encoder.fit(tissue_meta['AGE'])
            encoded_Y = encoder.transform(tissue_meta['AGE'])
            dummy_y = tf.keras.utils.to_categorical(encoded_Y)
            
            #Training the model for the current tissue type
            acc, fileName = k_model.model_train(numpy_matrix, dummy_y)
            if acc != -1:
                try:
                    acc = acc * 100
                    tissue_type.append(TISSUE)
                    tissue_model_persist.append(fileName)
                    tissue_model_accuracy.append(acc)
                    #print("Final accuracy:", acc)
                except:
                    print("Error occurred for tissue type: ", TISSUE)
            #print("\n")
            break
            
keras_model = pd.DataFrame(
    data = {
        'tissue_type': tissue_type,
        'model_file': tissue_model_persist,
        'model_accuracy': tissue_model_accuracy
    }
)

keras_model = keras_model.sort_values('tissue_type')

keras_model.to_csv("keras_model_results.tsv",sep="\t",index=False)

mpl.style.use('default')
fig = plt.figure(figsize=(5,3))
ax = fig.add_subplot(111)

model_plot = sns.barplot(x = keras_model['tissue_type'],
            y=keras_model['model_accuracy'],
            data = keras_model,
            orient='v',
            estimator=np.mean,
            capsize=0.1)

tick_labels = keras_model['tissue_type'].as_matrix()
ax.set_xticklabels(tick_labels,rotation=30)
ax.set(ylim=(0,50))
ax.set(title='Keras: accuracy across tissue types')
ax.set(xlabel='Tissue types')
ax.set(ylabel='Training accuracy')

plt.show()
test = model_plot.get_figure()
test.savefig('keras_accuracy.png', bbox_inches='tight')