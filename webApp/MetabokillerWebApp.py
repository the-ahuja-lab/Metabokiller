import os
from collections import Counter
import deepchem as dc
import matplotlib
import numpy as np
import pandas as pd
from boruta import BorutaPy
from deepchem.data import NumpyDataset
from flask import Flask, render_template, request, send_from_directory
from imblearn.over_sampling import SMOTE
from mordred import Calculator, descriptors
from PIL import Image
from rdkit import Chem
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier
import torch
matplotlib.use('Agg')
import time

import joblib
import matplotlib.pyplot as plt
import seaborn as sns
from fpdf import FPDF

# from webApp import app

app = Flask('__name__',template_folder='templates')

def extract_features(smiles):
  calc = Calculator(descriptors)
  mols = [Chem.MolFromSmiles(smi) for smi in smiles]  
  df = calc.pandas(mols)
  return df

def perform_column_pruning(data,th=75):
  data = data.replace(r'\s+', np.nan, regex=True)
  data[data == np.inf] = np.nan
  data = data.replace(r'^\s*$', np.nan, regex=True)  

  na_sum_series = data.isna().mean()
  org_data = data.copy()

  NAN_data = pd.DataFrame({0: na_sum_series.index, 1: na_sum_series.values})
  dropped = []
  for i in range(len(NAN_data)):
      if NAN_data.iloc[i][1] >= (th / 100):  # TODO check if sum or sum/length i.e. avg greater than threshold
          dropped.append(NAN_data.iloc[i][0])
  data = data.drop(dropped, axis=1)
  return data

def handle_missing_values(data):
  data = data.replace([np.inf, -np.inf, "", " "], np.nan)
  data = data.replace(["", " "], np.nan)
  data.fillna(data.mean(), inplace=True)
  return data

def smote_imb(xtrain, ytrain):
        print('Inside SMOTE, Original dataset shape %s' % Counter(ytrain))
        # cols = xtrain.columns
        sm = SMOTE(random_state=50)
        xtrain_new, ytrain_new = sm.fit_resample(xtrain, ytrain)
        
        # xtrain_df = pd.DataFrame(xtrain_new, columns = cols)
        print('Inside SMOTE, Resampled dataset shape %s' % Counter(ytrain_new))
        return xtrain_new, ytrain_new

def boruta_fs(xtrain, xtest, ytrain):
    print("Inside BorutaFS, Before Shape Train: ", xtrain.shape)
    print("Inside BorutaFS, Before Shape Test: ", xtest.shape)
    rfc = RandomForestClassifier(n_estimators=100, n_jobs=-1)
    boruta_selector = BorutaPy(rfc, n_estimators='auto', random_state=50)
    boruta_selector.fit(xtrain.values, ytrain)
    xtrain_sel = boruta_selector.transform(xtrain.values)
    xtest_sel = boruta_selector.transform(xtest.values)
#     print(boruta_selector.support_)
    # sel_cols = xtrain.columns[boruta_selector.support_]
    print("Inside BorutaFS, After Shape Train: ", xtrain_sel.shape)
    print("Inside BorutaFS,  After Shape Test: ", xtest_sel.shape)
    # print("Inside BorutaFS, IN FeatureSelector get_feature_names ", sel_cols)
    return xtrain_sel, xtest_sel

def get_predictions(probs,thresh):
    preds = []
    for prob in probs:
        if prob >= thresh:
            preds.append(1)    
        else:
            preds.append(0)
    return preds

class model_epigenetic:
    def __init__(self,test):                
        self.test = test
    def extract_feature(self,smiles):        
        labels = [0]*len(smiles)        
        featurizer = dc.feat.MolGraphConvFeaturizer()
        features = featurizer.featurize(smiles)            
        dataset = NumpyDataset(features,labels)  
        return dataset 
    def get_labels(self,pred_test): #Getting discrete labels from probability values    
        test_pred = []        
        for i in range(pred_test.shape[0]):
            if(pred_test[i][0]>pred_test[i][1]):
                test_pred.append(0)
            else:
                test_pred.append(1)
        return test_pred     
       
    def test_model(self):
        test = self.test
        test_dataset = self.extract_feature(test)
        model = dc.models.GCNModel(n_tasks=1,mode='classification',graph_conv_layers=[64,64], predictor_hidden_feats=64,
        predictor_dropout = 0.2, Learning_rate = 0.001, batch_size = 32,model_dir='models/model_epigenetic')
        model.restore()
        probs = model.predict(test_dataset)    
        preds = self.get_labels(probs)
        return probs,preds

class model_oxidative:    
    def __init__(self,test):                
        self.test = test
    def extract_feature(self,smiles):        
        labels = [0]*len(smiles)        
        featurizer = dc.feat.MolGraphConvFeaturizer(use_edges=True)
        features = featurizer.featurize(smiles)            
        dataset = NumpyDataset(features,labels)  
        return dataset 
    def get_labels(self,pred_test): #Getting discrete labels from probability values    
        test_pred = []        
        for i in range(pred_test.shape[0]):
            if(pred_test[i][0]>pred_test[i][1]):
                test_pred.append(0)
            else:
                test_pred.append(1)
        return test_pred     
       
    def test_model(self):
        test = self.test
        test_dataset = self.extract_feature(test)
        model = dc.models.AttentiveFPModel(n_tasks=1,mode='classification',num_layers=2,graph_feat_size=200,dropout=0.5,batch_size=32, learning_rate=0.001,model_dir='models/model_oxidative')
        model.restore()
        probs = model.predict(test_dataset)    
        preds = self.get_labels(probs)
        return probs,preds

class model_genotoxic:
    def __init__(self,test):                
        self.test = test
    def extract_feature(self,smiles):        
        labels = [0] * len(smiles)         
        featurizer = dc.feat.ConvMolFeaturizer()
        features = featurizer.featurize(smiles)            
        dataset = NumpyDataset(features,labels)  
        return dataset 
    def get_labels(self,pred_test): #Getting discrete labels from probability values    
        test_pred = []        
        for i in range(pred_test.shape[0]):
            if(pred_test[i][0][0]>pred_test[i][0][1]):
                test_pred.append(0)
            else:
                test_pred.append(1)
        return test_pred     
       
    def test_model(self):
        test = self.test
        test_dataset = self.extract_feature(test)
        model = dc.models.GraphConvModel(mode='classification',n_tasks=1,batch_size=32,learning_rate=0.001,model_dir='models/model_genotoxic')
        model.restore()
        probs = model.predict(test_dataset)    
        preds = self.get_labels(probs)
        return probs,preds

class model_apoptosis:
    def __init__(self,data,train_length):
        self.data = data
        self.train_length = train_length        

    def preprocess(self,labels): 
        length = self.train_length
        data = self.data                                   
        y_train = labels             
        data[data.columns] = data[data.columns].apply(pd.to_numeric,errors='coerce')        
        data = perform_column_pruning(data)
        train = data.iloc[:length,:]
        test = data.iloc[length:,:]
        train = handle_missing_values(train)        
        test = handle_missing_values(test) 
        x_train = train
        x_test = test        
        x_train,y_train = smote_imb(x_train,y_train)     
        x_train, x_test = boruta_fs(x_train, x_test, y_train)
        data = pd.DataFrame(x_test)  
        data.fillna(data.mean(), inplace=True) 
        x_test = np.array(data)          
        self.x_train,self.y_train,self.x_test = x_train,y_train,x_test
    
    def train_model(self):
        x_train,y_train = self.x_train,self.y_train
        rf = RandomForestClassifier(n_estimators=100,min_samples_split=2,min_samples_leaf=1,max_features='sqrt',max_depth=100,bootstrap=False)
        rf.fit(x_train,y_train)
        joblib.dump(rf, 'models/model_apoptosis.sav')

    def test_model(self):
        x_test = self.x_test 
        x_test = np.nan_to_num(x_test)                         
        rf = joblib.load('models/model_apoptosis.sav')        
        pred = rf.predict(x_test)
        prob = rf.predict_proba(x_test)
        return prob,pred

class model_proliferation:
    def __init__(self,data,train_length):
        self.data = data
        self.train_length = train_length  

    def preprocess(self,labels): 
        data = self.data                             
        y_train = labels            
        data[data.columns] = data[data.columns].apply(pd.to_numeric,errors='coerce')        
        data = perform_column_pruning(data)
        train = data.iloc[:self.train_length,:]
        test = data.iloc[self.train_length:,:]        
        train = handle_missing_values(train)        
        test = handle_missing_values(test)
        x_train = train
        x_test = test        
        x_train,y_train = smote_imb(x_train,y_train)     
        x_train, x_test = boruta_fs(x_train, x_test, y_train)
        self.x_train,self.y_train,self.x_test = x_train,y_train,x_test

    def train_model(self):
        x_train,y_train = self.x_train,self.y_train
        et = ExtraTreesClassifier(max_depth=15,n_estimators=90)
        et.fit(x_train,y_train)        
        joblib.dump(et, 'models/model_proliferations.sav')        

    def test_model(self):
        x_test = self.x_test.astype(np.float)
        x_test = np.nan_to_num(x_test,nan=0.0)         
        et = joblib.load('models/model_proliferations.sav')        
        pred = et.predict(x_test)
        prob = et.predict_proba(x_test)
        return prob,pred

class model_electrophiles:
    def __init__(self,test):                
        self.test = test
    def extract_feature(self,smiles):        
        labels = [0] * len(smiles)        
        featurizer = dc.feat.MolGraphConvFeaturizer(use_edges=True)
        features = featurizer.featurize(smiles)            
        dataset = NumpyDataset(features,labels)  
        return dataset 
    def get_labels(self,pred_test): #Getting discrete labels from probability values    
        test_pred = []   
        print(pred_test.shape)     
        for i in range(pred_test.shape[0]):
            if(pred_test[i][0]>pred_test[i][1]):
                test_pred.append(0)
            else:
                test_pred.append(1)
        return test_pred     

    def test_model(self):
        test = self.test
        test_dataset = self.extract_feature(test)
        model = dc.models.AttentiveFPModel(n_tasks=1,mode='classification',num_layers=2,graph_feat_size=200,dropout=0,batch_size=32, learning_rate=0.001,model_dir='models/model_electrophiles')
        model.restore()
        probs = model.predict(test_dataset)    
        preds = self.get_labels(probs)
        return probs,preds

@app.route("/")
def index():
    return render_template("index.html")

@app.route("/predict",methods=['GET','POST'])
def predict_():
    if request.method == 'GET':
        return render_template("predict.html",predicted=False)
    smiles = request.form['textSmiles'].replace("\r","").split('\n')
    while "" in smiles:
        smiles.remove("")
    smiles = smiles[:1]
    predictions = pd.DataFrame(columns=['smiles','epigenetic_0','epigenetic_1','epigenetic_preds',
    'oxidative_0','oxidative_1','oxidative_preds','genotoxic_0','genotoxic_1','genotoxic_preds',
    'apoptosis_0','apoptosis_1','apoptosis_preds','proliferations_0','proliferations_1','proliferations_preds',
    'electrophiles_0','electrophiles_1','electrophiles_preds','ensemble_0','ensemble_1','ensemble_preds'])    
    predictions['smiles'] = smiles
    print("B0", smiles)
    
    m1 = model_electrophiles(smiles)
    probs, preds = m1.test_model()
    print("B1", probs)
    predictions['electrophiles_0'] = probs[:,0]
    predictions['electrophiles_1'] = probs[:,1]
    predictions['electrophiles_preds'] = preds  
    
    m2 = model_genotoxic(smiles)
    probs, preds = m2.test_model()
    print("B2", probs)
    predictions['genotoxic_0'] = probs[:,0,0]
    predictions['genotoxic_1'] = probs[:,0,1]
    predictions['genotoxic_preds'] = preds
    
    m3 = model_oxidative(smiles)
    probs, preds = m3.test_model()
    print("B3", probs)
    predictions['oxidative_0'] = probs[:,0]
    predictions['oxidative_1'] = probs[:,1]
    predictions['oxidative_preds'] = preds
    
    m4 = model_epigenetic(smiles)
    probs, preds = m4.test_model()
    print("B4", probs)
    predictions['epigenetic_0'] = probs[:,0]
    predictions['epigenetic_1'] = probs[:,1]
    predictions['epigenetic_preds'] = preds 
    
    apop = pd.read_csv('Data/apoptosis.csv')
    smiles_apop = list(apop['smiles'])
    length = len(smiles_apop)
    # smiles_apop.extend(smiles)
    mord_apop_input = extract_features(smiles)
    mord_apop_train = pd.read_csv("Data/apoptosis_train_mordred.csv", low_memory=False)
    mord_apop = mord_apop_train.append(mord_apop_input, ignore_index=True)
    m5 = model_apoptosis(mord_apop,length)
    m5.preprocess(apop['status'])
    m5.train_model()
    probs, preds = m5.test_model()
    print("B5", probs)
    predictions['apoptosis_0'] = probs[:,0]
    predictions['apoptosis_1'] = probs[:,1]
    predictions['apoptosis_preds'] = preds  
    
    proliferations = pd.read_csv('Data/cell_proliferation.csv')
    smiles_prol = list(proliferations['smiles'])
    length = len(smiles_prol)
    # smiles_prol.extend(smiles)
    mord_prol_input = extract_features(smiles)
    mord_prol_train = pd.read_csv("Data/proliferation_train_mordred.csv", low_memory=False)
    mord_prol = mord_prol_train.append(mord_prol_input,ignore_index=True)
    m6 = model_proliferation(mord_prol,length)
    m6.preprocess(proliferations['status'])
    m6.train_model()
    probs, preds = m6.test_model()
    print("B6", probs)
    predictions['proliferations_0'] = probs[:,0]
    predictions['proliferations_1'] = probs[:,1]
    predictions['proliferations_preds'] = preds
    m1_1 = list(predictions['epigenetic_1'])    
    m2_1 = list(predictions['oxidative_1'])    
    m3_1 = list(predictions['genotoxic_1'])
    m4_1 = list(predictions['apoptosis_1'])
    m5_1 = list(predictions['proliferations_1'])
    m6_1 = list(predictions['electrophiles_1'])
    weights = {'epigenetic': 0.4640625, 'oxidative': 0.459375, 'genotoxic': 0.6296875, 'apoptosis': 1.0, 'proliferations': 0.853125, 'electrophiles': 0.5625}   
    probs = []
    sum = weights['epigenetic'] + weights['oxidative'] + weights['genotoxic'] + weights['apoptosis'] + weights['proliferations'] + weights['electrophiles']
    for i in range(len(m1_1)):
        prob_ensemble = (m1_1[i]*weights['epigenetic']+m2_1[i]*weights['oxidative']+m3_1[i]*weights['genotoxic']+m4_1[i]*weights['apoptosis']+m5_1[i]*weights['proliferations']+m6_1[i]*weights['electrophiles'])/sum
        probs.append(prob_ensemble)
    predictions['ensemble_1'] = probs
    probs_status0 = [1-prob for prob in probs]
    predictions['ensemble_0'] = probs_status0
    preds = get_predictions(probs,0.5)
    predictions['ensemble_preds'] = preds 
    # file =  'result_'+str(time.time())+'.csv' 
    file = 'result.csv'    
    # os.remove(file) 
    predictions.to_csv(file,index=False)
    torch.cuda.empty_cache()
    return render_template("predict.html",predicted=True)

@app.route('/download') 
def download_results():
    return send_from_directory('','result.csv',as_attachment=True, cache_timeout=-1)

@app.route('/visualize')
def visualize():
    result = pd.read_csv('result.csv')
    
    cols_probs = ['apoptosis_1', 'proliferations_1', 'oxidative_1', 'genotoxic_1', 
            'epigenetic_1', 'electrophiles_1', 'ensemble_1']
    column_names_bargraph = ['Anti-apoptotic Properties','Proliferative Properties','Oxidative Stress','Genotoxicity', 
                            'Epigenetic Modifications', 'Electrophilic Properties', 'Carcinogenicity']
    plt1 = plt.figure(1)
    bargraph_df = list(result.loc[0,cols_probs] )
    bar = sns.barplot(x = cols_probs, y = bargraph_df, color = 'black')
    plt.ylabel('Probability')
    bar.set_xticklabels(column_names_bargraph, rotation=90) 
    plt.tight_layout()
    plt.savefig('static/images/barplot.jpg')
    plt.close()

    ans_row = {'SMILE':result['smiles'][0]}
    final_ans = result['ensemble_preds'][0]
    if final_ans == 1:
        ans_row['Carcinogenicity'] = 'Carcinogen'
        ans_row['Probability'] = result['ensemble_1'][0]
    else:
        ans_row['Carcinogenicity'] = 'Non-carcinogen'
        ans_row['Probability'] = result['ensemble_0'][0]
    
    im1 = Image.open("static/images/barplot.jpg")
    im1 = im1.resize((1000,800))
    filename = "static/images/plots.pdf"
    im1.save(filename, "PDF" ,resolution=100.0, save_all=True)    
    return render_template('visualize.html',plots_generated=True, table = ans_row)

@app.route('/download_plots') 
def download_plots():
    return send_from_directory('','static/images/plots.pdf',as_attachment=True, cache_timeout=-1)

@app.route('/contact')
def contact():
    return render_template('contact.html')

import logging
logging.basicConfig(filename='/home/gaurav/MetabokillerWebApp/execution_log.log', filemode='a+', format=' [%(filename)s:%(lineno)s:%(funcName)s()]- %(asctime)s - %(levelname)s - %(message)s', level=logging.INFO)
 
gunicorn_logger = logging.getLogger('gunicorn.error')
app.logger.handlers = gunicorn_logger.handlers
app.logger.setLevel(gunicorn_logger.level)

if __name__ == '__main__':
    app.run()

