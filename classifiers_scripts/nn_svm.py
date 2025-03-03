import json
import torch.optim as optim
from sys import argv
import random
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from typing import Self, Any
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn import svm
from pydantic import BaseModel
import pickle
from pathlib import Path
import torch
from torch.utils.data import DataLoader, Dataset, TensorDataset
from torch import nn
import torch.nn.functional as F
from sklearn.svm import SVC
import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler()
logger.addHandler(console_handler)

class FeatureExtractor(nn.Module):
    def __init__(self, input_size, feature_dim):
        super(FeatureExtractor, self).__init__()
        self.input = nn.Dropout(0.2)
        self.fc1 = nn.Linear(input_size, 50)
        self.mid_layer = nn.BatchNorm1d(50)
        self.fc2 = nn.Linear(50, feature_dim)

    def forward(self, x):
        x = self.input(x)
        x = F.relu(self.fc1(x))
        x = self.mid_layer(x)
        x = self.fc2(x)
        return x  # This is passed to the SVM loss

class SVMLoss(nn.Module):
    def __init__(self, margin=1.0):
        super(SVMLoss, self).__init__()
        self.margin = margin  # SVM margin

    def forward(self, features, labels):
        labels = labels * 2 - 1  # Convert 0/1 labels to -1/+1 for hinge loss
        scores = torch.matmul(features, features.T)  # Compute similarity scores
        hinge_loss = torch.clamp(self.margin - labels * scores, min=0).mean()  # SVM hinge loss
        return hinge_loss

class NN_SVM_Classifier(BaseModel):
    """Pydantic_based class to load/retrain SVM classifiers

        Attributes:
                folder_path=path to JSON with RiPP entries
    """
    dataset: Any
    model_folder: Path = Path(__file__).parent.parent.joinpath('pretrained_models')
    feature_path: Path = Path(__file__).parent.parent.joinpath("feature_extraction_manager/feature_list.json")
    feature_list: list = []

    def load_feature_list(self:Self):
        with open(self.feature_path) as feature_json:
            self.feature_list=json.load(feature_json)["feature_list"]

    def set_device(self:Self):
        device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
        return device

    def svm_training_split(self:Self,dataset_file, feature_list):
        dataset = pd.read_csv(dataset_file)
        x = dataset[feature_list]
        y = dataset[["RiPP"]].to_numpy().ravel()
        x_train, x_test, y_train, y_test = train_test_split(x, y, stratify=y, test_size=0.2, random_state=42)
        return x_train, x_test, y_train, y_test

    def initiate_dataloader(self:Self):

        x_train, x_test, y_train, y_test = self.svm_training_split(self.dataset,self.feature_list)

        label_dict = {'RiPP': 1, 'No_RiPP': 0}

        x_train = torch.tensor(x_train.values)

        y_train = torch.tensor(pd.Series(y_train).map(label_dict).tolist())
        x_test = torch.tensor(x_test.values)

        y_test = torch.tensor(pd.Series(y_test).map(label_dict).tolist())

        training_data = TensorDataset(x_train.to(torch.float32), y_train)
        test_data = TensorDataset(x_test.to(torch.float32), y_test)

        train_dataloader = DataLoader(training_data, batch_size=256)
        test_dataloader = DataLoader(test_data, batch_size=256)
        return train_dataloader, test_dataloader

    def hybrid_model_predict(self,neurnet, svm_classifier, inputs,device):
        neurnet.eval()
        with torch.no_grad():
            features = neurnet(inputs.to(device).float()).cpu().numpy()
        return svm_classifier.predict(features)

    def init_weights(self:Self,layer):
        if type(layer) == nn.Linear or type(layer) == nn.Conv2d:
            nn.init.xavier_uniform_(layer.weight)

    def retrain_model(self:Self,model_type):
        self.load_feature_list()
        device = self.set_device()
        neurnet = FeatureExtractor(len(self.feature_list), 10)
        svm_loss_fn = SVMLoss()
        optimizer = optim.Adam(neurnet.parameters(), lr=0.001)
        neurnet.apply(self.init_weights)
        train_dataloader, test_dataloader = self.initiate_dataloader()

        for epoch in range(100):
            neurnet.train()
            for batch, data in enumerate(train_dataloader):
                inputs = data[0].to(device).float()
                labels = data[1].to(device).float()

                optimizer.zero_grad()

                # Extract features
                features = neurnet(inputs)

                # Compute SVM loss
                loss = svm_loss_fn(features, labels)

                # Backpropagate through the NN (this optimizes the NN based on SVM performance)
                loss.backward()
                optimizer.step()

        neurnet.eval()

        torch.save(neurnet, self.model_folder.joinpath('nn_svm.pth'))

        features_list, labels_list = [], []
        with torch.no_grad():
            for data in train_dataloader:
                inputs = data[0].to(device).float()
                labels = data[1].cpu().numpy()
                features = neurnet(inputs).cpu().numpy()
                features_list.append(features)
                labels_list.append(labels)
        X_train_svm = np.concatenate(features_list, axis=0)
        y_train_svm = np.concatenate(labels_list, axis=0)
        if model_type=='poly3' or model_type=='POLY3':
            svm_clf = SVC(C=100, kernel="poly",degree=3, random_state=0, class_weight="balanced")
            svm_clf.fit(X_train_svm, y_train_svm)
            svm_path = self.model_folder.joinpath('nn_svm_poly3.pkl')

        elif model_type=='rbf' or model_type=='RBF':
            svm_clf = SVC(C=1000, kernel="rbf", random_state=0, class_weight="balanced")
            svm_clf.fit(X_train_svm, y_train_svm)
            svm_path = self.model_folder.joinpath('nn_svm_rbf.pkl')

        else:
            logger.fatal('Error: invalid kernel')
            raise RuntimeError

        with open(svm_path, 'wb') as f:
            pickle.dump(svm_clf, f)

    def load_model(self:Self,model_type):
        neurnet= torch.load(self.model_folder.joinpath('nn_svm.pth'))
        if model_type=='poly3' or model_type=='POLY3':
            svm_path = self.model_folder.joinpath('nn_svm_poly3.pkl')
        elif model_type=='rbf' or model_type=='RBF':
            svm_path = self.model_folder.joinpath('nn_svm_rbf.pkl')
        else:
            logger.fatal('Error: invalid kernel')
            raise RuntimeError
        with open(svm_path, 'rb') as f:
            svm_clf = pickle.load(f)
        return neurnet,svm_clf