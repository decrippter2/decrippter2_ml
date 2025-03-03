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
import torch
import torchvision
from torch import nn
import torch.optim as optim
import torch.nn.functional as F
from torchvision import datasets
from torchvision.transforms import ToTensor
from torch.utils.data import DataLoader, Dataset, TensorDataset
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
import pandas as pd
import numpy as np
import json
from sklearn import svm
from matplotlib import pyplot as plt
from classifier_class import BaseClassifier

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler()
logger.addHandler(console_handler)

class FeatureExtractor(nn.Module):
    def __init__(self, input_size):
        super(FeatureExtractor, self).__init__()
        self.input = nn.Dropout(0.2)
        self.fc0 = nn.Linear(input_size, 50)
        # self.mid_layer = nn.BatchNorm1d(50)
        self.fc1 = nn.Linear(50, 10)
        self.fc2 = nn.Dropout(0.5)
        self.fc3 = nn.Linear(10, 2)

    def forward(self, x):
        x = self.input(x)
        x = self.fc0(x)
        x = F.relu(x)
        x = F.relu(self.fc1(x))
        #x = self.mid_layer(x)
        x = self.fc2(x)
        x = self.fc3(x)
        return x  # This is passed to the SVM loss

class NN_Classifier(BaseClassifier):
    """Pydantic_based class to load/retrain SVM classifiers

        Attributes:
                folder_path=path to JSON with RiPP entries
    """

    def retrain_model(self:Self):
        self.load_feature_list()
        device = self.set_device()
        neurnet = FeatureExtractor(len(self.feature_list))
        loss_fn = nn.CrossEntropyLoss()
        optimizer = optim.Adam(neurnet.parameters(), lr=0.0001)
        neurnet.apply(self.init_weights)
        train_dataloader, test_dataloader = self.initiate_dataloader()

        for epoch in range(500):
            neurnet.train()
            for batch, data in enumerate(train_dataloader):
                inputs = data[0].to(device).float()
                labels = data[1].to(device).float()

                optimizer.zero_grad()

                # Extract features
                logits = neurnet(inputs)

                # Compute SVM loss
                loss = loss_fn(logits, labels.long())

                # Backpropagate through the NN (this optimizes the NN based on SVM performance)
                loss.backward()
                optimizer.step()

        neurnet.eval()

        torch.save(neurnet, self.model_folder.joinpath('nn.pth'))

    def load_model(self:Self):
        neurnet = torch.load(self.model_folder.joinpath('nn.pth'))
        return neurnet
