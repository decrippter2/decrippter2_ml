from pydantic import BaseModel
from typing import Self, Any
from pathlib import Path
import json
import torch
import pandas as pd
from sklearn.model_selection import train_test_split
from torch.utils.data import DataLoader, TensorDataset
from torch import nn
class BaseClassifier(BaseModel):
    """Pydantic_based class to load/retrain all classifiers

        Attributes:
            dataset = Pandas dataframe used to retrain models or make predictions
            model folder = path to directory containing pretrained models
            feature_path = path to JSON file containing feature names
            feature_list = list containing names of features
    """
    dataset: Any = None
    model_folder: Path = Path(__file__).parent.parent.parent.joinpath('pretrained_models')
    feature_path: Path = Path(__file__).parent.parent.joinpath("feature_extraction_manager/feature_list.json")
    feature_list: list = []

    def load_feature_list(self:Self):
        """Loads feature list to self. variable"""
        with open(self.feature_path) as feature_json:
            self.feature_list=json.load(feature_json)["feature_list"]

    def set_device(self:Self):
        """Sets cuda or cpu as device"""
        device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
        return device

    def svm_training_split(self:Self,dataset_file, feature_list):
        """Performs 80/20 training/validation split. (Is it necessary?)

        Args:
            dataset_file: Pandas Dataframe containing all features
            feature_list: list of extractes features

        returns:
            4 dataframes containing x and y values for training and testing
        """
        dataset = pd.read_csv(dataset_file)
        x = dataset[feature_list]
        y = dataset[["RiPP"]].to_numpy().ravel()
        x_train, x_test, y_train, y_test = train_test_split(x, y, stratify=y, test_size=0.2, random_state=42)
        return x_train, x_test, y_train, y_test

    def initiate_train_dataloader(self:Self):
        """Initiates DataLoader object to be passed to NN when re-training

        returns:
            DataLoader object"""
        x = self.dataset[self.feature_list]
        y = self.dataset[["RiPP"]].to_numpy().ravel()

        label_dict = {'RiPP': 1, 'No_RiPP': 0}

        x_train = torch.tensor(x.values)

        y_train = torch.tensor(pd.Series(y).map(label_dict).tolist())

        training_data = TensorDataset(x_train.to(torch.float32), y_train)

        dataloader = DataLoader(training_data, batch_size=256)
        return dataloader
    def initiate_predict_dataloader(self:Self,dataset):
        """Initiates DataLoader object to be passed to NN for predicting

        Args:
            dataset: Pandas dataframe containing extracted features

        Returns:
            DataLoader object"""
        self.load_feature_list()
        x=dataset[self.feature_list]
        x['length']=x.length.astype(float)
        x=torch.tensor(x.values, dtype=torch.float32)
        x_data=TensorDataset(x.to(torch.float32))
        dataloader = DataLoader (x_data,batch_size=128)
        return dataloader

    def init_weights(self:Self,layer):
        """Function for weight initialization for NN re-training"""
        if type(layer) == nn.Linear or type(layer) == nn.Conv2d:
            nn.init.xavier_uniform_(layer.weight)