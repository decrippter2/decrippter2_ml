
from pydantic import BaseModel
from typing import Self, Any
from pathlib import Path
import json
import torch
import pandas as pd
from sklearn.model_selection import train_test_split
from torch.utils.data import DataLoader, Dataset, TensorDataset
from torch import nn
class BaseClassifier(BaseModel):
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