
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
        x = self.dataset[self.feature_list]
        y = self.dataset[["RiPP"]].to_numpy().ravel()

        label_dict = {'RiPP': 1, 'No_RiPP': 0}

        x_train = torch.tensor(x.values)

        y_train = torch.tensor(pd.Series(y).map(label_dict).tolist())

        training_data = TensorDataset(x_train.to(torch.float32), y_train)

        dataloader = DataLoader(training_data, batch_size=256)
        return dataloader
    def hybrid_model_predict(self,neurnet, svm_classifier, inputs,device):
        neurnet.eval()
        with torch.no_grad():
            features = neurnet(inputs.to(device).float()).cpu().numpy()
        return svm_classifier.predict(features)

    def init_weights(self:Self,layer):
        if type(layer) == nn.Linear or type(layer) == nn.Conv2d:
            nn.init.xavier_uniform_(layer.weight)