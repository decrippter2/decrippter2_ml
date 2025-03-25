from typing import Self
import logging
import torch
from torch import nn
import torch.optim as optim
import torch.nn.functional as F
from .classifier_class import BaseClassifier

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler()
logger.addHandler(console_handler)

class ExtractorNetwork(nn.Module):
    """Class containing NN architecture for binary classification"""
    def __init__(self, input_size):
        super(ExtractorNetwork, self).__init__()
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
        return x

class NN_Classifier(BaseClassifier):
    """Pydantic_based class to load/retrain NN classifier"""

    def retrain_model(self:Self):
        """Retrains NN model"""
        self.load_feature_list()
        device = self.set_device()
        neurnet = ExtractorNetwork(len(self.feature_list))
        loss_fn = nn.CrossEntropyLoss()
        optimizer = optim.Adam(neurnet.parameters(), lr=0.0001)
        neurnet.apply(self.init_weights)
        dataloader = self.initiate_train_dataloader()

        for epoch in range(500):
            neurnet.train()
            for batch, data in enumerate(dataloader):
                inputs = data[0].to(device).float()
                labels = data[1].to(device).float()

                optimizer.zero_grad()

                # Obtain logits
                logits = neurnet(inputs)

                # Compute loss
                loss = loss_fn(logits, labels.long())

                # Backpropagate through the NN
                loss.backward()
                optimizer.step()

        neurnet.eval()

        torch.save(neurnet, self.model_folder.joinpath('nn.pth'))

    def load_model(self:Self):
        """Loads pre-trained NN"""
        neurnet = torch.load(self.model_folder.joinpath('nn.pth'))
        return neurnet

    def predict(self:Self,input):
        """Runs prediction of RiPPs

        Args:
            input: Pandas Dataframe containing the sequence(s) for prediction and extracted features

        Returns:
            list of predictions as 0 (No_RiPP) and 1 (RiPP) corresponding to the indices in input pd"""
        neurnet=self.load_model()
        neurnet.eval()
        self.load_feature_list()
        input=input[self.feature_list]
        input_dataloader=self.initiate_predict_dataloader(input)
        device=self.set_device()
        pred_list=[]
        with torch.no_grad():
            for data in input_dataloader:
                logits=neurnet(data[0].to(device).float())
                preds=logits.softmax(dim=1).argmax(1)
                pred_list.extend(preds.tolist())
        return pred_list