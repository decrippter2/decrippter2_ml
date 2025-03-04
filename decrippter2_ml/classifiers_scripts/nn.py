from typing import Self
import logging
import torch
from torch import nn
import torch.optim as optim
import torch.nn.functional as F
from decrippter2_ml.classifiers_scripts.classifier_class import BaseClassifier
from torch.utils.data import DataLoader

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
        dataloader = self.initiate_dataloader()

        for epoch in range(500):
            neurnet.train()
            for batch, data in enumerate(dataloader):
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

    def predict(self:Self,input):
        neurnet=self.load_model()
        neurnet.eval()
        input=input[self.feature_list]
        input_dataloader=DataLoader(input,batch_size=128)
        device=self.set_device()
        pred_list=[]
        with torch.no_grad():
            for data in input_dataloader:
                logits=neurnet(data.to(device).float())
                preds=logits.softmax(dim=1).argmax(1)
                pred_list.extend(preds)
        return pred_list