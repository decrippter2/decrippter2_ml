import json
from sys import argv
import random
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from typing import Self, Any
import seaborn as sns
from sklearn import svm
from pydantic import BaseModel
import pickle
from pathlib import Path

class SVM_Classifier(BaseModel):
    """Pydantic_based class to load/retrain SVM classifiers_scripts

    Attributes:
            folder_path=path to JSON with RiPP entries
    """
    dataset: Any
    model_folder: Path = Path(__file__).parent.parent.joinpath('pretrained_models')
    feature_path: Path = Path(__file__).parent.parent.joinpath("feature_extraction_manager/feature_list.json")
    feature_list: list =[]

    def load_feature_list(self:Self):
        with open(self.feature_path) as feature_json:
            self.feature_list=json.load(feature_json)["feature_list"]
    def load_poly_svm(self:Self):
        poly_svm=self.model_folder.joinpath('svm_poly3.pkl')
        with open(poly_svm,'rb') as f:
            poly_clf=pickle.load(f)
        return poly_clf

    def load_rbf_svm(self:Self):
        rbf_svm=self.model_folder.joinpath('svm_rbf.pkl')
        with open(rbf_svm,'rb') as f:
            rbf_clf=pickle.load(f)
        return rbf_clf

    def retrain_poly_svm(self:Self):
        self.load_feature_list()
        x = self.dataset[self.feature_list]
        y = self.dataset[["RiPP"]].to_numpy().ravel()
        clf = svm.SVC(C=100, kernel="poly", degree=3, random_state=0, class_weight="balanced")
        clf.fit(x, y)
        poly_svm = self.model_folder.joinpath('svm_poly3.pkl')
        print(poly_svm)
        with open(poly_svm, 'wb') as f:
            pickle.dump(clf, f)

    def retrain_rbf_svm(self:Self):
        self.load_feature_list()
        x = self.dataset[self.feature_list]
        y = self.dataset[["RiPP"]].to_numpy().ravel()
        clf = svm.SVC(C=1000,kernel="rbf",random_state=0,class_weight="balanced")
        clf.fit(x, y)
        rbf_svm=self.model_folder.joinpath('svm_rbf.pkl')
        with open(rbf_svm, 'wb') as f:
            pickle.dump(clf, f)

