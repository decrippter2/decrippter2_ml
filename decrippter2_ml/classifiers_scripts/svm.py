from typing import Self
import pickle
from .classifier_class import BaseClassifier
import logging
from sklearn.svm import SVC

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler()
logger.addHandler(console_handler)

class SVM_Classifier(BaseClassifier):
    """Pydantic_based class to load/retrain SVM classifiers_scripts

    Attributes:
            model_type = str containing SVM kernel, 'poly3' or 'rbf'
    """
    model_type: str

    def load_model(self:Self):
        """Loads pre-trained SVM model"""
        if self.model_type=='poly3' or self.model_type=='POLY3':
            svm_path = self.model_folder.joinpath('svm_poly3.pkl')
        elif self.model_type=='rbf' or self.model_type=='RBF':
            svm_path = self.model_folder.joinpath('svm_rbf.pkl')
        else:
            logger.fatal('Error: invalid kernel')
            raise RuntimeError
        with open(svm_path, 'rb') as f:
            svm_clf = pickle.load(f)
        return svm_clf

    def retrain_model(self:Self):
        """Retrains SVM model"""
        self.load_feature_list()
        x = self.dataset[self.feature_list]
        y = self.dataset[["RiPP"]].to_numpy().ravel()
        if self.model_type == 'poly3' or self.model_type == 'POLY3':
            svm_clf = SVC(C=100, kernel="poly", degree=3, random_state=0, class_weight="balanced")
            svm_clf.fit(x, y)
            svm_path = self.model_folder.joinpath('svm_poly3.pkl')

        elif self.model_type == 'rbf' or self.model_type == 'RBF':
            svm_clf = SVC(C=1000, kernel="rbf", random_state=0, class_weight="balanced")
            svm_clf.fit(x, y)
            svm_path = self.model_folder.joinpath('svm_rbf.pkl')
        else:
            logger.fatal('Error: invalid kernel')
            raise RuntimeError

        with open(svm_path, 'wb') as f:
            pickle.dump(svm_clf, f)

    def predict(self:Self,input):
        """Runs prediction of RiPPs

            Args:
                input: Pandas Dataframe containing the sequence(s) for prediction and extracted features

            Returns:
                list of predictions as 0 (No_RiPP) and 1 (RiPP) corresponding to the indices in input pd"""
        svm_clf=self.load_model()
        input=input[self.feature_list]
        output=svm_clf.predict(input)
        return output