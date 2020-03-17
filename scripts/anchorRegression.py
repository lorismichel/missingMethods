import numpy as np
import math
import random
import matplotlib.pyplot as plt
import copy
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import LinearRegression
from itertools import compress
import sys




class anchorRegression:
    def __init__(self,tradeoff_param):
        self.tradeoff_param = tradeoff_param
        
    def fit(self,X,Y,A):
        PiA = np.dot(np.dot(A, np.linalg.inv(np.dot(np.transpose(A),A))), np.transpose(A))
        H = np.diag(np.ones(np.shape(PiA)[0])) + (np.sqrt(self.tradeoff_param)-1)*PiA
        Ytilda = np.dot(H,np.transpose(Y))
        Xtilda = np.dot(H, X)
        fit = LinearRegression(fit_intercept="false").fit(Xtilda,Ytilda)
        self.coefs = fit.coef_
        
    def predict(self,X):
        Y_hat = np.squeeze(np.asarray(np.dot(X, self.coefs.T)))
        return Y_hat 
    
    
