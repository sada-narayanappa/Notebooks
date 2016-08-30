import matplotlib.pyplot as plt
from numpy import *
from collections import Counter
import numpy as np
import pylab as pl
from matplotlib.colors import ListedColormap
from sklearn import neighbors, datasets, cluster, preprocessing, decomposition, svm, datasets
from sklearn.decomposition import PCA
import pandas as pd
from pandas.tools.plotting import scatter_matrix
import numpy.random as random
from mpl_toolkits.mplot3d import Axes3D
import glob
import os
import time
from IPython.display import display
from IPython.display import Image

import scipy;
from scipy	import stats;
import sklearn;
import sklearn.ensemble;
import sklearn.neighbors
from sklearn.preprocessing import StandardScaler
from sklearn import metrics
from sklearn.metrics import *
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import label_binarize
from sklearn.cross_validation import train_test_split
from sklearn.multiclass import OneVsRestClassifier

import os
import subprocess
from sklearn.tree import DecisionTreeClassifier, export_graphviz
from sklearn import tree
from sklearn import datasets
from IPython.display import Image
#import seaborn as sns
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pydotplus
from sklearn.externals.six import StringIO

# -*- coding: utf-8 -*-
def run_cv(X,y,clf_class,printDebug = False , clf=None):
    # Construct a kfolds object
    kf = sklearn.cross_validation.KFold(len(y),n_folds=5,shuffle=True)
    y_pred = y.copy()
    
    # Iterate through folds\
    i = 0;
    for train_index, test_index in kf:
        X_train, X_test = X[train_index], X[test_index]
        y_train = y[train_index]
        # Initialize a classifier with key word arguments
        clf = clf_class(**kwargs) if (clf is None)  else clf;
        if (printDebug): print ("*",i, end ="");
        clf.fit(X_train,y_train)
        y_pred[test_index] = clf.predict(X_test)
        i = i +1;
    if (printDebug): print ("*");
    return y_pred, clf

# -*- coding: utf-8 -*-
def run_cvTT(X,y,clf_class,printDebug = True , clf=None):
    X_train, X_test, y_train, y_test = cross_validation.train_test_split( X, y, test_size=0.2, random_state=0);
            
    
def accuracy(y_true,y_pred):
    # NumPy interprets True and False as 1. and 0.
    return np.mean(y_true == y_pred)

#
# Call:
# cms = [("Decision Tree", [[25,24],[23,22]])]
#
def draw_confusion_matrices(confusion_matricies,class_names):
    class_names = class_names.tolist() if type(class_names) != list else class_names
    fig = plt.figure(figsize = (20,5))
    for i,cm in enumerate(confusion_matricies):
        classifierName, matrix = cm[0], cm[1]
        #cmstr = str(cm)

        ax = fig.add_subplot(1,8,i+1)
        plt.subplots_adjust(wspace = .4);
        cax = ax.matshow(matrix, cmap='seismic', interpolation='nearest')
        #plt.title('CM: %s' % classifier + "\n" + cmstr)
        plt.title(classifierName)
        plt.grid(None)
        if (i ==0 ):
            i=0;
            #fig.colorbar(cax);
            
            
        ax.set_xticklabels([''] + class_names)
        ax.set_yticklabels([''] + class_names)
        plt.xlabel('Predicted')
        plt.ylabel('True')
        for (ii, jj), z in np.ndenumerate(matrix):
            ax.text(jj, ii, '{:0.1f}'.format(z), ha='center', va='center',
                    bbox=dict(facecolor='white', edgecolor='0.3'))
    plt.show()


##
# Draw most 15 significant 
#
def DrawFeatureImportance(dft,clf, title="", ax =None, m=15):
    if ( not hasattr(clf,'feature_importances_') ):
        print ("No Feature Importance matrix for this classifier:", clf)
        return;
        
    # Get Feature Importance from the classifier
    feature_importance = clf.feature_importances_
    # Normalize The Features
    feature_importance = 100.0 * (feature_importance / feature_importance.max())
    sorted_idx = np.argsort(feature_importance)
    sorted_idx10 = sorted_idx[-m:]; #[0:5]  # TODO: Print Top 10 only ??
    pos = np.arange(sorted_idx10.shape[0]) + .5
    fc10=np.asanyarray(dft.columns.tolist())[sorted_idx10];
    
    if ( ax == None):
        plt.figure(figsize=(5, 3));
    plt.barh(pos, feature_importance[sorted_idx10], align='center', color='#7A68A6')
    plt.yticks(pos, fc10)
    plt.xlabel('Relative: '+ title)
    plt.title('Variable Importance')
    if ( ax == None): plt.show()

def DrawFeatureImportanceMatrix(df, clfs):
    fig = plt.figure(figsize = (25,3))
    plt.subplots_adjust(wspace = 1.2);
    for i in range(int (len(clfs)/2)) :
        classifierName, clf = clfs[i*2], clfs[i*2+1]
        if ( not hasattr(clf,'feature_importances_') ):
            continue;
        #cmstr = str(cm)
        #print classifierName;
        plt.subplots_adjust(wspace = 1.4);
        ax = fig.add_subplot(1,8,i+1)
        DrawFeatureImportance(df, clf, classifierName, ax)
        #plt.title(classifierName)
            
    plt.show()
    
# This function will take the categorical value or any column and
# either replaces inline or create a new column with enumeration of 
# Values for a left column will be changed as shown to the right array:
# For ex. [a,a,a,b,b,c,c,c,c,d, d] => [0,0,0,1,1,2,2,2,2,3, 3]
#
def encodeCategorical(df, columnName, newColumnName = None, makeCopy = False):
    df_mod = df.copy() if makeCopy else df;
    targets = df_mod[columnName].unique()
    mapToInt = {name: n for n, name in enumerate(targets)}
    newColumnName = newColumnName if (newColumnName!=None) else columnName;
    df_mod[newColumnName] = df_mod[columnName].replace(mapToInt)

    return (df_mod, targets, mapToInt)
    
#Classification Problems 
# Usually in case of classifcation, it is best to draw scatter plot of 
# Target Varible using Scatter plot
# df, t,m = encodeCategorical(dfL, "FiveTile1", "Target" );
# scatter_matrix(dfL, alpha=1, figsize=(10,10), s=100, c=df.Target);
# print "categorical Plot {}".format(m)
#
#
# Df - Data Frame that does not have a Predict Column
# y  - Predict Column Series 
def Classify(df, y, 
             printDebug = True ,
             drawConfusionMatrix = True,
             classifiers = None,
             scale =True
             ):
    if ( df is None or y is None):
        raise Exception("No Data Given");

    t = df.select_dtypes(exclude=[np.number])
    if ( len(t.columns) > 0) :
        raise Exception("nonnumeric columns? "  + t.columns);

    l = preprocessing.LabelEncoder()
    class_names = y.unique()
    y=l.fit_transform(y);

    df.fillna(0, inplace=True)
    X = df.as_matrix().astype(np.float)

    if (scale):
        scaler = StandardScaler()
        X = scaler.fit_transform(X)
      
    print ("Feature space holds %d observations and %d features" % X.shape)
    print ("Unique target labels:", class_names)
         
    cls = [# Note: SVM takes long time to run - so commenting out
           #"SVM"               , sklearn.svm.SVC(), 
           "Random Forest"     , sklearn.ensemble.RandomForestClassifier(),
           #"K-NN"              , sklearn.neighbors.KNeighborsClassifier(),
           "DecisionTree Gini" , DecisionTreeClassifier(max_depth=4, criterion="gini"), 
           "DecisionTree Entr" , DecisionTreeClassifier(max_depth=4, criterion="entropy"), 
           "Gradient Boosting" , sklearn.ensemble.RandomForestClassifier(),
           #"Logit Regression"  , sklearn.linear_model.LogisticRegression()
           ];

    if (not classifiers is None ):
        cls = classifiers;

    y_preds = {}
    ret_accuracy = [];
    cms = [];
    clfs = {}
    for i in arange( int (len(cls)/2) ):
        nm = cls[i*2];
        cl = cls[i*2 +1]
        y_pred, clfi = run_cv(X,y, None, clf=cl, printDebug=printDebug)
        y_preds[nm] = y_pred
        clfs[nm] = clfi
        ac  = accuracy(y, y_pred);
        cm = confusion_matrix(y, y_pred )
        ret_accuracy.append( (nm, ac, cm) )
        if (printDebug): 
            print ("%20s accuracy: %03f "% (nm, ac) );
            #print('{}\n'.format(metrics.classification_report(y, y_pred)))
            #print("%20s r^2 score: %03f"% (nm,sklearn.metrics.r2_score(y, y_pred, sample_weight=None, multioutput=None)))
            print("%20s r^2 score: %03f"% (nm,sklearn.metrics.r2_score(y, y_pred, sample_weight=None)))
        cms.append( (nm, cm) );
    if (drawConfusionMatrix): 
        #print cms, class_names
        draw_confusion_matrices(cms, class_names);
        DrawFeatureImportanceMatrix(df, cls)

    return (X, y, ret_accuracy,cls, y_preds,clfs);

def visualizeTree(dcls, feature_names, class_names= None):
    dot_data = StringIO()  
    tree.export_graphviz(dcls, out_file=dot_data,  
                         feature_names= feature_names,  
                         class_names= class_names,  
                         filled=True, rounded=True,  
                         special_characters=True)  
    graph = pydotplus.graph_from_dot_data(dot_data.getvalue())
    display(Image(graph.create_png()))
    return "";

##======================= DRAW Decision Trees here 
def DrawDecisionTree(X,y, cls, class_names=None):
    imgs=[]
    
    if ( class_names is None):
        class_names = y.unique().astype(str);
        class_names.sort()
    if (str(type(cls)).find('DecisionTreeClassifier') > 0 ):
        visualizeTree(cls, X.columns, class_names=class_names)
    else:
        for k in range(int( len(cls)/2) ) :
            dcls = cls[k*2+1];
            if (str(type(dcls)).find('DecisionTreeClassifier') > 0):
                visualizeTree(dcls, X.columns, class_names=class_names)

    
#======================== Get COde For DecisionTree        
# Here is how to use this code
#if __name__ == '__main__':
#    print("\n-- get data:")
#    df = dfL;
#
#    print("\n-- df.head():")
#    print(df.head(), end="\n\n")
#
#    features = ["SepalLength", "SepalWidth", "PetalLength", "PetalWidth"]
#    df, targets,mi = encodeCategorical(df, "Name", "Target")
#    y = df["Target"]
#    X = df[features]
#
#    dt = DecisionTreeClassifier(min_samples_split=20, random_state=99)
#    dt.fit(X, y)
#
#    print("\n-- get_code:")
##    get_code(dt, features, targets)
#
#    print("\n-- look back at original data using pandas")
#    print("-- df[df['PetalLength'] <= 2.45]]['Name'].unique(): ",
#          df[df['PetalLength'] <= 2.45]['Name'].unique(), end="\n\n")
#
#    visualizeTree(dt, features)

def getCodeOfDecisionTree(tree, feature_names, target_names, spacer_base="    "):
    """Produce psuedo-code for decision tree.
    
    Args
    ----
    tree -- scikit-leant DescisionTree.
    feature_names -- list of feature names.
    target_names -- list of target (class) names.
    spacer_base -- used for spacing code (default: "    ").

    Notes
    -----
    based on http://stackoverflow.com/a/30104792.
    """
    left      = tree.tree_.children_left
    right     = tree.tree_.children_right
    threshold = tree.tree_.threshold
    features  = [feature_names[i] for i in tree.tree_.feature]
    value = tree.tree_.value
   
    def recurse(left, right, threshold, features, node, depth):
        spacer = spacer_base * depth
        if (threshold[node] != -2):
            print(spacer + "if ( " + features[node] + " <= " + \
                  str(threshold[node]) + " ) {")
            if left[node] != -1:
                    recurse (left, right, threshold, features, left[node],
                            depth+1)
            print(spacer + "}\n" + spacer +"else {")
            if right[node] != -1:
                    recurse (left, right, threshold, features, right[node],
                             depth+1)
            print(spacer + "}")
        else:
            target = value[node]
            for i, v in zip(np.nonzero(target)[1], target[np.nonzero(target)]):
                target_name = target_names[i]
                target_count = int(v)
                print(spacer + "return " + str(target_name) + " ( " + \
                      str(target_count) + " examples )")
    
    recurse(left, right, threshold, features, 0, 0)


