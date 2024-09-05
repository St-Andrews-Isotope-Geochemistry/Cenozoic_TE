from pathlib import Path
import os
import pandas as pd
import numpy as np
import seaborn as sns
from sklearn.feature_selection import RFE
from sklearn.linear_model import LogisticRegression  # Replace with your model of choice
import matplotlib.pyplot as plt

## import data
data_path=Path(os.getcwd())/"data"


df=pd.read_csv(data_path/"foram_dataframe_240625_matched_noU.csv")
isotopes=['B11', 'Mg24', 'Al27', 'Na23', 'Li7', 'Ba138',  'U238', 'Cd111', 'Sr88', 'Mn55']




#principal component analysis
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# Standardize the data
scaler = StandardScaler()
X_scaled = scaler.fit_transform(df[isotopes])

# Apply PCA
pca = PCA(n_components=3)
pc_forams = pca.fit_transform(X_scaled)

pc_forams_df=pd.DataFrame(pc_forams, columns=['PC1', 'PC2', 'PC3'])
pc_forams_df['species']=df['species_simple']

#require at least 10 samples per species
species_counts=pc_forams_df['species'].value_counts()
species_counts=species_counts[species_counts>10]
pc_forams_df=pc_forams_df[pc_forams_df['species'].isin(species_counts.index)]

fig, ax = plt.subplots(figsize=(12, 12))   
sns.scatterplot(data=pc_forams_df, x='PC1', y='PC2', hue='species')
fig, ax = plt.subplots(figsize=(12, 12))   
sns.scatterplot(data=pc_forams_df, x='PC2', y='PC3', hue='species')

fig, ax = plt.subplots(figsize=(12, 12))   
sns.scatterplot(data=pc_forams_df, x='PC1', y='PC3', hue='species')


## train test split

from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(df[isotopes], df['species'], test_size=0.2, random_state=42)


## scaling
from sklearn.preprocessing import StandardScaler, MinMaxScaler


# Standardize features (zero mean, unit variance)
scaler = StandardScaler()
X_train_standardized = scaler.fit_transform(X_train)
X_test_standardized = scaler.transform(X_test)


# Apply Min-Max scaling to ensure no negative values
min_max_scaler = MinMaxScaler()
X_train_scaled = min_max_scaler.fit_transform(X_train_standardized)
X_test_scaled = min_max_scaler.transform(X_test_standardized)

# Sample data
X, y = X_train_scaled, y_train

# Initialize your machine learning model
model = LogisticRegression()  # Example with logistic regression

# Initialize RFE with your model and the desired number of features
rfe_selector = RFE(estimator=model, n_features_to_select=7)

# Fit and transform the data
X_new = rfe_selector.fit_transform(X, y)

# Selected features
selected_features_bool = rfe_selector.support_
selected_features=[i for i in range(len(isotopes)) if selected_features_bool[i]]



from sklearn.feature_selection import SelectKBest, chi2


# Initialize SelectKBest with the chi-squared test
k_best_selector = SelectKBest(score_func=chi2, k=7)

# Fit and transform the data
X_new = k_best_selector.fit_transform(X, y)

# Selected features
selected_features = k_best_selector.get_support(indices=True)





X_train_select=X_train_scaled[:, selected_features]
X_test_select=X_test_scaled[:, selected_features]



X_train, X_test, y_train, y_test = train_test_split(pc_forams_df[['PC1', 'PC2', 'PC3']], pc_forams_df['species'], test_size=0.2, random_state=42)
X_train=X_train_standardized
X_test=X_test_standardized

from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import StackingClassifier
from sklearn.metrics import accuracy_score



# Bagging: Random Forest
rf_clf = RandomForestClassifier(n_estimators=100, random_state=42)
rf_clf.fit(X_train, y_train)
rf_pred = rf_clf.predict(X_test)
rf_accuracy = accuracy_score(y_test, rf_pred)
print(f"Accuracy of Random Forest: {rf_accuracy:.2f}")

# AdaBoost
ada_clf = AdaBoostClassifier(n_estimators=100, random_state=42)
ada_clf.fit(X_train, y_train)
ada_pred = ada_clf.predict(X_test)
ada_accuracy = accuracy_score(y_test, ada_pred)
print(f"Accuracy of AdaBoost: {ada_accuracy:.2f}")

# Gradient Boosting
gb_clf = GradientBoostingClassifier(n_estimators=100, random_state=42)
gb_clf.fit(X_train, y_train)
gb_pred = gb_clf.predict(X_test)
gb_accuracy = accuracy_score(y_test, gb_pred)
print(f"Accuracy of Gradient Boosting: {gb_accuracy:.2f}")





#grid search random forest
from sklearn.model_selection import GridSearchCV
params_dt = {'max_depth': [3, 4, 5, 6], 
             'min_samples_leaf':[0.04, 0.06, 0.08], 
             'max_features':[0.2, 0.4, 0.6, 0.8], 
             'n_estimators':[100, 200, 300]}

grid_rf=GridSearchCV(estimator=RandomForestClassifier(random_state=42), 
                     param_grid=params_dt, cv=3, n_jobs=-1, verbose=1)



grid_rf.fit(X_train, y_train)

best_hyperparams = grid_rf.best_params_

best_CV_score = grid_rf.best_score_




## stacking

from sklearn.svm import SVC


# Define base models
base_models = [
    ('rf', RandomForestClassifier(n_estimators=100, random_state=42)),
    ('gb', GradientBoostingClassifier(n_estimators=100, random_state=42)),
    ('svc', SVC(probability=True, random_state=42))
]

# Define meta-model
meta_model = LogisticRegression()

# Create stacking classifier
stacking_clf = StackingClassifier(estimators=base_models, final_estimator=meta_model)

# Train stacking classifier
stacking_clf.fit(X_train, y_train)

# Make predictions
y_pred = stacking_clf.predict(X_test)

# Evaluate the model
accuracy = accuracy_score(y_test, y_pred)
print(f"Accuracy of stacking classifier: {accuracy:.2f}")
