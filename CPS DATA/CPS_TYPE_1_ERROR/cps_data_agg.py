import pandas as pd
import statsmodels.api as sm
from cps_data_prep import cps_data

def process_cps_data(file_path):

    df = cps_data(file_path)

    X = df[['High School', "Master's Degree", 'Up to Grade 10', 'AGE']] # the covariates used 1st stage of data aggregation
    y = df['INCWAGE']

    X = sm.add_constant(X)

    model = sm.OLS(y, X).fit()

    y_pred = model.predict(X)  # obtaining the predicted values of the model

    residuals = y - y_pred

    df['Residuals'] = residuals

    residuals_mean_by_state_year = df.groupby(['STATEFIP', 'YEAR'])['Residuals'].mean().reset_index()

    residuals_mean_by_state_year 

    return residuals_mean_by_state_year
