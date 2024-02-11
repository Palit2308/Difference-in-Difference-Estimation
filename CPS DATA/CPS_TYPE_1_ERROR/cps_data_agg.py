import pandas as pd
import statsmodels.api as sm
from cps_data_prep import (cps_data)

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

    dummy_df = pd.get_dummies(residuals_mean_by_state_year['STATEFIP'], prefix='STATEFIP', drop_first=True)

    residuals_mean_by_state_year1 = pd.concat([residuals_mean_by_state_year, dummy_df], axis=1)   # concatenating the dummy variables with the original DataFrame

    dummy_df2 = pd.get_dummies(residuals_mean_by_state_year1['YEAR'], prefix='YEAR', drop_first=True)


    residuals_mean_by_state_year1 = pd.concat([residuals_mean_by_state_year1, dummy_df2], axis=1)  # concatenating the dummy variables with the original DataFrame

    boolean_columns = [ 'STATEFIP_2',
        'STATEFIP_4', 'STATEFIP_5', 'STATEFIP_6', 'STATEFIP_8', 'STATEFIP_9',
        'STATEFIP_10', 'STATEFIP_12', 'STATEFIP_13', 'STATEFIP_15',
        'STATEFIP_16', 'STATEFIP_17', 'STATEFIP_18', 'STATEFIP_19',
        'STATEFIP_20', 'STATEFIP_21', 'STATEFIP_22', 'STATEFIP_23',
        'STATEFIP_24', 'STATEFIP_25', 'STATEFIP_26', 'STATEFIP_27',
        'STATEFIP_28', 'STATEFIP_29', 'STATEFIP_30', 'STATEFIP_31',
        'STATEFIP_32', 'STATEFIP_33', 'STATEFIP_34', 'STATEFIP_35',
        'STATEFIP_36', 'STATEFIP_37', 'STATEFIP_38', 'STATEFIP_39',
        'STATEFIP_40', 'STATEFIP_41', 'STATEFIP_42', 'STATEFIP_44',
        'STATEFIP_45', 'STATEFIP_46', 'STATEFIP_47', 'STATEFIP_48',
        'STATEFIP_49', 'STATEFIP_50', 'STATEFIP_51', 'STATEFIP_53',
        'STATEFIP_54', 'STATEFIP_55', 'STATEFIP_56', 'YEAR_1981', 'YEAR_1982',
        'YEAR_1983', 'YEAR_1984', 'YEAR_1985', 'YEAR_1986', 'YEAR_1987',
        'YEAR_1988', 'YEAR_1989', 'YEAR_1990', 'YEAR_1991', 'YEAR_1992',
        'YEAR_1993', 'YEAR_1994', 'YEAR_1995', 'YEAR_1996', 'YEAR_1997',
        'YEAR_1998', 'YEAR_1999', 'YEAR_2000']

    
    residuals_mean_by_state_year1[boolean_columns] = residuals_mean_by_state_year1[boolean_columns].astype(int)  # converting True and False to 1 and 0 in the specified columns

    return residuals_mean_by_state_year1
