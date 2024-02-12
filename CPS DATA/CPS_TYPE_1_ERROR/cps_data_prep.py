import pandas as pd
import numpy as np
import gzip

file_path = r'C:\Users\Biswajit Palit\Downloads\cps_00006.csv.gz' # path for the data

def cps_data(file_path):
  
    df = pd.read_csv(file_path, compression='gzip', header=0) # importing the data

    
    df = df[(df['INCWAGE'] != 99999999) & (df['INCWAGE'] != 0) & (df['INCWAGE'] != 999)]  # dropping the rows containing invalid INCWAGE

    df['INCWAGE'] = np.log(df['INCWAGE'])  # taking log of weekly earnings which will be our dependent variable in the regressions


    df = df[(df['EDUC'] != 0) & (df['EDUC'] != 1)]  # dropping education levels 0 and 1 as they are invalid entries as per the labels of the dataset

    df = df[(df['YEAR'] >= 1980) & (df['YEAR'] <= 2000)]  # taking the time frame from 1980 to 2000

    def categorize_education(educ_code): # creating a fucntion to categorize education levels
        if educ_code <= 10:
            return 'Up to Grade 10'
        elif 10 < educ_code <= 70:
            return 'High School'
        elif 70 < educ_code <= 123:
            return "Master's Degree"
        else:
            return 'Doctorate Degree'

    
    df['Education_Category'] = df['EDUC'].apply(categorize_education) # applying the function to create a new 'Education_Category' column

    df = pd.get_dummies(df, columns=['Education_Category'], prefix='', prefix_sep='', drop_first=True)
    boolean = ['Up to Grade 10', 'High School', "Master's Degree"]
    df[boolean] = df[boolean].astype(int)


    df = df[~((df['STATEFIP'] > 56) | (df['STATEFIP'] == 11))]  # taking only the 50 states of the States and exclusing the regions as per the labels of the dataset

    df = df[(df['AGE'] >= 25) & (df['AGE'] <= 50)]  # taking the age group from 25 to 50

    df = df[df['SEX'] == 2] # taking only female respondents

    return df