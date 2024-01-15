import pandas as pd
import numpy as np
import gzip

# Specify the path to the compressed CSV file
file_path = r'C:\Users\Biswajit Palit\Downloads\cps_00006.csv.gz'

def cps_data(file_path):
    # Use Pandas to read the compressed CSV file directly
    # The compression parameter is set to 'gzip'
    df = pd.read_csv(file_path, compression='gzip', header=0)

    # Drop rows where INCWAGE is 99999999
    df = df[(df['INCWAGE'] != 99999999) & (df['INCWAGE'] != 0) & (df['INCWAGE'] != 999)]
    df['INCWAGE'] = np.log(df['INCWAGE'])

    df = df[(df['EDUC'] != 0) & (df['EDUC'] != 1)]

    df = df[(df['YEAR'] >= 1980) & (df['YEAR'] <= 2000)]

    dummy_df = pd.get_dummies(df['YEAR'], prefix='YEAR', drop_first=True)
    df = pd.concat([df, dummy_df], axis=1)

    def categorize_education(educ_code):
        if educ_code <= 10:
            return 'Up to Grade 10'
        elif 10 < educ_code <= 70:
            return 'High School'
        elif 70 < educ_code <= 123:
            return "Master's Degree"
        else:
            return 'Doctorate Degree'

    # Apply the function to create a new 'Education_Category' column
    df['Education_Category'] = df['EDUC'].apply(categorize_education)
    df = pd.get_dummies(df, columns=['Education_Category'], prefix='', prefix_sep='', drop_first=True)

    df = df[~((df['STATEFIP'] > 56) | (df['STATEFIP'] == 11))]

    dummy_df = pd.get_dummies(df['STATEFIP'], prefix='STATEFIP', drop_first=True)
    df = pd.concat([df, dummy_df], axis=1)

    df = df[(df['AGE'] >= 25) & (df['AGE'] <= 50)]

    df = df[df['SEX'] == 2]

    boolean_columns = ['YEAR_1981', 'YEAR_1982',
                       'YEAR_1983', 'YEAR_1984', 'YEAR_1985', 'YEAR_1986', 'YEAR_1987',
                       'YEAR_1988', 'YEAR_1989', 'YEAR_1990', 'YEAR_1991', 'YEAR_1992',
                       'YEAR_1993', 'YEAR_1994', 'YEAR_1995', 'YEAR_1996', 'YEAR_1997',
                       'YEAR_1998', 'YEAR_1999', 'YEAR_2000', 'High School', "Master's Degree",
                       'Up to Grade 10', 'STATEFIP_2', 'STATEFIP_4', 'STATEFIP_5',
                       'STATEFIP_6', 'STATEFIP_8', 'STATEFIP_9', 'STATEFIP_10', 'STATEFIP_12',
                       'STATEFIP_13', 'STATEFIP_15', 'STATEFIP_16', 'STATEFIP_17',
                       'STATEFIP_18', 'STATEFIP_19', 'STATEFIP_20', 'STATEFIP_21',
                       'STATEFIP_22', 'STATEFIP_23', 'STATEFIP_24', 'STATEFIP_25',
                       'STATEFIP_26', 'STATEFIP_27', 'STATEFIP_28', 'STATEFIP_29',
                       'STATEFIP_30', 'STATEFIP_31', 'STATEFIP_32', 'STATEFIP_33',
                       'STATEFIP_34', 'STATEFIP_35', 'STATEFIP_36', 'STATEFIP_37',
                       'STATEFIP_38', 'STATEFIP_39', 'STATEFIP_40', 'STATEFIP_41',
                       'STATEFIP_42', 'STATEFIP_44', 'STATEFIP_45', 'STATEFIP_46',
                       'STATEFIP_47', 'STATEFIP_48', 'STATEFIP_49', 'STATEFIP_50',
                       'STATEFIP_51', 'STATEFIP_53', 'STATEFIP_54', 'STATEFIP_55',
                       'STATEFIP_56']

    # Convert True and False to 1 and 0 in the specified columns
    df[boolean_columns] = df[boolean_columns].astype(int)

    
    
    return df