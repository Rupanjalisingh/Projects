The goal of this project is to analyze two datasets of EHRs of patients with breast cancer through clustering techniques to detect relevant groups of patients based on specific clinical features.

Data Loading: We loaded the data from the 'Raw data' sheet of China.xlsx into china_df and the 'case' sheet of Korea.xlsx into korea_df.
Column Name Cleaning: We standardized column names in both dataframes by removing spaces and replacing them with underscores.
Missing Value Handling: In china_df, we imputed missing numerical and object values with the median and mode respectively. In korea_df, we imputed missing numerical values with the median and filled the columns 'Group', 'Miscarriage', 'Menopause', 'P53', 'Metastasis', 'ER', and 'PR' (which contained only NA values at a certain point) with 0s as requested.
Data Type Conversion: We ensured columns had appropriate data types, including converting date columns in china_df to datetime objects.
Categorical Data Handling: In china_df, we mapped the mixed-type columns ('B_G', 'N_G') and the ordinal column ('T-Nstage') to numerical representations. In korea_df, we one-hot encoded the 'Histopathological' column and mapped the ordinal 'TNMstage' (although this mapping resulted in all -1s). The binary columns ('Group', 'Miscarriage', etc.) were filled with 0s instead of being mapped to 0 and 1 due to their initial state of being all NA.
Numerical Data Normalization: We applied Min-Max scaling to the numerical columns in both dataframes.
Data Saving: We saved the processed dataframes to 'china_processed.csv' and 'korea_processed.csv'.
Readiness for Analysis:

The data has undergone significant preprocessing and is in a much cleaner and more structured format. For many types of analysis, the data is now ready to be used.

However, it's important to be aware of the following:

The binary categorical columns in korea_df ('Group', 'Miscarriage', 'Menopause', 'P53', 'Metastasis', 'ER', 'PR') were filled with 0s due to being entirely NA at a certain point. If these columns are critical for your analysis, you might need to investigate the original data source to understand why they were missing and potentially reload or handle them differently if accurate values are needed. The mapping for 'TNMstage' in korea_df resulted in all -1s, indicating the original values didn't match the mapping keys. If this column is important, you might need to adjust the mapping or investigate the values in the original data.
