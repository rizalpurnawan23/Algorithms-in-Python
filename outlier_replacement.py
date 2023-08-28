# REPLACING OUTLIERS:
class ReplaceOutliers:
    """
    Title       : ReplaceOutliers
    Type        : Python Class
    Author      : Rizal Purnawan
    Date        : 08/28/2023

    About the class
    ---------------
    This class is used to replace outliers from a dataset in a form
    of pandas dataframe.
    """

    # ---------------------------------------------------------------
    # Initialisation:
    def __init__(self, df):
        self.df = df.copy()

        # Necessary library:
        import pandas as pd
        import numpy as np
        import matplotlib.pyplot as plt
        import seaborn as sns

        self.__pd = pd
        self.__np = np
        self.__plt = plt
        self.__sns = sns

    # ---------------------------------------------------------------
    # Function for identifying outliers:
    def is_outlier(self, col, x):
        """
        Notes:
        ------
        This function takes a column of the dataframe and a numeric
        value then returns a boolean value.

        col         : The corresponding column in the dataframe.
        x           : A numeric value.
        """
        df = self.df.copy()
        df_col = df[col]
        Q1, Q3 = df_col.describe()["25%"], df_col.describe()["75%"]
        IQR = Q3 - Q1
        upper_bound = Q3 + 1.5 *IQR
        lower_bound = Q1 - 1.5 *IQR
        if not (lower_bound <= x <= upper_bound):
            return True
        else:
            return False

    # ---------------------------------------------------------------
    # Function for replacing outliers:
    def replace(self, col, replace_with= "mean", cleaned= True):
        """
        Notes:
        ------
        This function takes a column of the dataframe, a parametric
        string and a boolean value then returns a process pandas
        dataframe consisting of the replaced outliers.

        replace_with    : A string of either 'mean', 'mode',
                          'median' or 'delete' which indicates the
                          method used to replace the outliers.
        cleaned         : A boolean value. If True, then the function
                          returns the dataframe cleaned from
                          outliers. Otherwise, the function returs
                          a list of outliers.
        """
        df = self.df.copy()
        df_col = df[col]
        cleaned_list = list()
        outliers = list()
        if replace_with != 'delete':
            if replace_with == "mean":
                subs = df_col.mean()
            elif replace_with == "mode":
                subs = df_col.mode().min()
            elif replace_with == "median":
                subs = df_col.describe()["50%"]
            else:
                raise ValueError
            for x in df_col:
                if self.is_outlier(col, x):
                    cleaned_list.append(subs)
                    outliers.append(x)
                else:
                    cleaned_list.append(x)
            df[col] = cleaned_list
        else:
            ndf = df.copy()
            ndf["outlier"] = [
                True if self.is_outlier(col, x) else False
                for x in df_col
                ]
            out_df = ndf[(ndf["outlier"] == True)]
            ndf = ndf.drop(index= out_df.index)
            ndf = ndf.drop(columns= ["outlier"])
            df = ndf
            outliers = [x for x in df_col if self.is_outlier(col, x)]
        if cleaned == True:
            return df
        else:
            return outliers
