# Bismillahirrahmanirrahim
# Date of creation: 06/15/2023

"""
    Title       : cpt_forecast
    Type        : Python module
    Author      : Rizal Purnawan
    Date        : 06/15/2023

    About the module
    ----------------

    This module contains a Python class for forecasting Cone
    Penetration Test (CPT) data.
"""

class CPTForecast:
    # INITIALISATION-------------------------------------------------
    def __init__(self, cpt):
        self.cpt = cpt
        
        # Python standard statistics and math module
        from statistics import mean
        from math import ceil, sqrt
        self.__mean = mean
        self.__ceil = ceil
        self.__sqrt = sqrt

        # Additional libraries:
        import numpy as np
        import pandas as pd
        import matplotlib.pyplot as plt
        import seaborn as sns

        from sklearn.linear_model import LinearRegression
        from sklearn.metrics import mean_absolute_error, mean_squared_error
        from xgboost import XGBRegressor
        from statsmodels.tsa.deterministic import DeterministicProcess
        from statsmodels.graphics.tsaplots import plot_pacf

        self.__np = np
        self.__pd = pd
        self.__plt = plt
        self.__sns = sns
        self.__linear_regression = LinearRegression
        self.__xgbregressor = XGBRegressor
        self.__mae = mean_absolute_error
        self.__mse = mean_squared_error
        self.__dp = DeterministicProcess
        self.__plot_pacf = plot_pacf
        
    # SETTING MOVING AVERAGE TREND-----------------------------------
    def moving_average(self, window):
        """
        Title       : moving_average
        Type        : Class Function
        Author      : Rizal Purnawan
        Date        : 06/15/2023

        About the function
        ------------------
        
        The argument 'window' is a positive integer. It represents
        the window of the moving average trends.
        """
        pd = self.__pd
        # np = self.__np
        # plt = self.__plt
        # sns = self.__sns
        # LinearRegression = self.__linear_regression
        # XGBRegressor = self.__xgbregressor
        # mean_absolute_error = self.__mae
        # mean_squared_error = self.__mse
        # DeterministicProcess = self.__dp
        # plot_pacf = self.__plot_pacf

        cpt_df = self.cpt.copy()
        ma_trends_dict = dict()
        ceil = self.__ceil
        for col in cpt_df.columns:
            y = cpt_df[col].squeeze()
            y_ma = y.rolling(window= window,
                             min_periods= ceil(window / 2),
                             center= True).mean()
            ma_trends_dict[col] = y_ma
        ma_trends = pd.DataFrame(ma_trends_dict, index= cpt_df.index)
        # Instantiating the result
        self.ma_trends = ma_trends
        print(">>> Moving average trends are set ;)")
        
    # SETTING MAJOR TRENDS-------------------------------------------
    def set_major_trends(self, max_ord= 10):
        """
        Title       : set_major_trends
        Type        : Class Function
        Author      : Rizal Purnawan
        Date        : 06/15/2023
        
        About the function
        ------------------
        
        With this method we will approximate the major trends by
        making use of polynomial functions. The argument 'max_ord'
        refers to the maximum possible order of the polynomial. The
        chosen order of the polynomial will be determined by the
        algorithm.
        """
        pd = self.__pd
        # np = self.__np
        # plt = self.__plt
        # sns = self.__sns
        LinearRegression = self.__linear_regression
        # XGBRegressor = self.__xgbregressor
        # mean_absolute_error = self.__mae
        # mean_squared_error = self.__mse
        DeterministicProcess = self.__dp
        # plot_pacf = self.__plot_pacf

        cpt = self.cpt.copy()
        ma_trends = self.ma_trends.copy()
        orders = list()
        dp_of_major_trends = list()
        models_of_major_trends = list()
        major_dict = dict()
        for col in ma_trends.columns:
            mae_ord = list()
            for ord in range(max_ord):
                dp = DeterministicProcess(
                    index= cpt.index,
                    constant= True,
                    order= ord
                )
                X = dp.in_sample()
                y = cpt[col].squeeze()
                model = LinearRegression(fit_intercept= False)
                _ = model.fit(X, y)
                y_maj = model.predict(X)
                mean = self.__mean
                mae = mean(
                    [abs(a - b) for a, b in zip(list(ma_trends[col]),
                                                list(y_maj)
                                               )
                    ]
                )
                mae_ord.append((ord, mae))
            ord_ = [a for a, b in mae_ord
                    if b == min([m[1] for m in mae_ord])][0]
            orders.append(ord_)
            
            dp = DeterministicProcess(
                index= range(len(cpt)), constant= True, order= ord_
            )
            dp_of_major_trends.append(dp)
            X = dp.in_sample()
            X.index = cpt.index
            y = cpt[col].squeeze()
            model = LinearRegression(fit_intercept= False)
            _ = model.fit(X, y)
            y_maj = model.predict(X)
            major_dict[col] = list(y_maj)
            models_of_major_trends.append(model)
        major_trends = pd.DataFrame(major_dict, index= cpt.index)
        # Instantiating the results
        self.orders = orders
        self.__dp_of_major_trends = dp_of_major_trends
        self.__models_of_major_trends = models_of_major_trends
        self.major_trends = major_trends
        print(">>> Major trends are set ;)")
        
    # SETTING NOISES-------------------------------------------------
    def set_noises(self):
        """
        Title       : set_noises
        Type        : Class Function
        Author      : Rizal Purnawan
        Date        : 06/15/2023
        
        About the function
        ------------------
        
        Noises are set by detrending the CPT data with the major
        trends that we set earlier. That is, we subtract the CPT data
        with the major trends.
        """
        pd = self.__pd
        # np = self.__np
        # plt = self.__plt
        # sns = self.__sns
        # LinearRegression = self.__linear_regression
        # XGBRegressor = self.__xgbregressor
        # mean_absolute_error = self.__mae
        # mean_squared_error = self.__mse
        # DeterministicProcess = self.__dp
        # plot_pacf = self.__plot_pacf

        cpt = self.cpt.copy()
        major_trends = self.major_trends.copy()
        noises_dict = dict()
        for col in cpt.columns:
            noises_dict[col] = list(cpt[col] - major_trends[col])
        noises = pd.DataFrame(noises_dict, index= cpt.index)
        # Instantiating the result
        self.noises = noises
        print(">>> Noises are set")
        
    # SHIFTING THE NOISES FOR FORECAST-------------------------------
    def shift_noises(self, num_fore= 5):
        """
        Title       : shift_noises
        Type        : Class Function
        Author      : Rizal Purnawan
        Date        : 06/15/2023
        
        About the function
        ------------------
        
        With this method, we shift the noises for forecast. The
        argument 'num_fore' is the number of datapoints in the
        forecast horizon.
        """
        shifted_noises = self.noises.copy()
        max_id = shifted_noises.index[-1]
        index = list(shifted_noises.index)[num_fore:] \
                + [max_id + 0.2*(n + 1) for n in range(num_fore)]
        shifted_noises.index = index
        # Instantiating the results
        self.__num_fore = num_fore
        self.shifted_noises = shifted_noises
        print(">>> Noises are shifted ;)")
        
    # PACF OF SHIFTED NOISES-----------------------------------------
    def shifted_noises_pacfs(self, max_lags= 15):
        """
        Title       : shifted_noises_pacfs
        Type        : Class Function
        Author      : Rizal Purnawan
        Date        : 06/15/2023
        
        About the function
        ------------------
        
        This method helps us determine the partial autocorrelation
        functions (PACFs) of the shifted noises that we have set
        earlier.
        """
        plot_pacf = self.__plot_pacf
        shifted_noises = self.shifted_noises.copy()
        for col in shifted_noises.columns:
            plot_pacf(
                shifted_noises[col],
                lags= max_lags,
                title= f"PACF of {col}"
            );
            
    # SETTING SHIFTED NOISES LAGS------------------------------------
    def set_shifted_noises_lags(self, lags_list):
        """
        Title       : set_shifted_noises_lags
        Type        : Class Function
        Author      : Rizal Purnawan
        Date        : 06/15/2023
        
        About the function
        ------------------
        
        This method helps us determine the lagged features for the
        noises forecast. The argument 'lags_list' is a list of lists
        of integers. Each list in the bigger list consists of the
        numbers of lags for each column which are determined from the
        PACFs of the shifted noises.
        """
        shifted_noises = self.shifted_noises.copy()
        lagged_shifted_noises = list()
        for col, lags in zip(shifted_noises.columns, lags_list):
            df = shifted_noises.copy()[col].to_frame()
            df.columns = ["no_lag"]
            for lag in lags:
                df[str(lag)] = df["no_lag"].shift(lag)
                df = df.fillna(0.0)
            lagged_shifted_noises.append(df)
        # Instantiating the result
        self.lagged_shifted_noises = lagged_shifted_noises
        print(">>> Shifted noises have been lagged ;)")
        
    # THE FORECAST ALGORITHM-----------------------------------------
    def forecast(self, noises_models= None, show_validation= False):
        """
        Title       : forecast
        Type        : Class Function
        Author      : Rizal Purnawan
        Date        : 06/15/2023
        
        About the function
        ------------------
        
        This method provides us with the forecast algorithm. The
        argument 'noises_models', if set, has to be a list of chosen
        regression algorithm (such as 'LinearRegression' or
        'XGBRegressor'). By default, we set it 'None', and then we
        reset it into a list of all 'LinearRegression's. The argument
        'show_validation' is by default set 'False'. If it is set
        'True' then the method returns a data frame consisting of the
        mean absolute error values of the actual values compared to
        the predicted values. The data frame gives as a confidence on
        how the forecast works.
        """
        pd = self.__pd
        LinearRegression = self.__linear_regression

        dp_of_major_trends = self.__dp_of_major_trends
        models_of_major_trends = self.__models_of_major_trends
        
        cpt = self.cpt.copy()
        noises = self.noises.copy()
        shifted_noises = self.shifted_noises.copy()
        lagged_shifted_noises = [
            lsn.copy() for lsn in self.lagged_shifted_noises
        ]
        
        num_fore = self.__num_fore
        max_id = cpt.index[-1]
        columns = cpt.copy()
        pred_index = list(cpt.index)[num_fore:]
        fore_index = [max_id + (n + 1)*0.2 for n in range(num_fore)]
        pred_dict = dict()
        fore_dict = dict()
        
        if noises_models is None:
            noises_models = [LinearRegression] * len(noises.columns)
        
        for k, col in zip(range(len(columns)), columns):
            # Major trends
            major_dp = dp_of_major_trends[k]
            major_model = models_of_major_trends[k]
            X_major = major_dp.in_sample()
            X_major.index = range(len(X_major))
            X_major = X_major[num_fore:]
            X_major.index = pred_index
            X_major_fore = major_dp.out_of_sample(steps= num_fore)
            X_major_fore.index = fore_index
            y_major_pred = major_model.predict(X_major)
            y_major_fore = major_model.predict(X_major_fore)
            
            # Noises
            X_noise = lagged_shifted_noises[k]
            X_noise.index = range(len(X_noise))
            X_noise_fore = X_noise[-num_fore:]
            X_noise = X_noise[:-num_fore]
            X_noise.index = pred_index
            X_noise_fore.index = fore_index
            y_noise = list(noises[col])[num_fore:]
            noise_model = noises_models[k]()
            _ = noise_model.fit(X_noise, y_noise)
            y_noise_pred = noise_model.predict(X_noise)
            y_noise_fore = noise_model.predict(X_noise_fore)
            
            # Combined
            y_pred = [a + b for a, b in zip(list(y_major_pred),
                                            list(y_noise_pred)
                                           )
                     ]
            y_fore = [a + b for a, b in zip(list(y_major_fore),
                                            list(y_noise_fore)
                                           )
                     ]
            pred_dict[col] = y_pred
            fore_dict[col] = y_fore
            
        prediction_df = pd.DataFrame(pred_dict, index= pred_index)
        forecast_df = pd.DataFrame(fore_dict, index= fore_index)
        # Instantiating the result
        self.forecast_df = forecast_df
        
        if show_validation == False:
            print(">>> CPT has been forecasted ;)")
        elif show_validation == True:
            valid_df = cpt.reset_index()[num_fore:]
            mae_dict = dict()
            mean = self.__mean
            for col in columns:
                mae_val = mean(
                    [abs(a - b) for a, b in zip(
                        list(valid_df[col]),
                        list(prediction_df[col])
                    )]
                )
                extreme_val = max(
                    abs(valid_df[col].min()),
                    abs(valid_df[col].max())
                )
                normalised = abs(mae_val) / extreme_val
                mae_dict[col] = [mae_val, normalised]
            mae_df = pd.DataFrame(mae_dict,
                                  index= ["mean_absolute_error",
                                          "normalised_mae"]
                                 )
            return mae_df
        
    # VISUALISATION METHOD-------------------------------------------
    def visualise(self, data= "cpt",
                  titles= None, xlabs= None, ylabs= None,
                  figsize= (18, 6)):
        """
        Title       : visualise
        Type        : Class Function
        Author      : Rizal Purnawan
        Date        : 06/15/2023
        
        About the function
        ------------------
        
        The argument 'data' has to be one of the following strings:
        - "cpt"            : visualise the CPT data
        - "ma_trends"      : visualise the moving average trends
        - "major_trends"   : visualise the major trends
        - "noises"         : visualise the noises
        - "forecast"       : visualise the forecast result
        - "cpt-forecast"   : visualise the combined cpt with forecast
        """
        plt = self.__plt

        if data == "cpt":
            df = self.cpt.copy()
            if titles is None:
                titles = [f"Plot {col}" for col in df.columns]
        elif data == "ma_trends":
            df = self.ma_trends.copy()
            if titles is None:
                titles = [f"Moving Average: {col}"
                          for col in df.columns]
        elif data == "major_trends":
            df = self.major_trends.copy()
            if titles is None:
                titles = [f"Major Trend: {col}"
                          for col in df.columns]
        elif data == "noises":
            df = self.noises.copy()
            if titles is None:
                titles = [f"Noise: {col}"
                          for col in df.columns]
        elif data == "forecast":
            df = self.forecast_df.copy()
            if titles is None:
                titles = [f"Forecast: {col}"
                          for col in df.columns]
        elif data == "cpt-forecast":
            df_1 = self.cpt.copy()
            df_2 = self.forecast_df.copy()
            if titles is None:
                titles = [f"Investigated and Forecast: {col}"
                          for col in df_1]
        
        columns = list(df_1.columns) if data == "cpt-forecast" \
                    else list(df.columns)
        n_row = len(columns)
        if xlabs is None:
            xlabs = ["Kedalaman (m)"] *n_row
        if ylabs is None:
            ylabs = ["kg/sq-cm"] *n_row
        
        fig = plt.figure(figsize= (figsize[0], n_row *figsize[1]))
        gs = fig.add_gridspec(n_row, 1)
        ax = [
                [fig.add_subplot(gs[i, j]) for j in range(1)]
                for i in range(n_row)
        ]
        axs = [None] *n_row
        
        for i, ax_, col, tit, xlab, ylab in zip(
                range(n_row), axs, columns, titles, xlabs, ylabs
                ):
            if data == "cpt-forecast":
                ax_ = df_1[col].plot(marker= "o", ax= ax[i][0])
                _ = df_2[col].plot(marker= "o", ax= ax[i][0])
            else:
                ax_ = df[col].plot(marker= "o", ax= ax[i][0])
            ax_.set_title(tit, pad= 11, fontsize= 13, weight= "bold")
            ax_.set(xlabel= xlab, ylabel= ylab)
            plt.subplots_adjust(hspace= 0.3);