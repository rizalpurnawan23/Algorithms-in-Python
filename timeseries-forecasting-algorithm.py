# Bismillahirrahmanirrahim
# 16 August 2022
# Author: Rizal Purnawan

class Forecast:
    def __init__(self, ts):
        # 'ts' stands for timeseries
        self.ts = ts
        
        from statistics import mean
        from math import ceil, sqrt
        self.__mean = mean
        self.__ceil = ceil
        self.__sqrt = sqrt
    
    def read_me(self):
        return (
        """
        ----------------------------------------------------------------
        ___About the class
        
            This class contains the machinary to deal with the
        forecasting tasks. The class takes a pandas data frame as an
        input which is then instantiated as a class attribute 'ts'.
        The data frame is generally a time series data, which can have
        multiple columns, in the sense that, there are finite
        multi-dimensional observations at a given time.
            The outline of the forecasting tasks offered by this class
        is described as follows:
            1. First, we determine the major trend of each observation
               (column) in the time series by making use of the moving
               average trend.
            2. Then we create a model that approximates the moving
               average trends. These approximated trends are then
               considered as the major trends.
            3. We determine the noises of the data by detrending the
               data. That is, the noises are defined as the subtraction
               of the original data from the major trends.
            4. Then we create a model to approximate the noises.
            5. Finally, we can perform a forecasting by making use of
               the model we created to approximate the major trends and
               the model we created to approximate the noises.
               
        ___Notice
        
            This algorithm uses python libraries namely 'Pandas', 'Numpy',
        'Scikit Learn', 'Statsmodels', 'Matplotlib', 'Seaborn' and
        'Plotly'. If one uses a machine that does not have these
        libraries installed, we suggest to install them first.
        
            If one works in a notebook, the following scripts can be
        copied and run in a cell for installation:
        
            !pip install pandas
            !pip install numpy
            !pip install sklearn
            !pip install statsmodels
            !pip install matplotlib
            !pip install seaborn
            !pip install plotly
        
            Once the installation finishes, please copy the following
        scripts and run them in another cell if one works in a notebook,
        or copy them at the top of the page if one works in an IDE:
            
           import pandas as pd
           import numpy as np
           import matplotlib.pyplot as plt
           import plotly.graph_objects as go
           sns.set_style("dark")
           from sklearn.linear_model import LinearRegression
           from statsmodels.tsa.deterministic import CalendarFourier, \
                                                     DeterministicProcess
           from statsmodels.graphics.tsaplots import plot_pacf
        
                                                    Author,
                                                    Rizal Purnawan
        ----------------------------------------------------------------
        """
        )
        
    # Some useful probability theoretic methods:
    def __expected_value(self, X):
        # Argument 'X' has to be a list that represents the
        # numerical values of the observation (discrete random variable).
        # We will use definition 1 to build the algorithm --however,
        # in case that our algo is only compatible for discrete
        # random variables with fair probability distributions.
        uniq_vals = list(set(X))
        count_vals = [X.count(x) for x in uniq_vals]
        prob_mass = [n / len(X) for n in count_vals]
        return sum([p * x for p, x in zip(uniq_vals, prob_mass)])
    
    def __covariance(self, X, Y):
        # Arguments "X" and "Y" have to be lists that represents the
        # numerical values of the observations (discrete random
        # variables).
        ev_X, ev_Y = self.__expected_value(X), self.__expected_value(Y)
        XY = [x * y for x, y in zip(X, Y)]
        ev_XY = self.__expected_value(XY)
        return ev_XY - (ev_X *ev_Y)
    
    def __variance(self, X):
        # Argument "X" has to be a list that represent the numerical
        # values of the observation (discrete random variable).
        return self.__covariance(X, X)
    
    def __standard_dev(self, X):
        # Argument "X" has to be a list that represent the numerical
        # values of the observation (discrete random variable).
        sqrt = self.__sqrt
        return sqrt(self.__variance(X))
    
    def __correlation(self, X, Y):
        # Arguments "X" and "Y" have to be lists that represents the
        # numerical values of the observations (discrete random
        # variables).
        sqrt = self.__sqrt
        return (
            self.__covariance(X, Y)
            / sqrt( self.__variance(X) * self.__variance(Y) )
        )
        
    def moving_avg_trends(self, window):
        """
        ___Additional Info
        
            The argument 'window' is the window for rolling method in
        pandas.
        """
        ceil = self.__ceil
        ts = self.ts.copy()
        cols = ts.columns
        ma_trends = dict()
        for col in cols:
            y = ts[col].squeeze()
            y_ma = y.rolling(window= window,
                             min_periods= ceil(window/2),
                             center= True).mean()
            ma_trends[col] = y_ma
        ma_trends = pd.DataFrame(ma_trends, index= ts.index)
        self.ma_trends = ma_trends
        print(">>> Moving average trends are set")
        
    def __det_polynom_ord(self, max_ord= 10):
        """
        ___Additional Info
        
            This method is created to determine the order of the
        polynomial used to approximate the major trends based on the
        moving average trends that have been computed. In doing so,
        we will iterate all the possible polynomial orders such that it
        forms trends that have the smallest mean absolute error to the
        moving average trends. However we set the maximum possible
        order of the polynomial (argument 'max_ord') to 10 in default,
        since we do not want to obtain a misleading forecasting with
        considerably higher polynomial orders.
        """
        ts = self.ts.copy()
        ma_trends = self.ma_trends.copy()
        cols = ts.columns
        dp_ords = list()
        mean = self.__mean
        for col in cols:
            mae_list = list()
            for ord in range(1, max_ord + 1):
                dp = DeterministicProcess(index= ts.index,
                                          constant= True,
                                          order= ord)
                X_maj = dp.in_sample()
                y = ts[col].squeeze()
                model = LinearRegression(fit_intercept= False)
                model.fit(X_maj, y)
                y_maj = model.predict(X_maj)
                mae_list.append(
                    [ord, mean([abs(ma_trends[col][k] - y_maj[k])
                                for k in range(len(y_maj))])
                    ]
                )
            ord_ = [a for a, b in mae_list
                      if b == min([m[1] for m in mae_list])][0]
            dp_ords.append(ord_)
        self.dp_ords = dp_ords
        
    def set_major_trends(self, dp_params= None):
        """
        ___Additional Info
        
            Argument 'dp_params' has to be a list of positive
        integres. It contains the polynomial order for each column
        in the time series. By default we set it 'None', and it takes
        the result from method '__det_polynom_ord'. Otherwise, one may
        determine the value as desired.
        """
        if dp_params is None:
            self.__det_polynom_ord()
            dp_params = self.dp_ords
        ts = self.ts.copy()
        ts_cols = ts.columns
        X_maj_fts = [None] *len(ts_cols)
        dps = [None] *len(ts_cols)
        models = [None] *len(ts_cols)
        maj_trends = list()
        for k, col in zip(range(len(ts_cols)), ts_cols):
            dps[k] = DeterministicProcess(index= ts.index,
                                          constant= True,
                                          order= dp_params[k])
            X_maj = dps[k].in_sample()
            X_maj_fts[k] = X_maj
            y = ts[col].squeeze()
            model = LinearRegression(fit_intercept= False)
            model.fit(X_maj, y)
            models[k] = model
            maj_trends.append(list(model.predict(X_maj)))
        mt_cols = [f"maj_trend_{col}" for col in ts_cols]
        maj_trends = pd.DataFrame({mt_cols[k]: maj_trends[k]
                                   for k in range(len(mt_cols))},
                                  index= ts.index)
        self.maj_trends = maj_trends
        self.maj_models = models
        self.dps = dps
        print(">>> Major trends are set")
        
    def set_noises(self):
        """
        ___Additional Info
        
            This function helps us determine the noises from the
        original time series by subtracting the original time series
        whith the major trends determined earlier.
        """
        ts = self.ts.copy()
        ts_fitted = self.maj_trends
        noises = pd.DataFrame(
            {f"noise_{col_ts}": ts[col_ts] - ts_fitted[col_fit]
             for col_ts, col_fit in zip(ts.columns, ts_fitted.columns)},
            index= ts.index)
        self.noises = noises
        print(">>> Noises are set")
        
    def __forecast_corr(self, n_fore= 1):
        """
        ___Additional Info
        
            This function will tell the correlation of the available
        time series data with the forecasting result for some future
        points forecasting. In fact it checks the correlation of the
        original noise and the expected future noise. The correlation
        values indicate the accuracy of the forecast given the future
        data points to be forecasted (forecast horizon).
        """
        noise = self.noises
        cols = noise.columns
        fore_core = pd.DataFrame(
            {f"{col}_{n_fore}_shifted":
                 noise[col].shift(n_fore).fillna(0.0)
                 for col in cols},
            index= noise.index
        )
        corr_pair = list(zip(cols, fore_core.columns))
        corr_df = pd.DataFrame(
            {str(cp[0]): [self.__correlation(list(noise[cp[0]]),
                                         list(fore_core[cp[1]])
                                        )
                      ] for cp in corr_pair},
            index= ["forecast_correlation"]
        )
        return corr_df
    
    def set_noises_shift(self, n_fore):
        """
        ___Additional Info
        
            Argument 'n_fore' is the number of datapoints to be
        forecasted.
        """
        self.n_fore = n_fore
        df_cols_shifted = list()
        noises = self.noises.copy()
        cols = noises.columns
        dp = DeterministicProcess(index= self.ts.index,
                                  constant= False,
                                  order= 1)
        id_time = dp.in_sample()[n_fore:].index
        id_fore = dp.out_of_sample(steps= n_fore).index
        id_all = id_time.append(id_fore)
        noises_shifted = pd.DataFrame(
            {col: list(noises[col]) for col in cols},
            index= id_all)
        self.noises_shifted = noises_shifted
        print(">>> Shifted noises are set")
            
    def pacf_shifted_noises(self, max_lags= 30):
        """
        ___Additional Info
        
            Argument 'max_lags' has to be a positive integer
        representing the desired maximum number of lags for the
        shifted noises. We set it to be '30' for the default value.
        """
        noises_shifted = self.noises_shifted.copy()
        cols = noises_shifted.columns
        for col in cols:
            plot_pacf(noises_shifted[col],
                      lags= max_lags,
                      title= f"PACF of {col}");
            
    def set_shifted_noises_lags(self, lags):
        """
        ___Additional Info
        
            Argument 'lags' has to be a list consisting the lags
        determined from the observation of partial autocorrelation
        function.
        """
        noises = self.noises_shifted.copy()
        cols = noises.columns
        noises_shifted_lagged = list()
        for l, col in zip(range(len(cols)), cols):
            df = pd.DataFrame(
                {f"{col}_lag_{lags[l][k]}":
                 noises[col].shift(lags[l][k]).fillna(0.0)
                 for k in range(len(lags[l]))},
                index= noises.index
            )
            noises_shifted_lagged.append(df)
        self.noises_shifted_lagged = noises_shifted_lagged
        print(">>> Lagged shifted noises are set")
    
    def forecast(self, noise_mods= None,
                 show_validation= False):
        """
        ___Additional Info
        
            Argument 'noise_mods' must be a list of available machine
        learning algorithms for regression. We set it 'None' and then it
        is redefined as a list whose members are all 'LinearRegression'
        from sklearn. One may change it with the desired machine learning
        algorithms.
        """
        n_fore = self.n_fore
        ts = self.ts.copy()
        ts_columns = ts.columns
        col_nums = len(ts_columns)
        
        maj_trends = self.maj_trends.copy()
        maj_models = self.maj_models
        
        shifted_noises = self.noises_shifted_lagged.copy()
        noises = self.noises.copy()
        forecast_list = list()
        noises_models = list()
        
        if noise_mods is None:
            noise_mods = [LinearRegression] *col_nums
        for k in range(col_nums):
            # Forecasting the major trend
            dp_maj = self.dps[k]
            Xm_time = dp_maj.in_sample()
            Xm_fore = dp_maj.out_of_sample(steps= n_fore)
            maj_model = maj_models[k]
            ym_fore = maj_model.predict(Xm_fore)
            
            # Forecasting the noise
            Xn_time = shifted_noises[k][0:-n_fore]
            Xn_fore = shifted_noises[k][-n_fore:]
            yn = noises[noises.columns[k]][n_fore:].squeeze()
            n_model = noise_mods[k]()
            n_model.fit(Xn_time, yn)
            yn_fore = n_model.predict(Xn_fore)
            noises_models.append(n_model)
            
            # Combined
            y_fore = np.array([a + b for a, b in zip(ym_fore, yn_fore)])
            forecast_list.append(y_fore)
        fore_index = shifted_noises[0][-n_fore:].index
        forecast_df = pd.DataFrame(
            {ts_columns[k]: forecast_list[k]
                            for k in range(col_nums)},
            index= fore_index)
        
        # Instantiating useful results:
        self.forecast_df = forecast_df
        self.hist_fore = pd.concat([ts, forecast_df], axis= 0)
        self.noises_models = noises_models
        
        # Validation with Mean Absolute Error:
        if show_validation == True:
            mean = self.__mean
            mae_list = list()
            for k in range(col_nums):
                X_maj = self.dps[k].in_sample()
                apx_maj_trend = maj_models[k].predict(X_maj)

                X_noi = shifted_noises[k][0:-n_fore]
                apx_noise = noises_models[k].predict(X_noi)

                apx_ts = [a + b for a, b in
                          zip(apx_maj_trend, apx_noise)]
                ori_ts = list(ts[ts.columns[k]])

                abs_err = [abs(a - b) for a, b in zip(ori_ts, apx_ts)]
                mae = mean(abs_err)
                mae_list.append(mae)
            mae_df = pd.DataFrame(
                {ts_columns[k]: [mae_list[k]] for k in range(col_nums)},
                index= ["mean_abs_error"]
            )
            display(mae_df)
        print(">>> Forecasting has been done")        
    
    def visualise(self, data= "ts",
                  titles= None,
                  figsize= (15, 5),
                  xlabs= None, ylabs= None,
                  save= False):
        """
        ___Additional Info
        
            Argument 'data' should be one of the following:
        - 'ts'        : visualise the original timeseries
        - 'ma_trends' : visualise the moving average trends
        - 'maj_trends': visualise the preset major trends
        - 'noises'    : visualise the noises
        - 'forecast'  : visualise the forecast result
        - 'hist_fore' : visualise the historical data as well as
                        the forecast.
        """
        df, cols = None, None
        if data == "ts":
            df = self.ts.copy()
        elif data == "ma_trends":
            df = self.ma_trends.copy()
        elif data == "maj_trends":
            df = self.maj_trends.copy()
        elif data == "noises":
            df = self.noises
        elif data == "forecast":
            df = self.forecast_df
        elif data == "hist_fore":
            df = self.hist_fore
            df_1 = self.ts
            df_2 = self.forecast_df
            
        cols = df.columns
        if titles is None:
            titles = list(cols)
        if xlabs is None:
            xlabs = ["Time"] *len(cols)
        if ylabs is None:
            ylabs = ["Value"] *len(cols)
        n_row = len(cols)
        fig = plt.figure(figsize= (figsize[0], n_row *figsize[1]))
        gs = fig.add_gridspec(n_row, 1)
        ax = [[fig.add_subplot(gs[i, j]) for j in range(1)]
              for i in range(n_row)]
        axs = [None] *n_row
        for i, ax_, col, tit, xl, yl in zip(range(n_row),
                                            axs,
                                            cols,
                                            titles,
                                            xlabs,
                                            ylabs):
            if data != "hist_fore":
                if data == "forecast":
                    ax_ = df[col].plot(ax= ax[i][0],
                                       marker= "o",
                                       color= "orange")
                else:
                    ax_ = df[col].plot(ax= ax[i][0])
            else:
                ax_ = df_1[col].plot(ax= ax[i][0])
                _ = df_2[col].plot(ax= ax[i][0])
            ax_.set_title(tit,
                          pad= 12,
                          fontsize= 15,
                          fontweight= "bold")
            ax_.set(xlabel= xl, ylabel= yl)
            plt.subplots_adjust(hspace= 0.5);
            
    def plot_candlestick(self, data= "ts",
                         title= "Candlestick Plot",
                         yax_title = "Value"):
        """
        ___Additional Info
        
            This method provides us with the candlestick plot of
        a given timeseries data. Argument 'data' has to be one
        of the following:
        - 'ts'        : visualise the candlestick of the original
                        timeseries.
        - 'forecast'  : visualise the candlestick of the forecast
                        result.
        - 'hist_fore' : visualise the candlestick of the original
                        timeseries combined with the forecast result.
        """
        if data == "ts":
            ts = self.ts.copy()
        elif data == "forecast":
            ts = self.forecast_df.copy()
        elif data == "hist_fore":
            ts = self.hist_fore.copy()
        index = ts.index.strftime("%Y-%m-%d")
        used_cols = ts.columns[0:4]
        figori_cs = go.Figure(
            data= [go.Candlestick(x= index,
                                  open= ts[used_cols[0]],
                                  high= ts[used_cols[1]],
                                  low= ts[used_cols[2]],
                                  close= ts[used_cols[3]],
                                  increasing_line_color= "cyan",
                                  decreasing_line_color= "gray"
                                 )
                  ]
        )
        figori_cs.update_layout(
            title= title,
            yaxis_title= yax_title)
        figori_cs.show()
