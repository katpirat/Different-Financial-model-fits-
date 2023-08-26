##### Nessesary packages to run the script #####

import pandas as pd
from pandas_datareader import data as pdr
import yfinance as yf
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt

yf.pdr_override() ### Needed for the yfinance download of data to work. 

class VaR: 
    def __init__(self, ticker_list, lookback = 300, mc_sims = 100, T = 10, initialPortfolio = 10000) -> None:
        self.s = [i + '.AX' for i in ticker_list]

        self.ticker_list = ticker_list
        self.endDate = dt.datetime.now()
        self.startDate = self.endDate-dt.timedelta(days = lookback)
        self.meanReturns, self.covMatrix = self.get_data(self.s, self.startDate, self.endDate)
        self.mc_sims = mc_sims # Number of simulations of paths; 
        self.T = T
        self.initialPortfolio = initialPortfolio
        self.simulated_stock = self.simulate_stocks()
        
    
    def get_data(self, stocks, start, end):
            stockData = pdr.get_data_yahoo(stocks, start, end)
            stockData = stockData['Close']
            returns = stockData.pct_change()
            meanReturns = returns.mean()
            covMatrix = returns.cov()
            return meanReturns, covMatrix

    def simulate_stocks(self): 
        weights = np.random.random(len(self.meanReturns))
        weights /= np.sum(weights)
        Ini = self.initialPortfolio

        meanM = np.full(shape=(self.T, len(weights)), fill_value=self.meanReturns)
        meanM = meanM.T
        portfolio_sims = np.full(shape=(self.T, self.mc_sims), fill_value=self.initialPortfolio)
        for m in range(0, self.mc_sims):
            Z = np.random.normal(size=(self.T, len(weights))) # uncorrelated normal RV's

            L = np.linalg.cholesky(self.covMatrix) # Cholesky decomposition to Lower Triangular Matrix
            dailyReturns = meanM + np.inner(L, Z) # Correlated daily returns for individual stocks 
            portfolio_sims[:,m] = np.cumprod(np.inner(weights, dailyReturns.T)+1)*self.initialPortfolio
        portfolio_sims[0] = self.initialPortfolio
        return portfolio_sims
    
    def plot_simulation(self): 

        print(self.simulated_stock)
        print(self.simulated_stock[0])
        print(self.simulated_stock[1])

        plt.plot(self.simulated_stock)
        plt.ylabel('Portfolio Value ($)')
        plt.xlabel('Days')
        plt.title('MC simulation of a stock portfolio')
        plt.show()
    
    def plot_cont_Vars(self, alpha = 5): 
        res_VaR = np.zeros(shape = self.T)
        res_CVaR = np.zeros(shape = self.T)
        for i in range(self.T): 
            res_VaR[i] = self.mcVaR(pd.Series(self.simulated_stock[i]), alpha)
            res_CVaR[i] = self.mcCVaR(pd.Series(self.simulated_stock[i]), alpha)
        fig, ax = plt.subplots()

        ax.plot(res_VaR, 'r-', label = 'VaR')
        ax.plot(res_CVaR, 'b-', label = 'CVaR')
        legend = ax.legend(loc='upper right', shadow=True, fontsize='medium')
        plt.show()
    
    def mcVaR(self, returns, alpha = 5):
        if isinstance(returns, pd.Series):
            return np.percentile(returns, alpha)
        else:
            raise TypeError('Expected a pandas data series. ')

    def mcCVaR(self, returns, alpha = 5):

        if isinstance(returns, pd.Series):
            belowVaR = (returns <= self.mcVaR(returns, alpha=alpha))
            return returns[belowVaR].mean()
        else:
            raise TypeError('Expected a pandas data series. ')
    
    def V(self, alpha = 5): 
        pfResults = pd.Series(self.simulated_stock[-1, :])

        VaR = self.initialPortfolio - self.mcVaR(pfResults, alpha=alpha)
        CVaR = self.initialPortfolio - self.mcCVaR(pfResults, alpha=alpha)

        print(f'VaR of the portfolio with initial value{self.initialPortfolio} is : {round(VaR, 2)}')
        print(f'CVaR of the portfolio with initial value{self.initialPortfolio} is : {round(CVaR, 2)}')


stockList = ['CBA', 'BHP', 'TLS', 'NAB', 'WBC', 'STO'] # Example stock list: 
obj1 = VaR(stockList, lookback=300, mc_sims=2000, T = 30, initialPortfolio=10000) # Create the object and inset a list of stock tickers as arguments.  
obj1.V() # Calulates the value at risk and conditional value at risk. The portfolio weights are randomly assigned.  
obj1.plot_cont_Vars(1)
