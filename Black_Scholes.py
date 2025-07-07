# The above code calculates the Black-Scholes option price and its Greeks for a given set of parameters.
# It uses the py_vollib library for the calculations and scipy.stats for statistical functions.
from py_vollib.black_scholes import black_scholes as bs
from py_vollib.black_scholes.greeks.analytical import delta, gamma, vega, theta, rho
import numpy as np
from scipy.stats import norm

# Function to calculate the Black-Scholes price
r = 0.01  # Risk-free interest rate
sigma = 0.20  # Volatility 
T = 260/365  # Time to expiration
K = 100  # Strike price
S = 90  # stock price
option_type = 'c'  # Option type: 'c' for call, 'p' for put

#define valuables and calculate
def calculate_black_scholes_price(S, K, T, r, sigma, option_type='c'):
    "Calculate BS price of call/put"
    d1 = (np.log(S / K) + (r + 0.5 * sigma**2) * T) / (sigma * np.sqrt(T))
    d2 = d1 - sigma * np.sqrt(T)
    try:
        if option_type == 'c':
            price = S * norm.cdf(d1, 0, 1) - K * np.exp(-r * T) * norm.cdf(d2, 0, 1)
        elif option_type == 'p':
            price = K * np.exp(-r * T) * norm.cdf(-d2, 0, 1) - S * norm.cdf(-d1, 0, 1)
        else:
            raise ValueError("Invalid option type. Use 'c' or 'p'.")
        return price, bs(option_type, S, K, T, r, sigma)
    except:
        print("Error! Please check the input option type. Only from 'c' or 'p'.")
        return None

def delta_calc(S, K, T, r, sigma, option_type='c'):
    "Calculate BS delta of call/put"
    d1 = (np.log(S / K) + (r + 0.5 * sigma**2) * T) / (sigma * np.sqrt(T))
    try:
        if option_type == 'c':
            delta_calc = norm.cdf(d1, 0, 1)
        elif option_type == 'p':
            delta_calc = -1*norm.cdf(-d1, 0, 1)
        else:
            raise ValueError("Invalid option type. Use 'c' or 'p'.")
        return delta_calc, delta(option_type, S, K, T, r, sigma)
    except:
        print("Error! Please check the input option type. Only from 'c' or 'p'.")
        return None

def gamma_calc(S, K, T, r, sigma, option_type='c'):
    "Calculate BS gamma of call/put"
    d1 = (np.log(S / K) + (r + sigma**2 /2) * T )/ (sigma*np.sqrt(T))
    try:
        gamma_calc = norm.pdf(d1, 0, 1) / (S * sigma * np.sqrt(T))
        return gamma_calc, gamma(option_type, S, K, T, r, sigma)
    except:
        print("Error! Please check the input option type. Only from 'c' or 'p'.")
        return None

def vega_calc(S, K, T, r, sigma, option_type='c'):
    "Calculate BS vega of call/put"
    d1 = (np.log(S / K) + (r + sigma**2 /2) * T )/ (sigma*np.sqrt(T))
    d2 = d1 - sigma * np.sqrt(T)
    try:
        vega_calc = S * norm.pdf(d1, 0, 1) * np.sqrt(T)
        return vega_calc*0.01, vega(option_type, S, K, T, r, sigma)
    except:
        print("Error! Please check the input option type. Only from 'c' or 'p'.")
        return None

def theta_calc(S, K, T, r, sigma, option_type='c'):
    "Calculate BS theta of call/put"
    d1 = (np.log(S / K) + (r + sigma**2 /2) * T )/ (sigma*np.sqrt(T))
    d2 = d1 - sigma * np.sqrt(T)
    try:
        if option_type == 'c':
            theta_calc = (-S * norm.pdf(d1, 0, 1) * sigma / (2 * np.sqrt(T)) - r * K * np.exp(-r * T) * norm.cdf(d2, 0, 1)) / 365
        elif option_type == 'p':
            theta_calc = (-S * norm.pdf(d1, 0, 1) * sigma / (2 * np.sqrt(T)) + r * K * np.exp(-r * T) * norm.cdf(-d2, 0, 1)) / 365
        else:
            raise ValueError("Invalid option type. Use 'c' or 'p'.")
        return theta_calc, theta(option_type, S, K, T, r, sigma)
    except:
        print("Error! Please check the input option type. Only from 'c' or 'p'.")  
        return None

def rho_calc(S, K, T, r, sigma, option_type='c'):
    "Calculate BS rho of call/put"
    d1 = (np.log(S / K) + (r + sigma**2 /2) * T )/ (sigma*np.sqrt(T))
    d2 = d1 - sigma * np.sqrt(T)
    try:
        if option_type == 'c':
            rho_calc = K * T * np.exp(-r * T) * norm.cdf(d2, 0, 1) / 100
        elif option_type == 'p':
            rho_calc = -K * T * np.exp(-r * T) * norm.cdf(-d2, 0, 1) / 100
        else:
            raise ValueError("Invalid option type. Use 'c' or 'p'.")
        return rho_calc, rho(option_type, S, K, T, r, sigma)
    except:
        print("Error! Please check the input option type. Only from 'c' or 'p'.")
        return None
        
# Example usage
try:
    print("Option Price:", [round(float(x), 3) for x in calculate_black_scholes_price(S, K, T, r, sigma, option_type)])
except:
    print("Option Price: Calculation failed.")
try:
    print("       Delta:", [round(float(x), 3) for x in delta_calc(S, K, T, r, sigma, option_type)])
except:
    print("Delta: Calculation failed.")
try:
    print("       Gamma:", [round(float(x), 3) for x in gamma_calc(S, K, T, r, sigma, option_type)])
except:
    print("Gamma: Calculation failed.")
try:
    print("        Vega:", [round(float(x), 3) for x in vega_calc(S, K, T, r, sigma, option_type)])
except:
    print("Vega: Calculation failed.")
try:
    print("       Theta:", [round(float(x), 3) for x in theta_calc(S, K, T, r, sigma, option_type)])
except:
    print("Theta: Calculation failed.")
try:
    print("         Rho:", [round(float(x), 3) for x in rho_calc(S, K, T, r, sigma, option_type)])
except:
    print("Rho: Calculation failed.")

#draw Vega Sensitivity Graph
import matplotlib.pyplot as plt

sigma_range = np.linspace(0.05, 1.0, 100)
vega_values = [vega_calc(S, K, T, r, s, option_type)[0] for s in sigma_range]

plt.figure(figsize=(12, 6))
plt.plot(sigma_range, vega_values, linewidth=2.5, color='#2c7bb6')

current_vega = vega_calc(S, K, T, r, sigma, option_type)[0]
plt.scatter(sigma, current_vega, color='red', s=100, 
           label=f'Current (σ={sigma}, Vega={current_vega:.2f})')

plt.title(f'Vega Sensitivity (S={S}, K={K}, T={T:.2f}y)')
plt.xlabel('Volatility (σ)')
plt.ylabel('Vega Value')
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()

plt.savefig('vega_analysis.png', dpi=150)
plt.show()

#draw Theta Decay Graph
T_days = np.linspace(5, 365, 100)  # Days to expiration
theta_values = [theta_calc(S, K, d/365, r, sigma, option_type)[0] for d in T_days]

plt.figure(figsize=(10, 5))
plt.plot(T_days, theta_values, 'r-', linewidth=2, label='Theta (Daily Time Decay)')
# Mark current position
current_days = T * 365
current_theta = theta_calc(S, K, T, r, sigma, option_type)[0]
plt.scatter(current_days, current_theta, color='black', s=80,
           label=f'Current Position ({current_days:.0f} DTE)')

plt.title(f'Option Time Decay Analysis\n(Strike K={K}, Vol σ={sigma:.0%})', pad=20)
plt.xlabel('Days to Expiration')
plt.ylabel('Theta Value (Daily Loss)')
plt.axhline(0, color='gray', ls='--', alpha=0.5)
plt.grid(axis='y', alpha=0.3)
plt.legend()

plt.tight_layout()
plt.savefig('theta_analysis_en.png', dpi=120, bbox_inches='tight')
plt.show()