# Black-Scholes Option Price & Greeks Calculation and Picture Show

This project use Black-Scholes model to calculate European option price and its five important Greeks number.

## Main Function
- Can calculate option price by Black-Scholes formula way
- Can get Greeks number (Delta, Gamma, Vega, Theta, Rho) two method: 
  - Use `py_vollib` library ready-made function
  - Use `scipy.stats` manual calculate
- Make Vega number change picture when volatility change
- Make Theta number change picture when time pass

## Visualizations
- `vega_sensitivity.png`: Vega vs. Volatility curve
- `theta_decay.png`: Theta vs. Time to expiration curve

## Requirements
- Python 3.x
- `py_vollib`
- `numpy`, `scipy`, `matplotlib`
  
## How to Run
Run the script in terminal or Jupyter:
```bash
python Black_Scholes.py
