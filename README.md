# ME494 HW3: Qube System Identification  

## Description  
This repository contains the third homework assignment for **ME494**, focusing on system identification using data from the Quanser Qube motor with encoder. The homework involves ordinary least squares (OLS) estimation, residual analysis, collinearity checks, and cross-correlation analysis. The repository includes MATLAB scripts, datasets, and a PDF containing problem descriptions.  

## Files Included  

### **Part 1: Ordinary Least Squares Estimation**  
- **File:** SID_HW3.m  
- **File:** sweep_5v_300s.mat  
- **File:** qube_data_multistep.mat  
- **Topics Covered:**  
  - OLS parameter estimation for Qube motion  
  - Model structure identification  
  - Comparison of model accuracy for different input types  
  - Confidence interval calculation for estimated parameters  

### **Part 2: Residual Analysis and Model Fit Assessment**  
- **File:** deriv.m  
- **Topics Covered:**  
  - Numerical differentiation for angular acceleration  
  - Residual plots for model evaluation  
  - Histogram and scatter plot visualization of residuals  

### **Part 3: Normalized Regressor and Collinearity Analysis**  
- **File:** quadroll.mat  
- **Topics Covered:**  
  - Normalization of regressor matrices  
  - Correlation coefficient matrix analysis  
  - Variance inflation factor (VIF) calculations  
  - Condition number estimation using SVD  

### **Part 4: Cross-Correlation and Time Lag Analysis**  
- **File:** SID_HW3_Dowell.pdf  
- **Topics Covered:**  
  - Cross-correlation to quantify time lag  
  - Time shift adjustments for better model fit  
  - Comparison of time lag between different input types  

## Installation  
Ensure MATLAB is installed before running the scripts. No additional toolboxes are required.  

## Usage  

### **Running the Qube Motion OLS Estimation**  
1. Open MATLAB.  
2. Load `sweep_5v_300s.mat` and `qube_data_multistep.mat` into the workspace.  
3. Run the script:  
   ```SID_HW3```  
4. View the estimated parameters and plotted model comparisons.  

### **Performing Residual Analysis**  
1. Open MATLAB.  
2. Run the script:  
   ```SID_HW3```  
3. Observe residuals over time and against predicted acceleration.  

### **Checking Collinearity in System Identification**  
1. Open MATLAB.  
2. Load `quadroll.mat` into the workspace.  
3. Run the collinearity check functions inside `SID_HW3.m`.  
4. View the correlation coefficient matrix, VIF values, and condition number results.  

### **Cross-Correlation Analysis for Time Lag Quantification**  
1. Open MATLAB.  
2. Run the script:  
   ```SID_HW3```  
3. View the cross-correlation results and time-shifted plots.  

## Example Output  

- **Qube Motion Model Estimation**  
  - Estimated Parameters: `[A, B, C, D]` values computed using least squares  
  - Comparison plot of measured vs. modeled angular acceleration  

- **Residual Analysis**  
  - Time series residual plot showing deviations  
  - Histogram of residuals highlighting distribution  

- **Collinearity Assessment**  
  - Correlation coefficient matrix visualization  
  - VIF values and condition number analysis  

- **Cross-Correlation and Time Lag Analysis**  
  - Time-shifted plots for frequency sweep and step-input data  
  - Observed time lag for different input types  

## Contributions  
This repository is intended for academic research and educational use. Contributions and modifications are welcome.  

## License  
This project is open for research and educational purposes.  

---  
**Author:** Alexander Dowell  

