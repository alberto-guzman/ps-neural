# Deep Neural Network for Propensity Score Analysis

---

## Simulation Plan

---

### LowDim (p = 10)

- 5 normal
- 3 binary
  # 1 ordinal
  1 poison

# Setoguchi et al. 2008 epi (flipped)

## HighDim (p = 100)

- 50 normal

  # 30 binary

  # 20 ordinal

## Vary mixed data composition?

# Confounders (Y,Tx)

# Affect (Y only)

# Affect (Tx only)

# Data Generating

## Propensity Score Model

    # Main Effects Logistic Regression (simple)
    # Complex with interactions and ^order poly

## Outcome Model

    # Main Effects Outcome model
    # Complex with interactions and ^order poly

## Should I vary sample size?

# (n = 50, 500, 2,000)

# Compare:

    # DNN (Keras)
    # Logistic
    # CART (Boosted)
    # Neural Networks
    # Random Forest

#simulate to application example

# look into simulated data

# Analyze Covar Balance

# Optimal Matching to access ATE

# Contributions

# Realistic/Effecient simulation plan for SocialScience != epi

# DNN

# Testing mispecification of not just PS generating model, but outcome model
