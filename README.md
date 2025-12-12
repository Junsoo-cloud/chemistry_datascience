# 🧪 데이터과학을 위한 통계방법론(AI20045)
### QSAR Modeling for Log K<sub>ow</sub> Prediction

---

## 📘 Overview
This project aims to build a **Quantitative Structure–Activity Relationship (QSAR)** model to predict the **Log K<sub>ow</sub> (octanol–water partition coefficient)** of Opera Model dataset(2018).

The Log K<sub>ow</sub> value is a critical physicochemical property widely used to estimate the environmental behavior and bioaccumulation potential of organic chemicals.  
Our goal is to construct an interpretable and reliable predictive model for these compounds.

---

## 🎯 Objectives
- Perform **feature selection** among more than 2,000 molecular descriptors.
- - Develop a **QSAR model** for Log K<sub>ow</sub> prediction.  
---

## ⚙️ Workflow

1. **EDA**
   - Create Boxplot, Histogram to see data.
   - Data cleaning : remove redundant rows.

3. **Feature Selection**
   - OLS Linear Regression to see p-value of coefficient
   - LASSO regularization
   - SHAP

4. **Modeling**
   - QSAR model development using:
     - Random Forest, XGBoost, or other ML algorithms.

---

## 👥 Team
**Project Members:**  
- 서울시립대학교 컴퓨터과학부 김준수
- 서울시립대학교 환경원예학과 윤채영
- 서울시립대학교 응용화학과 석서과정 심현수
