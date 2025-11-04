# 🧪 Organophosphorus Compounds Data Analysis Project  
### _QSAR Modeling for Log K<sub>ow</sub> Prediction_

---

## 📘 Overview
This project aims to build a **Quantitative Structure–Activity Relationship (QSAR)** model to predict the **Log K<sub>ow</sub> (octanol–water partition coefficient)** of **organophosphorus compounds**.

The Log K<sub>ow</sub> value is a critical physicochemical property widely used to estimate the environmental behavior and bioaccumulation potential of organic chemicals.  
Our goal is to construct an interpretable and reliable predictive model for these compounds.

---

## 🎯 Objectives
- Develop a **QSAR model** for Log K<sub>ow</sub> prediction.  
- Perform **feature selection** among more than 2,000 molecular descriptors.  
- Apply **explainable AI (XAI)** techniques to interpret model predictions.  

---

## ⚙️ Workflow

1. **Data Collection**
   - Dataset of organophosphorus compounds with known Log K<sub>ow</sub> values.
   - Molecular descriptors generated using cheminformatics tools (e.g., RDKit, Mordred).

2. **Exploratory Data Analysis (EDA)**
   - Descriptive statistics, distribution visualization, and correlation analysis.
   - Identification of redundant or non-informative descriptors.

3. **Feature Selection**
   - EDA
   - LASSO regularization

4. **Modeling**
   - QSAR model development using:
     - Random Forest, XGBoost, or other ML algorithms.
   - If we have some time and more data, then we'll apply DL method(GCN)

5. **Explainability**
   - Apply **XAI (Explainable AI)** methods such as:
     - **LIME** – Local Interpretable Model-agnostic Explanations  
     - **SHAP** – SHapley Additive exPlanations  
   - Visualize and interpret feature contributions to Log K<sub>ow</sub> predictions.

---

## 👥 Team
**Project Members:**  
- Each of our team member is Computer Science, Chemistry, Environmental Horticulture
