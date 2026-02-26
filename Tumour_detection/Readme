Tumor Classification using Random Forest

This project aims to build a machine learning model that can classify breast tumors as **malignant** or **benign** based on various medical features. Using a Random Forest Classifier,
the model is trained and evaluated on real-world diagnostic data.

## Dataset Description

The dataset used in this project is derived from breast cancer diagnostic records. It consists of **569 instances** and **32 columns**, including:

* A `diagnosis` column indicating the tumor type: **Malignant (M)** or **Benign (B)**
* **30 numerical features** describing characteristics of the cell nuclei present in digitized images of fine needle aspirates (FNAs) of breast masses
* An `id` column which is dropped during preprocessing

## Workflow Overview

1. **Data Loading and Cleaning**

   * The dataset is loaded from a CSV file.
   * The `id` column is removed as it is not useful for classification.
   * The `diagnosis` column is converted from categorical to numerical form (Malignant → 1, Benign → 0).

2. **Exploratory Data Analysis (EDA)**

   * A count plot is generated to observe the class distribution.
   * Correlation analysis is conducted using a heatmap to identify which features are most correlated with the diagnosis.

3. **Feature Engineering**

   * Features are divided into three types:

     * Mean values
     * Standard error
     * Worst-case values
   * Summary statistics are computed to understand data distribution.

4. **Data Preparation**

   * Features (`X`) and target (`y`) are separated.
   * Data is split into training and test sets (70% training, 30% testing).
   * Standardization is applied to scale the feature values.

5. **Model Building**

   * A Random Forest Classifier is trained on the training data.
   * Predictions are made on the test data.

6. **Model Evaluation**

   * The model is evaluated using accuracy score.
   * In this case, the classifier achieved approximately **95.9% accuracy**.

## Key Insights

* The dataset is well-balanced between malignant and benign cases.
* Certain features like **radius\_mean**, **concave points\_mean**, and **area\_worst** show high correlation with the tumor diagnosis.
* Random Forest performs well even without extensive tuning.

## Tools and Libraries

* **Pandas** and **NumPy** for data manipulation
* **Seaborn** and **Matplotlib** for visualization
* **Scikit-learn** for preprocessing, modeling, and evaluation

## Future Enhancements

* Apply additional classifiers for comparison (e.g., Support Vector Machine, Gradient Boosting)
* Perform hyperparameter tuning to optimize the Random Forest model
* Use techniques like cross-validation to further validate performance
* Deploy the model using a web interface for interactive predictions


