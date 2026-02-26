# COVID-19 vs World Happiness Report: Data Analysis

This project explores the relationship between the COVID-19 infection rates and various socioeconomic indicators from the World Happiness Report 2019. It combines COVID-19 time-series case
data with happiness metrics such as GDP per capita, social support, life expectancy, and freedom to make life choices. The analysis includes visualizations and correlation insights.

## 🔧 Steps Performed

### 1. Data Loading and Preprocessing
- Imported and cleaned COVID-19 data.
- Aggregated by country and removed unnecessary columns.
- Calculated **maximum daily infection rate** for each country using the first derivative of case counts.

### 2. World Happiness Data Cleaning
- Dropped unused columns: `Overall rank`, `Score`, `Generosity`, `Perceptions of corruption`.
- Set `Country or region` as the index for easier merging.

### 3. Merging Datasets
- Performed an inner join on country names.
- Created a new DataFrame with relevant columns:
  - `Maximum_Infection_Rates`
  - `GDP per capita`
  - `Social support`
  - `Healthy life expectancy`
  - `Freedom to make life choices`

### 4. Exploratory Data Analysis
- Calculated and displayed the correlation matrix.
- Visualized the relationships between `Maximum_Infection_Rates` and:
  - GDP per capita
  - Social support
  - Healthy life expectancy
  - Freedom to make life choices

### 5. Visualization
- Used `Seaborn` and `Matplotlib` for:
  - Scatter plots (with log scale for infection rates)
  - Regression plots to visualize trend lines

---

## 📊 Sample Visualizations

- **Scatterplot:** GDP per capita vs log(Max Infection Rate)
- **Regression Plot:** Social support vs log(Max Infection Rate)

> All plots help uncover trends between national well-being and vulnerability to infection spread.

---

## 💡 Findings

- Moderate positive correlation between GDP and infection rate.
- Life expectancy and social support show similar trends.
- Interpretation: Countries with better well-being may have higher reported infections due to better tracking/testing or urbanization.

---

## 📦 Libraries Used

- `pandas`
- `numpy`
- `matplotlib`
- `seaborn`

---

## 📝 How to Run

1. Clone this repository or copy the notebook/scripts.
2. Ensure the following files are in the same directory:
   - `covid19_Confirmed_dataset.csv`
   - `worldwide_happiness_report.csv`
3. Run the script or Jupyter Notebook.

---
