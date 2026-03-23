## Project Title: **Austin Weather Analysis & Precipitation Prediction Using Linear Regression**

### 📌 Description

This project involves analyzing historical weather data from Austin, Texas, to explore patterns and predict precipitation levels. The main goal is to understand how different weather 
attributes relate to precipitation and visualize these relationships using plots. A linear regression model is trained to predict precipitation based on multiple weather factors.

### 📁 Dataset Information

* **File Used**: `austin_weather.csv`
* **Rows**: 1320 entries of daily weather data
* **Key Columns**:

  * Temperature (High, Average, Low)
  * Dew Point
  * Humidity
  * Sea Level Pressure
  * Visibility
  * Wind Speed
  * Precipitation
  * Events (Rain, Thunderstorm, etc.)

### 🧹 Data Preprocessing

* Dropped unnecessary columns: `Date`, `Events`, `SeaLevelPressureLowInches`
* Replaced missing/special values like `'T'` and `'-'` with `0.0`
* Exported cleaned data to `austin_weather_final.csv`

### 📈 Model Used

* **Algorithm**: Linear Regression
* **Library**: `scikit-learn`
* **Target Variable**: `PrecipitationSumInches`
* **Input Features**: All other numerical weather parameters

### 📊 Visualizations

1. **Precipitation Trend Over Days**:

   * A scatterplot showing the change in precipitation over time
   * A specific day (index 798) is highlighted in red

2. **Precipitation vs. Weather Attributes**:

   * Subplots for each selected weather attribute against precipitation
   * Helps understand which factors might influence rainfall

### ⚙️ Requirements

* Python 3.x
* Libraries:

  * `pandas`
  * `numpy`
  * `matplotlib`
  * `scikit-learn`

### 📂 Files Generated

* `austin_weather_final.csv` – Cleaned dataset used for analysis
* Scatterplot visualizations – Shown using `matplotlib`



