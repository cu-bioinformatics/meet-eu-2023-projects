# Meet-EU 2023
## Sorbonne 2
![Meet-EU 2023](https://cu-bioinformatics.github.io/meet-eu-2023/assets/img/4eu.png) 
1. **Data Preparation and Feature Extraction**:
   Converting SMILES strings to RDKit molecule objects and then to MACCS fingerprints, which serve as features for our models.

2. **Clustering**:
    We perform K-means clustering based on the fingerprints, which is a common approach for grouping molecules with similar structural features.
   
4. **Dimensionality Reduction**:
   PCA is attempted on a similarity matrix, which is a valid approach for visualizing high-dimensional data. 

5. **Random Forest Regression**:
   We use Random forest for predicting molecular properties based on fingerprints. We ensure that our training and test data are properly aligned and that the model's hyperparameters are tuned for optimal performance.

6. **Neural Network**:
   We define a simple neural network architecture for regression. While the architecture is a good starting point, we did try experimenting with different layers, neurons, activation functions, and epochs to improve model performance. 

7. **Model Evaluation**:
       Implementing metrics such as Mean Squared Error (MSE), R-squared, or Mean Absolute Error (MAE) would be beneficial for assessing model performance.

8. **Stability and Reproducibility**:
   We attempt to increase the stability of predictions by averaging results over multiple iterations. This approach can help mitigate variability in predictions but also highlights the potential for model uncertainty.
   
9. **Chemical Data Manipulation**:
   Adding known inhibitors and merging them with clusters for prediction.

10. **Exporting Results**:
    Saving results to CSV files for documentation and further analysis.
