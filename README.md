# Relational persistent homology (PH) for tumor-immune point clouds

We use Dowker PH and Multispecies Witness PH to vectorise point clouds for input to ML algorithms
# Dowker PH
In the `dowker analysis` folder:
- `compute_persistence_diagrams.ipynb` computes Vietoris-Rips (via [Ripser](https://mtsch.github.io/Ripserer.jl/dev/)) and Dowker persistence diagrams (via `extension_method.jl`  - see https://github.com/irishryoon/Dowker_persistence) on the point cloud data 
- `compute_persistence_images.ipynb` computes vectorisations of persistence diagrams using persistence images
- `compute_non_topological_stats.ipynb` computes non topological features of the point cloud data for comparison to topological methods
- `compute_labels.ipynb` computes ground truth binary labels for the point cloud data based on the dominant macrophage phenotype
- `component_svm.ipynb` performs SVM classification on individual Vietoris-Rips and Dowker features as well as non-topological features and saves results to the  `analysis_results` folder
- `combined_svm.ipynb` performs SVM classification on combined Vietoris-Rips and combined Dowker features and saves results to the  `analysis_results` folder

`dimensionality_reduction.ipynb` peforms dimenionality reduction for combined Vietoris-Rips and combined Dowker features.
`plot_dowker_results.ipynb` plots the accuracies of the SVM classifiers from the `analysis_results` folder.

# Multispecies Witness PH
- `MultiplexTrivialTest.py` performs k-means clustering on non topological features
- `MultiplexWitnessPH.py` performs k-means clustering on the topological distance vectors computed via Multispecies PH of the point cloud data without knowledge of macrophage phenotype
- `MultiplexWitnessPH.py` performs k-means clustering on the topological distance vectors computed via Multispecies PH of the point cloud data with knowledge of macrophage phenotype
