import os
import numpy as np
import pandas as pd
import requests
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from tensorflow.keras import layers, models
import seaborn as sns
import matplotlib.pyplot as plt
from typing import Tuple, List, Set
from pathlib import Path
import logging

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class ProteomicsAnalyzer:
    def __init__(self, input_dim: int, n_clusters: int = 3):
        """
        Initialize the ProteomicsAnalyzer with specified dimensions and number of clusters.

        Args:
            input_dim: Dimension of input features
            n_clusters: Number of clusters for KMeans (default: 3)
        """
        self.input_dim = input_dim
        self.n_clusters = n_clusters
        self.scaler = StandardScaler()
        self.pca = PCA(n_components=2)
        self.kmeans = KMeans(n_clusters=n_clusters, random_state=42)
        self.autoencoder, self.encoder = self._build_autoencoder()
        self.year_markers = {
            0: 'o',  # circle
            1: '^',  # triangle up
            2: 's',  # square
        }
        self.cluster_colors = {
            2013: '#440154',  # dark purple
            2015: '#21918c',  # teal
            2016: '#fde725',  # yellow
            2017: '#f46d43',  # orange
            2018: '#3b528b',  # dark blue
            2019: '#fa9fb5',  # pink
            2020: '#4a4a4a',  # gray

        }

    def _build_autoencoder(self) -> Tuple[models.Model, models.Model]:
        """Build and compile the autoencoder model."""
        # Encoder
        inputs = layers.Input(shape=(self.input_dim,))
        x = layers.Dense(128, activation='relu')(inputs)
        x = layers.Dropout(0.2)(x)  # Add dropout for regularization
        x = layers.Dense(64, activation='relu')(x)
        x = layers.Dropout(0.2)(x)
        encoded = layers.Dense(32, activation='relu')(x)

        # Decoder
        x = layers.Dense(64, activation='relu')(encoded)
        x = layers.Dense(128, activation='relu')(x)
        decoded = layers.Dense(self.input_dim, activation='sigmoid')(x)

        # Create models
        autoencoder = models.Model(inputs, decoded)
        encoder = models.Model(inputs, encoded)

        # Compile with better optimizer settings
        autoencoder.compile(
            optimizer='adam',
            loss='mse',
            metrics=['mae']  # Add metrics for better monitoring
        )

        return autoencoder, encoder

    def download_data(self, url: str, filename: str) -> None:
        """
        Download data if not already present.

        Args:
            url: URL to download from
            filename: Name of file to save
        """
        if not Path(filename).exists():
            logger.info(f"Downloading {filename} from {url}")
            response = requests.get(url, stream=True)
            response.raise_for_status()  # Check for download errors

            with open(filename, "wb") as handle:
                for chunk in response.iter_content(chunk_size=8192):
                    handle.write(chunk)
            logger.info(f"Download completed: {filename}")

    def prepare_features(self, data: pd.DataFrame, properties: List[str]) -> np.ndarray:
        """
        Prepare and scale features for analysis.

        Args:
            data: Input DataFrame
            properties: List of properties to use as features

        Returns:
            Scaled features array
        """
        features = data[properties].fillna(0)
        return self.scaler.fit_transform(features)

    def analyze(self, data: pd.DataFrame, properties: List[str]) -> pd.DataFrame:
        """
        Perform complete analysis pipeline.

        Args:
            data: Input DataFrame
            properties: List of properties to analyze

        Returns:
            DataFrame with analysis results
        """
        # Prepare features
        features_scaled = self.prepare_features(data, properties)

        # Train autoencoder
        logger.info("Training autoencoder...")
        history = self.autoencoder.fit(
            features_scaled,
            features_scaled,
            epochs=50,
            batch_size=256,
            validation_split=0.2,  # Add validation
            shuffle=True
        )

        # Get encoded features and cluster
        encoded_features = self.encoder.predict(features_scaled)
        data['Cluster'] = self.kmeans.fit_predict(encoded_features)

        # Add PCA components
        pca_components = self.pca.fit_transform(features_scaled)
        data['PCA1'] = pca_components[:, 0]
        data['PCA2'] = pca_components[:, 1]

        return data

    def save_sample_ids(self, data: pd.DataFrame, output_file: str = "best_samples.txt") -> Set[str]:
        """
        Save sample IDs from the best performing cluster to a file and return them as a set.

        Args:
            data: DataFrame with analysis results
            output_file: Name of the output file

        Returns:
            Set of sample IDs from the best cluster
        """
        # Identify best cluster (highest average peptide sequences)
        best_cluster = data.groupby('Cluster')['Peptide Sequences Identified'].mean().idxmax()

        # Get samples from best cluster
        best_samples = data[data['Cluster'] == best_cluster]

        # Sort by number of peptide sequences identified (descending)
        best_samples_sorted = best_samples.sort_values('Peptide Sequences Identified', ascending=False)

        # Store sample IDs in instance variable
        self.best_samples = set(best_samples_sorted.index)

        # Prepare output content
        output_content = [
            "# Best Performing Samples",
            f"# Cluster: {best_cluster}",
            f"# Total Samples: {len(best_samples_sorted)}",
            "# Format: Sample_ID, Peptide_Sequences_Identified",
            "",
            *[f"{idx}, {count}" for idx, count in
              zip(best_samples_sorted.index, best_samples_sorted['Peptide Sequences Identified'])]
        ]

        # Write to file
        with open(output_file, 'w') as f:
            f.write('\n'.join(output_content))

        logger.info(f"Sample IDs saved to {output_file}")
        return self.best_samples

    def filter_sdrf(self, sdrf_data: pd.DataFrame) -> pd.DataFrame:
        """
        Filter SDRF data to keep only rows where the raw file name matches best samples.

        Args:
            sdrf_data: DataFrame containing SDRF data

        Returns:
            Filtered DataFrame containing only rows matching best samples
        """
        if not self.best_samples:
            raise ValueError("No best samples available. Run analysis first.")

        # Extract raw file names from the comment[data file] column
        def extract_filename(filepath: str) -> str:
            return Path(filepath).stem

        # Create a set of raw file names from best samples
        best_samples_raw = {extract_filename(sample) for sample in self.best_samples}

        # Filter SDRF data
        filtered_data = sdrf_data[
            sdrf_data['comment[data file]'].apply(extract_filename).isin(best_samples_raw)
        ]

        logger.info(f"Filtered SDRF data from {len(sdrf_data)} to {len(filtered_data)} rows")
        return filtered_data

    def plot_results(self, data: pd.DataFrame) -> None:
        """Generate all analysis plots."""
        # Create a figure with subplots
        fig, axes = plt.subplots(2, 2, figsize=(15, 15))

        # Plot 1: Distribution of Peptide Sequences
        sns.histplot(data['Peptide Sequences Identified'], bins=20, kde=True, ax=axes[0, 0])
        axes[0, 0].set_title('Distribution of Peptide Sequences Identified')

        # Plot 2: Clusters in PCA Space
        sns.scatterplot(
            x='PCA1',
            y='PCA2',
            hue='Cluster',
            data=data,
            palette='viridis',
            s=100,
            ax=axes[0, 1]
        )
        axes[0, 1].set_title('Clusters in PCA Space')

        # Plot 3: Peptides vs Clusters
        sns.boxplot(
            x='Cluster',
            y='Peptide Sequences Identified',
            data=data,
            palette='Set2',
            ax=axes[1, 0]
        )
        axes[1, 0].set_title('Peptide Sequences vs Clusters')

        # Plot 4: Best Cluster Distribution
        best_cluster = data.groupby('Cluster')['Peptide Sequences Identified'].mean().idxmax()
        selected_cluster = data[data['Cluster'] == best_cluster]
        sns.histplot(
            selected_cluster['Peptide Sequences Identified'],
            bins=20,
            kde=True,
            ax=axes[1, 1]
        )
        axes[1, 1].set_title(f'Peptide Distribution in Best Cluster ({best_cluster})')

        plt.tight_layout()
        plt.show()

        data['Year'] = pd.to_datetime(data['Content Creation Date']).dt.year

        plt.figure(figsize=(15, 15))

        data['Year'] = pd.to_datetime(data['Content Creation Date']).dt.year

        for year in sorted(data['Year'].unique()):
            year_data = data[data['Year'] == year]
            for cluster in range(3):
                cluster_data = year_data[year_data['Cluster'] == cluster]
                plt.scatter(
                    cluster_data['PCA1'],
                    cluster_data['PCA2'],
                    c=self.cluster_colors[year],
                    marker=self.year_markers[cluster],
                    label=f'{year} - Cluster {cluster}',
                    s=100,
                    alpha=0.7
                )

        plt.title('PCA Analysis by Year and Cluster', pad=20)
        plt.xlabel('First Principal Component')
        plt.ylabel('Second Principal Component')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.grid(True, alpha=0.3)
        plt.show()





# Usage example
if __name__ == "__main__":
    # Define properties for analysis
    properties = [
        'Number of MS1 spectra',
        'Number of MS2 spectra',
        "MS min RT",
        "MS max RT",
        "MS min MZ",
        "MS max MZ",
        "Number of scans",
        "MS/MS Submitted",
        "Mass Standard Deviation [ppm]"
    ]

    # Initialize analyzer
    analyzer = ProteomicsAnalyzer(input_dim=len(properties))

    # Download and load data
    url = 'https://ftp.pride.ebi.ac.uk/pride/data/archive/2023/12/PXD042233/pride_metadata.csv'
    analyzer.download_data(url, "pride_metadata.csv")
    data = pd.read_csv('pride_metadata.csv', index_col=0)

    # Perform analysis
    results = analyzer.analyze(data, properties)

    # Generate plots
    analyzer.plot_results(results)

    # Save sample IDs from best cluster
    analyzer.save_sample_ids(results, "best_samples.txt")

    # Print summary statistics
    print("\nCluster Statistics:")
    print(results.groupby('Cluster')['Peptide Sequences Identified'].describe())

    # Load and filter SDRF data
    sdrf_data = pd.read_csv('PXD042233.sdrf.tsv', sep='\t')  # Adjust filename as needed
    filtered_sdrf = analyzer.filter_sdrf(sdrf_data)

    # Save filtered SDRF data
    filtered_sdrf.to_csv('PXD042233-filtered.sdrf.tsv', sep='\t', index=False)

    # Print summary
    print("\nOriginal SDRF rows:", len(sdrf_data))
    print("Filtered SDRF rows:", len(filtered_sdrf))
    print("\nSample columns in filtered data:")
    for col in filtered_sdrf.columns:
        print(f"- {col}")