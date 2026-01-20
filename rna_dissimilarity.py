#!/usr/bin/env python3
"""
16S RNA Dissimilarity Analysis Tool

This module performs dissimilarity analysis on 16S RNA sequences based on k-mer frequencies.
It includes functionality to:
- Extract k-mers from RNA sequences
- Calculate dissimilarity matrices
- Generate heatmaps
- Visualize histograms of dissimilarities
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
from itertools import combinations
import warnings
# Filter only matplotlib-specific warnings to avoid clutter in visualization output
warnings.filterwarnings('ignore', category=UserWarning, module='matplotlib')


class RNADissimilarityAnalyzer:
    """
    Analyzer for computing dissimilarity between 16S RNA sequences using k-mer based methods.
    """
    
    def __init__(self, sequences, labels=None, k=3):
        """
        Initialize the RNA Dissimilarity Analyzer.
        
        Parameters:
        -----------
        sequences : list of str
            List of RNA sequences (strings containing A, C, G, U/T)
        labels : list of str, optional
            Labels for each sequence. If None, auto-generates labels.
        k : int, default=3
            Size of k-mers to use for analysis
        """
        self.sequences = [seq.upper().replace('T', 'U') for seq in sequences]
        self.k = k
        self.labels = labels if labels else [f"Seq_{i+1}" for i in range(len(sequences))]
        self.kmer_profiles = None
        self.dissimilarity_matrix = None
        
    def extract_kmers(self, sequence):
        """
        Extract all k-mers from a sequence.
        
        Parameters:
        -----------
        sequence : str
            RNA sequence
            
        Returns:
        --------
        Counter
            Dictionary with k-mer counts
        """
        kmers = []
        for i in range(len(sequence) - self.k + 1):
            kmers.append(sequence[i:i+self.k])
        return Counter(kmers)
    
    def build_kmer_profiles(self):
        """
        Build k-mer frequency profiles for all sequences.
        
        Returns:
        --------
        list of Counter
            K-mer frequency profiles for each sequence
        """
        self.kmer_profiles = [self.extract_kmers(seq) for seq in self.sequences]
        return self.kmer_profiles
    
    def calculate_dissimilarity(self, profile1, profile2, method='bray_curtis'):
        """
        Calculate dissimilarity between two k-mer profiles.
        
        Parameters:
        -----------
        profile1, profile2 : Counter
            K-mer frequency profiles
        method : str, default='bray_curtis'
            Dissimilarity metric to use. Options: 'bray_curtis', 'jaccard', 'euclidean'
            
        Returns:
        --------
        float
            Dissimilarity value (0 = identical, 1 = completely different)
        """
        # Get all unique k-mers
        all_kmers = set(profile1.keys()) | set(profile2.keys())
        
        if method == 'bray_curtis':
            # Bray-Curtis dissimilarity
            numerator = sum(abs(profile1.get(k, 0) - profile2.get(k, 0)) for k in all_kmers)
            denominator = sum(profile1.get(k, 0) + profile2.get(k, 0) for k in all_kmers)
            return numerator / denominator if denominator > 0 else 0
            
        elif method == 'jaccard':
            # Jaccard dissimilarity
            intersection = len(set(profile1.keys()) & set(profile2.keys()))
            union = len(set(profile1.keys()) | set(profile2.keys()))
            return 1 - (intersection / union) if union > 0 else 0
            
        elif method == 'euclidean':
            # Normalized Euclidean distance
            vec1 = np.array([profile1.get(k, 0) for k in all_kmers])
            vec2 = np.array([profile2.get(k, 0) for k in all_kmers])
            distance = np.sqrt(np.sum((vec1 - vec2) ** 2))
            # Normalize by maximum possible distance
            max_distance = np.sqrt(np.sum(vec1 ** 2) + np.sum(vec2 ** 2))
            return distance / max_distance if max_distance > 0 else 0
        
        else:
            raise ValueError(f"Unknown method: {method}")
    
    def compute_dissimilarity_matrix(self, method='bray_curtis'):
        """
        Compute pairwise dissimilarity matrix for all sequences.
        
        Parameters:
        -----------
        method : str, default='bray_curtis'
            Dissimilarity metric to use
            
        Returns:
        --------
        numpy.ndarray
            n x n dissimilarity matrix
        """
        if self.kmer_profiles is None:
            self.build_kmer_profiles()
        
        n = len(self.sequences)
        self.dissimilarity_matrix = np.zeros((n, n))
        
        for i in range(n):
            for j in range(i+1, n):
                dissim = self.calculate_dissimilarity(
                    self.kmer_profiles[i], 
                    self.kmer_profiles[j], 
                    method
                )
                self.dissimilarity_matrix[i, j] = dissim
                self.dissimilarity_matrix[j, i] = dissim
        
        return self.dissimilarity_matrix
    
    def plot_heatmap(self, method='bray_curtis', figsize=(10, 8), cmap='YlOrRd', 
                     save_path=None, title=None):
        """
        Generate and display a heatmap of the dissimilarity matrix.
        
        Parameters:
        -----------
        method : str, default='bray_curtis'
            Dissimilarity metric to use
        figsize : tuple, default=(10, 8)
            Figure size
        cmap : str, default='YlOrRd'
            Colormap to use
        save_path : str, optional
            Path to save the figure. If None, displays only.
        title : str, optional
            Custom title for the plot
        """
        if self.dissimilarity_matrix is None:
            self.compute_dissimilarity_matrix(method)
        
        plt.figure(figsize=figsize)
        
        # Create heatmap
        sns.heatmap(
            self.dissimilarity_matrix,
            xticklabels=self.labels,
            yticklabels=self.labels,
            cmap=cmap,
            annot=True,
            fmt='.3f',
            square=True,
            cbar_kws={'label': 'Dissimilarity'}
        )
        
        plot_title = title if title else f'16S RNA Dissimilarity Heatmap (k={self.k}, {method})'
        plt.title(plot_title, fontsize=14, fontweight='bold')
        plt.xlabel('Sequence', fontsize=12)
        plt.ylabel('Sequence', fontsize=12)
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Heatmap saved to {save_path}")
        
        plt.show()
    
    def plot_histogram(self, method='bray_curtis', bins=20, figsize=(10, 6), 
                       save_path=None, title=None):
        """
        Generate and display a histogram of dissimilarity values.
        
        Parameters:
        -----------
        method : str, default='bray_curtis'
            Dissimilarity metric to use
        bins : int, default=20
            Number of histogram bins
        figsize : tuple, default=(10, 6)
            Figure size
        save_path : str, optional
            Path to save the figure. If None, displays only.
        title : str, optional
            Custom title for the plot
        """
        if self.dissimilarity_matrix is None:
            self.compute_dissimilarity_matrix(method)
        
        # Extract upper triangle values (excluding diagonal)
        n = self.dissimilarity_matrix.shape[0]
        dissimilarities = []
        for i in range(n):
            for j in range(i+1, n):
                dissimilarities.append(self.dissimilarity_matrix[i, j])
        
        plt.figure(figsize=figsize)
        
        # Create histogram
        plt.hist(dissimilarities, bins=bins, color='steelblue', alpha=0.7, edgecolor='black')
        
        # Add statistics
        mean_dissim = np.mean(dissimilarities)
        median_dissim = np.median(dissimilarities)
        plt.axvline(mean_dissim, color='red', linestyle='--', linewidth=2, 
                   label=f'Mean: {mean_dissim:.3f}')
        plt.axvline(median_dissim, color='green', linestyle='--', linewidth=2, 
                   label=f'Median: {median_dissim:.3f}')
        
        plot_title = title if title else f'Distribution of Dissimilarities (k={self.k}, {method})'
        plt.title(plot_title, fontsize=14, fontweight='bold')
        plt.xlabel('Dissimilarity', fontsize=12)
        plt.ylabel('Frequency', fontsize=12)
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Histogram saved to {save_path}")
        
        plt.show()
        
        return dissimilarities
    
    def analyze_kmer_sizes(self, k_values, method='bray_curtis'):
        """
        Compare dissimilarity analysis across different k-mer sizes.
        
        Parameters:
        -----------
        k_values : list of int
            List of k-mer sizes to analyze
        method : str, default='bray_curtis'
            Dissimilarity metric to use
            
        Returns:
        --------
        dict
            Dictionary mapping k values to dissimilarity matrices
        """
        results = {}
        original_k = self.k
        
        for k in k_values:
            self.k = k
            self.kmer_profiles = None
            self.dissimilarity_matrix = None
            matrix = self.compute_dissimilarity_matrix(method)
            results[k] = matrix.copy()
        
        # Restore original k
        self.k = original_k
        self.kmer_profiles = None
        self.dissimilarity_matrix = None
        
        return results
    
    def compare_kmer_sizes(self, k_values, method='bray_curtis', figsize=(15, 5)):
        """
        Visualize comparison of dissimilarity matrices for different k-mer sizes.
        
        Parameters:
        -----------
        k_values : list of int
            List of k-mer sizes to compare
        method : str, default='bray_curtis'
            Dissimilarity metric to use
        figsize : tuple, default=(15, 5)
            Figure size
        """
        results = self.analyze_kmer_sizes(k_values, method)
        
        fig, axes = plt.subplots(1, len(k_values), figsize=figsize)
        if len(k_values) == 1:
            axes = [axes]
        
        for idx, k in enumerate(k_values):
            sns.heatmap(
                results[k],
                xticklabels=self.labels,
                yticklabels=self.labels,
                cmap='YlOrRd',
                annot=True,
                fmt='.2f',
                square=True,
                ax=axes[idx],
                cbar_kws={'label': 'Dissimilarity'}
            )
            axes[idx].set_title(f'k={k}', fontsize=12, fontweight='bold')
        
        plt.suptitle(f'Comparison of K-mer Sizes ({method})', 
                    fontsize=14, fontweight='bold', y=1.02)
        plt.tight_layout()
        plt.show()


def main():
    """
    Example usage demonstrating the RNA Dissimilarity Analysis tool.
    """
    print("16S RNA Dissimilarity Analysis Tool")
    print("=" * 50)
    
    # Example 16S RNA sequences (simplified for demonstration)
    sequences = [
        "AGAGUUUGAUCCUGGGCUCAGAUUGAACGCUGGCGGCAGGCCUAACACAUGCAAGUCGAACGGUAACAGGAAGAAGCUUGCUUCUUUGCUGACGAGUG",
        "AGAGAUUGAUCAUGGCUCAGAUUGAACGCUGGCGGCAGGCCUAACACAUGCAAGUCGAGCGGCAGCACAGAGGAACUUGCUUCUUUGCUGACGAGCG",
        "AGAGUUUGAUCCUGGGCUCAGAUUGAACGCUGGCGGCAGGCCUAACACAUGCAAGUCGAGCGGCAGCAUAGGGAGCUUGCUUCUUUGCUGACGAGCG",
        "UUGGCGGCAGGCCUAACACAUGCAAGUCGAACGAUACCUUUGCUGACGAGUAGCUAAUGAGGAAGAAGCUUGCUUCUUUGCUGACGAGCGGCGGACG",
        "AGAGUUUGAUCAUGGCUCAGAUUGAACGCUGGCGGCAUGCCUAACACAUGCAAGUCGAGCGGCAGCACAGAGGAACUUGCCAGAGUUUUGACGAGCG"
    ]
    
    labels = ["E.coli_1", "E.coli_2", "Salmonella", "Shigella", "Klebsiella"]
    
    # Analyze with k=3
    print("\nAnalyzing sequences with k=3...")
    analyzer = RNADissimilarityAnalyzer(sequences, labels, k=3)
    
    # Compute and display dissimilarity matrix
    print("\nComputing dissimilarity matrix...")
    matrix = analyzer.compute_dissimilarity_matrix(method='bray_curtis')
    print("\nDissimilarity Matrix:")
    print(matrix)
    
    # Generate heatmap
    print("\nGenerating heatmap...")
    analyzer.plot_heatmap(save_path='dissimilarity_heatmap.png')
    
    # Generate histogram
    print("\nGenerating histogram...")
    dissimilarities = analyzer.plot_histogram(save_path='dissimilarity_histogram.png')
    
    # Compare different k-mer sizes
    print("\nComparing different k-mer sizes (k=2, 3, 4)...")
    analyzer.compare_kmer_sizes([2, 3, 4])
    
    print("\n" + "=" * 50)
    print("Analysis complete!")
    print(f"Mean dissimilarity: {np.mean(dissimilarities):.3f}")
    print(f"Median dissimilarity: {np.median(dissimilarities):.3f}")
    print(f"Min dissimilarity: {np.min(dissimilarities):.3f}")
    print(f"Max dissimilarity: {np.max(dissimilarities):.3f}")


if __name__ == "__main__":
    main()
