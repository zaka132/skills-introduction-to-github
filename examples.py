#!/usr/bin/env python3
"""
Example script demonstrating various use cases of the RNA Dissimilarity Analysis tool.
"""

from rna_dissimilarity import RNADissimilarityAnalyzer
import numpy as np


def example_1_basic_analysis():
    """
    Example 1: Basic dissimilarity analysis with sample sequences
    """
    print("\n" + "="*60)
    print("EXAMPLE 1: Basic Dissimilarity Analysis")
    print("="*60)
    
    # Sample 16S RNA sequences (partial sequences for demonstration)
    sequences = [
        "AGAGUUUGAUCMUGGGCUCAGAUUGAACGCUGGCGGCAGGCCUAACACAUGCAAGUCGAACGGUAACAGGAAGAAGCUUGCUUCUUUGCUGACGAGUG",
        "AGAGAUUGAUCAUGGCUCAGAUUGAACGCUGGCGGCAGGCCUAACACAUGCAAGUCGAGCGGCAGCACAGAGGAACUUGCUUCUUUGCUGACGAGCG",
        "AGAGUUUGAUCMUGGGCUCAGAUUGAACGCUGGCGGCAGGCCUAACACAUGCAAGUCGAGCGGCAGCAUAGGGAGCUUGCUUCUUUGCUGACGAGCG",
    ]
    
    labels = ["Organism_A", "Organism_B", "Organism_C"]
    
    # Create analyzer
    analyzer = RNADissimilarityAnalyzer(sequences, labels=labels, k=3)
    
    # Compute dissimilarity matrix
    matrix = analyzer.compute_dissimilarity_matrix(method='bray_curtis')
    
    print("\nDissimilarity Matrix (Bray-Curtis):")
    print(matrix)
    
    # Generate visualizations
    analyzer.plot_heatmap(title="Example 1: Basic Analysis Heatmap")
    analyzer.plot_histogram(title="Example 1: Dissimilarity Distribution")


def example_2_compare_metrics():
    """
    Example 2: Compare different dissimilarity metrics
    """
    print("\n" + "="*60)
    print("EXAMPLE 2: Comparing Different Metrics")
    print("="*60)
    
    sequences = [
        "AGAGUUUGAUCMUGGGCUCAGAUUGAACGCUGGCGGCAGGCCUAACACAUGCAAGUCGAACGGUAACAGGAAGAAGCUUGCUUCUUUGCUGACGAGUG",
        "AGAGAUUGAUCAUGGCUCAGAUUGAACGCUGGCGGCAGGCCUAACACAUGCAAGUCGAGCGGCAGCACAGAGGAACUUGCUUCUUUGCUGACGAGCG",
        "AGAGUUUGAUCMUGGGCUCAGAUUGAACGCUGGCGGCAGGCCUAACACAUGCAAGUCGAGCGGCAGCAUAGGGAGCUUGCUUCUUUGCUGACGAGCG",
        "UUGGCGGCAGGCCUAACACAUGCAAGUCGAACGAUACCUUUGCUGACGAGUAGCUAAUGAGGAAGAAGCUUGCUUCUUUGCUGACGAGCGGCGGACG",
    ]
    
    labels = ["Sample_1", "Sample_2", "Sample_3", "Sample_4"]
    analyzer = RNADissimilarityAnalyzer(sequences, labels=labels, k=3)
    
    metrics = ['bray_curtis', 'jaccard', 'euclidean']
    
    for metric in metrics:
        print(f"\n{metric.upper()} Dissimilarity:")
        matrix = analyzer.compute_dissimilarity_matrix(method=metric)
        print(matrix)
        print(f"Mean: {np.mean(matrix[np.triu_indices_from(matrix, k=1)]):.3f}")


def example_3_kmer_optimization():
    """
    Example 3: Find optimal k-mer size
    """
    print("\n" + "="*60)
    print("EXAMPLE 3: K-mer Size Optimization")
    print("="*60)
    
    sequences = [
        "AGAGUUUGAUCMUGGGCUCAGAUUGAACGCUGGCGGCAGGCCUAACACAUGCAAGUCGAACGGUAACAGGAAGAAGCUUGCUUCUUUGCUGACGAGUG",
        "AGAGAUUGAUCAUGGCUCAGAUUGAACGCUGGCGGCAGGCCUAACACAUGCAAGUCGAGCGGCAGCACAGAGGAACUUGCUUCUUUGCUGACGAGCG",
        "AGAGUUUGAUCMUGGGCUCAGAUUGAACGCUGGCGGCAGGCCUAACACAUGCAAGUCGAGCGGCAGCAUAGGGAGCUUGCUUCUUUGCUGACGAGCG",
        "UUGGCGGCAGGCCUAACACAUGCAAGUCGAACGAUACCUUUGCUGACGAGUAGCUAAUGAGGAAGAAGCUUGCUUCUUUGCUGACGAGCGGCGGACG",
        "AGAGUUUGAUCAUGGCUCAGAUUGAACGCUGGCGGCAUGCCUAACACAUGCAAGUCGAGCGGCAGCACAGAGGAACUUGCCAGAGUUUUGACGAGCG"
    ]
    
    labels = ["Seq_1", "Seq_2", "Seq_3", "Seq_4", "Seq_5"]
    analyzer = RNADissimilarityAnalyzer(sequences, labels=labels, k=3)
    
    # Test different k values
    k_values = [2, 3, 4, 5]
    results = analyzer.analyze_kmer_sizes(k_values)
    
    print("\nMean dissimilarity for different k values:")
    for k, matrix in results.items():
        mean_dissim = np.mean(matrix[np.triu_indices_from(matrix, k=1)])
        print(f"  k={k}: {mean_dissim:.3f}")
    
    # Visualize comparison
    analyzer.compare_kmer_sizes(k_values)


def example_4_large_dataset():
    """
    Example 4: Analysis with more sequences
    """
    print("\n" + "="*60)
    print("EXAMPLE 4: Larger Dataset Analysis")
    print("="*60)
    
    # Simulate variations of a base sequence
    base_seq = "AGAGUUUGAUCMUGGGCUCAGAUUGAACGCUGGCGGCAGGCCUAACACAUGCAAGUCGAACGGUAACAGGAAGAAGCUUGCUUCUUUGCUGACGAGUG"
    
    sequences = [
        base_seq,
        base_seq.replace("AACGG", "AUCGG"),  # Single variation
        base_seq.replace("CAGGCC", "CAAGCC"),  # Different variation
        base_seq.replace("UGACG", "UGCCG"),  # Another variation
        base_seq.replace("GAAGAAG", "GAAGCAG"),  # Yet another
        base_seq.replace("AACGCUGG", "AACGAUGG"),  # More variations
        base_seq[:80],  # Truncated sequence
        base_seq[20:],  # Different region
    ]
    
    labels = [f"Variant_{i+1}" for i in range(len(sequences))]
    
    analyzer = RNADissimilarityAnalyzer(sequences, labels=labels, k=3)
    
    # Compute and visualize
    matrix = analyzer.compute_dissimilarity_matrix()
    
    print(f"\nAnalyzing {len(sequences)} sequences...")
    print(f"Matrix shape: {matrix.shape}")
    
    analyzer.plot_heatmap(
        figsize=(12, 10),
        title="Example 4: Large Dataset Heatmap"
    )
    
    dissimilarities = analyzer.plot_histogram(
        bins=30,
        title="Example 4: Distribution of Dissimilarities"
    )
    
    print(f"\nStatistics:")
    print(f"  Total pairwise comparisons: {len(dissimilarities)}")
    print(f"  Mean dissimilarity: {np.mean(dissimilarities):.3f}")
    print(f"  Std deviation: {np.std(dissimilarities):.3f}")
    print(f"  Range: [{np.min(dissimilarities):.3f}, {np.max(dissimilarities):.3f}]")


def example_5_custom_analysis():
    """
    Example 5: Custom analysis workflow
    """
    print("\n" + "="*60)
    print("EXAMPLE 5: Custom Analysis Workflow")
    print("="*60)
    
    # Define sequences representing different bacterial groups
    sequences = [
        # Enterobacteriaceae family
        "AGAGUUUGAUCMUGGGCUCAGAUUGAACGCUGGCGGCAGGCCUAACACAUGCAAGUCGAACGGUAACAGGAAGAAGCUUGCUUCUUUGCUGACGAGUG",
        "AGAGAUUGAUCAUGGCUCAGAUUGAACGCUGGCGGCAGGCCUAACACAUGCAAGUCGAGCGGCAGCACAGAGGAACUUGCUUCUUUGCUGACGAGCG",
        # Different family
        "UUGGCGGCAGGCCUAACACAUGCAAGUCGAACGAUACCUUUGCUGACGAGUAGCUAAUGAGGAAGAAGCUUGCUUCUUUGCUGACGAGCGGCGGACG",
        "AGAGUUUGAUCMUGGGCUCAGAUUGAACGCUGGCGGCAGGCCUAACACAUGCAAGUCGAGCGGCAGCAUAGGGAGCUUGCUUCUUUGCUGACGAGCG",
    ]
    
    labels = ["E.coli", "Salmonella", "Bacillus", "Shigella"]
    
    # Analyze with k=4 for more specificity
    analyzer = RNADissimilarityAnalyzer(sequences, labels=labels, k=4)
    
    # Build k-mer profiles
    profiles = analyzer.build_kmer_profiles()
    
    print("\nK-mer Profile Summary:")
    for i, (label, profile) in enumerate(zip(labels, profiles)):
        print(f"  {label}: {len(profile)} unique 4-mers")
    
    # Compute dissimilarity with Bray-Curtis
    matrix = analyzer.compute_dissimilarity_matrix(method='bray_curtis')
    
    # Find most similar and most dissimilar pairs
    n = len(sequences)
    max_dissim = 0
    min_dissim = 1
    max_pair = None
    min_pair = None
    
    for i in range(n):
        for j in range(i+1, n):
            dissim = matrix[i, j]
            if dissim > max_dissim:
                max_dissim = dissim
                max_pair = (labels[i], labels[j])
            if dissim < min_dissim:
                min_dissim = dissim
                min_pair = (labels[i], labels[j])
    
    print(f"\nMost similar pair: {min_pair[0]} - {min_pair[1]} (dissimilarity: {min_dissim:.3f})")
    print(f"Most dissimilar pair: {max_pair[0]} - {max_pair[1]} (dissimilarity: {max_dissim:.3f})")
    
    # Generate final visualizations
    analyzer.plot_heatmap(
        cmap='coolwarm',
        title="Example 5: Custom Analysis"
    )


def main():
    """
    Run all examples
    """
    print("\n" + "#"*60)
    print("# RNA DISSIMILARITY ANALYSIS - EXAMPLES")
    print("#"*60)
    
    examples = [
        ("Basic Analysis", example_1_basic_analysis),
        ("Compare Metrics", example_2_compare_metrics),
        ("K-mer Optimization", example_3_kmer_optimization),
        ("Large Dataset", example_4_large_dataset),
        ("Custom Workflow", example_5_custom_analysis),
    ]
    
    print("\nAvailable examples:")
    for i, (name, _) in enumerate(examples, 1):
        print(f"  {i}. {name}")
    
    print("\nRunning all examples...\n")
    
    for name, example_func in examples:
        try:
            example_func()
        except Exception as e:
            print(f"\nError in {name}: {str(e)}")
    
    print("\n" + "#"*60)
    print("# All examples completed!")
    print("#"*60)


if __name__ == "__main__":
    main()
