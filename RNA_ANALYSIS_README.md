# 16S RNA Dissimilarity Analysis Tool

A Python-based tool for performing dissimilarity analysis on 16S ribosomal RNA sequences using k-mer based methods. This tool enables researchers to compare bacterial sequences, generate heatmaps, and visualize dissimilarity distributions.

## Features

- **K-mer Based Analysis**: Extract and analyze k-mers of various sizes from 16S RNA sequences
- **Multiple Dissimilarity Metrics**: Support for Bray-Curtis, Jaccard, and Euclidean dissimilarity measures
- **Heatmap Generation**: Create publication-quality heatmaps of dissimilarity matrices
- **Histogram Visualization**: Visualize the distribution of pairwise dissimilarities
- **K-mer Size Comparison**: Compare results across different k-mer sizes
- **Flexible Input**: Accept any number of RNA sequences with custom labels

## Installation

### Prerequisites

- Python 3.7 or higher
- pip package manager

### Install Dependencies

```bash
pip install -r requirements.txt
```

The following packages will be installed:
- numpy (>=1.21.0)
- matplotlib (>=3.4.0)
- seaborn (>=0.11.0)

## Usage

### Basic Usage

```python
from rna_dissimilarity import RNADissimilarityAnalyzer

# Define your 16S RNA sequences
sequences = [
    "AGAGUUUGAUCMUGGGCUCAGAUUGAACGCUGGCGGCAGGCCUAACACAUGCAAGUCGAACGGUAACAGGAAGAAGCUUGCUUCUUUGCUGACGAGUG",
    "AGAGAUUGAUCAUGGCUCAGAUUGAACGCUGGCGGCAGGCCUAACACAUGCAAGUCGAGCGGCAGCACAGAGGAACUUGCUUCUUUGCUGACGAGCG",
    # Add more sequences...
]

# Create analyzer with k=3
analyzer = RNADissimilarityAnalyzer(sequences, k=3)

# Compute dissimilarity matrix
matrix = analyzer.compute_dissimilarity_matrix(method='bray_curtis')

# Generate visualizations
analyzer.plot_heatmap(save_path='heatmap.png')
analyzer.plot_histogram(save_path='histogram.png')
```

### Running the Example

The module includes a built-in example with sample bacterial sequences:

```bash
python rna_dissimilarity.py
```

This will:
1. Load example 16S RNA sequences from common bacteria
2. Compute dissimilarity matrices
3. Generate heatmaps and histograms
4. Compare different k-mer sizes (k=2, 3, 4)
5. Save output images

### Advanced Usage

#### Using Different Dissimilarity Metrics

```python
# Bray-Curtis dissimilarity (default)
analyzer.compute_dissimilarity_matrix(method='bray_curtis')

# Jaccard dissimilarity
analyzer.compute_dissimilarity_matrix(method='jaccard')

# Euclidean distance
analyzer.compute_dissimilarity_matrix(method='euclidean')
```

#### Comparing Multiple K-mer Sizes

```python
# Compare k=2, 3, 4, 5
analyzer.compare_kmer_sizes([2, 3, 4, 5])

# Get raw results for each k value
results = analyzer.analyze_kmer_sizes([2, 3, 4, 5])
```

#### Custom Sequence Labels

```python
sequences = ["AGAGUUUGAUCM...", "AGAGAUUGAUCA..."]
labels = ["Sample_A", "Sample_B"]

analyzer = RNADissimilarityAnalyzer(sequences, labels=labels, k=3)
```

## Dissimilarity Metrics

### Bray-Curtis Dissimilarity
- **Range**: 0 (identical) to 1 (completely different)
- **Best for**: Comparing abundance-based data
- **Formula**: BC = Σ|x_i - y_i| / Σ(x_i + y_i)

### Jaccard Dissimilarity
- **Range**: 0 (identical) to 1 (no overlap)
- **Best for**: Presence/absence comparisons
- **Formula**: J = 1 - (intersection / union)

### Euclidean Distance
- **Range**: 0 (identical) to 1 (normalized)
- **Best for**: Continuous numerical comparisons
- **Formula**: Normalized sqrt(Σ(x_i - y_i)²)

## Output Files

The tool can generate the following outputs:

1. **Heatmap** (`dissimilarity_heatmap.png`):
   - Color-coded matrix showing pairwise dissimilarities
   - Annotated with numerical values
   - Labeled with sequence identifiers

2. **Histogram** (`dissimilarity_histogram.png`):
   - Distribution of all pairwise dissimilarities
   - Mean and median lines for reference
   - Statistical summary

3. **K-mer Comparison** (displayed interactively):
   - Side-by-side heatmaps for different k values
   - Enables selection of optimal k-mer size

## API Reference

### RNADissimilarityAnalyzer

#### Constructor
```python
RNADissimilarityAnalyzer(sequences, labels=None, k=3)
```
- `sequences`: List of RNA sequence strings
- `labels`: Optional list of sequence identifiers
- `k`: K-mer size (default: 3)

#### Methods

##### `extract_kmers(sequence)`
Extract all k-mers from a sequence.
- **Returns**: Counter object with k-mer frequencies

##### `build_kmer_profiles()`
Build k-mer frequency profiles for all sequences.
- **Returns**: List of Counter objects

##### `compute_dissimilarity_matrix(method='bray_curtis')`
Compute pairwise dissimilarity matrix.
- **Parameters**: 
  - `method`: 'bray_curtis', 'jaccard', or 'euclidean'
- **Returns**: NumPy array (n×n matrix)

##### `plot_heatmap(method='bray_curtis', figsize=(10,8), cmap='YlOrRd', save_path=None, title=None)`
Generate heatmap visualization.
- **Parameters**:
  - `method`: Dissimilarity metric
  - `figsize`: Figure dimensions
  - `cmap`: Matplotlib colormap
  - `save_path`: File path to save image
  - `title`: Custom plot title

##### `plot_histogram(method='bray_curtis', bins=20, figsize=(10,6), save_path=None, title=None)`
Generate histogram of dissimilarities.
- **Parameters**:
  - `bins`: Number of histogram bins
  - `save_path`: File path to save image
- **Returns**: List of dissimilarity values

##### `compare_kmer_sizes(k_values, method='bray_curtis', figsize=(15,5))`
Visualize comparison across k-mer sizes.
- **Parameters**:
  - `k_values`: List of k values to compare

## Example Output

When running the example, you should see output similar to:

```
16S RNA Dissimilarity Analysis Tool
==================================================

Analyzing sequences with k=3...

Computing dissimilarity matrix...

Dissimilarity Matrix:
[[0.000 0.123 0.156 0.234 0.145]
 [0.123 0.000 0.089 0.267 0.078]
 [0.156 0.089 0.000 0.289 0.134]
 [0.234 0.267 0.289 0.000 0.298]
 [0.145 0.078 0.134 0.298 0.000]]

==================================================
Analysis complete!
Mean dissimilarity: 0.185
Median dissimilarity: 0.156
Min dissimilarity: 0.078
Max dissimilarity: 0.298
```

## Use Cases

- **Bacterial Classification**: Compare 16S sequences to identify taxonomic relationships
- **Microbiome Analysis**: Assess diversity within microbial communities
- **Quality Control**: Verify sequence similarity in cloning experiments
- **Phylogenetic Studies**: Preliminary analysis before tree construction
- **Method Optimization**: Determine optimal k-mer size for your dataset

## Technical Details

### K-mer Selection
- **Smaller k (2-3)**: More sensitive to local variations, higher resolution
- **Larger k (4-6)**: More specific, better for distant relationships
- **Recommendation**: Try k=3 or k=4 for 16S rRNA sequences

### Input Format
- Sequences can contain A, C, G, U (RNA) or T (DNA - automatically converted to U)
- Both uppercase and lowercase letters are accepted
- Sequences should be aligned to the same region for meaningful comparison

## Troubleshooting

**Issue**: ImportError for matplotlib or seaborn
- **Solution**: Ensure all dependencies are installed: `pip install -r requirements.txt`

**Issue**: All dissimilarity values are 0 or 1
- **Solution**: Check that sequences are different and k-mer size is appropriate

**Issue**: Memory error with large datasets
- **Solution**: Process sequences in batches or use a smaller k value

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use this tool in your research, please cite:

```
16S RNA Dissimilarity Analysis Tool
GitHub: https://github.com/zaka132/skills-introduction-to-github
```

## Contact

For questions or support, please open an issue on GitHub.

## Acknowledgments

This tool is designed for educational and research purposes in bioinformatics and microbial ecology.
