# Quick Start Guide - 16S RNA Dissimilarity Analysis

## Installation

```bash
# Install required dependencies
pip install -r requirements.txt
```

## Run the Demo

```bash
# Run the built-in example with sample bacterial sequences
python rna_dissimilarity.py
```

This will:
- Analyze 5 sample 16S RNA sequences
- Generate a dissimilarity heatmap (`dissimilarity_heatmap.png`)
- Generate a histogram of dissimilarities (`dissimilarity_histogram.png`)
- Compare k-mer sizes (k=2, 3, 4)
- Print statistical summary

## Quick Example

```python
from rna_dissimilarity import RNADissimilarityAnalyzer

# Your RNA sequences
sequences = [
    "AGAGUUUGAUCCUGGGCUCAGAUUGAACGCUGGCGGC...",
    "AGAGAUUGAUCAUGGCUCAGAUUGAACGCUGGCGGC...",
    "UUGGCGGCAGGCCUAACACAUGCAAGUCGAACGAU..."
]

# Create analyzer
analyzer = RNADissimilarityAnalyzer(sequences, k=3)

# Compute and visualize
analyzer.compute_dissimilarity_matrix()
analyzer.plot_heatmap(save_path='my_heatmap.png')
analyzer.plot_histogram(save_path='my_histogram.png')
```

## Run All Examples

```bash
# See 5 different examples demonstrating various use cases
python examples.py
```

## Documentation

See `RNA_ANALYSIS_README.md` for complete documentation including:
- Detailed API reference
- Advanced usage examples
- Explanation of dissimilarity metrics
- Troubleshooting tips

## Features at a Glance

✅ K-mer based dissimilarity analysis  
✅ Multiple dissimilarity metrics (Bray-Curtis, Jaccard, Euclidean)  
✅ Heatmap generation  
✅ Histogram visualization  
✅ K-mer size comparison  
✅ Comprehensive examples and documentation  

## Output Files

- `dissimilarity_heatmap.png` - Color-coded matrix of pairwise dissimilarities
- `dissimilarity_histogram.png` - Distribution of dissimilarity values

## Next Steps

1. Replace sample sequences with your own 16S RNA data
2. Experiment with different k-mer sizes (2-6 recommended)
3. Try different dissimilarity metrics
4. Customize visualization parameters

For questions and support, see the main README.
