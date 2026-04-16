#scripts for developing and benchmarking CMDTarget

## Background energy correction

See `bg-correction`

- `scripts/train-dimer-background.py`: learn interaction energy background with an dense neural network

## Hfq binding prediction

See `hfq`

- `hfq/RIL-seq.targets.and.CLIP.peak.wo.KP.fa`: the training data


## Performance evaluation and visualization
- Performance comparison between different settings: `plots-CMDTarget.ipynb`

- Benchmarking across 9 species: `other-clade-benchmark.ipynb`

- Visluzation of turnover in different clade: `plot-clade-variations.ipynb`
