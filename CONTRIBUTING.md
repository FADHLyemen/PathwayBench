# Contributing to PathwayBench

Thank you for your interest in contributing to PathwayBench.

## How to Contribute

### Reporting Issues

- Use [GitHub Issues](https://github.com/fadhlyemen/PathwayBench/issues) to report bugs or request features.
- Include the dataset, scoring method, and error log when reporting bugs.

### Pull Requests

1. Fork the repository.
2. Create a feature branch: `git checkout -b feature/your-feature`.
3. Make your changes and ensure the Snakemake pipeline still runs: `snakemake -n`.
4. Commit with clear messages.
5. Open a pull request against `main`.

### Adding a New Scoring Method

1. Implement the scoring function in `workflow/scripts/03_run_scoring.R`.
2. Add the method name to `config/config.yaml` under `methods`.
3. Update the evaluation script if the new method requires special handling.
4. Add tests with the example data in `tool/example_data/`.

### Adding a New Dataset

1. Add the dataset entry to `config/datasets.yaml` with the CellxGene query parameters.
2. Add ground-truth pathway directions to `config/pathways.yaml`.
3. Run the pipeline for the new dataset: `bash quickstart.sh your_dataset_id`.

## Code Style

- **R scripts**: Follow tidyverse style. Use `library()` calls at the top.
- **Python scripts**: Follow PEP 8. Use type hints where practical.
- **YAML configs**: Use 2-space indentation.

## License

By contributing, you agree that your contributions will be licensed under the MIT License.
