# gene2dis/mgap: Contributing Guidelines

Thank you for your interest in contributing to gene2dis/mgap!

## Contribution workflow

1. Fork the repository
2. Create a new branch from `dev` for your changes
3. Make your changes
4. Run tests locally: `nextflow run . -profile test,docker --outdir ./results`
5. Run linting: `nf-core pipelines lint`
6. Submit a pull request to the `dev` branch

## Code style

- Follow [nf-core coding guidelines](https://nf-co.re/docs/contributing/code_formatting)
- Use meaningful variable and process names
- Add comments for complex logic
- Include docstrings in custom functions

## Adding new modules

If you want to add a new tool to the pipeline:

1. Check if the module already exists in [nf-core/modules](https://github.com/nf-core/modules)
2. If it exists, install it: `nf-core modules install <module_name>`
3. If not, create a local module following the nf-core module template
4. Add the module to the appropriate workflow file
5. Update documentation and tests

## Reporting bugs

Please use the GitHub issue tracker and provide:

- A clear description of the bug
- The command you used
- Relevant log files
- System information (Nextflow version, executor, container engine)

## Feature requests

Open a GitHub issue with:

- A clear description of the feature
- The scientific use case
- Any relevant references or examples
