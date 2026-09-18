# Bioconda recipe handoff

`meta.yaml.template` is intentionally a template, not a submission-ready
recipe. Bioconda recipes pin the checksum of an immutable release archive,
which does not exist before MAGE-Bin's first release.

After publishing the matching version to PyPI:

1. download the PyPI source distribution and calculate its SHA-256 checksum;
2. replace all `REPLACE_WITH_...` fields and rename the file to `meta.yaml`;
3. add it as `recipes/magebin/meta.yaml` in a `bioconda-recipes` fork;
4. run `bioconda-utils lint recipes config.yml --packages magebin` and submit
   the recipe pull request.

The dependency name mapping is deliberate: PyPI's `igraph` is
`python-igraph` in conda, and PyPI's `torch` is `pytorch` in conda.
