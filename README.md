# ASC-ODE
A package for solving ordinary differential equations

Read the [documentation](https://tuwien-asc.github.io/ASC-ODE/)

Find theory behind here: https://jschoeberl.github.io/IntroSC/ODEs/ODEs.html

## Documentation

The documentation in `docs/` is a [Jupyter Book](https://jupyterbook.org) (version 2).
To preview it locally:

    pip install "jupyter-book>=2"
    cd docs
    jupyter-book start

On every push to `main`, the GitHub workflow `.github/workflows/docs.yml` builds it and publishes
it on GitHub Pages. Once per repository, set Settings -> Pages -> Source to "GitHub Actions";
in a fork, also enable workflows in the Actions tab.
