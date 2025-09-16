## Preliminar description of workflow

We will use this repository for all 3 different setups and will keep the key files here. We will try to use a bit advanced workflow utilising forking and pull requests.
Below are the steps required - after the initial cloning it is important to keep the main repository up to date via pull request and also keep up to date in your own work area with regular git pulls.

  1. Fork the *mews* repository
  2. Work in your private repository
  3. When ready - create a *pull request - PR*

The best way of working inside a fork and simultaneously updating from the central repo, is by setting up a remote.

On your PC, once:
1. git clone https://github.com/USERNAME/mews.git
2. git remote add upstream https://github.com/BoldingBruggeman/mews.git

On your PC, iteratively, in order to get the latest updates from the central repository on your fork:
1. git fetch upstream
2. git merge upstream/main
3. git push
(this should be done regularly, or at least before making PRs, to keep up with other contributions)

A full description of how to do this, for different OS and different programs to interact with github, can be found here: https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/fork-a-repo

## On the proper naming of data describing the depth of lakes

From https://en.wikipedia.org/wiki/Bathymetry:

Bathymetry (/bəˈθɪmətri/; from Ancient Greek βαθύς (bathús) 'deep', and μέτρον (métron) 'measure')[1][2] is the study of underwater depth of ocean floors (seabed topography), lake floors, or river floors. In other words, bathymetry is the underwater equivalent to hypsometry or topography. The first recorded evidence of water depth measurements are from Ancient Egypt over 3000 years ago.
