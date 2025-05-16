name: py-mosaic
channels:
  - conda-forge
dependencies:
  - python=3.10
  - pip
  - pandas
  - numpy
  - seaborn
  - matplotlib
  - pip:
      - missionbio.mosaic
