# scMyeloids_ICH
Characterizing the immune cell landscape in intracerebral hemorrhage stroke

# Contents

```
├── Figures
│   ├── Figure1
│   │   ├── 1.Figure1B_CCAmetadata.ipynb
│   │   ├── 1.Figure1C_UMAP.ipynb
│   │   ├── 1.Figure1D_CT_Markers.ipynb
│   │   ├── 1.Figure1E_CrumblR.ipynb
│   │   └── 1.Figure1F_Dreamlet.ipynb
│   ├── Figure2
│   │   ├── 2.Figure2_A_UMAP.ipynb
│   │   ├── 2.Figure2_B_scRDS.ipynb
│   │   ├── 2.Figure2_C_CrumblR.ipynb
│   │   ├── 2.Figure2_D_CL_DEGs.ipynb
│   │   └── 2.Figure2_E_GRN.R
│   ├── Figure3
│   │   ├── 3.Figure3_A_Liana.ipynb
│   │   ├── 3.Figure3_B_CellChat.ipynb
│   │   └── 3.Figure3_C_IREA.ipynb
│   ├── SFigure1
│   │   └── SFigure1.R
│   ├── SFigure2
│   │   ├── SFigure2_A_Co_FreshMG_ICH.ipynb
│   │   └── SFigure2_BCD.ipynb
│   ├── SFigure3
│   │   ├── SFigure_3A.ipynb
│   │   ├── SFigure_3B.ipynb
│   │   └── SFigure_3C.ipynb
│   └── SFigure4
│       └── SFigure4.R
├── README.md
├── Tables
└── scripts
    ├── 0.Functions.R
    ├── 0.Preprocess.R
    ├── 1.Merge_datasets.R
    ├── 2.Pegasus_Functions.py
    ├── 2.Pegasus_run.py
    ├── 3.Manual_Major_Annotation.py
    ├── 4.Transfer_MAnnot_To_Seurat.R
    ├── 5.Ref_SCANVI_Annot.py
    ├── 6.Metacells.R
    ├── 7.DEG_Metacells.R
    ├── 8.scDRS.ipynb
    ├── 9.Liana.ipynb
    └── utils.R
```
    

# System requirements
```
crumblr_0.99.6  
metafor_4.2-0  
numDeriv_2016.8-1.1  
metadat_1.2-0  
Matrix_1.5-4.1  
broom_1.0.5  
muscat_1.14.0  
ComplexHeatmap_2.16.0  
circlize_0.4.15  
aplot_0.1.10  
ggtree_3.8.0  
doParallel_1.0.17  
iterators_1.0.14  
foreach_1.5.2  
scater_1.28.0  
scuttle_1.10.1  
SingleCellExperiment_1.22.0  
SummarizedExperiment_1.30.2  
GenomicRanges_1.52.0  
GenomeInfoDb_1.36.3  
MatrixGenerics_1.12.3  
matrixStats_1.0.0  
dreamlet_0.99.16  
variancePartition_1.31.9  
BiocParallel_1.34.2  
limma_3.56.2  
SCopeLoomR_0.13.0  
scCustomize_1.1.1  
hdWGCNA_0.2.19  
WGCNA_1.72-1  
fastcluster_1.2.3  
dynamicTreeCut_1.63-1  
clusterProfiler_4.8.1  
lubridate_1.9.2  
forcats_1.0.0  
stringr_1.5.0  
purrr_1.0.1  
tidyr_1.3.0  
tibble_3.2.1  
tidyverse_2.0.0  
org.Hs.eg.db_3.17.0  
AnnotationDbi_1.62.2  
IRanges_2.34.1  
S4Vectors_0.38.1  
Biobase_2.60.0  
BiocGenerics_0.46.0  
RColorBrewer_1.1-3  
dittoSeq_1.12.0  
harmony_0.1.1  
Rcpp_1.0.11  
igraph_1.5.0  
ggplot2_3.5.0  
patchwork_1.2.0  
cowplot_1.1.1  
readr_2.1.4  
Seurat_4.9.9.9049  
SeuratObject_4.9.9.9086  
sp_2.0-0  
dplyr_1.1.2
# ==== Python
diff-match-patch==20200713
xarray==2024.1.0
tzdata==2023.4
certifi==2022.12.7
fsspec==2022.11.0
regex==2022.7.9
pytz==2022.7
dask==2022.7.0
distributed==2022.7.0
imagecodecs==2021.8.26
tifffile==2021.7.2
setuptools==65.6.3
cryptography==39.0.1
keyring==23.4.0
conda==23.3.1
pyzmq==23.2.0
pyopenssl==23.0.0
boltons==23.0.0
black==22.6.0
pip==22.3.1
twisted==22.2.0
attrs==22.1.0
packaging==22.0
contextlib2==21.6.0
incremental==21.3.0
argon2-cffi==21.3.0
argon2-cffi-bindings==21.2.0
hyperlink==21.0.0
virtualenv==20.25.0
automat==20.2.0
service-identity==18.1.0
constantly==15.1.0
rich==13.7.0
websockets==12.0
pillow==9.4.0
pyobjc-framework-coreservices==9.0
pyobjc-framework-cocoa==9.0
pyobjc-core==9.0
pyobjc-framework-fsevents==9.0
ipython==8.10.0
natsort==8.4.0
click==8.0.4
tenacity==8.0.1
jupyter-client==7.3.4
pytest==7.1.2
inflect==7.0.0
ipykernel==6.19.2
deepdiff==6.7.1
sip==6.6.2
nbconvert==6.5.4
notebook==6.5.2
pydocstyle==6.3.0
importlib-resources==6.1.1
tornado==6.1
multidict==6.0.4
flake8==6.0.0
pyyaml==6.0
isort==5.9.3
plotly==5.9.0
psutil==5.9.0
traitlets==5.7.1
nbformat==5.7.0
spyder==5.4.1
qtconsole==5.4.0
ujson==5.4.0
zope.interface==5.4.0
smart-open==5.2.1
jupyter-core==5.2.0
tzlocal==5.2
decorator==5.1.1
astropy==5.1
python-slugify==5.0.2
sphinx==5.0.2
gdown==5.0.0
tqdm==4.64.1
protobuf==4.25.2
fonttools==4.25.0
transformers==4.24.0
jsonschema==4.17.3
importlib-metadata==4.11.3
beautifulsoup4==4.11.1
lxml==4.9.1
typing-extensions==4.9.0
pexpect==4.8.0
hdf5plugin==4.3.0
gensim==4.3.0
textdistance==4.2.1
platformdirs==4.1.0
ordered-set==4.1.0
bleach==4.1.0
pyodbc==4.0.34
readchar==4.0.5
mock==4.0.3
async-timeout==4.0.3
chardet==4.0.0
conda-build==3.24.0
filelock==3.13.1
ply==3.11
zipp==3.11.0
aiohttp==3.9.1
tables==3.7.0
matplotlib==3.7.0
nltk==3.7
h5py==3.7.0
flit-core==3.6.0
pre-commit==3.6.0
rpy2==3.5.14
jupyterlab==3.5.3
anyio==3.5.0
conda-verify==3.4.2
markdown==3.4.1
idna==3.4
cfgv==3.4.0
opt-einsum==3.3.0
inquirer==3.2.1
tldextract==3.2.0
bcrypt==3.2.0
lz4==3.1.3
jinja2==3.1.2
intervaltree==3.1.0
prompt-toolkit==3.0.36
openpyxl==3.0.10
pyparsing==3.0.9
cython==3.0.8
wurlitzer==3.0.2
qdarkstyle==3.0.2
colorcet==3.0.1
pyflakes==3.0.1
markdown-it-py==3.0.0
requests==2.28.1
imageio==2.26.0
pycparser==2.21
jupyterlab-server==2.19.0
pygments==2.17.2
pylint==2.16.2
fastjsonschema==2.16.2
astroid==2.14.2
babel==2.11.0
pycodestyle==2.10.0
libarchive-c==2.9
numexpr==2.8.4
networkx==2.8.4
python-dateutil==2.8.2
scrapy==2.8.0
pybind11==2.6.1
identify==2.5.33
bokeh==2.4.3
spyder-kernels==2.4.1
anaconda-navigator==2.4.0
pydantic-core==2.4.0
pyjwt==2.4.0
sortedcontainers==2.4.0
soupsieve==2.3.2.post1
pylint-venv==2.3.0
termcolor==2.3.0
werkzeug==2.2.2
flask==2.2.2
backoff==2.2.1
qtpy==2.2.0
threadpoolctl==2.2.0
pandas==2.2.0
snowballstemmer==2.2.0
watchdog==2.1.6
pytorch-lightning==2.1.3
gmpy2==2.1.2
pydantic==2.1.1
markupsafe==2.1.1
zict==2.1.0
jsonpointer==2.1
absl-py==2.1.0
lightning==2.0.9.post0
pydispatcher==2.0.5
asttokens==2.0.5
charset-normalizer==2.0.4
pyviz-comms==2.0.2
pyhamcrest==2.0.2
conda-package-handling==2.0.2
greenlet==2.0.1
xlrd==2.0.1
tomli==2.0.1
itsdangerous==2.0.1
pyerfa==2.0.0
sphinxcontrib-htmlhelp==2.0.0
cloudpickle==2.0.0
biopython==1.81
botocore==1.34.23
boto3==1.34.23
jsonpatch==1.32
urllib3==1.26.14
numpy==1.23.5
jupyter-server==1.23.4
w3lib==1.21.0
blessed==1.20.0
annoy==1.17.3
six==1.16.0
holoviews==1.15.4
cffi==1.15.1
wrapt==1.14.1
param==1.12.3
torch==1.12.1
anaconda-client==1.11.2
sympy==1.11.1
py==1.11.0
scipy==1.10.0
scanpy==1.9.6
pkginfo==1.9.6
yarl==1.9.4
pyro-ppl==1.8.6
xmod==1.8.1
send2trash==1.8.0
nodeenv==1.8.0
cookiecutter==1.7.3
pep8==1.7.1
pysocks==1.7.1
python-lsp-server==1.7.1
rope==1.7.0
texttable==1.7.0
tblib==1.7.0
editor==1.6.5
backports.functools-lru-cache==1.6.4
celltypist==1.6.2
parsel==1.6.0
etils==1.6.0
lazy-object-proxy==1.6.0
autopep8==1.6.0
nest-asyncio==1.5.6
debugpy==1.5.1
requests-file==1.5.1
queuelib==1.5.0
pandocfilters==1.5.0
decoupler==1.5.0
numpydoc==1.5.0
sqlalchemy==1.4.39
appdirs==1.4.4
kiwisolver==1.4.4
pywavelets==1.4.1
imagesize==1.4.1
frozenlist==1.4.1
croniter==1.4.1
array-api-compat==1.4
pooch==1.4.0
atomicwrites==1.4.0
bottleneck==1.3.5
aiosignal==1.3.1
anndata2ri==1.3.1
torchmetrics==1.3.0.post0
text-unidecode==1.3
starsessions==1.3.0
pytoolconfig==1.2.5
arrow==1.2.3
clyent==1.2.2
qtawesome==1.2.2
tinycss2==1.2.1
python-lsp-black==1.2.1
scikit-learn==1.2.1
unidecode==1.2.0
partd==1.2.0
sniffio==1.2.0
runs==1.2.0
exceptiongroup==1.2.0
sphinxcontrib-serializinghtml==1.1.5
munkres==1.1.4
appscript==1.1.2
iniconfig==1.1.1
joblib==1.1.1
et-xmlfile==1.1.0
cssselect==1.1.0
conda-repo-cli==1.0.41
omnipath==1.0.8
contourpy==1.0.5
schpl==1.0.5
itemloaders==1.0.4
scvi-tools==1.0.4
liana==1.0.3
scdrs==1.0.3
msgpack==1.0.3
sphinxcontrib-qthelp==1.0.3
gseapy==1.0.3
sphinxcontrib-devhelp==1.0.2
sphinxcontrib-applehelp==1.0.2
whatthepatch==1.0.2
rtree==1.0.1
pathlib==1.0.1
sphinxcontrib-jsmath==1.0.1
heapdict==1.0.1
backports.weakref==1.0.post1
python-lsp-jsonrpc==1.0.0
backports.tempfile==1.0
locket==1.0.0
pluggy==1.0.0
newick==1.0.0
session-info==1.0.0
fastapi==0.109.0
websocket-client==0.58.0
numba==0.56.4
llvmlite==0.39.1
wheel==0.38.4
starlette==0.35.1
yapf==0.31.0
xlwings==0.29.1
uvicorn==0.26.0
scikit-image==0.19.3
zstandard==0.19.0
future==0.18.3
jedi==0.18.1
docutils==0.18.1
pyrsistent==0.18.0
ruamel-yaml-conda==0.17.21
ruamel.yaml==0.17.21
terminado==0.17.1
sparse==0.15.1
datashader==0.14.4
panel==0.14.3
prometheus-client==0.14.1
statsmodels==0.14.0
h11==0.14.0
numpyro==0.13.2
seaborn==0.13.2
plotnine==0.12.4
cytoolz==0.12.0
toolz==0.12.0
tokenizers==0.11.4
tomlkit==0.11.1
anaconda-project==0.11.1
docstring-to-markdown==0.11
cycler==0.11.0
igraph==0.10.8
pathspec==0.10.3
anndata==0.10.3
toml==0.10.2
huggingface-hub==0.10.1
imbalanced-learn==0.10.1
lightning-utilities==0.10.1
leidenalg==0.10.1
jmespath==0.10.0
stdlib-list==0.10.0
s3transfer==0.10.0
json5==0.9.6
mizani==0.9.3
requests-toolbelt==0.9.1
corneto==0.9.1a5
jellyfish==0.9.0
tabulate==0.8.10
kneed==0.8.5
mistune==0.8.4
executing==0.8.3
parso==0.8.3
hvplot==0.8.2
tensorly==0.8.1
alabaster==0.7.12
pickleshare==0.7.5
flax==0.7.5
cell2cell==0.7.3
defusedxml==0.7.1
brotlipy==0.7.0
ptyprocess==0.7.0
conda-package-streaming==0.7.0
glob2==0.7
mccabe==0.7.0
dateutils==0.6.12
intake==0.6.7
pycosat==0.6.4
python-snappy==0.6.1
multipledispatch==0.6.0
annotated-types==0.6.0
conda-pack==0.6.0
statannotations==0.6.0
lightning-cloud==0.5.61
nbclient==0.5.13
pynndescent==0.5.11
scarches==0.5.10
umap-learn==0.5.5
datashape==0.5.4
patsy==0.5.3
nbclassic==0.5.2
webencodings==0.5.1
inflection==0.5.1
pycirclize==0.5.1
orbax-checkpoint==0.5.0
pyct==0.5.0
poyo==0.5.0
fire==0.5.0
jax==0.4.23
jaxlib==0.4.23
pyasn1==0.4.8
colorama==0.4.6
pydeseq2==0.4.4
binaryornot==0.4.4
mypy-extensions==0.4.3
pyls-spyder==0.4.0
entrypoints==0.4
conda-token==0.4.0
distlib==0.3.7
dill==0.3.6
docrep==0.3.2
ml-dtypes==0.3.2
navigator-updater==0.3.0
scikit-misc==0.3.0
itemadapter==0.3.0
applaunchservices==0.3.0
pyasn1-modules==0.2.8
ruamel.yaml.clib==0.2.6
wcwidth==0.2.5
scrublet==0.2.3
mudata==0.2.3
pure-eval==0.2.2
qstylizer==0.2.2
notebook-shim==0.2.2
tbb==0.2
stack-data==0.2.0
backcall==0.2.0
ipython-genutils==0.2.0
jinja2-time==0.2.0
tensorstore==0.1.52
protego==0.1.16
optax==0.1.8
dm-tree==0.1.8
chex==0.1.7
matplotlib-inline==0.1.6
muon==0.1.5
conda-content-trust==0.1.3
appnope==0.1.2
mdurl==0.1.2
comm==0.1.2
jupyterlab-pygments==0.1.2
pyro-api==0.1.2
ml-collections==0.1.1
three-merge==0.1.1
python-multipart==0.0.6
pyqt5-sip==12.11.0
mpmath==1.2.1
pycurl==7.45.1
```

# License
MIT License