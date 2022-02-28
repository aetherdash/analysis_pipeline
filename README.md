## Analytics Database Interface Code with Additional Shared Analytics Utilities ##

Primarily for interfacing with the Analytics DB (specpipe) and use by the automated pipelines and dashboards.  Generic data analysis routines which could be used by notebooks, dashboards, or other repos/images may also be found in this repo.

Put code in this repo if:

* it modifies the structure (DDL commands) of the analytics database (postgres 'specpipe' DB)
* it makes a connection directly to the DB
* it inserts or deletes rows from any of the DB tables (DML commands)
* it provide generic methods for data querying/retrieval (DQL commands)
* it implements a GraphQL or SQLAlchemy layer for the DB
* it is to be shared across other repos or Docker images

Do not put code in this repo if:  

* it is large/complex enough to warrant its own repository (eg a purpose-specific docker module or large ML model)
* it is a generic machine learning (ML) tool or utility (put it in ml_tools)
* it is a generic bioinformatics tool or utility (put it in bio_scripts)
* it is a dashboard or a more complex front-end UI

Please consider putting an example of how to use your code in examples/.  That is a good place for notebooks

Tests policy: none yet

can be installed as a package using:
```
pip install --extra-index-url http://pypi.aetherbio.com/simple --trusted-host pypi.aetherbio.com aether-analytics-utils
```

and used as a normal module:
```
from analytics_tools.utils import retrieve_table
```

Maintained and originally populated by Charmaine February 2022.  Originally created by Louis Clark 2/18/2022.  

## Automated Pipeline Functionality ##

The associated BitBucket pipeline performs the following steps: 

1. automatically increments the aether-ml-tools module based on build number
2. performs a pipenv install
3. performs a make deploy
   * builds the module 
   * uploads module to http://pypi.aetherbio.com for consumption
