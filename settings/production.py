from settings import *

DEBUG = False
TEMPLATES[0]['OPTIONS']['debug'] = DEBUG

DATABASES = {
    'default': {
        'ENGINE':'django.db.backends.postgresql_psycopg2',
        'NAME': 'specpipe',
        'USER': 'postgres',
        'PASSWORD': 'postgres',
        'HOST': 'ml-db-dev.cvsjiwtmbmak.us-east-1.rds.amazonaws.com',
        'PORT': '5432',
    }
}
