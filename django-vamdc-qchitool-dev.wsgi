import os
import sys

os.environ['DJANGO_SETTINGS_MODULE'] = 'oacagliari.settings'

import django.core.handlers.wsgi
application = django.core.handlers.wsgi.WSGIHandler()
path = '/var/www/wsgi/qchitool-dev'
if path not in sys.path:
    sys.path.append(path)
path = '/var/www/wsgi/qchitool-dev/oacagliari'
if path not in sys.path:
    sys.path.append(path)


