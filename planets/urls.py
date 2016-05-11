from django.conf.urls import url

from . import views

urlpatterns = [
	url(r'^$', views.index, name='index'),
	url(r'position/(?P<latitude>\-?\d+(.\d+)?)/(?P<longitude>\-?\d+(.\d+)?)/$', 
		views.position, name='position'),
	url(r'riseset/(?P<latitude>\-?\d+(.\d+)?)/(?P<longitude>\-?\d+(.\d+)?)/(?P<timezone>\-?\d+(.\d+)?)/$', 
		views.riseset, name='position'),
	url(r'moon/(?P<latitude>\-?\d+(.\d+)?)/(?P<longitude>\-?\d+(.\d+)?)/(?P<timezone>\-?\d+(.\d+)?)/$', 
		views.moon, name='position'),
]