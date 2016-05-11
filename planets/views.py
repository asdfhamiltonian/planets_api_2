from django.http import HttpResponse
from math import copysign, floor, pi, atan2, sqrt, cos, sin, acos, asin
from functools import partial
import datetime
from collections import OrderedDict
from planets.planets import *
from planets.luna import *
import json

def index(request):
	return HttpResponse("Hello, world. You're at the polls index.")

def position(request, latitude, longitude):
	latitude = float(latitude)
	longitude = float(longitude)
	location = Planets(latitude, longitude)
	now = OrderedDict()
	now["mercury"] = location.calcMercuryNow()
	now["venus"] = location.calcVenusNow()
	now["mars"] = location.calcMarsNow()
	now["jupiter"] = location.calcJupiterNow()
	now["saturn"] = location.calcSaturnNow()
	now["uranus"] = location.calcUranusNow()
	now["neptune"] = location.calcNeptuneNow()
	now["pluto"] = location.calcPlutoThatIsNotAPlanetNow()
	return HttpResponse(json.dumps(now))

def moon(self, latitude, longitude, timezone):
	latitude = float(latitude)
	longitude = float(longitude)
	timezone = float(timezone)
	location = Luna(latitude, longitude)
	moonList = location.now()
	moonDict = OrderedDict()
	moonDict["azimuth"] = moonList[0]
	moonDict["altitude"] = moonList[1]
	moonDict["right ascension"] = moonList[2]
	moonDict["declination"] = moonList[3]
	moonDict["phase"] = moonList[4]
	moonDict["topocentric right ascension"] = moonList[5]
	moonDict["topocentric declination"] = moonList[6]
	moonDict["mpar"] = moonList[7]
	moonDict["msd"] = moonList[8]
	moonDict["longitude degrees"] = moonList[9]
	moonDict["latitude degrees"] = moonList[10]
	moonDict["gmsto"] = moonList[11]
	moonDict["epoch day"] = moonList[12]
	moonDict["riseset1"] = location.risesetlist(1, timezone)[0:3]
	moonDict["riseset2"] = location.risesetlist(2, timezone)[0:3]
	moonDict["riseset3"] = location.risesetlist(3, timezone)[0:3]
	return HttpResponse(json.dumps(moonDict))

def riseset(self, latitude, longitude, timezone):
	latitude = float(latitude)
	longitude = float(longitude)
	timezone = float(timezone)
	location = Planets(latitude, longitude)
	riseSetDict = location.dictRiseSet(-7)
	return HttpResponse(json.dumps(riseSetDict))