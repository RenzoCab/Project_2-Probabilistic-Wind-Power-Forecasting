import json
import os

class Test(object):
    def __init__(self, data):
	    self.__dict__ = json.load(data)
